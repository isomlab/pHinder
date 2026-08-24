"""Run a pHinder calculation and report progress as it goes.

Lifted out of ``phinder_main_gui.main()``, which built the window and ran the
science in one 300-line closure. The calls and their order are unchanged; what
is new is that each calculation announces itself, output is routed to the
caller instead of to a globally-replaced ``sys.stdout``, and the run can be
asked to stop between steps.

One behavioural change, deliberate: the three shared prerequisites
(triangulation, surface, sidechain classification) are memoised. The original
guarded them with ``if not topologyCalculation:`` inside *each* dependent
calculation, so asking for both interface classification and virtual screening
without topology ran the whole triangulation twice.
"""

import os
import sys
import time
from os import sep

# Residue-set shorthands pHinder recognises; anything else becomes a customSet.
AA_SETS = {
    ('Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'His', 'Ile', 'Leu',
     'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val'): "allSet",
    ('Ala', 'Ile', 'Leu', 'Met', 'Phe', 'Pro', 'Trp', 'Val'): "apolarSet",
    ('Asn', 'Cys', 'Gln', 'Ser', 'Thr', 'Tyr'): "polarSet",
    ('Arg', 'Asp', 'Cys', 'Glu', 'His', 'Lys'): "ionizableSet",
    ('Arg', 'Asp', 'Glu', 'His', 'Lys'): "ionizableSetNoCys",
    ('Asp', 'Glu'): "acidicSet",
    ('Arg', 'His', 'Lys'): "basicSet",
}

# "Group Chains" rides along in the same dict as the chain checkboxes; it is an
# option, not a chain, and must never reach pHinderInstance.chains.
GROUP_CHAINS = "Group Chains"


# Each pHinder call the GUI makes, with what to show while it runs. Anything not
# listed still reports, just under its method name.
STEP_LABELS = {
    "setQuerySet": "Preparing the query set",
    "openPDBs": "Reading the structure",
    "hetLigand4D": "Locating ligands",
    "hydrogens": "Handling hydrogens",
    "makeAtomCollections": "Collecting atoms",
    "makeVertices4D": "Building vertices",
    "selectTscTriangulationAtoms": "Selecting triangulation atoms",
    "triangulateTscAtoms": "Triangulating sidechain atoms",
    "writeTriangulation": "Writing the triangulation",
    "pruneTriangulation": "Pruning the triangulation",
    "minimizePrunedTriangulation": "Minimising the pruned network",
    "identifyTightBonds": "Identifying tight bonds",
    "calculateNetworkParity": "Calculating network parity",
    "surface": "Computing the molecular surface",
    "writeSurface": "Writing the surface",
    "surfaceLigands": "Surfacing ligands",
    "writeLigandSurfaces": "Writing ligand surfaces",
    "selectTscClassificationAtoms": "Selecting atoms to classify",
    "classifySidechainLocation": "Classifying sidechain locations",
    "identifyMissingTscAtoms": "Checking for missing sidechain atoms",
    "writeSidechainClassificationResults": "Writing classification results",
    "classifyInterfaceSidechains": "Classifying interface sidechains",
    "makeSamplingGridUsingProteinSurface": "Building the sampling grid",
    "filterSamplingPointsUsingClashes": "Removing clashing grid points",
    "triangulateRemainingGridPoints": "Triangulating grid points",
    "identifyAndParseIndividualSamplingVoids": "Parsing sampling voids",
    "calculateSamplingVoidSurfaces": "Surfacing sampling voids",
}

# The calls each block performs, in order -- used both to run them and to size
# the progress bar before anything starts.
PREP_STEPS = ["setQuerySet", "openPDBs", "hetLigand4D", "hydrogens",
              "makeAtomCollections", "makeVertices4D"]
TRIANGULATE_STEPS = ["selectTscTriangulationAtoms", "triangulateTscAtoms",
                     "writeTriangulation", "pruneTriangulation",
                     "minimizePrunedTriangulation", "identifyTightBonds",
                     "calculateNetworkParity"]
SURFACE_STEPS = ["surface", "writeSurface", "surfaceLigands", "writeLigandSurfaces"]
CLASSIFY_STEPS = ["selectTscClassificationAtoms", "classifySidechainLocation",
                  "identifyMissingTscAtoms", "writeSidechainClassificationResults"]
SCREEN_STEPS = ["makeSamplingGridUsingProteinSurface", "filterSamplingPointsUsingClashes",
                "triangulateRemainingGridPoints", "identifyAndParseIndividualSamplingVoids",
                "calculateSamplingVoidSurfaces"]

# What each calculation needs, so the total can be counted without running it.
NEEDS = {
    "topologyCalculation": ["triangulate"],
    "surfaceCalculation": ["surface"],
    "sidechainClassification": ["triangulate", "surface", "classify"],
    "interfaceClassification": ["triangulate", "surface", "classify", "interface"],
    "virtualScreenSurfacesCalculation": ["triangulate", "surface", "screen"],
}
BLOCK_STEPS = {
    "triangulate": TRIANGULATE_STEPS, "surface": SURFACE_STEPS,
    "classify": CLASSIFY_STEPS, "interface": ["classifyInterfaceSidechains"],
    "screen": SCREEN_STEPS,
}


def count_steps(calc):
    """How many calls a run will make, counting each shared block once."""
    blocks = []
    for key, needed in NEEDS.items():
        if calc.get(key):
            for b in needed:
                if b not in blocks:
                    blocks.append(b)
    return len(PREP_STEPS) + sum(len(BLOCK_STEPS[b]) for b in blocks)


class Cancelled(Exception):
    """Raised to unwind when the user asks the run to stop."""


class _Out:
    """Forward whatever the calculation prints to the progress panel.

    pHinder prints its own narration, so the panel shows the real thing rather
    than a paraphrase. The original replaced ``sys.stdout`` globally and never
    put it back; this restores it in a finally block.
    """

    def __init__(self, report, original):
        self._report, self._original, self._buf = report, original, ""

    def write(self, text):
        if self._original is not None:
            try:
                self._original.write(text)
            except Exception:
                pass
        self._buf += text
        while "\n" in self._buf:
            line, self._buf = self._buf.split("\n", 1)
            self._report.write(line)

    def flush(self):
        if self._buf:
            self._report.write(self._buf)
            self._buf = ""


class _Steps:
    """The shared prerequisites, each performed at most once per run."""

    def __init__(self, inst, report, surface_kwargs):
        self._inst, self._report, self._surface_kwargs = inst, report, surface_kwargs
        self._done = set()

    def _check(self):
        if self._report.cancelled:
            raise Cancelled()

    def do(self, name, *args, **kwargs):
        """Announce a call, honour a stop request, then make it.

        Reporting per call is what makes the bar move inside a long stage --
        stage boundaries alone can be minutes apart.
        """
        self._check()
        self._report.substep(STEP_LABELS.get(name, name))
        return getattr(self._inst, name)(*args, **kwargs)

    def run_block(self, names):
        for name in names:
            self.do(name)

    def triangulate(self):
        if "triangulate" in self._done:
            return
        self.run_block(TRIANGULATE_STEPS)
        self._done.add("triangulate")

    def surface(self):
        if "surface" in self._done:
            return
        self.do("surface", **self._surface_kwargs)
        for name in SURFACE_STEPS[1:]:
            self.do(name)
        self._done.add("surface")

    def classify(self):
        if "classify" in self._done:
            return
        self.triangulate()
        self.surface()
        self.run_block(CLASSIFY_STEPS)
        self._done.add("classify")


def _residue_set(selections, report):
    chosen = tuple(sorted(k for k, v in selections.items() if v))
    if chosen in AA_SETS:
        report.write(f"Residue set: {AA_SETS[chosen]}")
        return AA_SETS[chosen]
    custom = "customSet:" + ",".join(a.upper() for a in chosen)
    report.write(f"Residue set: {custom}")
    return custom


def _configure(inst, results, report):
    """Copy the collected parameters onto the pHinder instance."""
    file_path = results["file_path"]
    pdb_file_name = file_path.split(sep)[-1]

    inst.gui = True
    try:
        import multiprocessing as mp
        inst.processes = max(1, mp.cpu_count() - 1)
    except Exception:
        inst.processes = 1

    inst.pdbFilePath = sep.join(file_path.split(sep)[:-1]) + sep
    inst.pdbFileName = pdb_file_name
    inst.outPath = (results["save_path"] or sep.join(file_path.split(sep)[:-1])) + sep + "pHinderResults/"
    inst.pdbFormat = "mmCIF" if ".cif" in pdb_file_name else "pdb"
    inst.zip = 1 if ".gz" in pdb_file_name else 0

    chains = results["chains"]
    inst.chains = sorted(c for c, on in chains.items() if on and c != GROUP_CHAINS)
    inst.group_chains = chains.get(GROUP_CHAINS, 0)

    inst.residueSet = _residue_set(results["amino_acid_selections"], report)

    net = results["network_options"]
    inst.maxNetworkEdgeLength = net['MAX_NETWORK_EDGE_LENGTH']
    inst.minNetworkSize = net['MIN_NETWORK_SIZE']
    inst.saveNetworkTriangulation = net['SAVE_NETWORK_TRIANGULATION']

    surf = results["surface_options"]
    inst.highResolutionSurface = surf['HIGH_RESOLUTION_SURFACE']
    inst.saveSurface = surf['SAVE_SURFACE']
    inst.allowSmallSurfaces = surf['ALLOW_SMALL_SURFACES']
    inst.saveLigandSurfaces = surf['SAVE_LIGAND_SURFACES']
    inst.writeSurfaceCreationAnimation = surf['WRITE_SURFACE_CREATION_ANIMATION']

    sc = results["sidechain_classification_options"]
    inst.coreCutoff = sc['CORE_CUTOFF']
    inst.marginCutoff = sc['MARGIN_CUTOFF']

    inst.interface_distance_filter = results["interface_options"]["INTERFACE_DISTANCE_FILTER"]

    vs = results["virtual_screening_options"]
    inst.virtualClashCutoff = vs['VIRTUAL_CLASH_CUTOFF']
    inst.inIterations = vs['IN_ITERATIONS']
    inst.inIterationsStepSize = vs['IN_ITERATIONS_STEP_SIZE']
    inst.outIterations = vs['OUT_ITERATIONS']
    inst.outIterationsStepSize = vs['OUT_ITERATIONS_STEP_SIZE']

    adv = results["advanced_options"]
    inst.includeHydrogens = adv['INCLUDE_HYDROGENS']
    inst.includeWater = adv['INCLUDE_WATER']
    inst.includeIons = adv['INCLUDE_IONS']

    return {"circumSphereRadiusLimit": surf['circumSphereRadiusLimit'],
            "minArea": surf['minArea']}


def run(results, report):
    """Perform the selected calculations, reporting each one as it goes."""
    from pHinder.pHinder_7_0 import pHinder
    from pHinder.gui.phinder_main_gui import write_phinder_log

    original_stdout, original_stderr = sys.stdout, sys.stderr
    sys.stdout = sys.stderr = _Out(report, original_stdout)
    # pHinder is CPU-bound pure Python. On the default 5 ms switch interval it
    # holds the GIL long enough to starve the Tk main loop, so the window goes
    # sluggish and progress cannot paint until the run ends. Yielding more often
    # costs a little throughput and buys back a responsive UI.
    original_switch = sys.getswitchinterval()
    sys.setswitchinterval(0.001)
    started = time.time()
    calc = results["calculation_options"]

    try:
        name = results["file_path"].split(sep)[-1]
        report.write(f"Running pHinder for {name}")
        report.write("-" * 60)

        inst = pHinder()
        surface_kwargs = _configure(inst, results, report)

        limit = int(results["advanced_options"]['PYTHON_RECURSION_LIMIT'])
        report.write(f"Recursion limit {sys.getrecursionlimit()} -> {limit}")
        sys.setrecursionlimit(limit)

        report.set_total_steps(count_steps(calc))
        report.status("Reading structure…", name)
        steps = _Steps(inst, report, surface_kwargs)
        steps.do("setQuerySet")
        steps.do("openPDBs", inst.pdbFilePath, inst.pdbFileName, zip_status=inst.zip)
        for _name in PREP_STEPS[2:]:
            steps.do(_name)

        if results["advanced_options"].get('SAVE_LOG_FILE', 1):
            # The parameter log is the first thing written, so it is the first
            # thing to fail if the results directory is not there yet.
            os.makedirs(inst.outPath, exist_ok=True)
            write_phinder_log(inst, inst.outPath + "pHinder_parameters.log")

        # (key, label, what it needs, what it then does) -- in run order.
        plan = [
            ("topologyCalculation", "Residue network topology", steps.triangulate, None),
            ("surfaceCalculation", "Molecular surface", steps.surface, None),
            ("sidechainClassification", "Sidechain classification", steps.classify, None),
            ("interfaceClassification", "Interface classification", steps.classify,
             lambda: steps.do("classifyInterfaceSidechains")),
            ("virtualScreenSurfacesCalculation", "Virtual screening surfaces",
             lambda: (steps.triangulate(), steps.surface()),
             lambda: _virtual_screen(steps, results)),
        ]

        current = None
        try:
            for key, label, prereq, extra in plan:
                if not calc.get(key):
                    continue
                current = key
                report.begin(key)
                report.status(f"{label}…", name)
                elapsed = time.time()
                prereq()
                if extra:
                    extra()
                report.done(key, f"{time.time() - elapsed:.1f}s")
                current = None
        except Cancelled:
            # Mark the calculation that was in flight, so the checklist shows
            # where the run actually stopped rather than going quiet.
            if current:
                report.failed(current, "stopped")
            report.write("Stopped before the next step.", error=True)
            report.status("Stopped", "Partial results may have been written.")
            return
        except Exception as exc:
            if current:
                report.failed(current, "failed")
            report.write(f"{type(exc).__name__}: {exc}", error=True)
            raise

        report.write(f"\nRuntime: {time.time() - started:.1f}s")
        report.write("pHinder completed successfully.")
        report.status("Finished", f"Results in {inst.outPath}")

    finally:
        try:
            sys.stdout.flush()
        except Exception:
            pass
        sys.setswitchinterval(original_switch)
        sys.stdout, sys.stderr = original_stdout, original_stderr


def _virtual_screen(steps, results):
    vs = results["virtual_screening_options"]
    steps.do("makeSamplingGridUsingProteinSurface")
    steps.do("filterSamplingPointsUsingClashes")
    steps.do("triangulateRemainingGridPoints")
    steps.do("identifyAndParseIndividualSamplingVoids",
             maxVoidNetworkEdgeLength=vs['MAX_VOID_NETWORK_EDGE_LENGTH'],
             minVoidNetworkEdgeLength=0.0,
             minVoidNetworkSize=vs['MIN_VOID_NETWORK_SIZE'],
             psa=1)
    steps.do("calculateSamplingVoidSurfaces", extend_sampling=True)
