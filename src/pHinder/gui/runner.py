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

    def triangulate(self):
        if "triangulate" in self._done:
            return
        self._check()
        inst = self._inst
        inst.selectTscTriangulationAtoms()
        inst.triangulateTscAtoms()
        inst.writeTriangulation()
        self._check()
        inst.pruneTriangulation()
        inst.minimizePrunedTriangulation()
        self._check()
        inst.identifyTightBonds()
        inst.calculateNetworkParity()
        self._done.add("triangulate")

    def surface(self):
        if "surface" in self._done:
            return
        self._check()
        inst = self._inst
        inst.surface(**self._surface_kwargs)
        inst.writeSurface()
        self._check()
        inst.surfaceLigands()
        inst.writeLigandSurfaces()
        self._done.add("surface")

    def classify(self):
        if "classify" in self._done:
            return
        self.triangulate()
        self.surface()
        self._check()
        inst = self._inst
        inst.selectTscClassificationAtoms()
        inst.classifySidechainLocation()
        self._check()
        inst.identifyMissingTscAtoms()
        inst.writeSidechainClassificationResults()
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
    inst.reducedNetworkRepresentation = net['REDUCED_NETWORK_REPRESENTATION']
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
    inst.marginCutoffCoreNetwork = sc['MARGIN_CUTOFF_CORE_NETWORK']

    inst.interface_distance_filter = results["interface_options"]["INTERFACE_DISTANCE_FILTER"]

    vs = results["virtual_screening_options"]
    inst.virtualClashCutoff = vs['VIRTUAL_CLASH_CUTOFF']
    inst.inIterations = vs['IN_ITERATIONS']
    inst.inIterationsStepSize = vs['IN_ITERATIONS_STEP_SIZE']
    inst.outIterations = vs['OUT_ITERATIONS']
    inst.outIterationsStepSize = vs['OUT_ITERATIONS_STEP_SIZE']

    adv = results["advanced_options"]
    inst.allowCysCoreSeeding = adv['ALLOW_CYS_CORE_SEEDING']
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

        report.status("Reading structure…", name)
        inst.setQuerySet()
        inst.openPDBs(inst.pdbFilePath, inst.pdbFileName, zip_status=inst.zip)
        inst.hetLigand4D()
        inst.hydrogens()
        inst.makeAtomCollections()
        inst.makeVertices4D()

        if results["advanced_options"].get('SAVE_LOG_FILE', 1):
            # The parameter log is the first thing written, so it is the first
            # thing to fail if the results directory is not there yet.
            os.makedirs(inst.outPath, exist_ok=True)
            write_phinder_log(inst, inst.outPath + "pHinder_parameters.log")

        steps = _Steps(inst, report, surface_kwargs)

        # (key, label, what it needs, what it then does) -- in run order.
        plan = [
            ("topologyCalculation", "Residue network topology", steps.triangulate, None),
            ("surfaceCalculation", "Molecular surface", steps.surface, None),
            ("sidechainClassification", "Sidechain classification", steps.classify, None),
            ("interfaceClassification", "Interface classification", steps.classify,
             lambda: inst.classifyInterfaceSidechains()),
            ("virtualScreenSurfacesCalculation", "Virtual screening surfaces",
             lambda: (steps.triangulate(), steps.surface()),
             lambda: _virtual_screen(inst, results, report)),
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
        sys.stdout, sys.stderr = original_stdout, original_stderr


def _virtual_screen(inst, results, report):
    vs = results["virtual_screening_options"]
    inst.makeSamplingGridUsingProteinSurface()
    inst.filterSamplingPointsUsingClashes()
    inst.triangulateRemainingGridPoints()
    inst.identifyAndParseIndividualSamplingVoids(
        maxVoidNetworkEdgeLength=vs['MAX_VOID_NETWORK_EDGE_LENGTH'],
        minVoidNetworkEdgeLength=0.0,
        minVoidNetworkSize=vs['MIN_VOID_NETWORK_SIZE'],
        psa=1,
    )
    inst.calculateSamplingVoidSurfaces(extend_sampling=True)
