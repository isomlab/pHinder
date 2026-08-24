"""Tests for the GUI runner.

A stand-in pHinder class records which methods get called, so the ordering,
the memoised prerequisites, cancellation and parameter marshalling can all be
checked without the science dependencies or a structure file.
"""

import sys
import types

import pytest


@pytest.fixture
def fake_phinder(monkeypatch):
    """Install a recording stand-in for pHinder.pHinder_7_0.pHinder."""
    calls = []

    class Fake:
        def __init__(self):
            object.__setattr__(self, "attrs", {})

        def __setattr__(self, key, value):
            self.attrs[key] = value

        def __getattr__(self, name):
            if name.startswith("_") or name == "attrs":
                raise AttributeError(name)
            if name in self.attrs:            # a value that was set, not a method
                return self.attrs[name]

            def record(*args, **kwargs):
                calls.append(name)

            return record

    holder = {"instance": None}

    def factory():
        holder["instance"] = Fake()
        return holder["instance"]

    module = types.ModuleType("pHinder.pHinder_7_0")
    module.pHinder = factory
    monkeypatch.setitem(sys.modules, "pHinder.pHinder_7_0", module)
    return calls, holder


class Report:
    """Minimal stand-in for the progress panel's reporting surface."""

    def __init__(self, stop_after=None):
        self.stages, self.lines, self.substeps = [], [], []
        self.total = None
        self._stop_after, self._begun = stop_after, 0

    @property
    def cancelled(self):
        return self._stop_after is not None and self._begun > self._stop_after

    def begin(self, key):
        self._begun += 1
        self.stages.append(("begin", key))

    def done(self, key, note=""):
        self.stages.append(("done", key))

    def failed(self, key, note=""):
        self.stages.append(("failed", key))

    def skipped(self, key, note=""):
        self.stages.append(("skipped", key))

    def status(self, text, detail=""):
        pass

    def set_total_steps(self, n):
        self.total = n

    def substep(self, text):
        self.substeps.append(text)

    def write(self, text, error=False):
        self.lines.append(text)

    def begun(self):
        return [k for kind, k in self.stages if kind == "begin"]


CALC_KEYS = [
    "topologyCalculation",
    "surfaceCalculation",
    "sidechainClassification",
    "interfaceClassification",
    "virtualScreenSurfacesCalculation",
]


def make_results(**calc):
    options = dict.fromkeys(CALC_KEYS, 0)
    options.update(calc)
    return {
        "file_path": "/data/1UBQ.pdb",
        "save_path": "/data/out",
        "chains": {"A": 1, "B": 0, "Group Chains": 1},
        "amino_acid_selections": {"Asp": 1, "Glu": 1, "Lys": 1, "Arg": 1, "His": 1, "Ala": 0},
        "calculation_options": options,
        "sidechain_classification_options": {
            "CORE_CUTOFF": -3.0, "MARGIN_CUTOFF": -2.0},
        "network_options": {
            "MAX_NETWORK_EDGE_LENGTH": 10.0, "MIN_NETWORK_SIZE": 1,
            "SAVE_NETWORK_TRIANGULATION": 1},
        "surface_options": {
            "HIGH_RESOLUTION_SURFACE": 1, "SAVE_SURFACE": 1, "ALLOW_SMALL_SURFACES": 0,
            "SAVE_LIGAND_SURFACES": 0, "WRITE_SURFACE_CREATION_ANIMATION": 0,
            "circumSphereRadiusLimit": 10.0, "minArea": 1.0},
        "interface_options": {"INTERFACE_DISTANCE_FILTER": 8.0},
        "virtual_screening_options": {
            "MAX_VOID_NETWORK_EDGE_LENGTH": 2.0, "MIN_VOID_NETWORK_SIZE": 10,
            "VIRTUAL_CLASH_CUTOFF": 2.5, "IN_ITERATIONS": 1, "IN_ITERATIONS_STEP_SIZE": 2.0,
            "OUT_ITERATIONS": 1, "OUT_ITERATIONS_STEP_SIZE": 2.0},
        "advanced_options": {
            "INCLUDE_HYDROGENS": 0, "INCLUDE_WATER": 0,
            "INCLUDE_IONS": 0, "SAVE_LOG_FILE": 0, "PYTHON_RECURSION_LIMIT": 10000},
    }


def run(results, report):
    from pHinder.gui import runner
    runner.run(results, report)


def test_shared_prerequisites_run_once(fake_phinder):
    """Interface + screening without topology must not triangulate twice.

    The pane-based GUI guarded the prerequisite inside each dependent
    calculation, so this combination did the whole triangulation and surface
    twice over.
    """
    calls, _ = fake_phinder
    report = Report()
    run(make_results(interfaceClassification=1, virtualScreenSurfacesCalculation=1), report)
    assert calls.count("triangulateTscAtoms") == 1
    assert calls.count("writeSurface") == 1
    assert "classifyInterfaceSidechains" in calls
    assert "calculateSamplingVoidSurfaces" in calls


def test_stages_reported_in_run_order(fake_phinder):
    _, _ = fake_phinder
    report = Report()
    run(make_results(topologyCalculation=1, surfaceCalculation=1, sidechainClassification=1), report)
    assert report.begun() == ["topologyCalculation", "surfaceCalculation", "sidechainClassification"]
    assert sum(1 for kind, _ in report.stages if kind == "done") == 3


def test_unticked_calculations_do_not_run(fake_phinder):
    calls, _ = fake_phinder
    report = Report()
    run(make_results(surfaceCalculation=1), report)
    assert report.begun() == ["surfaceCalculation"]
    assert "triangulateTscAtoms" not in calls


def test_group_chains_is_an_option_not_a_chain(fake_phinder):
    """'Group Chains' shares the chain dict but must never become a chain id."""
    _, holder = fake_phinder
    run(make_results(surfaceCalculation=1), Report())
    attrs = holder["instance"].attrs
    assert attrs["chains"] == ["A"]
    assert attrs["group_chains"] == 1


def test_residue_set_maps_to_a_named_set(fake_phinder):
    _, holder = fake_phinder
    run(make_results(surfaceCalculation=1), Report())
    assert holder["instance"].attrs["residueSet"] == "ionizableSetNoCys"


def test_custom_residue_selection_becomes_a_custom_set(fake_phinder):
    _, holder = fake_phinder
    results = make_results(surfaceCalculation=1)
    results["amino_acid_selections"] = {"Asp": 1, "Trp": 1}
    run(results, Report())
    assert holder["instance"].attrs["residueSet"] == "customSet:ASP,TRP"


def test_cancel_stops_and_marks_the_stage_in_flight(fake_phinder):
    _, _ = fake_phinder
    report = Report(stop_after=1)
    run(make_results(topologyCalculation=1, surfaceCalculation=1, sidechainClassification=1), report)
    assert any(kind == "failed" for kind, _ in report.stages)
    assert sum(1 for kind, _ in report.stages if kind == "done") < 3


def test_stdout_is_restored(fake_phinder):
    """The pane GUI replaced sys.stdout globally and never put it back."""
    _, _ = fake_phinder
    before = sys.stdout
    run(make_results(surfaceCalculation=1), Report())
    assert sys.stdout is before


def test_progress_reports_each_call_not_just_each_stage(fake_phinder):
    """The bar must move inside a stage; stage boundaries can be minutes apart."""
    _, _ = fake_phinder
    report = Report()
    run(make_results(topologyCalculation=1), report)
    assert "Reading the structure" in report.substeps
    assert "Triangulating sidechain atoms" in report.substeps
    # One sub-step per pHinder call: 6 prep + 7 triangulation.
    assert len(report.substeps) == 13


def test_total_steps_counts_shared_blocks_once(fake_phinder):
    """The bar is sized before the run, so it must predict the memoisation."""
    from pHinder.gui import runner

    _, _ = fake_phinder
    report = Report()
    calc = make_results(interfaceClassification=1,
                        virtualScreenSurfacesCalculation=1)["calculation_options"]
    predicted = runner.count_steps(calc)
    run(make_results(interfaceClassification=1, virtualScreenSurfacesCalculation=1), report)
    assert report.total == predicted
    assert len(report.substeps) == predicted


def test_switch_interval_is_restored(fake_phinder):
    """The run lowers it to keep the UI responsive; it must not leak."""
    _, _ = fake_phinder
    before = sys.getswitchinterval()
    run(make_results(surfaceCalculation=1), Report())
    assert sys.getswitchinterval() == before


def test_every_parameter_and_calculation_has_help():
    """Hover help must cover the whole surface, or it is not help."""
    from pHinder.gui import help_text
    from pHinder.gui import phinder_main_gui as legacy
    from pHinder.gui.app import CALCULATIONS

    groups = ["sidechain_classification_options", "network_options", "surface_options",
              "interface_options", "virtual_screening_options", "advanced_options"]
    undocumented = []
    for group in groups:
        for key in getattr(legacy, f"default_{group}"):
            if not help_text.for_option(key)[1]:
                undocumented.append(f"{group}.{key}")
    for key, _, _ in CALCULATIONS:
        if key not in help_text.CALCULATIONS:
            undocumented.append(f"calculation.{key}")
    assert not undocumented, f"no hover help for: {undocumented}"


def test_counts_are_entries_not_checkboxes():
    """A default of 1 does not make a parameter boolean.

    IN_ITERATIONS, OUT_ITERATIONS and MIN_NETWORK_SIZE all default to 1 but are
    counts; treating them as flags capped them at 1 in the GUI.
    """
    from pHinder.gui import help_text

    for key in ("IN_ITERATIONS", "OUT_ITERATIONS", "MIN_NETWORK_SIZE",
                "PYTHON_RECURSION_LIMIT", "MIN_VOID_NETWORK_SIZE"):
        assert not help_text.is_boolean(key), f"{key} should be a numeric field"
    for key in ("SAVE_SURFACE", "INCLUDE_WATER", "HIGH_RESOLUTION_SURFACE"):
        assert help_text.is_boolean(key), f"{key} should be a checkbox"


INERT = ("MARGIN_CUTOFF_CORE_NETWORK", "REDUCED_NETWORK_REPRESENTATION",
         "ALLOW_CYS_CORE_SEEDING")


def test_inert_parameters_are_gone_from_the_gui():
    """They were offered as controls but nothing read them.

    The algorithms that consumed them were rewritten for 7.0; see
    docs/inert_parameters.md for what each one did.
    """
    from pHinder.gui import phinder_main_gui as legacy

    groups = ["sidechain_classification_options", "network_options", "surface_options",
              "interface_options", "virtual_screening_options", "advanced_options"]
    offered = {k for g in groups for k in getattr(legacy, f"default_{g}")}
    assert not (offered & set(INERT))


def test_runner_does_not_assign_inert_parameters(fake_phinder):
    """Assigning a value nothing reads only invites the confusion back."""
    _, holder = fake_phinder
    run(make_results(surfaceCalculation=1), Report())
    attrs = holder["instance"].attrs
    for name in ("marginCutoffCoreNetwork", "reducedNetworkRepresentation",
                 "allowCysCoreSeeding"):
        assert name not in attrs


def test_defaults_agree_across_entry_points():
    """The GUI, the CLI and the class must start from the same values.

    They had drifted: the class began with marginCutoff 1.05 against -2.0
    everywhere else, and virtualClashCutoff 3.0 against 2.5. Anyone
    constructing pHinder() directly got different classification bands than
    the GUI and CLI produce.
    """
    import re
    from pathlib import Path

    from pHinder.gui import phinder_main_gui as legacy

    src = Path(__file__).resolve().parents[1] / "src" / "pHinder"

    def class_default(attr):
        text = (src / "pHinder_7_0.py").read_text()
        m = re.search(rf"self\.{attr}\s*=\s*(-?[\d.]+)", text)
        assert m, f"no class default found for {attr}"
        return float(m.group(1))

    def cli_default(flag):
        text = (src / "command_line.py").read_text()
        m = re.search(rf'"--{flag}".*?default=(-?[\d.]+)', text)
        assert m, f"no CLI default found for --{flag}"
        return float(m.group(1))

    checks = [
        ("marginCutoff", "margin-cutoff",
         legacy.default_sidechain_classification_options["MARGIN_CUTOFF"]),
        ("coreCutoff", "core-cutoff",
         legacy.default_sidechain_classification_options["CORE_CUTOFF"]),
        ("virtualClashCutoff", "virtual-clash-cutoff",
         legacy.default_virtual_screening_options["VIRTUAL_CLASH_CUTOFF"]),
        ("maxNetworkEdgeLength", "max-network-edge-length",
         legacy.default_network_options["MAX_NETWORK_EDGE_LENGTH"]),
        ("interface_distance_filter", "interface-distance-filter",
         legacy.default_interface_options["INTERFACE_DISTANCE_FILTER"]),
    ]
    mismatches = []
    for attr, flag, gui_value in checks:
        cls, cli = class_default(attr), cli_default(flag)
        if not (cls == cli == float(gui_value)):
            mismatches.append(f"{attr}: class={cls} cli={cli} gui={gui_value}")
    assert not mismatches, "defaults disagree -> " + "; ".join(mismatches)


def test_margin_cutoff_sits_between_core_and_the_surface():
    """The margin band only exists if coreCutoff < marginCutoff < 0.

    inLocalSurface tests `aveDistance <= coreCutoff` for core, then
    `aveDistance <= marginCutoff` for margin. A positive marginCutoff would
    sweep exposed residues into the margin class, which is what 1.05 did.
    """
    from pHinder.gui import phinder_main_gui as legacy

    sc = legacy.default_sidechain_classification_options
    assert sc["CORE_CUTOFF"] < sc["MARGIN_CUTOFF"] < 0
