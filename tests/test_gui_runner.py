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
        self.stages, self.lines = [], []
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
            "CORE_CUTOFF": -3.0, "MARGIN_CUTOFF": -2.0, "MARGIN_CUTOFF_CORE_NETWORK": -2.0},
        "network_options": {
            "MAX_NETWORK_EDGE_LENGTH": 10.0, "MIN_NETWORK_SIZE": 1,
            "REDUCED_NETWORK_REPRESENTATION": 1, "SAVE_NETWORK_TRIANGULATION": 1},
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
            "ALLOW_CYS_CORE_SEEDING": 0, "INCLUDE_HYDROGENS": 0, "INCLUDE_WATER": 0,
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
