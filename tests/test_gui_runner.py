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


def test_tooltip_does_not_survive_a_tab_change():
    """A hidden tab unmaps its children without sending <Leave>.

    The tip is a topmost window, so one left behind sits over the newly shown
    tab until the pointer happens to cross another widget -- which looked like
    the tab failing to draw its own contents.

    Drives the real window, so it exercises the app's own binding rather than
    one the test installed.
    """
    import tkinter as tk

    import pytest

    from pHinder.gui import phinder_main_gui as legacy
    from pHinder.gui.app import PHinderApp

    groups = ["calculation_options", "sidechain_classification_options",
              "network_options", "surface_options", "interface_options",
              "virtual_screening_options", "advanced_options"]
    try:
        app = PHinderApp({g: getattr(legacy, f"default_{g}") for g in groups})
    except tk.TclError:                       # no display
        pytest.skip("no display")
    try:
        for var in app.vars["calculation_options"].values():
            var.set(1)
        app.update_idletasks()

        names = [app._nb.tab(i, "text") for i in range(len(app._nb.tabs()))]
        first, second = names.index("Screening"), names.index("Networks")
        app._nb.select(first)
        app.update_idletasks()

        label = None
        def find(w):
            nonlocal label
            if label is None and "<Enter>" in w.bind() and w.winfo_class() == "TLabel":
                label = w
            for child in w.winfo_children():
                find(child)
        find(app.nametowidget(app._nb.tabs()[first]))
        assert label is not None, "no widget with hover help on the Screening tab"

        app._tip = getattr(app, "_tip", None)
        from pHinder.gui.tooltip import _Tip
        if app._tip is None:
            app._tip = _Tip(app)
        app._tip.show(label, "Title", "some help")
        app.update_idletasks()
        assert app._tip._win is not None, "tip should be showing"

        app._nb.select(second)
        app.update_idletasks()
        app.update()
        assert app._tip._win is None, "tip must not outlive the tab that raised it"
    finally:
        app.destroy()


def test_scrollhost_refresh_repairs_a_stale_canvas():
    """A tab mapped with stale canvas state draws as an empty pane.

    A canvas only learns its width from a <Configure>, so one hidden across a
    resize can come back with item width 0, an oversized scrollregion and a
    scrolled yview -- which draws as nothing at all. ScrollHost.refresh() is
    what a newly shown tab calls to put all three right.

    Exercised directly: driving it through a real tab change depends on Tk
    servicing after_idle, which is not reliable under pytest.
    """
    import tkinter as tk

    import pytest

    from pHinder.gui import theme

    try:
        root = tk.Tk()
    except tk.TclError:                       # no display
        pytest.skip("no display")
    try:
        root.geometry("600x400")
        root.fonts = theme.Fonts()
        theme.apply_style(root.fonts)
        host = theme.ScrollHost(root)
        holder = tk.Frame(root, width=600, height=400)
        holder.pack(fill="both", expand=True)
        body = host.body(holder)
        tk.Label(body, text="content").pack()
        root.update_idletasks()

        canvas = host._canvases[-1]
        item = host._items[canvas]

        canvas.itemconfigure(item, width=0)
        canvas.configure(scrollregion=(0, 0, 600, 4000))
        canvas.yview_moveto(0.55)
        root.update_idletasks()
        assert int(canvas.itemcget(item, "width")) == 0
        assert canvas.yview()[0] > 0

        host.refresh(holder)
        root.update_idletasks()

        assert int(canvas.itemcget(item, "width")) > 1, "body never got the canvas width"
        assert canvas.yview()[0] == 0.0, "short content left scrolled out of view"
        assert canvas.cget("scrollregion") != "0 0 600 4000", "stale scrollregion kept"
    finally:
        root.destroy()


def test_scrollhost_refresh_does_not_touch_sibling_tabs():
    """Tk paths are dot-separated: ".!frame" must not match ".!frame3"."""
    import tkinter as tk

    import pytest

    from pHinder.gui import theme

    try:
        root = tk.Tk()
    except tk.TclError:
        pytest.skip("no display")
    try:
        root.fonts = theme.Fonts()
        theme.apply_style(root.fonts)
        host = theme.ScrollHost(root)
        first = tk.Frame(root, name="frame", width=300, height=200)
        second = tk.Frame(root, name="frame3", width=300, height=200)
        for f in (first, second):
            f.pack()
            host.body(f)
        root.update_idletasks()

        sibling = host._canvases[-1]
        sibling_item = host._items[sibling]
        sibling.itemconfigure(sibling_item, width=0)

        host.refresh(first)                    # refresh only the first tab
        root.update_idletasks()

        assert int(sibling.itemcget(sibling_item, "width")) == 0, \
            "refreshing one tab reached into its sibling"
    finally:
        root.destroy()


def test_closing_during_a_run_exits_cleanly(fake_phinder):
    """The red X must work while a calculation is in flight.

    Tearing the widget tree down with a queued after() callback and a tooltip
    Toplevel still live raised TclError from inside destroy(), which left the
    window half-dismantled and the process alive.
    """
    import tkinter as tk

    import pytest

    from pHinder.gui import phinder_main_gui as legacy
    from pHinder.gui.app import PHinderApp

    _, _ = fake_phinder
    groups = ["calculation_options", "sidechain_classification_options",
              "network_options", "surface_options", "interface_options",
              "virtual_screening_options", "advanced_options"]
    try:
        app = PHinderApp({g: getattr(legacy, f"default_{g}") for g in groups})
    except tk.TclError:
        pytest.skip("no display")

    closed = False
    try:
        app.vars["calculation_options"]["surfaceCalculation"].set(1)
        app.file_widget.file_path.set("/data/1UBQ.pdb")
        app.file_widget.show_file_specific_options("Chains:", ["A"])
        for name, var in app.file_widget.option_vars.items():
            var.set(1 if name == "A" else 0)
        app.update_idletasks()

        app.progress._run_clicked()
        app.update_idletasks()

        app.on_close()                     # exactly what the window button calls
        closed = True

        # A destroyed root cannot be queried at all -- that is the proof it went.
        with pytest.raises(tk.TclError):
            app.winfo_exists()
    finally:
        if not closed:
            app.destroy()


def test_progress_panel_drain_loop_can_be_stopped():
    """A pending drain after() firing mid-destroy is what raised from destroy()."""
    import tkinter as tk

    import pytest

    from pHinder.gui import theme
    from pHinder.gui.progress import ProgressPanel

    try:
        root = tk.Tk()
    except tk.TclError:
        pytest.skip("no display")
    try:
        root.fonts = theme.Fonts()
        theme.apply_style(root.fonts)
        panel = ProgressPanel(root, root.fonts, lambda: None, lambda: None)
        panel.pack()
        root.update_idletasks()

        assert panel._drain_after is not None, "drain loop should be scheduled"
        panel.shutdown()
        assert panel._drain_after is None, "shutdown must cancel the pending drain"
        assert panel._closing is True
    finally:
        root.destroy()


def test_file_and_save_path_have_distinct_help():
    """Both fields once shared one tooltip attached to the whole card."""
    from pHinder.gui import help_text

    file_title, file_body = help_text.INPUT["file_path"]
    save_title, save_body = help_text.INPUT["save_path"]
    assert file_body != save_body
    assert file_title != save_title


def test_residue_help_frames_the_selection_as_network_nodes():
    from pHinder.gui import help_text

    body = help_text.INPUT["residues"][1].lower()
    assert "network" in body


def test_surface_help_says_ligands_are_optional():
    """surfaceLigands() and writeLigandSurfaces() both no-op unless
    saveLigandSurfaces is set, so the help must not promise ligand surfaces."""
    from pHinder.gui import help_text

    body = help_text.CALCULATIONS["surfaceCalculation"][1].lower()
    assert "only if" in body
    assert "save ligand surfaces" in body


def test_dependent_calculations_pull_in_their_prerequisites():
    """Ticking classification or interfaces must not require ticking the rest."""
    from pHinder.gui import runner

    assert runner.NEEDS["sidechainClassification"] == ["triangulate", "surface", "classify"]
    assert runner.NEEDS["interfaceClassification"] == [
        "triangulate", "surface", "classify", "interface"]
    for key in ("sidechainClassification", "interfaceClassification"):
        body = runner and __import__(
            "pHinder.gui.help_text", fromlist=["x"]).CALCULATIONS[key][1].lower()
        assert "runs the networks" in body, f"{key} help should say what it runs for you"


def test_interface_chain_check_covers_every_selection():
    """Interface classification compares chains, so one chain gives it nothing.

    Either offer every chain the structure has, or say plainly that it has only
    one -- never run it silently against a single chain.
    """
    from pHinder.gui.app import interface_chain_check
    from pHinder.gui.runner import GROUP_CHAINS

    two_available = {"A": 1, "B": 0, GROUP_CHAINS: 0}
    both_selected = {"A": 1, "B": 1, GROUP_CHAINS: 0}
    one_only = {"A": 1, GROUP_CHAINS: 1}

    # Not asking for interfaces: never interfere.
    assert interface_chain_check(two_available, False) is None
    assert interface_chain_check(one_only, False) is None

    # Two or more chains chosen: fine.
    assert interface_chain_check(both_selected, True) is None

    # One chosen but more exist: offer them.
    assert interface_chain_check(two_available, True) == "offer_all"

    # The structure genuinely has one chain: nothing to offer.
    assert interface_chain_check(one_only, True) == "too_few"

    # "Group Chains" is an option, not a chain, and must not count as the second.
    assert interface_chain_check({"A": 1, GROUP_CHAINS: 1}, True) == "too_few"

    # Nothing selected at all still counts as too few to compare.
    assert interface_chain_check({"A": 0, "B": 0, GROUP_CHAINS: 0}, True) == "offer_all"


def test_runner_copes_when_group_chains_is_absent(fake_phinder):
    """A single-chain structure no longer offers the Group Chains box at all,
    so the key is missing from the chain dict rather than present and zero."""
    _, holder = fake_phinder
    results = make_results(surfaceCalculation=1)
    results["chains"] = {"A": 1}                 # no GROUP_CHAINS key
    run(results, Report())
    attrs = holder["instance"].attrs
    assert attrs["chains"] == ["A"]
    assert attrs["group_chains"] == 0


def test_interface_check_without_a_group_chains_key():
    from pHinder.gui.app import interface_chain_check

    assert interface_chain_check({"A": 1}, True) == "too_few"
    assert interface_chain_check({"A": 1, "B": 1}, True) is None


def test_group_chains_hidden_for_a_single_chain_structure():
    """Grouping is meaningless with one chain, and the box invites a click that
    cannot mean anything."""
    import tkinter as tk

    import pytest

    from pHinder.gui import phinder_main_gui as legacy
    from pHinder.gui.app import PHinderApp
    from pHinder.gui.runner import GROUP_CHAINS

    groups = ["calculation_options", "sidechain_classification_options",
              "network_options", "surface_options", "interface_options",
              "virtual_screening_options", "advanced_options"]
    try:
        app = PHinderApp({g: getattr(legacy, f"default_{g}") for g in groups})
    except tk.TclError:
        pytest.skip("no display")
    try:
        app.file_widget.show_file_specific_options("Chains:", ["A", "B"])
        app.update_idletasks()
        assert GROUP_CHAINS in app.file_widget.option_widgets

        app.file_widget.show_file_specific_options("Chains:", ["A"])
        app.update_idletasks()
        assert GROUP_CHAINS not in app.file_widget.option_widgets
        assert GROUP_CHAINS not in app.file_widget.get_values()
    finally:
        app.destroy()


def test_chain_boxes_get_hover_help_when_they_are_built():
    """The chain row does not exist until a file is read, so its help has to be
    attached at that point rather than when the form is laid out."""
    import tkinter as tk

    import pytest

    from pHinder.gui import phinder_main_gui as legacy
    from pHinder.gui.app import PHinderApp

    groups = ["calculation_options", "sidechain_classification_options",
              "network_options", "surface_options", "interface_options",
              "virtual_screening_options", "advanced_options"]
    try:
        app = PHinderApp({g: getattr(legacy, f"default_{g}") for g in groups})
    except tk.TclError:
        pytest.skip("no display")
    try:
        app.file_widget.show_file_specific_options("Chains:", ["A", "B"])
        app.update_idletasks()
        for name, box in app.file_widget.option_widgets.items():
            assert "<Enter>" in box.bind(), f"no hover help on the {name} box"
        assert "<Enter>" in app.file_widget.options_label_widget.bind()
    finally:
        app.destroy()


def test_the_defaults_note_appears_on_every_parameter_tab():
    """Every tab holding lab defaults carries the note -- the conditional five
    plus Advanced -- while Input and Calculations do not."""
    import tkinter as tk

    import pytest

    from pHinder.gui import phinder_main_gui as legacy
    from pHinder.gui.app import PHinderApp, DEFAULTS_NOTE, TAB_REQUIRES

    groups = ["calculation_options", "sidechain_classification_options",
              "network_options", "surface_options", "interface_options",
              "virtual_screening_options", "advanced_options"]
    try:
        app = PHinderApp({g: getattr(legacy, f"default_{g}") for g in groups})
    except tk.TclError:
        pytest.skip("no display")
    try:
        for var in app.vars["calculation_options"].values():
            var.set(1)
        app.update_idletasks()

        fragment = DEFAULTS_NOTE[:30]

        def has_note(frame):
            found = []

            def walk(widget):
                if isinstance(widget, tk.Label) and fragment in str(widget.cget("text")):
                    found.append(widget)
                for child in widget.winfo_children():
                    walk(child)

            walk(frame)
            return bool(found)

        carrying = {app._nb.tab(i, "text")
                    for i in range(len(app._nb.tabs()))
                    if has_note(app.nametowidget(app._nb.tabs()[i]))}
        assert carrying == set(TAB_REQUIRES) | {"Advanced"}
        # Input and Calculations are not parameter tabs.
        assert "Input" not in carrying
        assert "Calculations" not in carrying
    finally:
        app.destroy()


def test_min_area_help_does_not_describe_a_floor():
    """facetAreas() collects facets with area > minArea and dividePatch()
    subdivides them until none exceed it, so minArea is an upper bound. The
    help previously said the opposite -- that smaller facets were discarded."""
    from pHinder.gui import help_text

    body = help_text.for_option("minArea")[1].lower()
    assert "larger" in body
    assert "not a floor" in body or "maximum" in body


def test_probe_help_warns_about_fragmentation():
    """Too small a radius fragments the surface, and goFoSurface() then drops
    components under the size threshold -- the user needs to know why parts of
    the surface can vanish."""
    from pHinder.gui import help_text

    body = help_text.for_option("circumSphereRadiusLimit")[1].lower()
    assert "probe" in body
    assert "fragment" in body
    assert "allow small surfaces" in body


CIF_TWO_CHAINS = """data_model
loop_
_entity.id
_entity.type
1 polymer
2 polymer
#
loop_
_struct_asym.entity_id
_struct_asym.id
1 A
2 B
#
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.auth_asym_id
ATOM 1 N MET A 1 1 1.0 2.0 3.0 A
ATOM 2 CA MET A 1 1 1.0 2.0 3.0 A
ATOM 3 N GLY B 2 1 4.0 5.0 6.0 B
#
"""

PDB_TWO_CHAINS = (
    "ATOM      1  N   MET A   1      -1.000   2.000   3.000  1.00  0.00           N\n"
    "ATOM      2  CA  MET A   1      -1.000   2.000   3.000  1.00  0.00           C\n"
    "ATOM      3  N   GLY B   1       4.000   5.000   6.000  1.00  0.00           N\n"
)


def test_chain_ids_from_mmcif():
    """The reader used to slice PDB column 21, which finds nothing in mmCIF."""
    from pHinder.gui.app import chain_ids

    assert chain_ids(CIF_TWO_CHAINS) == ["A", "B"]


def test_entity_ids_are_not_mistaken_for_chains():
    """A two-chain model carries entity ids 1 and 2 alongside chains A and B.
    Reading the wrong column is what makes it look like four chains."""
    from pHinder.gui.app import chain_ids

    ids = chain_ids(CIF_TWO_CHAINS)
    assert "1" not in ids and "2" not in ids
    assert len(ids) == 2


def test_author_chain_id_wins_over_the_label():
    """auth_asym_id is the chain a viewer shows and what people mean by 'chain A'."""
    from pHinder.gui.app import chain_ids

    text = CIF_TWO_CHAINS.replace("3.0 A", "3.0 X").replace("6.0 B", "6.0 Y")
    assert chain_ids(text) == ["X", "Y"]


def test_chain_ids_from_pdb_still_work():
    from pHinder.gui.app import chain_ids

    assert chain_ids(PDB_TWO_CHAINS) == ["A", "B"]


def test_gzipped_structures_are_read(tmp_path):
    """The runner sets zip=1 for .gz, so the chain picker must read them too."""
    import gzip

    from pHinder.gui.app import _read_structure_text, chain_ids

    path = tmp_path / "model.cif.gz"
    with gzip.open(path, "wt") as fh:
        fh.write(CIF_TWO_CHAINS)
    assert chain_ids(_read_structure_text(str(path))) == ["A", "B"]
