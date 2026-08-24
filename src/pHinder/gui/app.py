"""pHinder's main window: tabbed parameters on the left, progress on the right.

Replaces the single scrolling column of collapsible panes. The parameter groups
that were panes are now tabs, and the console that was a bare black text box is
now a progress panel that stays visible whichever tab you are on -- so you can
watch a run without losing your place in the settings.

The calculation itself is untouched: this module only collects parameters and
hands them to the existing runner.
"""

import threading
import tkinter as tk
from tkinter import ttk, messagebox

from pHinder.gui import theme
from pHinder.gui.progress import ProgressPanel, ACTIVE, DONE, FAILED, SKIPPED
from pHinder.gui.runner import GROUP_CHAINS, NEEDS
from pHinder.gui.file_open import FilePathWidget
from pHinder.gui.dynamic_option_widget_amino_acid_selection import AminoAcidSelectionWidget

try:
    from pHinder import __version__ as VERSION
except Exception:                                    # running from a source tree
    VERSION = ""

# The pHinder set of ionizable groups.
DEFAULT_AA = ["Aspartic Acid (D)", "Glutamic Acid (E)", "Lysine (K)",
              "Arginine (R)", "Histidine (H)"]

# The five calculations, in the order the runner performs them. The keys are the
# ones the runner reads out of results["calculation_options"].
CALCULATIONS = [
    ("topologyCalculation", "Residue network topology",
     "Triangulate the selected residues, prune to a network, and analyse its topology."),
    ("surfaceCalculation", "Molecular surface",
     "Compute the surface used for depth and exposure measurements."),
    ("sidechainClassification", "Sidechain classification",
     "Classify each selected sidechain as core, margin, or surface. Needs the "
     "triangulation and the surface, and will compute them if not ticked above."),
    ("interfaceClassification", "Interface classification",
     "Identify interface sidechains between the selected chains. Needs the "
     "classification, and will compute it if not ticked above."),
    ("virtualScreenSurfacesCalculation", "Virtual screening surfaces",
     "Grid the surface, remove clashes, and parse void volumes for screening."),
]

# A parameter tab is shown when the block that reads it will actually run.
# Keyed on the block rather than the checkbox because the blocks are shared:
# ticking only sidechain classification still runs the triangulation and the
# surface, and those tabs hold the parameters being applied.
TAB_REQUIRES = {
    "Classification": "classify",
    "Networks": "triangulate",
    "Surfaces": "surface",
    "Interfaces": "interface",
    "Screening": "screen",
}

# Parameter groups -> the tab each belongs on.
TAB_FOR_GROUP = {
    "sidechain_classification_options": "Classification",
    "network_options": "Networks",
    "surface_options": "Surfaces",
    "interface_options": "Interfaces",
    "virtual_screening_options": "Screening",
    "advanced_options": "Advanced",
}


def prettify(key):
    """CONSTANT_CASE / camelCase -> a sentence-case label.

    Purely mechanical: no attempt is made to guess units or meaning, so a label
    never claims something the constant does not say.
    """
    if key.isupper() or "_" in key:
        words = key.replace("_", " ").lower().split()
    else:
        out, cur = [], ""
        for ch in key:
            if ch.isupper() and cur:
                out.append(cur); cur = ch.lower()
            else:
                cur += ch
        out.append(cur)
        words = out
    return " ".join(words).capitalize()


class PHinderApp(tk.Tk):
    def __init__(self, defaults, runner=None):
        """`defaults` maps group name -> {option: default}. `runner` is called on
        a worker thread as runner(results, report) and defaults to the real one."""
        super().__init__()
        self.title("pHinder electroinformatics for understanding how protons regulate protein structure-function relationships")
        self.configure(bg=theme.BG)
        self.minsize(1180, 760)
        sw, sh = self.winfo_screenwidth(), self.winfo_screenheight()
        self.geometry(f"{min(1520, sw - 80)}x{min(940, sh - 120)}")

        self.defaults = defaults
        self._runner = runner
        self._cancel = threading.Event()
        self.vars = {group: {} for group in defaults}

        self.fonts = theme.Fonts()
        theme.apply_style(self.fonts)
        self.scroll = theme.ScrollHost(self)

        theme.header(self, "pHinder",
                     "electroinformatics for understanding how protons regulate "
                     "protein structure-function relationships",
                     f"v{VERSION}" if VERSION else "")

        split = ttk.PanedWindow(self, orient="horizontal")
        split.pack(fill="both", expand=True)

        left = ttk.Frame(split, style="TFrame")
        self.progress = ProgressPanel(split, self.fonts, self._start_run, self._cancel_run)
        split.add(left, weight=3)
        split.add(self.progress, weight=2)

        self._split = split
        self._build_tabs(left)
        self._refresh_stages()
        # Size the split from what the tabs actually need. ttk places the sash
        # from the panes' requested widths, and a Canvas does not report the
        # width of the frame scrolling inside it -- so the left pane asked for
        # 301px while the file row alone needs 661, and Browse was clipped.
        self.after_idle(self._fit_layout)

    def _fit_layout(self):
        """Open wide enough for the widest tab, with the sash placed to match."""
        self.update_idletasks()
        SCROLLBAR, MARGIN = 18, 10
        needed = max((b.winfo_reqwidth() for b in self.scroll.bodies), default=600)
        needed += SCROLLBAR + MARGIN

        progress_min = max(self.progress.winfo_reqwidth(), 430)
        wanted = needed + progress_min
        screen = self.winfo_screenwidth() - 80
        if wanted > self.winfo_width():
            width = min(wanted, screen)
            self.geometry(f"{width}x{self.winfo_height()}")
            self.update_idletasks()

        # Leave the progress panel its minimum if the screen could not fit both.
        sash = min(needed, max(320, self.winfo_width() - progress_min))
        try:
            self._split.sashpos(0, int(sash))
        except tk.TclError:
            pass

    # --- tabs ---------------------------------------------------------------
    def _build_tabs(self, parent):
        nb = ttk.Notebook(parent)
        nb.pack(fill="both", expand=True)
        self._nb = nb

        self._tabs = {}
        for title, builder in (
            ("Input", self._build_input_tab),
            ("Calculations", self._build_calculations_tab),
            ("Classification", lambda b: self._build_group_tab(b, "sidechain_classification_options")),
            ("Networks", lambda b: self._build_group_tab(b, "network_options")),
            ("Surfaces", lambda b: self._build_group_tab(b, "surface_options")),
            ("Interfaces", lambda b: self._build_group_tab(b, "interface_options")),
            ("Screening", lambda b: self._build_group_tab(b, "virtual_screening_options")),
            ("Advanced", lambda b: self._build_group_tab(b, "advanced_options")),
        ):
            frame = ttk.Frame(nb, style="TFrame")
            nb.add(frame, text=title)
            self._tabs[title] = frame
            builder(self.scroll.body(frame))

    def _build_input_tab(self, body):
        card = theme.section(body, "Structure",
                             "Choose a PDB file, then the chains to include and where results go.")
        self.file_widget = FilePathWidget(card, process_file=self._read_chains,
                                          options_label="Chains:",
                                          extra_options={"Group Chains": 0})
        self.file_widget.frame.pack(fill="x", anchor="w")

        card = theme.section(body, "Residues",
                             "pHinder defaults to the ionizable set: D, E, K, R and H.")
        self.aa_widget = AminoAcidSelectionWidget(card, "Amino acids",
                                                  default_selections=DEFAULT_AA)
        self.aa_widget.frame.pack(fill="x", anchor="w")

    def _build_calculations_tab(self, body):
        card = theme.section(body, "Calculations to run",
                             "Each one you tick becomes a step in the progress list on the right.")
        self.vars["calculation_options"] = {}
        for row, (key, label, blurb) in enumerate(CALCULATIONS):
            var = tk.IntVar(value=self.defaults["calculation_options"].get(key, 0))
            var.trace_add("write", lambda *_: self._refresh_stages())
            self.vars["calculation_options"][key] = var
            ttk.Checkbutton(card, text=label, variable=var,
                            style="Card.TCheckbutton").grid(row=row * 2, column=0,
                                                            sticky="w", pady=(6, 0))
            ttk.Label(card, text=blurb, style="Help.TLabel", wraplength=520,
                      justify="left").grid(row=row * 2 + 1, column=0, sticky="w",
                                           padx=(22, 0), pady=(0, 4))

    def _build_group_tab(self, body, group):
        defaults = self.defaults[group]
        card = theme.section(body, prettify(group.replace("_options", "")) + " parameters")
        row = 0
        for key, default in defaults.items():
            var, widget = self._field(card, key, default, row)
            self.vars[group][key] = var
            row += 1
        if not defaults:
            ttk.Label(card, text="No parameters in this group.",
                      style="Help.TLabel").grid(row=0, column=0, sticky="w")
        actions = ttk.Frame(body, style="TFrame")
        actions.pack(fill="x", pady=(12, 0))
        ttk.Button(actions, text="Reset to defaults",
                   command=lambda g=group: self._reset(g)).pack(side="left")

    def _field(self, card, key, default, row):
        """Checkbox for 0/1 flags, entry for numbers -- matching how the panes did it."""
        label = prettify(key)
        if isinstance(default, int) and default in (0, 1) and not isinstance(default, bool):
            var = tk.IntVar(value=default)
            w = ttk.Checkbutton(card, text=label, variable=var, style="Card.TCheckbutton")
            w.grid(row=row, column=0, columnspan=2, sticky="w", pady=3)
            return var, w
        var = tk.DoubleVar(value=default) if isinstance(default, float) else tk.IntVar(value=default)
        ttk.Label(card, text=label, style="Field.TLabel").grid(row=row, column=0,
                                                               sticky="w", padx=(0, 12), pady=3)
        w = ttk.Entry(card, textvariable=var, width=14)
        w.grid(row=row, column=1, sticky="w", pady=3)
        return var, w

    def _reset(self, group):
        for key, var in self.vars[group].items():
            var.set(self.defaults[group][key])

    # --- input helpers ------------------------------------------------------
    def _read_chains(self, file_path):
        """Chain ids from a PDB, for the chain checkboxes."""
        chains = []
        try:
            with open(file_path) as fh:
                for line in fh:
                    if line.startswith(("ATOM", "HETATM")) and len(line) > 21:
                        cid = line[21]
                        if cid.strip() and cid not in chains:
                            chains.append(cid)
        except OSError as exc:
            messagebox.showerror("Could not read file", str(exc))
        return chains

    # --- progress wiring ----------------------------------------------------
    def _selected_calculations(self):
        return [(k, label) for k, label, _ in CALCULATIONS
                if self.vars.get("calculation_options", {}).get(k)
                and self.vars["calculation_options"][k].get()]

    def _active_blocks(self):
        """Blocks that will run, given what is ticked."""
        blocks = set()
        for key, needed in NEEDS.items():
            var = self.vars.get("calculation_options", {}).get(key)
            if var and var.get():
                blocks.update(needed)
        return blocks

    def _refresh_tabs(self):
        """Hide parameter tabs whose parameters nothing will read."""
        if not getattr(self, "_tabs", None):
            return
        active = self._active_blocks()
        for title, block in TAB_REQUIRES.items():
            frame = self._tabs.get(title)
            if frame is None:
                continue
            state = "normal" if block in active else "hidden"
            try:
                self._nb.tab(frame, state=state)
            except tk.TclError:
                pass

    def _refresh_stages(self):
        self._refresh_tabs()
        chosen = self._selected_calculations()
        self.progress.set_stages(chosen)
        if chosen:
            self.progress.status("Idle", f"{len(chosen)} calculation(s) selected.")
        else:
            self.progress.status("Idle", "No calculations selected — tick at least one.")

    def collect(self):
        """Everything the runner needs, in the shape it already expects."""
        results = {group: {k: v.get() for k, v in vars_.items()}
                   for group, vars_ in self.vars.items()}
        results["file_path"] = self.file_widget.get_file_path()
        results["save_path"] = self.file_widget.get_save_path() or ""
        results["chains"] = self.file_widget.get_values()
        results["amino_acid_selections"] = self.aa_widget.get_values()
        return results

    def _validate(self, results):
        if not results["file_path"]:
            return "Choose a PDB file on the Input tab."
        if not any(on for c, on in results["chains"].items() if c != GROUP_CHAINS):
            return "Select at least one chain on the Input tab."
        if not any(bool(v) for v in results["amino_acid_selections"].values()):
            return "Select at least one amino acid on the Input tab."
        if not self._selected_calculations():
            return "Tick at least one calculation on the Calculations tab."
        return None

    def _start_run(self):
        results = self.collect()
        problem = self._validate(results)
        if problem:
            messagebox.showwarning("Nothing to run", problem)
            self.progress.finished(ok=False, note=problem)
            return

        self._cancel.clear()
        self.progress.clear_log()
        self.progress.status("Running…", "Starting pHinder.")
        for key, _ in self._selected_calculations():
            self.progress.stage(key, "pending")

        runner = self._runner or _default_runner
        threading.Thread(target=self._run, args=(runner, results), daemon=True).start()

    def _run(self, runner, results):
        report = _Report(self.progress, self._cancel)
        try:
            runner(results, report)
            self.progress.finished(ok=not self._cancel.is_set())
        except Exception as exc:                       # surfaced, never swallowed
            self.progress.write(f"\nError: {exc}\n", error=True)
            self.progress.finished(ok=False, note=str(exc))

    def _cancel_run(self):
        self._cancel.set()
        self.progress.write("\nStop requested — finishing the current step.\n")


class _Report:
    """What a runner is handed: stage transitions, log lines, and a stop flag."""

    def __init__(self, panel, cancel):
        self._panel, self._cancel = panel, cancel

    def begin(self, key):
        self._panel.stage(key, ACTIVE)

    def done(self, key, note=""):
        self._panel.stage(key, DONE, note)

    def failed(self, key, note=""):
        self._panel.stage(key, FAILED, note)

    def skipped(self, key, note=""):
        self._panel.stage(key, SKIPPED, note)

    def status(self, text, detail=""):
        self._panel.status(text, detail)

    def set_total_steps(self, n):
        self._panel.set_total_steps(n)

    def substep(self, text):
        self._panel.substep(text)

    def write(self, text, error=False):
        self._panel.write(text if text.endswith("\n") else text + "\n", error)

    @property
    def cancelled(self):
        return self._cancel.is_set()


def _default_runner(results, report):
    from pHinder.gui.runner import run
    run(results, report)


def main():
    from pHinder.gui import phinder_main_gui as legacy
    defaults = {
        "calculation_options": legacy.default_calculation_options,
        "sidechain_classification_options": legacy.default_sidechain_classification_options,
        "network_options": legacy.default_network_options,
        "surface_options": legacy.default_surface_options,
        "interface_options": legacy.default_interface_options,
        "virtual_screening_options": legacy.default_virtual_screening_options,
        "advanced_options": legacy.default_advanced_options,
    }
    PHinderApp(defaults).mainloop()


if __name__ == "__main__":
    main()
