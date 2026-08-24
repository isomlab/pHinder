"""The right-hand progress panel: what pHinder is doing, and how far in.

A pHinder run is a sequence of independent calculations the user opted into, so
progress is reported as a checklist of those stages plus a live console. The
panel is deliberately passive -- it renders whatever the runner reports and
never reaches back into the calculation.

Every public method is safe to call from the worker thread: each one marshals
onto the Tk thread with ``after``, because Tk is not thread-safe and pHinder
runs off the main thread.
"""

import queue
import tkinter as tk
from tkinter import ttk

from pHinder.gui import theme

PENDING, ACTIVE, DONE, FAILED, SKIPPED = "pending", "active", "done", "failed", "skipped"

_MARK = {
    PENDING: ("○", theme.MUTED),
    ACTIVE:  ("▶", theme.RUNNING),
    DONE:    ("✓", theme.OK),
    FAILED:  ("✕", theme.FAIL),
    SKIPPED: ("–", theme.MUTED),
}


class ProgressPanel(ttk.Frame):
    def __init__(self, parent, fonts, on_run, on_cancel):
        super().__init__(parent, style="TFrame", padding=(16, 12, 16, 12))
        self.fonts = fonts
        self._on_run, self._on_cancel = on_run, on_cancel
        self._rows = {}                  # stage key -> (mark label, text label)
        self._q = queue.Queue()          # worker thread -> Tk thread
        self._running = False

        ttk.Label(self, text="Progress", style="Section.TLabel").pack(anchor="w")

        self._status = ttk.Label(self, text="Idle", style="Status.TLabel")
        self._status.pack(anchor="w", pady=(6, 2))

        self._bar = ttk.Progressbar(self, mode="determinate", maximum=1)
        self._bar.pack(fill="x", pady=(0, 2))

        self._detail = ttk.Label(self, text="Choose a structure and at least one calculation.",
                                 style="StatusMuted.TLabel", wraplength=380, justify="left")
        self._detail.pack(anchor="w", pady=(0, 10))

        self._stages = ttk.Frame(self, style="Card.TFrame", padding=10)
        self._stages.pack(fill="x")
        self._stages.columnconfigure(1, weight=1)

        ttk.Label(self, text="Output", style="Section.TLabel").pack(anchor="w", pady=(14, 4))
        wrap = ttk.Frame(self, style="TFrame")
        wrap.pack(fill="both", expand=True)
        self._log = tk.Text(wrap, wrap="word", state="disabled", height=10,
                            background=theme.CONSOLE_BG, foreground=theme.CONSOLE_FG,
                            insertbackground=theme.CONSOLE_FG, font=fonts.mono,
                            relief="flat", padx=10, pady=8)
        sb = ttk.Scrollbar(wrap, orient="vertical", command=self._log.yview)
        self._log.configure(yscrollcommand=sb.set)
        self._log.pack(side="left", fill="both", expand=True)
        sb.pack(side="right", fill="y")
        self._log.tag_configure("err", foreground="#ff9b8a")

        buttons = ttk.Frame(self, style="TFrame")
        buttons.pack(fill="x", pady=(12, 0))
        self.run_button = ttk.Button(buttons, text="Run pHinder", style="Accent.TButton",
                                     command=self._run_clicked)
        self.run_button.pack(side="left")
        self.cancel_button = ttk.Button(buttons, text="Stop", style="Stop.TButton",
                                        command=self._cancel_clicked, state="disabled")
        self.cancel_button.pack(side="left", padx=(8, 0))
        ttk.Button(buttons, text="Clear", command=self.clear_log).pack(side="right")

        self.after(30, self._drain)

    # --- stages -------------------------------------------------------------
    def set_stages(self, stages):
        """Render the checklist. `stages` is [(key, label), ...] in run order."""
        for child in self._stages.winfo_children():
            child.destroy()
        self._rows.clear()
        if not stages:
            ttk.Label(self._stages, text="No calculations selected.",
                      style="Help.TLabel").grid(row=0, column=0, sticky="w")
            self._bar.configure(maximum=1, value=0)
            return
        for row, (key, label) in enumerate(stages):
            mark = tk.Label(self._stages, text=_MARK[PENDING][0], fg=_MARK[PENDING][1],
                            bg=theme.CARD, font=self.fonts.label, width=2)
            mark.grid(row=row, column=0, sticky="w")
            text = tk.Label(self._stages, text=label, fg=theme.MUTED, bg=theme.CARD,
                            font=self.fonts.label, anchor="w")
            text.grid(row=row, column=1, sticky="w")
            self._rows[key] = (mark, text)
        self._bar.configure(maximum=len(stages), value=0)

    def set_total_steps(self, n):
        """Size the bar in sub-steps so it advances during a long stage."""
        self._q.put(("total", (max(1, n),)))

    # --- thread-safe reporting ---------------------------------------------
    def stage(self, key, state, note=""):
        self._q.put(("stage", (key, state, note)))

    def status(self, text, detail=""):
        self._q.put(("status", (text, detail)))

    def substep(self, text):
        self._q.put(("substep", (text,)))

    def write(self, text, error=False):
        self._q.put(("log", (text, error)))

    def finished(self, ok=True, note=""):
        self._q.put(("finished", (ok, note)))

    # --- Tk-thread handlers -------------------------------------------------
    def _drain(self):
        try:
            while True:
                kind, payload = self._q.get_nowait()
                getattr(self, f"_apply_{kind}")(*payload)
        except queue.Empty:
            pass
        self.after(30, self._drain)

    def _apply_stage(self, key, state, note):
        row = self._rows.get(key)
        if not row:
            return
        mark, text = row
        glyph, colour = _MARK[state]
        mark.configure(text=glyph, fg=colour)
        text.configure(fg=theme.TEXT if state in (ACTIVE, DONE) else theme.MUTED)
        if note:
            text.configure(text=f"{text.cget('text').split('  —')[0]}  — {note}")

    def _apply_total(self, n):
        self._bar.configure(maximum=n, value=0)

    def _apply_substep(self, text):
        self._bar.configure(value=min(self._bar["value"] + 1, self._bar["maximum"]))
        self._detail.configure(text=text)

    def _apply_status(self, text, detail):
        self._status.configure(text=text)
        if detail:
            self._detail.configure(text=detail)

    def _apply_log(self, text, error):
        self._log.configure(state="normal")
        self._log.insert("end", text, ("err",) if error else ())
        self._log.see("end")
        self._log.configure(state="disabled")

    def _apply_finished(self, ok, note):
        self._running = False
        self.run_button.configure(state="normal")
        self.cancel_button.configure(state="disabled")
        self._status.configure(text="Finished" if ok else "Stopped")
        self._detail.configure(text=note or ("Run complete." if ok else "Run did not complete."))
        if ok:
            self._bar.configure(value=self._bar["maximum"])

    # --- buttons ------------------------------------------------------------
    def _run_clicked(self):
        if self._running:
            return
        self._running = True
        self.run_button.configure(state="disabled")
        self.cancel_button.configure(state="normal")
        self._status.configure(text="Running…")
        self._on_run()

    def _cancel_clicked(self):
        self._detail.configure(text="Stopping after the current step…")
        self._on_cancel()

    def clear_log(self):
        self._log.configure(state="normal")
        self._log.delete("1.0", "end")
        self._log.configure(state="disabled")
