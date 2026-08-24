"""Hover help.

A tooltip appears after a short delay, stays while the pointer is over the
widget or the tip itself, and is clamped to the screen. Text is supplied by
help_text.py so the wording lives with the content, not the plumbing.
"""

import tkinter as tk
from tkinter import ttk

from pHinder.gui import theme

DELAY_MS = 450
MAX_WIDTH = 380
# Offset from the pointer, far enough that the tip never lands under the cursor
# (which would fire <Leave> on the widget and flicker the tip on and off).
CURSOR_DX, CURSOR_DY = 16, 22


class _Tip:
    """One shared popup, reused by every widget that registers for help."""

    def __init__(self, root):
        self.root = root
        self._win = None
        self._after = None
        self._title = None
        self._body = None

    def schedule(self, widget, title, body):
        self.cancel()
        self._after = widget.after(DELAY_MS, lambda: self.show(widget, title, body))

    def cancel(self):
        if self._after is not None:
            try:
                self.root.after_cancel(self._after)
            except Exception:
                pass
            self._after = None

    def show(self, widget, title, body):
        self.hide()
        if not widget.winfo_exists():
            return
        win = tk.Toplevel(self.root)
        win.wm_overrideredirect(True)
        try:
            win.wm_attributes("-topmost", True)
        except tk.TclError:
            pass
        frame = tk.Frame(win, background=theme.CARD,
                         highlightbackground=theme.RULE_STRONG,
                         highlightthickness=1, padx=12, pady=9)
        frame.pack()
        if title:
            tk.Label(frame, text=title, background=theme.CARD, foreground=theme.HEADER_BG,
                     font=self.root.fonts.label, justify="left", anchor="w").pack(anchor="w")
        tk.Label(frame, text=body, background=theme.CARD, foreground=theme.TEXT,
                 font=self.root.fonts.help, justify="left", anchor="w",
                 wraplength=MAX_WIDTH).pack(anchor="w", pady=(3 if title else 0, 0))

        win.update_idletasks()
        w, h = win.winfo_reqwidth(), win.winfo_reqheight()

        # Anchor to the pointer, not to the widget's box. Help is attached to
        # whole section cards as well as to single fields, and the Residues card
        # is 500px tall -- anchoring to its bottom edge put the tip half a
        # screen from whatever the user was actually pointing at.
        px, py = widget.winfo_pointerxy()
        x, y = px + CURSOR_DX, py + CURSOR_DY

        sw, sh = widget.winfo_screenwidth(), widget.winfo_screenheight()
        if x + w > sw - 8:
            x = max(8, px - w - 8)          # flip to the left of the pointer
        if y + h > sh - 8:
            y = max(8, py - h - 10)         # flip above the pointer
        win.wm_geometry(f"+{int(x)}+{int(y)}")
        win.lift()
        win.update_idletasks()
        self._win = win

    def hide(self):
        self.cancel()
        if self._win is not None:
            try:
                self._win.destroy()
            except Exception:
                pass
            self._win = None


def attach(root, widget, body, title=""):
    """Give `widget` hover help. `root` must carry a `_tip` and `fonts`."""
    if not body:
        return
    tip = getattr(root, "_tip", None)
    if tip is None:
        tip = root._tip = _Tip(root)
    widget.bind("<Enter>", lambda e: tip.schedule(widget, title, body), add="+")
    widget.bind("<Leave>", lambda e: tip.hide(), add="+")
    widget.bind("<ButtonPress>", lambda e: tip.hide(), add="+")
