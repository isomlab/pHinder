"""The shared Isom-lab tool look, applied to pHinder.

Palette, fonts and ttk styles here match cellog / scopelog / probelog /
plasmidlog so the tools read as one family. Kept in its own module because it is
generic: nothing in here knows anything about pHinder, and it is the natural
thing to lift into ``isomlab.tkwidgets`` once a second tool wants it.
"""

import sys
import tkinter as tk
from tkinter import ttk
import tkinter.font as tkfont

# Palette (shared look with the other lab tools) ------------------------------
BG = "#eef2f5"
CARD = "#ffffff"
HEADER_BG = "#1f3a5f"
HEADER_FG = "#ffffff"
SUBTITLE_FG = "#c3d0de"
TEXT = "#1f2a36"
MUTED = "#6b7a8d"
ACCENT = "#2e7d32"
ACCENT_ACTIVE = "#256628"

# pHinder runs a calculation rather than storing records, so it needs a few
# states the record tools never show.
STOP = "#b3261e"
STOP_ACTIVE = "#8f1d17"
OK = "#2e7d32"
FAIL = "#b3261e"
RUNNING = "#1f3a5f"
CONSOLE_BG = "#101820"
CONSOLE_FG = "#dbe4ee"


class Fonts:
    """Named fonts at the sizes the other tools use."""

    def __init__(self):
        fam = tkfont.nametofont("TkDefaultFont").actual("family")
        mono = "Menlo" if sys.platform == "darwin" else "Consolas"
        self.title = tkfont.Font(family=fam, size=19, weight="bold")
        self.sub = tkfont.Font(family=fam, size=11)
        self.section = tkfont.Font(family=fam, size=13, weight="bold")
        self.label = tkfont.Font(family=fam, size=11)
        self.help = tkfont.Font(family=fam, size=9)
        self.mono = tkfont.Font(family=mono, size=10)


def apply_style(fonts):
    """Configure the ttk styles every island screen is built from."""
    st = ttk.Style()
    try:
        st.theme_use("clam")
    except tk.TclError:
        pass
    st.configure("TFrame", background=BG)
    st.configure("Card.TFrame", background=CARD)
    st.configure("Header.TFrame", background=HEADER_BG)
    st.configure("Title.TLabel", background=HEADER_BG, foreground=HEADER_FG, font=fonts.title)
    st.configure("Sub.TLabel", background=HEADER_BG, foreground=SUBTITLE_FG, font=fonts.sub)
    st.configure("Version.TLabel", background=HEADER_BG, foreground=SUBTITLE_FG, font=fonts.help)
    st.configure("Section.TLabel", background=BG, foreground=HEADER_BG, font=fonts.section)
    st.configure("Field.TLabel", background=CARD, foreground=TEXT, font=fonts.label)
    st.configure("Help.TLabel", background=CARD, foreground=MUTED, font=fonts.help)
    st.configure("CardHelp.TLabel", background=BG, foreground=MUTED, font=fonts.help)
    st.configure("TLabelframe", background=CARD, bordercolor="#d7dee6", borderwidth=1)
    st.configure("TLabelframe.Label", background=CARD, foreground=HEADER_BG, font=fonts.label)
    st.configure("TCheckbutton", background=CARD, font=fonts.label)
    st.configure("TLabel", background=CARD, foreground=TEXT, font=fonts.label)
    st.configure("TEntry", fieldbackground="#ffffff")
    st.configure("Card.TCheckbutton", background=CARD, font=fonts.label)
    st.configure("Card.TRadiobutton", background=CARD, font=fonts.label)
    st.configure("Accent.TButton", font=fonts.label, foreground="#fff",
                 background=ACCENT, padding=(16, 8), borderwidth=0)
    st.map("Accent.TButton", background=[("active", ACCENT_ACTIVE),
                                         ("disabled", "#9bb39d")])
    st.configure("Stop.TButton", font=fonts.label, foreground="#fff",
                 background=STOP, padding=(16, 8), borderwidth=0)
    st.map("Stop.TButton", background=[("active", STOP_ACTIVE),
                                       ("disabled", "#c9a5a2")])
    st.configure("Status.TLabel", background=BG, foreground=TEXT, font=fonts.label)
    st.configure("StatusMuted.TLabel", background=BG, foreground=MUTED, font=fonts.help)
    # clam's stock notebook and progressbar are khaki; recolour both so the
    # panel reads as the same surface as the rest of the window.
    st.configure("TNotebook", background=BG, borderwidth=0, tabmargins=(2, 6, 2, 0))
    st.configure("TNotebook.Tab", font=fonts.label, padding=(14, 8),
                 background=BG, foreground=MUTED, borderwidth=0)
    st.map("TNotebook.Tab",
           background=[("selected", CARD), ("active", "#e3e9ef")],
           foreground=[("selected", HEADER_BG)],
           expand=[("selected", (0, 0, 0, 2))])
    st.configure("TProgressbar", troughcolor="#dfe6ec", bordercolor="#dfe6ec",
                 background=ACCENT, lightcolor=ACCENT, darkcolor=ACCENT,
                 thickness=10)
    st.configure("Vertical.TScrollbar", background="#cfd8e2", troughcolor=BG,
                 bordercolor=BG, arrowcolor=MUTED)
    return st


def header(parent, name, subtitle, version=""):
    """The navy title bar every island tool opens with."""
    bar = ttk.Frame(parent, style="Header.TFrame")
    bar.pack(fill="x")
    inner = ttk.Frame(bar, style="Header.TFrame")
    inner.pack(fill="x", padx=20, pady=12)
    ttk.Label(inner, text=name, style="Title.TLabel").pack(side="left")
    ttk.Label(inner, text=f"  {subtitle}", style="Sub.TLabel").pack(side="left", padx=(8, 0))
    if version:
        ttk.Label(inner, text=version, style="Version.TLabel").pack(side="right")
    return bar


def section(parent, title, help_text=""):
    """A section heading over a white card. Returns the card to fill."""
    ttk.Label(parent, text=title, style="Section.TLabel").pack(anchor="w", pady=(14, 4))
    if help_text:
        ttk.Label(parent, text=help_text, style="CardHelp.TLabel",
                  wraplength=560, justify="left").pack(anchor="w", pady=(0, 6))
    card = ttk.Frame(parent, style="Card.TFrame", padding=12)
    card.pack(fill="x")
    card.columnconfigure(1, weight=1)
    return card


class ScrollHost:
    """Scrollable bodies plus one central wheel binding that routes to whichever
    one the pointer is over.

    A ``bind_all`` per canvas would have each new scroll region clobber the last
    (the pam_scanning lesson), so the owner installs the bindings once and every
    body registers itself here.
    """

    def __init__(self, root):
        self.root = root
        self._canvases = []
        self.bodies = []      # inner frames, for measuring what the pane needs
        root.bind_all("<MouseWheel>", self._route)
        root.bind_all("<Button-4>", lambda e: self._route(e, -1))
        root.bind_all("<Button-5>", lambda e: self._route(e, 1))
        # ttk Comboboxes/Spinboxes change their VALUE on wheel by default, so a
        # user scrolling the page silently edits a parameter. Scroll instead.
        for cls in ("TCombobox", "TSpinbox"):
            root.bind_class(cls, "<MouseWheel>", self._guard)
            root.bind_class(cls, "<Button-4>", lambda e: self._guard(e, -1))
            root.bind_class(cls, "<Button-5>", lambda e: self._guard(e, 1))

    def body(self, parent, padx=16):
        canvas = tk.Canvas(parent, background=BG, highlightthickness=0)
        sb = ttk.Scrollbar(parent, orient="vertical", command=canvas.yview)
        inner = ttk.Frame(canvas, style="TFrame", padding=(padx, 4, padx, 16))
        inner.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        win = canvas.create_window((0, 0), window=inner, anchor="nw")
        canvas.bind("<Configure>", lambda e: canvas.itemconfigure(win, width=e.width))
        canvas.configure(yscrollcommand=sb.set)
        canvas.pack(side="left", fill="both", expand=True)
        sb.pack(side="right", fill="y")
        self._canvases.append(canvas)
        self.bodies.append(inner)
        return inner

    def _route(self, event, units=None):
        widget = self.root.winfo_containing(event.x_root, event.y_root)
        while widget is not None:
            if widget in self._canvases:
                widget.yview_scroll(units if units is not None else self._units(event), "units")
                return
            widget = getattr(widget, "master", None)

    def _guard(self, event, units=None):
        self._route(event, units)
        return "break"

    @staticmethod
    def _units(event):
        # macOS reports small per-notch deltas; Windows/X11 report multiples of 120.
        if sys.platform == "darwin":
            return -1 if event.delta > 0 else 1
        return int(-event.delta / 120) or (-1 if event.delta > 0 else 1)
