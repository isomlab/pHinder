# pHinder

Ionizable-residue network and surface analysis for protein structures.

This is pHinder packaged as an installable tool, **re-based onto the shared
[`isomlab`](https://github.com/isomlab/isomlab) core**: the computational-geometry
and PDB/mmCIF routines it used to bundle now come from `isomlab` as a dependency,
so only the pHinder algorithm, GUI, and CLI live in this repo.

## Install

```bash
pip install isomlab pHinder        # once published
# or from checkouts (install isomlab first):
pip install -e ../isomlab && pip install -e .
```

Requires Python ≥ 3.9, `isomlab`, and `openpyxl`.

## Use

```bash
phinder --help          # command-line interface
phinder-gui             # Tkinter GUI
```

## Layout

```
pHinder/
├── pHinder_7_0.py      the algorithm (imports geometry/structure from isomlab)
├── command_line.py     CLI entry point  ->  `phinder`
├── single_script.py    standalone batch script
└── gui/                Tkinter front-end ->  `phinder-gui`
```

## Provenance

Extracted from the legacy `pythonScripts/pHinder` source (active code only; the
2012–2013 `z_archives` were dropped). Algorithm bodies are unchanged — only import
statements were rewritten from bundled bare-module names to `isomlab.*` package
paths. PDB and mmCIF parsing are both handled by `isomlab.structure.pdbFile`.

## License

MIT © Daniel G. Isom (Isom Lab)
