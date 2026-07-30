# pHinder

Ionizable-residue network and surface analysis for protein structures.

This is pHinder packaged as an installable tool. It is **self-contained**: the
algorithm, GUI, CLI, and all the computational geometry it needs live in this
repository, so installing it requires no other lab repo.

## Install

**conda** (recommended — brings its own Python):

```bash
conda env create -f environment.yml     # from a checkout of this repo
conda activate pHinder
```

**pip**, from GitHub — pHinder is not on PyPI:

```bash
pip install git+https://github.com/isomlab/pHinder
```

Or from a local checkout:

```bash
pip install -e .
```

**No installing at all** — the double-click launchers in [`launchers/`](launchers) do the
conda setup for you (install [Miniforge](https://conda-forge.org/download/) once, then
double-click **Launch pHinder**: `.command` on Mac, `.bat` on Windows).

Requires Python ≥ 3.9, `numpy`, and `openpyxl`. No other lab repository is needed.

> **Do not run `pip install pHinder`.** That name belongs to an unrelated project on
> PyPI; it will not install this tool. Install from this repository instead.

## Use

```bash
phinder --help          # command-line interface
phinder-gui             # Tkinter GUI
```

## Layout

```
pHinder/
├── pHinder_7_0.py      the algorithm
├── command_line.py     CLI entry point  ->  `phinder`
├── single_script.py    standalone batch script
├── gui/                Tkinter front-end ->  `phinder-gui`
├── geometry/           convex hulls, spheres, network walks (pHinder's own)
├── structure/          surface construction and structure output
└── _vendor/            shared with isomlab: compGeometry, determinants, pdbFile
```

## Shared code

Three modules are shared with the rest of the lab and are **vendored** here rather
than taken as a dependency, so that pHinder installs and archives as one unit:
`compGeometry`, `determinants`, `pdbFile`.

[`isomlab`](https://github.com/isomlab/isomlab) is their upstream home and remains the
source of truth. Fix bugs there, then pull them in:

```bash
python tools/sync_vendor.py --from ../isomlab   # sync + restamp
python tools/sync_vendor.py --check             # report drift only
```

`src/pHinder/_vendor/VENDOR.json` records the upstream commit and a SHA-256 per file,
and `tests/test_vendor_integrity.py` fails if a vendored file is edited in place — the
copies cannot silently diverge.

## Provenance

Extracted from the legacy `pythonScripts/pHinder` source (active code only; the
2012–2013 `z_archives` were dropped). Algorithm bodies are unchanged — only import
statements were rewritten. PDB and mmCIF parsing are both handled by
`pHinder._vendor.pdbFile`.

## License

MIT © Daniel G. Isom (Isom Lab)
