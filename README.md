# pHinder

Ionizable-residue network and surface analysis for protein structures.

This is pHinder packaged as an installable tool. It is **self-contained**: the
algorithm, GUI, CLI, and all the computational geometry it needs live in this
repository, so installing it requires no other lab repo.

## Install & run

**Most people: just use the launcher — no typing.**

- **macOS:** in the `launchers` folder, double-click **`Launch pHinder.command`**
- **Windows:** double-click **`launchers\Launch pHinder.bat`**

The first launch sets everything up on its own (it needs [Miniforge](https://conda-forge.org/download/)
— a free, one-time install; brand new to this? see the
[install-from-scratch guide](docs/INSTALL.md)); after that it opens straight away.
Step-by-step, including how to download the code: **[getting started](docs/getting_started.md)**.

pHinder is a **public** repository, so nothing here needs a GitHub account.

<details>
<summary><b>Prefer the command line?</b> — conda or pip</summary>

```bash
# conda (its Python includes the GUI toolkit):
conda env create -f environment.yml
conda activate pHinder
phinder --help       # or: phinder-gui

# or plain pip, from a clone (needs a Python that already has tkinter):
pip install -e .

# or without cloning at all:
pip install git+https://github.com/isomlab/pHinder
```

Requires Python ≥ 3.9, `numpy`, and `openpyxl`. No other lab repository is needed.

> **Do not run `pip install pHinder`.** That name belongs to an unrelated project on
> PyPI; it will not install this tool. Install from this repository instead.

</details>

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
