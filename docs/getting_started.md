# Getting started with pHinder

pHinder analyses ionizable-residue networks and molecular surfaces in protein structures.

Everything the lab tools have in common — downloading, launching, updating, and what
to do when something goes wrong — is on one shared page: **[Getting started with a lab
tool](https://dangerisom.github.io/Isom-Lab/getting-started/)**. This guide covers
only what is specific to pHinder.

---

## Before you start

Your computer needs **conda** (Miniforge, Miniconda, or Anaconda).

> **First time on this computer?** Do the one-time **[Setting up your
> computer](https://dangerisom.github.io/Isom-Lab/setup/)** first, then come back here.
> pHinder's own install notes are in **[INSTALL.md](INSTALL.md)**.

---

## Get it and launch it

**1. Download it** from **[github.com/isomlab/pHinder](https://github.com/isomlab/pHinder)**.
   It is **public**, so no account or password is needed: **Download ZIP**,
   **GitHub Desktop**, or `git clone`. Step by step:
   **[Get the code](https://dangerisom.github.io/Isom-Lab/getting-started/#public-tools)**.

**2. Open the `launchers` folder inside it and double-click:**

- **Mac:** `Launch pHinder.command`
- **Windows:** `Launch pHinder.bat`

The **first** launch takes a few minutes while it builds a private, isolated conda
environment named `phinder` containing Python and everything the app needs. Every
launch after that opens straight away. You don't need to type anything.

If macOS blocks the file, or double-clicking does nothing, see **[Launch
it](https://dangerisom.github.io/Isom-Lab/getting-started/#launch-it)**.

---

## Use it

The launcher opens the graphical interface. From there you load a structure, pick a
residue set, and run the topology / surface / network calculations.

If you'd rather work in a Terminal — for batch jobs or scripting — activate the
environment the launcher built and use the command line:

```bash
conda activate phinder
phinder --help                       # see all options
phinder structure.pdb --chains A --sidechain-classification
```

`phinder-gui` opens the same window the launcher does.

---

## Updating later

Refresh the folder the way you got it — see **[Updating
later](https://dangerisom.github.io/Isom-Lab/getting-started/#updating-later)**. If a
release changes what the tool depends on, delete its environment and let the launcher
rebuild it on the next double-click:

```bash
conda env remove -n phinder
```

---

## If something goes wrong

The usual problems are on the shared page: **[If something goes
wrong](https://dangerisom.github.io/Isom-Lab/getting-started/#if-something-goes-wrong)**.

Stuck? Send Dan (<a href="&#109;&#97;&#105;&#108;&#116;&#111;&#58;&#100;&#105;&#115;&#111;&#109;&#64;&#109;&#105;&#97;&#109;&#105;&#46;&#101;&#100;&#117;">&#100;&#105;&#115;&#111;&#109;<span>&#64;</span>&#109;&#105;&#97;&#109;&#105;<span>&#46;</span>&#101;&#100;&#117;</a>) the exact command you ran and the message you got.

---

## Alternative: plain pip

The general recipe is on the shared page: **[Alternative: plain
pip](https://dangerisom.github.io/Isom-Lab/getting-started/#alternative-plain-pip)**.

**Do not run `pip install pHinder`.** That name belongs to an unrelated project on
PyPI and will not install this tool. Install from the GitHub address instead:

```bash
pip install git+https://github.com/isomlab/pHinder
```
