# Getting started with pHinder

This guide gets pHinder running on your computer. You do **not** need to know how to
code — follow the steps and copy-paste where asked. It takes about 10 minutes, once.
After that, launching it is a double-click.

> On Windows everything below works the same — use the **Miniforge Prompt** app
> wherever this says "Terminal", and see the Windows note in each step.

---

## Before you start

This guide assumes your computer already has **conda** (Miniforge, Miniconda, or
Anaconda) installed.

> **First time? Never installed conda?** Do the one-time
> **[install-from-scratch guide → INSTALL.md](INSTALL.md)** first, then come back here.

---

## Step 1 — Get the code

pHinder is a **public** repository, so no account or password is needed. Pick
whichever way you prefer.

**Option A — Download ZIP (fastest, nothing to install).**
1. Open **[github.com/isomlab/pHinder](https://github.com/isomlab/pHinder)**.
2. Click the green **`Code ▾`** button → **Download ZIP**.
3. Double-click the downloaded file to unzip it. You now have a folder called
   **`pHinder-main`** — move it somewhere easy, like your **Documents**.

*(To update later, download a fresh ZIP and replace the folder.)*

**Option B — GitHub Desktop (best if you'll update often).**
1. Open GitHub Desktop → **File ▸ Clone repository… ▸ URL**.
2. Paste `https://github.com/isomlab/pHinder` → pick a folder (e.g. Documents) →
   **Clone**. This makes a folder called **`pHinder`**.
3. To update later, open GitHub Desktop and click **Fetch/Pull origin**.

**Option C — `git clone` in Terminal.** Because the repo is public this just works,
with no password:

```bash
cd ~/Documents
git clone https://github.com/isomlab/pHinder.git
```

*(To update later: `cd ~/Documents/pHinder && git pull`.)*

Either way you now have a pHinder folder on your computer.

---

## Step 2 — Launch it

**Open the `launchers` folder inside your pHinder folder and double-click:**

- **Mac:** `Launch pHinder.command`
- **Windows:** `Launch pHinder.bat`

That's it. The **first** launch takes a few minutes: it builds a private, isolated
conda environment (named `pHinder`) containing Python and everything the app needs,
then opens the window. **Every launch after that opens straight away.**

You don't need to type anything, and it won't touch any other Python on your
computer.

> **Mac, first time only:** if macOS says *"cannot be opened because it is from an
> unidentified developer"*, right-click (or Control-click) the `.command` file and
> choose **Open**, then **Open** again in the dialog. You only do this once.

> **Mac, if double-clicking does nothing:** the file may have lost its executable
> flag in transit (common after unzipping). In Terminal, run
> `chmod +x ~/Documents/pHinder/launchers/"Launch pHinder.command"` once, then
> double-click again.

---

## Step 3 — Use it

The launcher opens the graphical interface. From there you load a structure, pick a
residue set, and run the topology / surface / network calculations.

If you'd rather work in a Terminal — for batch jobs or scripting — activate the
environment the launcher built and use the command line:

```bash
conda activate pHinder
phinder --help                       # see all options
phinder structure.pdb --chains A --sidechain-classification
```

`phinder-gui` opens the same window the launcher does.

---

## Updating later

- **GitHub Desktop (Option B):** open it and click **Fetch / Pull origin**.
- **`git clone` (Option C):** `cd ~/Documents/pHinder && git pull`.
- **Downloaded the ZIP (Option A):** download a fresh ZIP and replace the old
  folder's contents (keep the same folder name and location).

The environment installs the code in "editable" mode, so an update takes effect the
next time you launch — no reinstall. If a release changes the dependencies, delete
the environment and let the launcher rebuild it:

```bash
conda env remove -n pHinder
```

---

## If something goes wrong

- **"conda: command not found"** — close and reopen Terminal after installing conda
  (the installer needs a fresh window). On Mac, if it still isn't found, run
  `source ~/miniforge3/bin/activate` once.
- **The launcher says it can't find conda** — same cause. Install conda from the
  [install-from-scratch guide](INSTALL.md), then double-click the launcher again.
- **"phinder: command not found"** — you probably forgot `conda activate pHinder`
  first. Run it, then try again.
- **The window doesn't appear** — use the **conda** install above; its Python
  includes the Tk graphics toolkit the GUI needs. A plain system-Python `pip install`
  can be missing it.
- **Setup failed partway through** — remove the half-built environment with
  `conda env remove -n pHinder` and double-click the launcher again.

Stuck? Send Dan the exact command you ran and the message you got.

---

## Alternative: plain `pip` (if you don't use conda)

Conda is recommended because it guarantees the GUI toolkit is present. If you'd
rather use `pip`, from inside the pHinder folder:

```bash
python3 -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -e .
phinder --help
```

Or without cloning at all:

```bash
pip install git+https://github.com/isomlab/pHinder
```

This only works if your Python includes **tkinter**: macOS's built-in `python3`
does; Homebrew Python needs `brew install python-tk`; conda always does.

> **Do not run `pip install pHinder`.** That name belongs to an unrelated project on
> PyPI and will not install this tool.
