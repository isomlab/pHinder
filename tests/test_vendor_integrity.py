"""The vendored isomlab trio must not be edited in place.

pHinder ships its own copy of compGeometry, determinants and pdbFile so it installs
without a second repository. The risk that buys is silent divergence: someone patches
the copy, upstream never learns, and the two drift apart. These tests fail loudly if
the vendored bytes stop matching the hashes recorded when they were synced.

To land a real fix: change it in isomlab, then run `python tools/sync_vendor.py`.
"""

import hashlib
import json
import pathlib
import re

import pytest

VENDOR_DIR = pathlib.Path(__file__).resolve().parent.parent / "src" / "pHinder" / "_vendor"
STAMP = json.loads((VENDOR_DIR / "VENDOR.json").read_text())


def test_stamp_records_an_upstream_commit():
    assert len(STAMP["upstream_commit"]) == 40
    assert STAMP["upstream_repo"].endswith("isomlab/isomlab")


@pytest.mark.parametrize("name", sorted(STAMP["files"]))
def test_vendored_file_matches_recorded_hash(name):
    actual = hashlib.sha256((VENDOR_DIR / name).read_bytes()).hexdigest()
    assert actual == STAMP["files"][name], (
        f"{name} was edited in place. Fix it in isomlab and run tools/sync_vendor.py."
    )


IMPORTS_ISOMLAB = re.compile(r"^\s*(from\s+isomlab[.\s]|import\s+isomlab[.\s])", re.M)


def test_no_isomlab_imports_remain():
    """pHinder must not import isomlab at runtime -- that is the whole point.

    Matches import statements only: prose mentions of isomlab in docstrings, and the
    upstream URL in VENDOR.json, are expected and fine.
    """
    src = pathlib.Path(__file__).resolve().parent.parent / "src"
    offenders = [
        str(p.relative_to(src))
        for p in src.rglob("*.py")
        if ".egg-info" not in p.parts and IMPORTS_ISOMLAB.search(p.read_text())
    ]
    assert not offenders, f"stale isomlab imports in: {offenders}"


def test_install_files_need_no_sibling_checkout():
    """environment.yml and the launchers must not require a second lab repo.

    The vendoring only pays off if *installing* is self-contained too. This is the
    gap that shipped once: src/ was clean while environment.yml still had
    `-e ../isomlab`, so `conda env create` failed for anyone without a sibling
    checkout. README may mention ../isomlab -- that is the sync-tool instruction.
    """
    root = pathlib.Path(__file__).resolve().parent.parent
    targets = [root / "environment.yml", *(root / "launchers").iterdir()]
    offenders = [
        t.name for t in targets
        if t.is_file() and "isomlab" in t.read_text()
    ]
    assert not offenders, f"install files still reference isomlab: {offenders}"
