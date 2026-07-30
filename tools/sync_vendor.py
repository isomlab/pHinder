#!/usr/bin/env python3
"""Refresh pHinder's vendored copy of the shared isomlab trio.

pHinder vendors three modules that isomlab also ships -- compGeometry, determinants
and pdbFile -- so that pHinder installs and runs without a second repository. The
upstream copy in isomlab is the source of truth: fix bugs there, then run this to
pull the fix in and restamp VENDOR.json.

    python tools/sync_vendor.py --from ../isomlab      # local checkout
    python tools/sync_vendor.py --check                # report drift, change nothing

Import paths are rewritten on the way in (isomlab.* -> pHinder._vendor.*), so the
vendored files are drop-in.
"""

import argparse
import hashlib
import json
import pathlib
import subprocess
import sys

VENDOR_DIR = pathlib.Path(__file__).resolve().parent.parent / "src" / "pHinder" / "_vendor"
STAMP = VENDOR_DIR / "VENDOR.json"

# upstream relative path -> vendored filename
FILES = {
    "src/isomlab/geometry/compGeometry.py": "compGeometry.py",
    "src/isomlab/geometry/determinants.py": "determinants.py",
    "src/isomlab/structure/pdbFile.py": "pdbFile.py",
}

REWRITES = {
    "isomlab.geometry.compGeometry": "pHinder._vendor.compGeometry",
    "isomlab.geometry.determinants": "pHinder._vendor.determinants",
    "isomlab.structure.pdbFile": "pHinder._vendor.pdbFile",
}


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def rewrite(text: str) -> str:
    for old, new in sorted(REWRITES.items(), key=lambda kv: -len(kv[0])):
        text = text.replace(old, new)
    return text


def check() -> int:
    """Report whether the vendored files still match their recorded hashes."""
    stamp = json.loads(STAMP.read_text())
    drifted = [
        name
        for name, expected in stamp["files"].items()
        if sha256((VENDOR_DIR / name).read_bytes()) != expected
    ]
    if drifted:
        print("VENDOR DRIFT -- these files were edited in place:", ", ".join(drifted))
        print("Fix the bug in isomlab instead, then re-run this script without --check.")
        return 1
    print(f"vendored files match VENDOR.json (upstream {stamp['upstream_commit'][:10]})")
    return 0


def sync(source: pathlib.Path) -> int:
    if not (source / "src" / "isomlab").is_dir():
        sys.exit(f"not an isomlab checkout: {source}")
    commit = subprocess.run(
        ["git", "-C", str(source), "rev-parse", "HEAD"],
        capture_output=True, text=True,
    ).stdout.strip()

    hashes = {}
    for rel, name in FILES.items():
        text = rewrite((source / rel).read_text())
        (VENDOR_DIR / name).write_text(text)
        hashes[name] = sha256(text.encode())
        print("synced", name)

    STAMP.write_text(json.dumps({
        "upstream_repo": "https://github.com/isomlab/isomlab",
        "upstream_commit": commit,
        "note": "Shared trio vendored from isomlab. Fix bugs upstream, then run tools/sync_vendor.py.",
        "files": hashes,
    }, indent=2) + "\n")
    print("stamped upstream", commit[:10])
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--from", dest="source", default="../isomlab",
                    help="path to an isomlab checkout (default: ../isomlab)")
    ap.add_argument("--check", action="store_true",
                    help="only report drift; do not modify anything")
    args = ap.parse_args()
    return check() if args.check else sync(pathlib.Path(args.source).resolve())


if __name__ == "__main__":
    raise SystemExit(main())
