"""Smoke tests for the isomlab-rebased pHinder package.

The load-bearing check: importing pHinder.pHinder_7_0 pulls the whole
isomlab geometry/structure stack through the rewritten imports. If the
re-base is wrong, this fails at import time.
"""
import importlib

import pytest

MODULES = [
    "pHinder",
    "pHinder.pHinder_7_0",
    "pHinder.command_line",
    "pHinder.gui.input_widgets",
    "pHinder.gui.file_open",
    "pHinder.gui.collapsible_toggle",
    "pHinder.gui.dynamic_option_widget",
    "pHinder.gui.dynamic_option_widget_amino_acid_selection",
    "pHinder.gui.terminal_output_widget",
]


@pytest.mark.parametrize("name", MODULES)
def test_imports(name):
    assert importlib.import_module(name) is not None


def test_algorithm_pulls_isomlab():
    # pHinder_7_0 re-exports the isomlab geometry/structure symbols it uses.
    mod = importlib.import_module("pHinder.pHinder_7_0")
    assert hasattr(mod, "pHinder")          # main algorithm entry class/func
    assert hasattr(mod, "convexHull4D")     # came via isomlab.geometry.convexHull4D
    assert hasattr(mod, "PDBfile")          # came via isomlab.structure.pdbFile


def test_cli_entrypoint_exists():
    from pHinder.command_line import main
    assert callable(main)
