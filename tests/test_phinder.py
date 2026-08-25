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


def test_algorithm_pulls_geometry_and_structure():
    # pHinder_7_0 re-exports the geometry/structure symbols it uses.
    mod = importlib.import_module("pHinder.pHinder_7_0")
    assert hasattr(mod, "pHinder")          # main algorithm entry class/func
    assert hasattr(mod, "convexHull4D")     # came via pHinder.geometry.convexHull4D
    assert hasattr(mod, "PDBfile")          # came via pHinder._vendor.pdbFile


def test_cli_entrypoint_exists():
    from pHinder.command_line import main
    assert callable(main)


# mmCIF column alignment ------------------------------------------------------
# 8HS2 arrived as a bug report: "running phinder on this cif file fails". The
# reader matched atom fields to a fixed 26-column list instead of the columns the
# file declares, so any file omitting a column -- almost every PDB-deposited
# mmCIF omits the five *_esd fields -- had every field after the gap read from
# the wrong place. Chains came out as label_asym_id and the author chain the user
# selected "did not exist".

_CIF_NO_ESD = """data_TEST
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 N N . SER A 1 19 ? 144.929 112.907 89.226 1.00 176.74 ? 19 SER B N 1
ATOM 2 C CA . SER A 1 19 ? 145.000 113.000 90.000 1.00 100.00 ? 19 SER B CA 1
ATOM 3 N N . GLY C 1 5 ? 10.000 11.000 12.000 1.00 50.00 ? 305 GLY R N 1
"""


def _parse(tmp_path, text, name="test.cif"):
    from pHinder._vendor.pdbFile import PDBfile

    path = tmp_path / name
    path.write_text(text)
    return PDBfile(str(tmp_path) + "/", name)


def test_cif_columns_come_from_the_file_not_a_fixed_list(tmp_path):
    """Fields after an omitted column must still land in the right place."""
    protein = _parse(tmp_path, _CIF_NO_ESD)
    first = protein.atoms[1]
    assert first.residue_name == "SER"
    assert first.atom_name == "N"
    assert (first.x, first.y, first.z) == (144.929, 112.907, 89.226)


def test_cif_chains_are_the_author_ids(tmp_path):
    """B and R are what a viewer shows and what a person selects; A and C are
    the label ids and must not reach the user."""
    protein = _parse(tmp_path, _CIF_NO_ESD)
    assert sorted(protein.chains) == ["B", "R"]


def test_cif_residue_numbers_are_the_author_numbering(tmp_path):
    """Chain and residue number have to come from the same system."""
    protein = _parse(tmp_path, _CIF_NO_ESD)
    numbers = {a.chain_identifier: a.residue_sequence_number
               for a in protein.atoms.values()}
    assert numbers["R"] == 305, "label_seq_id is 5 here; the author number is 305"


def test_gzipped_cif_is_not_handed_to_the_pdb_parser(tmp_path):
    """A .cif.gz ends with .gz, so extension matching has to strip that first."""
    import gzip

    from pHinder._vendor.pdbFile import PDBfile

    path = tmp_path / "test.cif.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(_CIF_NO_ESD)
    protein = PDBfile(str(tmp_path) + "/", "test.cif.gz", zip_status=1)
    assert protein.format == "cif"
    assert sorted(protein.chains) == ["B", "R"]
