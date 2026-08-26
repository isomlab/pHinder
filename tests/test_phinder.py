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


# Pseudo-atoms in the triangulation -------------------------------------------
# convexHull4D seeds the hull with a template tetrahedron of PseudoAtom vertices.
# They have residue = None by construction. A large query set buries them, but a
# narrow residue selection leaves so few real vertices that they survive pruning
# -- which is why the CLI's own default, ionizableSetNoCys, crashed, along with
# acidicSet and basicSet, while allSet and the broader sets did not.


def test_min_sidechain_distance_tolerates_a_pseudo_atom():
    from pHinder._vendor.compGeometry import findMinSidechainDistance
    from pHinder._vendor.pdbFile import PseudoAtom

    psa = PseudoAtom()
    psa.x, psa.y, psa.z = 1.0, 0.0, 0.0
    psa.reinitialize()
    assert psa.residue is None, "a pseudo-atom stands for no residue"

    other = PseudoAtom()
    other.x, other.y, other.z = 0.0, 1.0, 0.0
    other.reinitialize()

    d = findMinSidechainDistance(psa, other)
    assert d == pytest.approx(2 ** 0.5), "it stands for itself, so plain distance"


def test_pruning_skips_pseudo_atom_nodes():
    """A node with no residue is scaffolding: drop its edges and move on."""
    from pHinder.geometry.goFo import pruneTriangulationGoFo

    class Data:
        # goFo identifies pseudo-atoms by residue_name; the crash came from
        # residue itself being None, so the guard keys on that.
        def __init__(self, residue, residue_name):
            self.residue = residue
            self.residue_name = residue_name

    class Residue:
        num, chn, name = 1, "A", "ASP"

    class Node:
        def __init__(self, residue, residue_name, s2s):
            self.s1 = type("V", (), {"data": Data(residue, residue_name)})()
            self.s2s = list(s2s)

    triangulation = {0: Node(None, "PSA", [(0, 1)]),
                     1: Node(Residue(), "ASP", [])}
    triangulation, networks = pruneTriangulationGoFo(triangulation, 5.0)

    assert triangulation[0].s2s == [], "the pseudo node keeps no edges"
    assert isinstance(networks, dict)


# Structure files are the user's, and residues have numbers ------------------
# Both bugs sat in the same place. A PDB with insertion codes triggered a
# renumbering pass that (a) rewrote the file on disk, four times over, and
# (b) keyed on atom_name == " N" -- PDB atom names are four characters, so the
# backbone nitrogen is " N  " and the test never fired. Every residue was
# therefore numbered 0, and residues collapsed by (0, name, chain): 281
# residues in 2PTC became 38.

_PDB_WITH_INSERTION_CODE = """\
ATOM      1  N   ILE E  16      18.871  65.715  12.731  1.00 21.86           N
ATOM      2  CA  ILE E  16      19.481  64.399  12.507  1.00 20.71           C
ATOM      3  N   VAL E  17      20.109  64.328  11.320  1.00 19.64           N
ATOM      4  CA  VAL E  17      20.759  63.086  10.939  1.00 18.86           C
ATOM      5  N   GLY E 184      21.000  62.000  10.000  1.00 18.00           N
ATOM      6  CA  GLY E 184      21.500  61.500   9.500  1.00 18.00           C
ATOM      7  N   GLY E 184A     22.000  61.000   9.000  1.00 18.00           N
ATOM      8  CA  GLY E 184A     22.500  60.500   8.500  1.00 18.00           C
END
"""


def test_a_structure_file_is_never_rewritten(tmp_path):
    import hashlib

    from pHinder._vendor.pdbFile import PDBfile

    path = tmp_path / "insertion.pdb"
    path.write_text(_PDB_WITH_INSERTION_CODE)
    before = hashlib.sha256(path.read_bytes()).hexdigest()

    PDBfile(str(tmp_path) + "/", "insertion.pdb")

    assert hashlib.sha256(path.read_bytes()).hexdigest() == before, \
        "pHinder must not rewrite the structure it was handed"


def test_renumbering_actually_numbers(tmp_path):
    """Four residues in, four residues out -- not one bucket numbered zero."""
    from pHinder._vendor.pdbFile import PDBfile

    path = tmp_path / "insertion.pdb"
    path.write_text(_PDB_WITH_INSERTION_CODE)
    protein = PDBfile(str(tmp_path) + "/", "insertion.pdb")

    numbers = sorted({r.num for r in protein.residues.values()})
    assert len(protein.residues) == 4, "ILE, VAL and both GLYs are distinct residues"
    assert 0 not in numbers, "residue 0 is the signature of the old broken test"
    assert len(numbers) == 4, f"every residue needs its own number, got {numbers}"


def test_a_file_without_insertion_codes_keeps_its_numbering(tmp_path):
    """The renumbering pass only fires when it has to."""
    from pHinder._vendor.pdbFile import PDBfile

    clean = "\n".join(l for l in _PDB_WITH_INSERTION_CODE.splitlines()
                      if "184A" not in l) + "\n"
    path = tmp_path / "clean.pdb"
    path.write_text(clean)
    protein = PDBfile(str(tmp_path) + "/", "clean.pdb")

    assert sorted({r.num for r in protein.residues.values()}) == [16, 17, 184]


def test_the_original_residue_label_is_carried(tmp_path):
    """Renumbering replaces the file's numbering; the label survives alongside
    it, or a result cannot be taken back to the structure it came from."""
    from pHinder._vendor.pdbFile import PDBfile

    path = tmp_path / "insertion.pdb"
    path.write_text(_PDB_WITH_INSERTION_CODE)
    protein = PDBfile(str(tmp_path) + "/", "insertion.pdb")

    labels = {r.num: r.num_original for r in protein.residues.values()}
    assert sorted(labels) == [1, 2, 3, 4], "renumbered 1..4"
    assert sorted(labels.values()) == ["16", "17", "184", "184A"], \
        f"the file's own labels, insertion code included: {labels}"


def test_a_file_that_is_not_renumbered_reports_its_own_numbers(tmp_path):
    from pHinder._vendor.pdbFile import PDBfile

    clean = "\n".join(l for l in _PDB_WITH_INSERTION_CODE.splitlines()
                      if "184A" not in l) + "\n"
    path = tmp_path / "clean.pdb"
    path.write_text(clean)
    protein = PDBfile(str(tmp_path) + "/", "clean.pdb")

    for r in protein.residues.values():
        assert r.num_original == str(r.num), "nothing was renumbered"


def test_the_reported_depth_is_the_one_the_class_was_decided_on():
    """The spreadsheet used to recover depth from a "%5.1f" rendering of it, so
    -3.04 and -2.96 both arrived as -3.0 and the number could not explain why one
    residue is core and the other margin. The writer reads the instance now."""
    import inspect

    from pHinder.pHinder_7_0 import pHinder

    src = inspect.getsource(pHinder.writeSidechainClassificationResults)
    assert "ci.depth" in src, "depth comes from the classification instance"
    assert "float(fullSplit" not in src, "not parsed back out of the printed string"


def test_the_class_letter_is_kept_on_the_instance():
    from pHinder.pHinder_7_0 import sidechainClassification

    import inspect
    src = inspect.getsource(sidechainClassification.setClassificationString)
    assert "self.locationKey = locationKey" in src


# Deterministic jostling -----------------------------------------------------
# Jostling is required: without it the orientation predicates are ambiguous at
# exact zero and facet rotation can cycle instead of terminating. None of that
# changes. What changed is only where the random numbers come from, so that the
# same structure gets the same perturbation. The switch below restores the old
# behaviour exactly, and these tests pin both directions.


def test_the_same_vertex_gets_the_same_nudge_every_time():
    from pHinder._vendor import compGeometry as cg

    def nudge():
        v = cg.Vertex((1.234, 5.678, 9.012))
        v.jostle()
        return (v.x, v.y, v.z)

    assert nudge() == nudge(), "same input, same perturbation"


def test_a_retry_still_moves_the_vertex_somewhere_new():
    """Escaping degeneracy depends on repeated jostles differing."""
    from pHinder._vendor import compGeometry as cg

    v = cg.Vertex((1.234, 5.678, 9.012))
    seen = []
    for _ in range(4):
        v.jostle()
        seen.append((v.x, v.y, v.z))
    assert len(set(seen)) == 4, "each retry has to land somewhere new"


def test_different_vertices_get_different_nudges():
    from pHinder._vendor import compGeometry as cg

    a = cg.Vertex((1.0, 2.0, 3.0)); a.jostle()
    b = cg.Vertex((1.5, 2.5, 3.5)); b.jostle()
    assert (a.x - 1.0, a.y - 2.0) != (b.x - 1.5, b.y - 2.5)


def test_the_revert_switch_is_wired_to_both_jostles():
    """PHINDER_JOSTLE=random must put every draw back on the global RNG."""
    import inspect

    from pHinder._vendor import compGeometry as cg

    assert cg.JOSTLE_DETERMINISTIC is True, "deterministic is the default"
    for method in (cg.Vertex.jostle, cg.Vertex4D.jostle):
        src = inspect.getsource(method)
        assert "JOSTLE_DETERMINISTIC" in src, f"{method.__qualname__} not switched"
    assert 'os.environ.get("PHINDER_JOSTLE"' in inspect.getsource(cg)


def test_the_perturbation_is_still_bounded_by_eps():
    """The revert has to be a change of source, not of magnitude."""
    from pHinder._vendor import compGeometry as cg

    eps = 1e-3
    for _ in range(20):
        v = cg.Vertex((0.0, 0.0, 0.0))
        v.jostle(eps=eps)
        assert abs(v.x) <= eps and abs(v.y) <= eps and abs(v.z) <= eps


# ---------------------------------------------------------------------------
# General position: the xy-only gp2D, and the unbounded jostle loops it hung.
# ---------------------------------------------------------------------------

def _v3(x, y, z):
    from pHinder._vendor.compGeometry import Vertex
    return Vertex((x, y, z))


def test_gp2D_is_a_3d_test_not_an_xy_projection():
    # gp2D compared only x and y. Both terms of that determinant carry a factor
    # of (x2-x1) or (y2-y1), so two vertices sharing those coordinates made it
    # report collinearity for EVERY third vertex -- and the callers respond to a
    # zero by jostling the third vertex, which cannot change either factor.
    from pHinder._vendor.compGeometry import gp2D
    assert gp2D(_v3(0.0, 0.0, 2.0), _v3(0.0, 0.0, 2.2), _v3(34.9, 13.2, 22.6)) == 1
    # Collinear in the xy projection, plainly not collinear in space.
    assert gp2D(_v3(0.0, 0.0, 0.0), _v3(1.0, 0.0, 1.0), _v3(2.0, 0.0, 0.0)) == 1


def test_gp2D_still_catches_real_collinearity():
    from pHinder._vendor.compGeometry import gp2D
    assert gp2D(_v3(0.0, 0.0, 0.0), _v3(0.0, 0.0, 1.0), _v3(0.0, 0.0, 2.0)) == 0
    assert gp2D(_v3(0.0, 0.0, 0.0), _v3(1.0, 2.0, 3.0), _v3(2.0, 4.0, 6.0)) == 0
    assert gp2D(_v3(0.0, 0.0, 0.0), _v3(1.0, 0.0, 0.0), _v3(0.0, 1.0, 0.0)) == 1


def test_jostling_the_third_vertex_can_now_reach_general_position():
    # The property every caller depends on: when the test fails, moving the
    # vertex the loop moves has to be able to fix it.
    from pHinder._vendor.compGeometry import gp2D
    v1, v2, v3 = _v3(0.0, 0.0, 2.0), _v3(0.0, 0.0, 2.2), _v3(0.0, 0.0, 2.4)
    attempts = 0
    while not gp2D(v1, v2, v3):
        v3.jostle()
        attempts += 1
        assert attempts < 100, "jostling never reached general position"


def test_jostle_budget_raises_instead_of_spinning_forever():
    from pHinder.geometry.general_position import JostleBudget, GeneralPositionError
    budget = JostleBudget("a seed simplex", max_jostles=5)
    v = _v3(0.0, 0.0, 0.0)
    for _ in range(5):
        budget.jostle(v)
    with pytest.raises(GeneralPositionError) as excinfo:
        budget.jostle(v)
    # The message has to name the loop and locate the vertex, or a bug report
    # from a lab member says only "it stopped".
    assert "a seed simplex" in str(excinfo.value)
    assert "5 jostles" in str(excinfo.value)


def test_a_degenerate_lattice_triangulates_instead_of_hanging():
    # 33 points on a 0.5 A cubic lattice, offset along z so that the two
    # lowest-u vertices -- the two convexHull4D sorts to the front -- share an x
    # and a y. That is the configuration the xy-only gp2D could never resolve:
    # before the fix this call did not return, jostling one vertex millions of
    # times while it wandered hundreds of angstroms away.
    import itertools
    from pHinder._vendor.compGeometry import Vertex4D
    from pHinder._vendor.pdbFile import PseudoAtom
    from pHinder.geometry.convexHull4D import convexHull4D

    n, h = 2, 0.5
    pts = [(i * h, j * h, k * h + 3.0)
           for i, j, k in itertools.product(range(-n, n + 1), repeat=3)
           if i * i + j * j + k * k <= n * n]
    pts.sort(key=lambda p: (p[0] ** 2 + p[1] ** 2 + p[2] ** 2, p[0], p[1], p[2]))
    assert (pts[0][0], pts[0][1]) == (pts[1][0], pts[1][1]), "trigger not reproduced"

    vertices = []
    for idx, (x, y, z) in enumerate(pts):
        psa = PseudoAtom()
        psa.x, psa.y, psa.z = x, y, z
        vertices.append(Vertex4D((x, y, z, x * x + y * y + z * z),
                                 data=psa, unique_id=idx))

    hull = convexHull4D(vertices)
    assert hull.hull4D, "hull built but empty"


class _JostleAtom:
    """Stand-in for the atom a vertex carries; Vertex4D.reinitialize writes to it."""
    x = y = z = 0.0
    v = None


def test_coincident_vertices_come_apart_when_jostled():
    # Seeding on position and jostle-count alone gave two vertices at the same
    # point the same nudge, so they moved together and stayed coincident however
    # many times they were jostled -- the one degeneracy the perturbation exists
    # to break, and the one case the old global RNG handled for free. Their ids
    # differ, so they have to separate.
    from pHinder._vendor.compGeometry import Vertex, Vertex4D

    a = Vertex((0.4, -0.2, 0.6), unique_id=7)
    b = Vertex((0.4, -0.2, 0.6), unique_id=8)
    for _ in range(5):
        a.jostle()
        b.jostle()
    assert (a.x, a.y, a.z) != (b.x, b.y, b.z)

    u = 0.4 ** 2 + 0.2 ** 2 + 0.6 ** 2
    c = Vertex4D((0.4, -0.2, 0.6, u), data=_JostleAtom(), unique_id=7)
    e = Vertex4D((0.4, -0.2, 0.6, u), data=_JostleAtom(), unique_id=8)
    for _ in range(5):
        c.jostle()
        e.jostle()
    assert (c.x, c.y, c.z) != (e.x, e.y, e.z)


def test_identity_does_not_cost_determinism():
    from pHinder._vendor.compGeometry import Vertex

    def nudge(uid):
        v = Vertex((1.234, 5.678, 9.012), unique_id=uid)
        v.jostle()
        return (v.x, v.y, v.z)

    # Rebuilt in a later run, the same vertex draws the same nudge.
    assert nudge(7) == nudge(7)
    assert nudge(7) != nudge(8)
    # A vertex with no id is seeded exactly as before identity was folded in.
    assert nudge("no id") == nudge("no id")
