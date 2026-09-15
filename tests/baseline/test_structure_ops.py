"""Structure processing: RMSD pruning, bonds, geometry measurement/manipulation, xyz/mol IO."""

import math

import pytest
from builders import (
    METHANE,
    RIGHT_ANGLE,
    rotated_permuted_methane,
    stretched_methane,
    system_with_atoms,
    write_xyz,
)

from accel import Box
from accel.base.atoms import BondType


def _xyz(c):
    return [a.xyz for a in c.atoms]


def _flat(rows):
    return [v for row in rows for v in row]


def test_rmsd_limit_removes_symmetry_equivalent_duplicate():
    box = Box(
        [
            system_with_atoms("c1", METHANE, energy=0.0, label="m"),
            system_with_atoms("c2", rotated_permuted_methane(), energy=1.0, label="m"),
            system_with_atoms("c3", stretched_methane(), energy=2.0, label="m"),
        ]
    )
    box.rmsd_limit()
    assert [c.state for c in box.contents] == [True, False, True]
    assert box.contents[1].history == "rmsd_limit"
    assert all(c.data["has_symm"] is True for c in box.contents)


def test_rmsd_limit_keeps_the_lower_energy_duplicate():
    box = Box(
        [
            system_with_atoms("c1", METHANE, energy=1.0, label="m"),
            system_with_atoms("c2", rotated_permuted_methane(), energy=0.0, label="m"),
        ]
    )
    box.rmsd_limit()
    assert [c.state for c in box.contents] == [False, True]


def test_rmsd_limit_never_compares_across_labels():
    box = Box(
        [
            system_with_atoms("c1", METHANE, energy=0.0, label="a"),
            system_with_atoms("c2", METHANE, energy=1.0, label="b"),
        ]
    )
    box.rmsd_limit()
    assert [c.state for c in box.contents] == [True, True]


def test_calc_bonds_methane():
    box = Box([system_with_atoms("m", METHANE)]).calc_bonds()
    bonds = box.contents[0].atoms.bonds.to_dict()
    assert bonds == {
        (1, 2): BondType.single,
        (1, 3): BondType.single,
        (1, 4): BondType.single,
        (1, 5): BondType.single,
    }


def test_geometry_measurements_store_data_with_default_keys():
    box = Box([system_with_atoms("t", RIGHT_ANGLE)])
    box.calc_length(1, 2).calc_angle(2, 1, 3).calc_length(2, 3, key="hh")
    data = box.contents[0].data
    assert data["length_O1-H2"] == pytest.approx(1.0)
    assert data["angle_H2-O1-H3"] == pytest.approx(90.0)
    assert data["hh"] == pytest.approx(math.sqrt(2.0))
    box.calc_angle(2, 1, 3, key="rad", radian=True)
    assert data["rad"] == pytest.approx(math.pi / 2)


@pytest.mark.parametrize(
    "d,expected",
    [((1.0, 1.0, 0.0), 0.0), ((1.0, -1.0, 0.0), 180.0), ((1.0, 0.0, 1.0), 90.0)],
)
def test_calc_dihedral_iupac_sign(d, expected):
    atoms = [("H", 0.0, 1.0, 0.0), ("C", 0.0, 0.0, 0.0), ("C", 1.0, 0.0, 0.0), ("H", *d)]
    box = Box([system_with_atoms("d", atoms)]).calc_dihedral(1, 2, 3, 4, key="phi")
    assert box.contents[0].data["phi"] == pytest.approx(expected)


def test_modify_length_moves_both_atoms_symmetrically():
    box = Box([system_with_atoms("t", RIGHT_ANGLE)]).modify_length(1, 2, 1.5)
    a, b = box.contents[0].atoms[0], box.contents[0].atoms[1]
    assert box.contents[0].atoms.get_length(1, 2) == pytest.approx(1.5)
    assert (a.x + b.x) / 2 == pytest.approx(0.5)


def test_modify_length_with_fixed_atom_and_followers():
    box = Box([system_with_atoms("t", RIGHT_ANGLE)])
    box.modify_length(1, 2, 1.5, fix_a=True, numbers_along_with_b=[3])
    atoms = box.contents[0].atoms
    assert atoms[0].xyz == [0.0, 0.0, 0.0]
    assert atoms[1].xyz == pytest.approx([1.5, 0.0, 0.0])
    assert atoms[2].xyz == pytest.approx([0.5, 1.0, 0.0])


def test_convert_to_mirror_without_centering_negates_coordinates():
    box = Box([system_with_atoms("t", stretched_methane())]).convert_to_mirror(centering=False)
    assert _xyz(box.contents[0]) == [[-x, -y, -z] for _, x, y, z in stretched_methane()]


def test_convert_to_mirror_with_centering_rounds_centre_to_input_precision():
    # the centre is rounded to the number of decimals in the input coordinates (here 1),
    # so the result is centred only approximately: mean(x) = 1/3 -> shift of 0.3
    box = Box([system_with_atoms("t", RIGHT_ANGLE)]).convert_to_mirror()
    assert _flat(_xyz(box.contents[0])) == pytest.approx([0.3, 0.3, 0.0, -0.7, 0.3, 0.0, 0.3, -0.7, 0.0])


def test_xyz_round_trip_without_centering(tmp_path):
    src = write_xyz(tmp_path / "m_1.xyz", stretched_methane())
    out = tmp_path / "out"
    box = Box(src).read_atoms().write_xyz(directory=out, centering=False)
    assert box.contents[0].path == (out / "m_1.xyz").absolute()
    again = Box(out / "*.xyz").read_atoms()
    assert _xyz(again.contents[0]) == _xyz(box.contents[0])
    assert [a.symbol for a in again.contents[0].atoms] == ["C", "H", "H", "H", "H"]


def test_write_xyz_centering_rounds_centre_to_input_precision(tmp_path):
    # same rounding rule as convert_to_mirror: centre = round(mean, input decimals)
    src = write_xyz(tmp_path / "t.xyz", RIGHT_ANGLE)
    Box(src).read_atoms().write_xyz(directory=tmp_path / "out")
    again = Box(tmp_path / "out" / "t.xyz").read_atoms()
    assert _flat(_xyz(again.contents[0])) == pytest.approx([-0.3, -0.3, 0.0, 0.7, -0.3, 0.0, -0.3, 0.7, 0.0])


def test_read_xyz_with_wrong_atom_count_deactivates(tmp_path):
    p = tmp_path / "bad.xyz"
    p.write_text("3\ncomment\nC 0 0 0\nH 1 0 0\n")
    box = Box(p).read_atoms()
    assert box.contents[0].state is False


def test_mol_round_trip(tmp_path):
    src = write_xyz(tmp_path / "m_1.xyz", METHANE)
    Box(src).read_atoms().write_mol(directory=tmp_path / "out")
    again = Box(tmp_path / "out" / "*.sdf").read_atoms()
    c = again.contents[0]
    assert c.filetype == "format/mol"
    assert [a.symbol for a in c.atoms] == ["C", "H", "H", "H", "H"]
    assert _flat(_xyz(c)) == pytest.approx(_flat([[x, y, z] for _, x, y, z in METHANE]), abs=1e-4)
    assert c.atoms.bonds.to_dict() == {(1, n): BondType.single for n in (2, 3, 4, 5)}


def test_mol_formal_charges_round_trip(tmp_path):
    src = write_xyz(tmp_path / "ion.xyz", [("N", 0.0, 0.0, 0.0), ("H", 1.01, 0.0, 0.0)])
    box = Box(src).read_atoms()
    box.contents[0].atoms[0].charge = 1
    box.write_mol(directory=tmp_path / "out")
    again = Box(tmp_path / "out" / "ion.sdf").read_atoms()
    assert [a.charge for a in again.contents[0].atoms] == [1, None]
    assert again.contents[0].charge == 1
