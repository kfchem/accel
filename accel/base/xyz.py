import math
from statistics import mean
from typing import Sequence, Tuple

import numpy as np
from accel.base.atoms import Atom
from accel.base.systems import System
from accel.base.tools import float_to_str


def edit_bond_length(
    system: System,
    atom_a: int,
    atom_b: int,
    target_length: float,
    fixed_atom: Tuple[bool, bool] = (False, False),
    move_along_with_a: Sequence[int] = (),
    move_along_with_b: Sequence[int] = (),
):
    vect = [0.0, 0.0, 0.0]
    for i in range(3):
        vect[i] = system.atoms.get(atom_b).xyz[i] - system.atoms.get(atom_a).xyz[i]
    distance = math.sqrt(sum(x**2 for x in vect))

    def move_atoms(atoms_list, vect_factor):
        for atom_no in atoms_list:
            system.atoms.get(atom_no).xyz = [
                _val + (vect_factor * vect[i] * (distance - target_length) / distance)
                for i, _val in enumerate(system.atoms.get(atom_no).xyz)
            ]

    if fixed_atom == (False, False):
        move_atoms([atom_a] + list(move_along_with_a), 0.5)
        move_atoms([atom_b] + list(move_along_with_b), -0.5)
    elif fixed_atom == (False, True):
        move_atoms([atom_a] + list(move_along_with_a), 1.0)
    elif fixed_atom == (True, False):
        move_atoms([atom_b] + list(move_along_with_b), -1.0)
    else:
        raise ValueError


def set_chirality(system: System, center_index: int, sub_index: list[int]):
    if len(sub_index) != 4:
        raise ValueError
    else:
        sorted_index = sorted(sub_index)

    _sub_xyzs = np.array([system.atoms.get(i).xyz for i in sorted_index[1:]]) - np.array(
        [system.atoms.get(sorted_index[0]).xyz for _ in range(3)]
    )
    _ret = np.linalg.det(_sub_xyzs)
    if _ret > 0:
        _ret = 1
    elif _ret < 0:
        _ret = -1
    else:
        _ret = 0
    system.data[f"chiral_{center_index}_to_{sub_index}"] = _ret


def calc_length(system: System, atom_index_a: int, atom_index_b: int, key: str = ""):
    a_ = system.atoms.get(atom_index_a).xyz
    b_ = system.atoms.get(atom_index_b).xyz
    d_ = [float(a_[i]) - float(b_[i]) for i in range(3)]
    distance = math.sqrt(sum(x**2 for x in d_))
    if key == "" or not isinstance(key, str):
        key = "distance_{}{}-{}{}".format(
            system.atoms.get(atom_index_a).symbol,
            str(atom_index_a),
            system.atoms.get(atom_index_b).symbol,
            str(atom_index_b),
        )
    system.data[key] = distance


def get_dihedral(atom_a: Atom, atom_b: Atom, atom_c: Atom, atom_d: Atom) -> float:
    va = np.array(atom_a.xyz)
    vb = np.array(atom_b.xyz)
    vc = np.array(atom_c.xyz)
    vd = np.array(atom_d.xyz)
    vab = va - vb
    vcb = vc - vb
    vdc = vd - vc
    pvac = np.cross(vab, vcb)
    pvbd = np.cross(vdc, vcb)
    dac = np.linalg.norm(pvac)
    dbd = np.linalg.norm(pvbd)
    angle = np.arccos(np.sum(pvac * pvbd) / (dac * dbd))
    if np.sum(pvac * np.cross(pvbd, vcb)) < 0:
        angle = -angle
    angle = float(np.rad2deg(angle))
    return angle


def calc_dihedral(
    system: System,
    atom_index_a: int,
    atom_index_b: int,
    atom_index_c: int,
    atom_index_d: int,
    key: str = "",
):
    if key == "" or not isinstance(key, str):
        key = "dihedral_{}{}-{}{}-{}{}-{}{}".format(
            system.atoms.get(atom_index_a).symbol,
            str(atom_index_a),
            system.atoms.get(atom_index_b).symbol,
            str(atom_index_b),
            system.atoms.get(atom_index_c).symbol,
            str(atom_index_c),
            system.atoms.get(atom_index_d).symbol,
            str(atom_index_d),
        )
    system.data[key] = get_dihedral(
        system.atoms.get(atom_index_a),
        system.atoms.get(atom_index_b),
        system.atoms.get(atom_index_c),
        system.atoms.get(atom_index_d),
    )


def get_angle(atom_a: Atom, atom_b: Atom, atom_c: Atom) -> float:
    va = np.array(atom_a.xyz)
    vb = np.array(atom_b.xyz)
    vc = np.array(atom_c.xyz)
    vba = vb - va
    vbc = vb - vc
    dba = np.linalg.norm(vba)
    dbc = np.linalg.norm(vbc)
    angle = np.arccos(np.sum(vba * vbc) / (dba * dbc))
    angle = float(np.rad2deg(angle))
    return angle


def calc_angle(
    system: System,
    atom_index_a: int,
    atom_index_b: int,
    atom_index_c: int,
    key: str = "",
):
    if key == "" or not isinstance(key, str):
        key = "angle_{}{}-{}{}-{}{}".format(
            system.atoms.get(atom_index_a).symbol,
            str(atom_index_a),
            system.atoms.get(atom_index_b).symbol,
            str(atom_index_b),
            system.atoms.get(atom_index_c).symbol,
            str(atom_index_c),
        )
    system.data[key] = get_angle(
        system.atoms.get(atom_index_a),
        system.atoms.get(atom_index_b),
        system.atoms.get(atom_index_c),
    )


def convert_to_mirror(system: System, centering=True):
    if centering:
        center = [mean([a.xyz[i] for a in system.atoms]) for i in range(3)]
        prec = max(max(len(str(_a.xyz[i]).split(".")[1]) for _a in system.atoms) for i in range(3))
        center = [round((-1) * _v, prec) for _v in center]
    for a in system.atoms:
        xyz = [(-1) * a.x, (-1) * a.y, (-1) * a.z]
        if centering:
            xyz = [float_to_str(round(_v - center[i], prec)) for i, _v in enumerate(xyz)]
        xyz = [float(float_to_str(_v)) for _v in xyz]
        a.x = xyz[0]
        a.y = xyz[1]
        a.z = xyz[2]
