import math
from statistics import mean
from typing import List, Sequence, Tuple

import numpy as np
from accel.base.atoms import Atom
from accel.base.systems import System
from accel.base.tools import float_to_str


def edit_bond_length(
    _c: System,
    atom_a: int,
    atom_b: int,
    target_length: float,
    fixed_atom: Tuple[bool, bool] = (False, False),
    move_along_with_a: Sequence[int] = (),
    move_along_with_b: Sequence[int] = (),
):
    _vect = [0.0, 0.0, 0.0]
    for i in range(3):
        _vect[i] = _c.atoms.get(atom_b).xyz[i] - _c.atoms.get(atom_a).xyz[i]
    _dist = math.sqrt(sum(x**2 for x in _vect))

    def move_atoms(atoms_list, vect_factor):
        for atom_no in atoms_list:
            _c.atoms.get(atom_no).xyz = [
                _val + (vect_factor * _vect[i] * (_dist - target_length) / _dist)
                for i, _val in enumerate(_c.atoms.get(atom_no).xyz)
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


def set_chirality(_c: System, center_index: int, sub_index: list[int]):
    if len(sub_index) != 4:
        raise ValueError
    else:
        sorted_index = sorted(sub_index)

    _sub_xyzs = np.array([_c.atoms.get(i).xyz for i in sorted_index[1:]]) - np.array(
        [_c.atoms.get(sorted_index[0]).xyz for _ in range(3)]
    )
    _ret = np.linalg.det(_sub_xyzs)
    if _ret > 0:
        _ret = 1
    elif _ret < 0:
        _ret = -1
    else:
        _ret = 0
    _c.data[f"chiral_{center_index}_to_{sub_index}"] = _ret


def calc_length(_c: System, atom_index_a: int, atom_index_b: int, key: str = ""):

    _a = _c.atoms.get(atom_index_a).xyz
    _b = _c.atoms.get(atom_index_b).xyz
    _d = [float(_a[i]) - float(_b[i]) for i in range(3)]
    _dist = math.sqrt(sum(x**2 for x in _d))
    if key == "" or not isinstance(key, str):
        key = "distance_{}{}-{}{}".format(
            _c.atoms.get(atom_index_a).symbol,
            str(atom_index_a),
            _c.atoms.get(atom_index_b).symbol,
            str(atom_index_b),
        )
    _c.data[key] = _dist


def get_dihedral(atom_a: Atom, atom_b: Atom, atom_c: Atom, atom_d: Atom) -> float:
    _va = np.array(atom_a.xyz)
    _vb = np.array(atom_b.xyz)
    _vc = np.array(atom_c.xyz)
    _vd = np.array(atom_d.xyz)
    _vab = _va - _vb
    _vcb = _vc - _vb
    _vdc = _vd - _vc
    _pvac = np.cross(_vab, _vcb)
    _pvbd = np.cross(_vdc, _vcb)
    _dac = np.linalg.norm(_pvac)
    _dbd = np.linalg.norm(_pvbd)
    _angle = np.arccos(np.sum(_pvac * _pvbd) / (_dac * _dbd))
    if np.sum(_pvac * np.cross(_pvbd, _vcb)) < 0:
        _angle = -_angle
    _angle = float(np.rad2deg(_angle))
    return _angle


def calc_dihedral(
    _c: System,
    atom_index_a: int,
    atom_index_b: int,
    atom_index_c: int,
    atom_index_d: int,
    key: str = "",
):
    if key == "" or not isinstance(key, str):
        key = "dihedral_{}{}-{}{}-{}{}-{}{}".format(
            _c.atoms.get(atom_index_a).symbol,
            str(atom_index_a),
            _c.atoms.get(atom_index_b).symbol,
            str(atom_index_b),
            _c.atoms.get(atom_index_c).symbol,
            str(atom_index_c),
            _c.atoms.get(atom_index_d).symbol,
            str(atom_index_d),
        )
    _c.data[key] = get_dihedral(
        _c.atoms.get(atom_index_a),
        _c.atoms.get(atom_index_b),
        _c.atoms.get(atom_index_c),
        _c.atoms.get(atom_index_d),
    )


def get_angle(atom_a: Atom, atom_b: Atom, atom_c: Atom) -> float:
    _va = np.array(atom_a.xyz)
    _vb = np.array(atom_b.xyz)
    _vc = np.array(atom_c.xyz)
    _vba = _vb - _va
    _vbc = _vb - _vc
    _dba = np.linalg.norm(_vba)
    _dbc = np.linalg.norm(_vbc)
    _angle = np.arccos(np.sum(_vba * _vbc) / (_dba * _dbc))
    _angle = float(np.rad2deg(_angle))
    return _angle


def calc_angle(
    _c: System,
    atom_index_a: int,
    atom_index_b: int,
    atom_index_c: int,
    key: str = "",
):
    if key == "" or not isinstance(key, str):
        key = "angle_{}{}-{}{}-{}{}".format(
            _c.atoms.get(atom_index_a).symbol,
            str(atom_index_a),
            _c.atoms.get(atom_index_b).symbol,
            str(atom_index_b),
            _c.atoms.get(atom_index_c).symbol,
            str(atom_index_c),
        )
    _c.data[key] = get_angle(
        _c.atoms.get(atom_index_a),
        _c.atoms.get(atom_index_b),
        _c.atoms.get(atom_index_c),
    )


def convert_to_mirror(_c: System, centering=True):
    if centering:
        _cnt = [mean([_a.xyz[i] for _a in _c.atoms]) for i in range(3)]
        _prec = max(max(len(str(_a.xyz[i]).split(".")[1]) for _a in _c.atoms) for i in range(3))
        _cnt = [round((-1) * _v, _prec) for _v in _cnt]
    for _a in _c.atoms:
        _xyz = [(-1) * _a.x, (-1) * _a.y, (-1) * _a.z]
        if centering:
            _xyz = [float_to_str(round(_v - _cnt[i], _prec)) for i, _v in enumerate(_xyz)]
        _xyz = [float(float_to_str(_v)) for _v in _xyz]
        _a.x = _xyz[0]
        _a.y = _xyz[1]
        _a.z = _xyz[2]
