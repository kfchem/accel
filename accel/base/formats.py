import datetime
from pathlib import Path
from statistics import mean

from accel.base.atoms import BondType
from accel.base.mols import Mol
from accel.base.tools import change_dir, float_to_str
from accel.util import FileType
from accel.util.constants import Elements
from accel.util.log import logger


def read_xyz(_c: Mol):
    with _c.path.open() as f:
        _ls = f.readlines()
    _axyz = [_l.split() for _l in _ls[2:] if len(_l.split()) == 4]
    if len(_axyz) != int(_ls[0]):
        _c.deactivate("read_atoms: orca from xyz file")
        return None
    _c.atoms.clear()
    for _l in _axyz:
        _c.atoms.append(_l)
    logger.debug(f"reading xyz of {_c.path.name}")


def read_xyz_coordinates(_c: Mol):
    with _c.path.open() as f:
        _ls = f.readlines()
    _axyz = [_l.split() for _l in _ls[2:] if len(_l.split()) == 4]
    for i, _l in enumerate(_axyz):
        _a = _c.atoms.get(i + 1)
        if _a.symbol != Elements.canonicalize(_l[0]):
            logger.error(
                "symbol did not matched on {}: original {}: xyz {}: coordinates did not updated".format(
                    _a, _a.symbol, _l[0]
                )
            )
            continue
        _a.xyz = _l[1:]
    logger.debug(f"reading xyz coordinate of {_c.path.name}")


def write_xyz(_c: Mol, output_dir=None, change_path=False, centering=True):
    _ls = [str(len(_c.atoms)) + "\n"]
    _ls.append(_c.name + "\n")
    if centering:
        _cnt = [mean([_a.xyz[i] for _a in _c.atoms]) for i in range(3)]
        try:
            _prec = max(max(len(float_to_str(_a.xyz[i]).split(".")[1]) for _a in _c.atoms) for i in range(3))
        except IndexError:
            logger.error(f"could not resolve the precision in converting to xyz file of {_c.name}")
            _prec = 10
        _cnt = [round(_v, _prec) for _v in _cnt]
    for _a in _c.atoms:
        _xyz = [_a.x, _a.y, _a.z]
        if centering:
            _xyz = [float_to_str(round(_v - _cnt[i], _prec)) for i, _v in enumerate(_xyz)]
        _xyz = [float_to_str(_v) for _v in _xyz]
        _ls.append(f"{_a.symbol:<2} {_xyz[0]:>15} {_xyz[1]:>15} {_xyz[2]:>15}\n")
    _p = change_dir(_c.path, output_dir, _c.name).with_suffix(".xyz")
    with _p.open("w", newline="\n") as f:
        f.writelines(_ls)
    logger.info(f"{_c.name} was exported to {_p.name}")
    if change_path:
        _c.path = _p


@FileType.add("format/xyz", 10)
def is_xyz_file(_p: Path):
    if _p.suffix != ".xyz":
        return False
    return True


conv_sdf_charge = [0, 3, 2, 1, 0, -1, -2, -3]


def read_mol(_c: Mol):
    mol_data = {_key: [] for _key in ["header", "counts", "atom", "bond", "middle", "properties"]}
    with _c.path.open("r") as f:
        mol_lines = f.readlines()
    mol_data["header"] = {}
    mol_data["header"]["name"] = mol_lines[0].replace("\n", "")
    mol_data["header"]["meta"] = mol_lines[1].replace("\n", "")
    mol_data["header"]["comments"] = mol_lines[2].replace("\n", "")
    for i in range(11):
        mol_data["counts"].append(mol_lines[3][i * 3 : (i + 1) * 3])
    mol_data["counts"].append(mol_lines[3][33:39])
    l_start = 4
    l_end = l_start + int(mol_data["counts"][0])
    for atom_line in mol_lines[l_start:l_end]:
        _a = []
        for i in range(3):
            _a.append(atom_line[i * 10 : (i + 1) * 10])
        _a.append(atom_line[31:34])
        _a.append(atom_line[34:36])
        for i in range(11):
            _a.append(atom_line[(i * 3) + 36 : ((i + 1) * 3) + 36])
        mol_data["atom"].append(_a)
    l_start = l_end
    l_end = l_start + int(mol_data["counts"][1])
    for bond_line in mol_lines[l_start:l_end]:
        _bond = []
        for i in range(7):
            _bond.append(bond_line[i * 3 : (i + 1) * 3].replace("\n", ""))
        mol_data["bond"].append(_bond)
    l_start = l_end
    l_end = l_start + int(mol_data["counts"][2]) + int(mol_data["counts"][5])
    mol_data["middle"] = mol_lines[l_start:l_end]
    l_start = l_end
    l_end = l_start + mol_lines.index("M  END\n") + 1
    mol_data["properties"] = mol_lines[l_start:l_end]
    _c.data["sdf"] = mol_data

    _c.atoms.clear()
    for _a in _c.data["sdf"]["atom"]:
        _new_a = _c.atoms.get_new_atom([_a[3], _a[0], _a[1], _a[2]])
        if 1 <= int(_a[5]) and int(_a[5]) <= 7:
            _new_a.charge = conv_sdf_charge[int(_a[5])]
    for _bond in _c.data["sdf"]["bond"]:
        _number_a = int(_bond[0])
        _number_b = int(_bond[1])
        if int(_bond[2]) == 1:
            _c.atoms.add_bond(_number_a, _number_b, BondType.single)
        elif int(_bond[2]) == 2:
            _c.atoms.add_bond(_number_a, _number_b, BondType.double)
        elif int(_bond[2]) == 3:
            _c.atoms.add_bond(_number_a, _number_b, BondType.triple)
        elif int(_bond[2]) == 4:
            _c.atoms.add_bond(_number_a, _number_b, BondType.aromatic)
        else:
            _c.atoms.add_bond(_number_a, _number_b, BondType.undefined)
    for _prop in _c.data["sdf"]["properties"]:
        _prop: str = _prop
        if _prop.startswith("M  CHG"):
            _prop = _prop.split()
            for i in range(int(_prop[2])):
                _c.atoms.get(int(_prop[(2 * i) + 3])).charge = int(_prop[(2 * i) + 4])
    _total_charge = 0
    for _a in _c.atoms:
        _total_charge += _a.charge
    _c.total_charge = _total_charge
    logger.debug(f"read {_c.name}")


def get_mol_str(_c: Mol, centering=True) -> str:
    _ls = []
    _ls.append(_c.name + "\n")
    _now = datetime.datetime.now()
    if _c.energy is None:
        _energy = 0.0
    else:
        _energy = _c.energy
    _energy = float_to_str(_energy)
    if len(_energy) > 12:
        _energy = _energy[:12]
    _ls.append("  ACCeL   {}3D            {:12}      \n".format(_now.strftime("%m%d%y%H%M"), _energy))
    _ls.append("\n")
    _ls.append(f"{len(_c.atoms):>3d}{len(_c.atoms.bonds):>3d}  0  0  1  0            999 V2000\n")
    _prec = 4
    if centering:
        _cnt = [mean([_a.xyz[i] for _a in _c.atoms]) for i in range(3)]
        _cnt = [round(_v, _prec) for _v in _cnt]
    for _a in _c.atoms:
        _xyz = [_a.x, _a.y, _a.z]
        if centering:
            _xyz = [float_to_str(round(_v - _cnt[i], _prec)) for i, _v in enumerate(_xyz)]
        else:
            _xyz = [float_to_str(round(_v, _prec)) for _v in _xyz]
        try:
            _chg = conv_sdf_charge.index(int(_a.charge))
        except ValueError:
            logger.error(f"could not recognize charge: {_a.charge}")
            _chg = 0
        _ls.append(
            "{:>10.4f}{:>10.4f}{:>10.4f} {:<3s} 0{:>3d}  0  0  0  0                  \n".format(
                float(_xyz[0]), float(_xyz[1]), float(_xyz[2]), _a.symbol, _chg
            )
        )
    for _from_to, _type in _c.atoms.bonds.items():
        if _type in [BondType.undefined, BondType.contact]:
            continue
        _ls.append(f"{_from_to[0]:>3d}{_from_to[1]:>3d}{_type:>3d}  0  0  0  0\n")
    _ls.append("M  END\n")
    return "".join(_ls)


def write_mol(_c: Mol, output_dir=None, change_path=False, centering=True) -> Path:
    mol_str = get_mol_str(_c, centering)
    _p = change_dir(_c.path, output_dir, _c.name).with_suffix(".sdf")
    with _p.open("w", newline="\n") as f:
        f.write(mol_str)
    logger.info(f"{_c.name} was exported to {_p.name}")
    if change_path:
        _c.path = _p
    return _p


@FileType.add("format/mol", 11)
def is_mol_file(_p: Path) -> bool:
    if _p.suffix not in (".sdf", ".mol", ".sd"):
        return False
    return True
