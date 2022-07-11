import gzip
import platform
import subprocess
import tempfile
from collections import OrderedDict
from pathlib import Path
from typing import Dict, List, Union

from accel.base.atoms import BondType
from accel.base.boxcore import BoxCore
from accel.base.formats import write_mol
from accel.base.mols import Mol, Mols
from accel.base.selector import Selectors
from accel.base.tools import change_dir
from accel.util import Execmd, FileType
from accel.util.constants import Elements, Unit, Units
from accel.util.log import logger


def _get_f_m_ct_list(lines: List[str]) -> List[Dict[str, Union[str, List]]]:
    ls = "".join(lines).replace("\n", " ").split()
    words: List[str] = []
    while True:
        try:
            words.append(ls.pop())
        except IndexError:
            break
        for ecs in ("#", '"'):
            while words[-1].count(ecs) == 1:
                words[-1] = " ".join([ls.pop(), words[-1]])
    f_m_ct_list = []
    while len(words) != 0:
        nw = words.pop()
        if "s_m_m2io_version" == nw:
            if words.pop() != ":::":
                raise ValueError
            nw = words.pop()
            if not ("2.0.0" == nw):
                logger.error("{}: s_m_m2io_version was not 2.0.0")
                raise ValueError
                break
        if "f_m_ct" == nw:
            if words.pop() != "{":
                raise ValueError
            dict_f_m_ct = OrderedDict()
            nw = words.pop()
            while nw != ":::":
                dict_f_m_ct[nw] = ""
                nw = words.pop()
            for key in dict_f_m_ct:
                dict_f_m_ct[key] = words.pop()
            nw = words.pop()
            while nw != "}":
                internal_list = []
                internal_dict_keys = []
                dict_f_m_ct[nw.split("[", 1)[0]] = internal_list
                num_of_content = int(nw.split("[", 1)[1].replace("]", ""))
                if words.pop() != "{":
                    raise ValueError
                nw = words.pop()
                while nw != ":::":
                    internal_dict_keys.append(nw)
                    nw = words.pop()
                for _ in range(num_of_content):
                    internal_dict = {}
                    for key in internal_dict_keys:
                        nw = words.pop()
                        internal_dict[key] = nw
                    internal_list.append(internal_dict)
                if words.pop() != ":::":
                    raise ValueError
                if words.pop() != "}":
                    raise ValueError
                nw = words.pop()
            f_m_ct_list.append(dict_f_m_ct)
            logger.debug("found {} No.{}".format(dict_f_m_ct.get("s_m_title", "unknown"), len(f_m_ct_list)))

    return f_m_ct_list


def _embed_f_m_ct(_c: Mol, f_m_ct: Dict[str, Union[str, List]]):
    for key in f_m_ct:
        if not isinstance(f_m_ct[key], str):
            continue
        ks = key.split("_", 1)
        if len(ks) == 1:
            _c.data[key] = f_m_ct[key]
        elif ks[0] == "s":
            _c.data[ks[1]] = f_m_ct[key]
        elif ks[0] == "r":
            _c.data[ks[1]] = float(f_m_ct[key])
        elif ks[0] == "b":
            _c.data[ks[1]] = bool(int(f_m_ct[key]))
        elif ks[0] == "i":
            _c.data[ks[1]] = int(f_m_ct[key])
        else:
            _c.data[key] = f_m_ct[key]
    for _m_atom in f_m_ct["m_atom"]:
        _c.atoms.append(
            [
                Elements.canonicalize(int(_m_atom["i_m_atomic_number"])),
                float(_m_atom["r_m_x_coord"]),
                float(_m_atom["r_m_y_coord"]),
                float(_m_atom["r_m_z_coord"]),
            ]
        )
        if "i_m_formal_charge" in _m_atom:
            _c.atoms[-1].charge = int(_m_atom["i_m_formal_charge"])
    _total_charge = 0
    for _a in _c.atoms:
        if _a.charge is None:
            continue
        _total_charge += _a.charge
    _c.total_charge = _total_charge
    for _m_bond in f_m_ct["m_bond"]:
        _number_a = int(_m_bond["i_m_from"])
        _number_b = int(_m_bond["i_m_to"])
        if int(_m_bond["i_m_order"]) == 1:
            _c.atoms.add_bond(_number_a, _number_b, BondType.single)
        elif int(_m_bond["i_m_order"]) == 2:
            _c.atoms.add_bond(_number_a, _number_b, BondType.double)
        elif int(_m_bond["i_m_order"]) == 3:
            _c.atoms.add_bond(_number_a, _number_b, BondType.triple)
        else:
            _c.atoms.add_bond(_number_a, _number_b, BondType.undefined)


def read_mae(_c: Mol):
    with _c.path.open() as f:
        ls = f.readlines()
    try:
        f_m_ct_list = _get_f_m_ct_list(ls)
    except ValueError:
        logger.error(f"error in parsing mae file {_c.path.name}")
        _c.deactivate("read_atoms: mae")
        return None
    if len(f_m_ct_list) != 1:
        logger.info(f"multiple molecule found in {_c.path.name}: used top one")
    f_m_ct = f_m_ct_list[0]
    _c.atoms.clear()
    _embed_f_m_ct(_c, f_m_ct)
    logger.debug(f"read data coordinates and bonds from {_c.path.name}")


@FileType.add("format/mae", 15)
def is_mae_format(_p: Path) -> bool:
    if _p.suffix not in (".mae"):
        return False
    return True


@FileType.add("format/maegz", 20)
def is_maegz_format(_p: Path) -> bool:
    if _p.suffix not in (".maegz"):
        return False
    return True


def read_potential_energy_from_mae(_c: Mol):
    with _c.path.open() as f:
        ls = f.readlines()
    try:
        f_m_ct_list = _get_f_m_ct_list(ls)
    except ValueError:
        logger.error(f"error in parsing mae file {_c.path.name}")
        return None
    if len(f_m_ct_list) != 1:
        logger.info(f"multiple molecule found in {_c.path.name}: used top one")
    for key in f_m_ct_list[0]:
        if key.startswith("r_mmod_Potential_Energy"):
            logger.info(f"{_c.name}: {key} was used for the potential energy")
            _c.energy = Units.kJ_mol(float(f_m_ct_list[0][key])).to_kcal_mol
            break
    else:
        _c.deactivate("no energy entry in mae file")
        return None


class MaeBox(BoxCore):
    @Selectors.read_atoms.add("format/mae")
    def read_mae(self):
        for _c in self.mols:
            read_mae(_c)
        logger.debug(f"done: {str(self)}")
        return self

    def read_maegz(self):
        for _c in self.mols:
            with gzip.open(_c.path, mode="rt") as f:
                ls = f.readlines()
            try:
                f_m_ct_list = _get_f_m_ct_list(ls)
            except ValueError:
                logger.error(f"error in parsing mae file {_c.path.name}")
                continue
            if len(f_m_ct_list) != 1:
                logger.info(f"multiple molecule found in {_c.path.name}: used top one")
            f_m_ct = f_m_ct_list[0]
            _c.atoms.clear()
            _embed_f_m_ct(_c, f_m_ct)
            logger.debug(f"read data coordinates and bonds from {_c.path.name}")
        logger.debug(f"done: {str(self)}")
        return self

    @Selectors.read_energy.add("format/mae")
    def read_energy(self):
        for _c in self.mols:
            read_potential_energy_from_mae(_c)
        logger.debug(f"done: {str(self)}")
        return self

    def get_unzip(self, zero_fill_digit=3, max_loading=999):
        boxls = []
        for _maegz in self.mols:
            logger.info(f"reading {_maegz.name}")
            with gzip.open(_maegz.path, mode="rt") as f:
                ls = f.readlines()
            try:
                f_m_ct_list = _get_f_m_ct_list(ls)
            except ValueError:
                logger.error(f"error in parsing mae file {_maegz.path.name}")
                continue
            if zero_fill_digit is None:
                max_zero = len(str(len(f_m_ct_list)))
            else:
                max_zero = zero_fill_digit
            for num, f_m_ct in enumerate(f_m_ct_list):
                if num >= max_loading:
                    logger.error(f"max loading count {max_loading} reached")
                    break
                new_c = Mol(file_path=_maegz.path)
                new_c.name = "{c}_{num:0={zero}}".format(c=f_m_ct["s_m_title"], num=num + 1, zero=max_zero)
                _embed_f_m_ct(new_c, f_m_ct)
                boxls.append(new_c)
                logger.info(f"read {new_c.name}")
        ret_box = Mols().bind(boxls)
        logger.debug(f"done: {str(self)}")
        return ret_box

    def write_mae(self, directory: Path = None):
        schrodinger_version = "2021-2"
        if platform.system() == "Windows":
            Execmd.add(
                "sdconvert", r"C:\Program Files\Schrodinger" + schrodinger_version + r"\utilities\sdconvert.exe"
            )
        elif platform.system() == "Linux":
            Execmd.add("sdconvert", r"/opt/schrodinger" + schrodinger_version + "/utilities/sdconvert")
        elif platform.system() == "Mac":
            pass
        with tempfile.TemporaryDirectory() as tmp_dir:
            for _c in self.mols:
                mol_path = write_mol(_c, output_dir=tmp_dir, centering=False)
                _p = change_dir(_c.path, directory, _c.name).with_suffix(".mae")
                _txt = [
                    Execmd.get("sdconvert"),
                    "-isd",
                    str(mol_path),
                    "-omae",
                    str(_p),
                ]
                if subprocess.run(_txt, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL).returncode != 0:
                    logger.info(f"error in conversion process of {_p}")
                    _c.deactivate("write_mae")
                else:
                    logger.info(f"successfully converted to {_p}")

        logger.debug(f"done: {str(self)}")
        return self

    def calc_energy(self, keys: List[str] = ["mmod_Potential_Energy-S-OPLS"], unit: Unit = Units.kJ_mol):
        return super().calc_energy(keys=keys, unit=unit)

    def is_mae_format(self):
        for _c in self.mols:
            if not is_mae_format(_c.path):
                _c.flag = False
        logger.debug(f"done: {str(self)}")
        return self

    def is_maegz_format(self):
        for _c in self.mols:
            if not is_maegz_format(_c.path):
                _c.flag = False
        logger.debug(f"done: {str(self)}")
        return self
