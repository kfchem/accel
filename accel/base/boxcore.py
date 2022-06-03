import csv
import math
import shutil
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence, Union

from accel.base import formats as formats
from accel.base import text as text
from accel.base import topology as topology
from accel.base import xyz as xyz
from accel.base.mols import Mol, Mols
from accel.base.selector import Selectors
from accel.base.tools import change_dir
from accel.util.constants import CONSTANTS, Unit, Units
from accel.util.datadict import Data
from accel.util.log import Log, logger
from accel.util.matrix import Matrix


def _add(_c, mulcos: "BoxCore"):
    if isinstance(_c, Iterable) and not isinstance(_c, str):
        for _ct in _c:
            _add(_ct, mulcos)
        return None
    if isinstance(_c, str) or isinstance(_c, Path):
        _pc = Path(_c)
        if "*" in _pc.name or "?" in _pc.name:
            return _add(sorted(_pc.parent.glob(_pc.name)), mulcos)
        if _pc.is_dir():
            return _add(_pc.iterdir(), mulcos)
    mulcos._mols.append(_c)
    logger.info(f"added: {str(mulcos._mols[-1].path.absolute())}")
    return None


class BoxCore:
    def __init__(self, files: Union[Iterable[Union[Mol, Path, str]], str, "BoxCore"] = None):
        if isinstance(files, BoxCore):
            self._mols = files._mols
            self.data = files.data
        else:
            self._mols: Mols = Mols()
            self.data = Data(self)
            if files is not None:
                _add(files, self)

    @property
    def mols(self) -> Mols:
        return self._mols.has_state(True)

    @property
    def allmols(self) -> Mols:
        return self._mols.has_state(None)

    def __len__(self) -> int:
        return len(self._mols)

    def __str__(self) -> str:
        return f"Box: {len(self.mols)}/{len(self._mols)} molecules"

    def add(self, contents: Any):
        if isinstance(contents, BoxCore):
            _add(contents.mols, self)
        else:
            _add(contents, self)
        logger.debug(f"done: {str(self)}")
        return self

    def bind(self, mols: Mols):
        if isinstance(mols, Mols):
            self._mols = mols
        else:
            logger.error("bind method accepts only Mols instance")
            raise TypeError
        return self

    def show(self):
        self._mols.show()
        return self

    def labeling(self, separator: str = "_", index_list: List[int] = [0, 1]):
        for _c in self._mols:
            sp_stem = _c.name.split(separator)
            _c.label = separator.join([sp_stem[i] for i in index_list if i < len(sp_stem)])
        logger.debug(f"done: {str(self)}")
        return self

    def set_data(self, key: str, value=None):
        for _c in self._mols:
            _c.data[key] = value
        logger.debug(f"done: {str(self)}")
        return self

    def set_state(self, flag: bool = True):
        for _c in self._mols:
            _c.flag = flag
        logger.debug(f"done: {str(self)}")
        return self

    def set_label(self, label: str = ""):
        for _c in self._mols:
            _c.label = label
        logger.debug(f"done: {str(self)}")
        return self

    def zero_fill(self, digit: int = 3, sepalator: str = "_", position: int = 3):
        if int(digit) <= 0:
            logger.error("invalid digit")
            return self
        for _c in self.mols:
            word_list = _c.name.split(sepalator)
            for i, _n in enumerate(word_list):
                if i == position - 1:
                    word_list[i] = _n.zfill(digit)
            _c.name = sepalator.join(word_list)
        logger.debug(f"done: {str(self)}")
        return self

    def count(self, comment: str = ""):
        labels = {}
        for label in self._mols.labels.keys():
            labels[label] = len(self.mols.labels.get(label, []))
        _l = ""
        for _label, _i in labels.items():
            _l += f"{_label}: {_i}, "
        _l += f"all: {sum(list(labels.values()))} conformers"
        if comment != "":
            _l = f"{comment}: {_l}"
        logger.info(_l)
        return self

    def duplicate(self):
        _n = self.__class__()
        _n._mols = self._mols.duplicate(_n)
        _n.data = self.data.duplicate(_n)
        return _n

    def export_data(self, filepath: Path):
        exported_rows = []
        column_keys = [
            "name",
            "path",
            "filetype",
            "state",
            "history",
            "label",
            "energy",
            "charge",
            "multiplicity",
        ]
        for _c in self._mols:
            data_dict = {}
            data_dict["name"] = _c.name
            data_dict["path"] = _c.path
            data_dict["filetype"] = _c.filetype
            data_dict["state"] = _c.flag
            data_dict["history"] = _c.history
            data_dict["label"] = _c.label
            data_dict["energy"] = _c.energy
            data_dict["charge"] = _c.charge
            data_dict["multiplicity"] = _c.multiplicity
            for _key, _val in _c.data._data.items():
                str_val = str(_val)
                if len(str_val) >= 256:
                    data_dict[_key] = "too long data"
                else:
                    data_dict[_key] = str_val
                if _key not in column_keys:
                    column_keys.append(_key)
            exported_rows.append(data_dict)
        _p = Path(filepath).with_suffix(".csv")
        with _p.open("w", newline="") as f:
            _w = csv.DictWriter(f, column_keys)
            _w.writeheader()
            _w.writerows(exported_rows)
        logger.info(f"{_p.name} was exported")
        Log.set_output_dir(_p.parent)
        return self

    def energy_limit(self, threshold=3.0, max_limit: int = None, in_label: bool = True):
        threshold = float(threshold)
        if max_limit is None:
            max_limit = 0
        max_limit = int(max_limit)
        if in_label:
            if len(self.mols) != len(self.mols.has_label()):
                logger.info("labeling called automatically")
                self.labeling()
            for confs in self.mols.labels.values():
                min_e = min(_c.energy for _c in confs)
                logger.info(f"the minimun energy: {min_e}")
                for _c in confs:
                    if (_c.energy - min_e) >= threshold:
                        _c.deactivate("energy_limit")
            if int(max_limit) > 0:
                for confs in self.mols.labels.values():
                    sorted_confs: List[Mol] = sorted(confs, key=lambda t: t.energy)
                    for rank, _c in enumerate(sorted_confs, start=1):
                        if rank > max_limit:
                            logger.info(f"{_c.name}: {rank}th stable (max_limit={max_limit})")
                            _c.deactivate("max limit in energy_limit")
        else:
            confs = self.mols
            min_e = min(_c.energy for _c in confs)
            logger.info(f"the minimun energy: {min_e}")
            for _c in confs:
                if (_c.energy - min_e) >= threshold:
                    _c.deactivate("energy_limit")
            if int(max_limit) > 0:
                sorted_confs: List[Mol] = sorted(confs, key=lambda t: t.energy)
                for rank, _c in enumerate(sorted_confs, start=1):
                    if rank > max_limit:
                        logger.info(f"{_c.name}: {rank}th stable (max_limit={max_limit})")
                        _c.deactivate("max limit in energy_limit")
        logger.debug(f"done: {str(self)}")
        return self

    def calc_rel_energy(self, in_label: bool = True):
        if in_label:
            if len(self.mols) != len(self.mols.has_label()):
                logger.info("labeling called automatically")
                self.labeling()
            for confs in self.mols.labels.values():
                min_e = min(_c.energy for _c in confs)
                logger.info(f"the minimun energy: {min_e}")
                for _c in confs:
                    _c.energy = _c.energy - min_e
        else:
            min_e = min(_c.energy for _c in self.mols)
            logger.info(f"the minimun energy: {min_e}")
            for _c in self.mols:
                _c.energy = _c.energy - min_e
        logger.debug(f"done: {str(self)}")
        return self

    def calc_distribution(self, in_label: bool = True, temperature: float = 298.15):
        def _calc_distr(confs: List[Mol]):
            if len(confs) == 0:
                logger.error("there are no conformers in the Confs")
                return False
            min_e = min(_c.energy for _c in confs)
            sum_bfac = 0.0
            for _c in confs:
                _c.cache["relative_energy"] = _c.energy - min_e
                _c.cache["distr_factor"] = math.exp(
                    -1
                    * (_c.cache["relative_energy"] * CONSTANTS.cal_to_J * 1000)
                    / (temperature * CONSTANTS.gas_constant_in_J_over_mol_K)
                )
                sum_bfac += _c.cache["distr_factor"]
            for _c in confs:
                _c.data["distribution"] = _c.cache["distr_factor"] / sum_bfac

        if in_label:
            if len(self.mols) != len(self.mols.has_label()):
                logger.info("labeling called automatically")
                self.labeling()
            for confs in self.mols.labels.values():
                _calc_distr(confs)
        else:
            _calc_distr(self.mols)
        logger.debug(f"done: {str(self)}")
        return self

    def calc_energy(self, keys: List[str] = [], unit: Unit = Units.kcal_mol):
        for _c in self.mols:
            _energy = 0.0
            for data_key in keys:
                _energy += _c.data[data_key]
            _c.energy = unit(_energy).to_kcal_mol
        logger.debug(f"done: {str(self)}")
        return self

    def copy_files(self, directory: Path, change_path: bool = False, suffix: str = None):
        for _c in self.mols:
            _p = change_dir(_c.path, directory)
            _p = _p.with_name(_c.name).with_suffix(_c.path.suffix)
            if suffix is not None:
                _p = _p.with_suffix(suffix)
            shutil.copy2(_c.path, _p)
            logger.info(f"{_c.path} was copied to {_p}")
            if change_path:
                _c.path = _p
        logger.debug(f"done: {str(self)}")
        return self

    def search(self, directory: Path = None, existing_check: bool = True, suffix: str = None):
        for _c in self.mols:
            if directory is not None and suffix is not None:
                _p = Path(directory).joinpath(_c.name).with_suffix(suffix)
            elif directory is None and suffix is not None:
                _p = _c.path.with_name(_c.name).with_suffix(suffix)
            elif directory is not None and suffix is None:
                _p = Path(directory).joinpath(_c.name).with_suffix(_c.path.suffix)
            if existing_check:
                if _p.exists():
                    _c.path = _p
                else:
                    logger.error(f"{_p.name} does not exist")
                    _c.deactivate("search_dir")
            else:
                _c.path = _p
        logger.debug(f"done: {str(self)}")
        return self

    @Selectors.read_atoms.add("format/xyz")
    def read_xyz(self):
        for _c in self.mols:
            formats.read_xyz(_c)
        logger.debug(f"done: {str(self)}")
        return self

    def write_xyz(self, directory: Path = None, link: bool = True, centering: bool = True):
        for _c in self.mols:
            formats.write_xyz(_c, directory, link, centering)
        Log.set_output_dir(directory)
        logger.debug(f"done: {str(self)}")
        return self

    @Selectors.read_atoms.add("format/mol")
    def read_mol(self):
        for _c in self.mols:
            formats.read_mol(_c)
        logger.debug(f"done: {str(self)}")
        return self

    def write_mol(self, directory: Path = None, link: bool = True, centering: bool = True):
        if len(self.mols) != len(self.mols.has_bonds()):
            self.calc_bonds()
        for _c in self.mols:
            formats.write_mol(_c, directory, link, centering)
        Log.set_output_dir(directory)
        logger.debug(f"done: {str(self)}")
        return self

    def write_input(self, template: Path, directory=None, link: bool = True, arg: Dict[str, str] = None):
        for _c in self.mols:
            text.write_input(_c, template, directory, link)
        Log.set_output_dir(directory)
        logger.debug(f"done: {str(self)}")
        return self

    def calc_bonds(self, cov_scaling: float = 1.1, vdw_scaling: float = 1.0, aromatize=True):
        for _c in self.mols:
            topology.embed_bonds(_c, cov_scaling, vdw_scaling)
            if aromatize:
                topology.aromatize(_c)
        logger.debug(f"done: {str(self)}")
        return self

    def calc_symm(self, calc_all: bool = False):
        if not calc_all:
            cfs = []
            if len(self.mols) != len(self.mols.has_label()):
                logger.info("labeling called automatically")
                self.labeling()
            for _confs in self.mols.labels.values():
                cfs.append(_confs.get())
                if len(_confs) != 1:
                    logger.info(
                        f"symmetry information of {_confs.get().name} will be duplicated from initial confomer to all"
                    )
            cfs = Mols().bind(cfs)
        else:
            cfs = self.mols

        if len(cfs) != len(cfs.has_bonds()):
            logger.info("embed_bonds called automatically")
            self.calc_bonds()

        for _c in cfs:
            topology.embed_symm(_c)
            _c.data["has_symm"] = True

        if not calc_all:
            rot_mats_dict: Dict[str, List[Matrix]] = {_c.label: _c.data["rotamer"] for _c in cfs}
            num_mats_dict: Dict[str, List[Matrix]] = {_c.label: _c.data["numisomer"] for _c in cfs}
            for _label, _confs in self.mols.labels.items():
                for _c in _confs:
                    atoms_list = _c.atoms.to_list()
                    _c.data["rotamer"] = [Matrix(atoms_list).bind(_m._matrix) for _m in rot_mats_dict[_label]]
                    _c.data["numisomer"] = [Matrix(atoms_list).bind(_m._matrix) for _m in num_mats_dict[_label]]
                    _c.data["has_symm"] = True
        logger.debug(f"done: {str(self)}")
        return self

    def rmsd_limit(
        self,
        threshold: float = 0.01,
        all_combinations_of_confs: bool = False,
        redundant_check: int = 3,
        all_perturbation_of_rotamers: bool = False,
    ):
        if len(self.mols) != len(self.mols.has_data("has_symm", True)):
            logger.info("embed_symm called automatically")
            self.calc_symm()
        topology.rmsd_pruning(
            self.mols,
            threshold,
            all_combinations_of_confs,
            redundant_check,
            all_perturbation_of_rotamers,
        )
        logger.debug(f"done: {str(self)}")
        return self

    def map_numbers(self, reference_box: "BoxCore" = None):
        if reference_box is None:
            reference_box = BoxCore([_confs.get() for _confs in self.mols.labels.values()])
        r_mulcos = reference_box.duplicate()
        r_mulcos.calc_bonds()
        r_mulcos.calc_symm()
        topology.map_numbers(self.mols, r_mulcos._mols)
        logger.debug(f"done: {str(self)}")
        return self

    def calc_length(self, number_a: int, number_b: int, key: str = ""):
        number_a = int(number_a)
        number_b = int(number_b)
        for _c in self.mols:
            xyz.calc_length(_c, number_a, number_b, key)
        logger.debug(f"done: {str(self)}")
        return self

    def calc_dihedral(self, number_a: int, number_b: int, number_c: int, number_d: int, key: str = ""):
        number_a = int(number_a)
        number_b = int(number_b)
        number_c = int(number_c)
        number_d = int(number_d)
        for _c in self.mols:
            xyz.calc_dihedral(_c, number_a, number_b, number_c, number_d, key)
        logger.debug(f"done: {str(self)}")
        return self

    def calc_angle(self, number_a: int, number_b: int, number_c: int, key: str = ""):
        number_a = int(number_a)
        number_b = int(number_b)
        number_c = int(number_c)
        for _c in self.mols:
            xyz.calc_angle(_c, number_a, number_b, number_c, key)
        logger.debug(f"done: {str(self)}")
        return self

    def get_average(self, keys: List[str] = [], keys_for_atoms: List[str] = []) -> Mols:
        new_ls: List[Mol] = []
        if len(self.mols) != len(self.mols.has_data("distribution")):
            logger.info("all of active conformers do not have distribution data")
            logger.info("calc_distr called automatically")
            self.calc_distribution()
        for label, _ref_confs in self.mols.labels.items():
            _c = Mol()
            _ref_one = _ref_confs.get()
            new_ls.append(_c)
            _c.name = label
            _c.energy = _ref_one.energy
            _c.label = label
            for _ref_a in _ref_one.atoms:
                _c.atoms.get_new_atom().axyz = _ref_a.axyz
            for _ref_c in _ref_confs:
                if len(_ref_c.atoms) != len(_c.atoms):
                    logger.error("atomic symbol does not matched")
                    raise Exception
                for _ref_a in _ref_c.atoms:
                    if _c.atoms.get(_ref_a.number).symbol != _ref_a.symbol:
                        logger.error("atomic symbol does not matched")
                        raise Exception
        for _key in keys:
            if len(self.mols) != len(self.mols.has_data(_key)):
                logger.info(f"all of active conformers do not have {_key} data")
                continue
            for _c in new_ls:
                for _oc in self.mols.labels.get(_c.name):
                    _c.data[_key] += _oc.data[_key] * _oc.data["distribution"]
        for _key in keys_for_atoms:
            for _c in new_ls:
                for _a in _c.atoms:
                    sum_val = 0.0
                    for _oc in self.mols.labels.get(_c.name):
                        tmp_val = _oc.atoms.get(_a.number).data.get(_key)
                        if tmp_val is None:
                            logger.error(f"atom {_oc.atoms.get(_a.number)} of {_oc} have no {_key} data")
                            break
                        sum_val += float(tmp_val) * _oc.data["distribution"]
                    else:
                        _a.data[_key] = sum_val
        return Mols().bind(new_ls)

    def modify_length(
        self,
        number_a: int,
        number_b: int,
        target: float,
        fix_a: bool = False,
        fix_b: bool = False,
        numbers_along_with_a: Sequence[int] = [],
        numbers_along_with_b: Sequence[int] = [],
    ):
        for _c in self.mols:
            xyz.edit_bond_length(
                _c, number_a, number_b, target, (bool(fix_a), bool(fix_b)), numbers_along_with_a, numbers_along_with_b
            )
        logger.debug(f"done: {str(self)}")
        return self

    def convert_to_mirror(self, centering: bool = True):
        for _c in self.mols:
            xyz.convert_to_mirror(_c, centering)
        logger.debug(f"done: {str(self)}")
        return self

    def only_minimum(self, in_label: bool = True):
        if in_label:
            if len(self.mols) != len(self.mols.has_label()):
                logger.info("labeling called automatically")
                self.labeling()
            for confs in self.mols.labels.values():
                _list = sorted(confs, key=lambda _c: _c.energy)
                for _c in _list[1:]:
                    _c.deactivate("only_minimum")
        else:
            _list = sorted(self.mols, key=lambda _c: _c.energy)
            for _c in _list[1:]:
                _c.deactivate("only_minimum")
        logger.debug(f"done: {str(self)}")
        return self
