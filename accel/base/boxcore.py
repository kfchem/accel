import csv
import math
import shutil
from pathlib import Path
from typing import Iterable, Sequence, Union

from accel.base import formats as formats
from accel.base import topology as topology
from accel.base.selector import Selectors
from accel.base.systems import System, Systems
from accel.base.tools import change_dir
from accel.util.constants import CONSTANTS, Unit, Units
from accel.util.datadict import Data
from accel.util.log import Log, logger
from accel.util.matrix import Matrix


def _add(content, boxcore: "BoxCore"):
    if isinstance(content, Iterable) and not isinstance(content, str):
        for _ct in content:
            _add(_ct, boxcore)
        return None
    if isinstance(content, str) or isinstance(content, Path):
        _pc = Path(content)
        if "*" in _pc.name or "?" in _pc.name:
            return _add(sorted(_pc.parent.glob(_pc.name)), boxcore)
        if _pc.is_dir():
            return _add(_pc.iterdir(), boxcore)
    boxcore.contents.append(content)
    if boxcore.contents[-1].path is None:
        logger.info(f"added: {boxcore.contents[-1].name}")
    else:
        logger.info(f"added: {str(boxcore.contents[-1].path.absolute())}")
    return None


class BoxCore:
    def __init__(self, contents: Union[Iterable[Union[System, Path, str]], str, "BoxCore"] = None):
        if isinstance(contents, BoxCore):
            self.contents = contents.contents
            self.data = contents.data
        else:
            self.contents: Systems = Systems()
            self.data = Data(self)
            if contents is not None:
                _add(contents, self)

    def get(self, arg: Union[bool, str] = True) -> Systems:
        if isinstance(arg, bool) or arg is None:
            return self.contents.has_state(arg)
        elif isinstance(arg, str):
            return self.contents.has_label(arg)
        else:
            logger.error(f"invalid arg: {arg}")
            raise ValueError

    # for backward compatibility
    @property
    def _mols(self) -> Systems:
        logger.error("_mols is deprecated: use contents")
        return self.contents

    # for backward compatibility
    @property
    def mols(self) -> Systems:
        logger.error("mols is deprecated: use get()")
        return self.contents.has_state(True)

    # for backward compatibility
    @property
    def allmols(self) -> Systems:
        logger.error("allmols is deprecated: use get(None)")
        return self.contents.has_state(None)

    def __len__(self) -> int:
        return len(self.contents)

    def __str__(self) -> str:
        return f"Box: active/all = {len(self.get())}/{len(self.contents)}"

    def add(self, files: Union[Iterable[Union[System, Path, str]], str, "BoxCore"]):
        if isinstance(files, BoxCore):
            _add(files.get(), self)
        else:
            _add(files, self)
        logger.debug(f"done: {str(self)}")
        return self

    def bind(self, contents: Systems):
        if isinstance(contents, Systems):
            self.contents = contents
        else:
            logger.error("bind method accepts only Systems instance")
            raise TypeError
        return self

    def show(self):
        self.contents.show()
        return self

    def labeling(self, separator: str = "_", index_list: list[int] = [0, 1]):
        for c in self.contents:
            sp_stem = c.name.split(separator)
            c.label = separator.join([sp_stem[i] for i in index_list if i < len(sp_stem)])
        logger.debug(f"done: {str(self)}")
        return self

    def set_data(self, key: str, value=None):
        for c in self.contents:
            c.data[key] = value
        logger.debug(f"done: {str(self)}")
        return self

    def set_state(self, state: bool = True):
        for c in self.contents:
            c.state = state
        logger.debug(f"done: {str(self)}")
        return self

    def set_label(self, label: str = ""):
        for c in self.contents:
            c.label = label
        logger.debug(f"done: {str(self)}")
        return self

    def zero_fill(self, digit: int = 3, sepalator: str = "_", position: int = 3):
        if int(digit) <= 0:
            logger.error("invalid digit")
            return self
        for c in self.get():
            word_list = c.name.split(sepalator)
            for i, _n in enumerate(word_list):
                if i == position - 1:
                    word_list[i] = _n.zfill(digit)
            c.name = sepalator.join(word_list)
        logger.debug(f"done: {str(self)}")
        return self

    def count(self, comment: str = ""):
        labels = {}
        for label in self.contents.labels.keys():
            labels[label] = len(self.get().labels.get(label, []))
        t = ""
        for label, i in labels.items():
            t += f"{label}: {i}, "
        t += f"all: {sum(list(labels.values()))} conformers"
        if comment != "":
            t = f"{comment}: {t}"
        logger.info(t)
        return self

    def duplicate(self):
        n = self.__class__()
        n.contents = self.contents.duplicate(n)
        n.data = self.data.duplicate(n)
        return n

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
            "distribution",
            "charge",
            "multiplicity",
        ]
        for _c in self.contents:
            data_dict = {}
            data_dict["name"] = _c.name
            data_dict["path"] = _c.path
            data_dict["filetype"] = _c.filetype
            data_dict["state"] = _c.state
            data_dict["history"] = _c.history
            data_dict["label"] = _c.label
            data_dict["energy"] = _c.energy
            data_dict["distribution"] = _c.distribution
            data_dict["charge"] = _c.charge
            data_dict["multiplicity"] = _c.multiplicity
            for _key, _val in _c.data._data.items():
                str_val = str(_val)
                if len(str_val) >= 256:
                    data_dict[_key] = "######"
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
            if len(self.get()) != len(self.get().has_label()):
                logger.info("labeling called automatically")
                self.labeling()
            for cs in self.get().labels.values():
                min_e = min(c.energy for c in cs)
                logger.info(f"the minimun energy: {min_e}")
                for c in cs:
                    if (c.energy - min_e) >= threshold:
                        c.deactivate("energy_limit")
            if int(max_limit) > 0:
                for cs in self.get().labels.values():
                    sorted_cs: list[System] = sorted(cs, key=lambda t: t.energy)
                    for rank, c in enumerate(sorted_cs, start=1):
                        if rank > max_limit:
                            logger.info(f"{c.name}: {rank}th stable (max_limit={max_limit})")
                            c.deactivate("max limit in energy_limit")
        else:
            cs = self.get()
            min_e = min(c.energy for c in cs)
            logger.info(f"the minimun energy: {min_e}")
            for c in cs:
                if (c.energy - min_e) >= threshold:
                    c.deactivate("energy_limit")
            if int(max_limit) > 0:
                sorted_cs: list[System] = sorted(cs, key=lambda t: t.energy)
                for rank, c in enumerate(sorted_cs, start=1):
                    if rank > max_limit:
                        logger.info(f"{c.name}: {rank}th stable (max_limit={max_limit})")
                        c.deactivate("max limit in energy_limit")
        logger.debug(f"done: {str(self)}")
        return self

    def calc_rel_energy(self, in_label: bool = True):
        if in_label:
            if len(self.get()) != len(self.get().has_label()):
                logger.info("labeling called automatically")
                self.labeling()
            for cs in self.get().labels.values():
                min_e = min(c.energy for c in cs)
                logger.info(f"the minimun energy: {min_e}")
                for c in cs:
                    c.energy = c.energy - min_e
        else:
            min_e = min(c.energy for c in self.get())
            logger.info(f"the minimun energy: {min_e}")
            for c in self.get():
                c.energy = c.energy - min_e
        logger.debug(f"done: {str(self)}")
        return self

    def calc_distribution(self, in_label: bool = True, temperature: float = 298.15):
        def _calc_distr(confs: list[System]):
            if len(confs) == 0:
                logger.error("there are no conformers in the Confs")
                return False
            min_e = min(c.energy for c in confs)
            sum_bfac = 0.0
            for c in confs:
                c.cache["relative_energy"] = c.energy - min_e
                c.cache["distr_factor"] = math.exp(
                    -1
                    * (c.cache["relative_energy"] * CONSTANTS.cal_to_J * 1000)
                    / (temperature * CONSTANTS.gas_constant_in_J_over_mol_K)
                )
                sum_bfac += c.cache["distr_factor"]
            for c in confs:
                c.distribution = c.cache["distr_factor"] / sum_bfac

        if in_label:
            if len(self.get()) != len(self.get().has_label()):
                logger.info("labeling called automatically")
                self.labeling()
            for cs in self.get().labels.values():
                _calc_distr(cs)
        else:
            _calc_distr(self.get())
        logger.debug(f"done: {str(self)}")
        return self

    def calc_energy(self, keys: list[str] = [], unit: Unit = Units.kcal_mol):
        for c in self.get():
            _energy = 0.0
            for data_key in keys:
                _energy += c.data[data_key]
            c.energy = unit(_energy).to_kcal_mol
        logger.debug(f"done: {str(self)}")
        return self

    def copy_files(self, directory: Path, change_path: bool = False, suffix: str = None):
        for c in self.get():
            p = change_dir(c.path, directory)
            p = p.with_name(c.name).with_suffix(c.path.suffix)
            if suffix is not None:
                p = p.with_suffix(suffix)
            shutil.copy2(c.path, p)
            logger.info(f"{c.path} was copied to {p}")
            if change_path:
                c.path = p
        logger.debug(f"done: {str(self)}")
        return self

    def search(self, directory: Path = None, existing_check: bool = True, suffix: str = None):
        for c in self.get():
            if directory is not None and suffix is not None:
                p = Path(directory).joinpath(c.name).with_suffix(suffix)
            elif directory is None and suffix is not None:
                p = c.path.with_name(c.name).with_suffix(suffix)
            elif directory is not None and suffix is None:
                p = Path(directory).joinpath(c.name).with_suffix(c.path.suffix)
            if existing_check:
                if p.exists():
                    c.path = p
                else:
                    logger.error(f"{p.name} does not exist")
                    c.deactivate("search_dir")
            else:
                c.path = p
        logger.debug(f"done: {str(self)}")
        return self

    @Selectors.read_atoms.add("format/xyz")
    def read_xyz(self):
        for c in self.get():
            formats.read_xyz(c)
        logger.debug(f"done: {str(self)}")
        return self

    def write_xyz(self, directory: Path = None, link: bool = True, centering: bool = True):
        for c in self.get():
            formats.write_xyz(c, directory, link, centering)
        Log.set_output_dir(directory)
        logger.debug(f"done: {str(self)}")
        return self

    @Selectors.read_atoms.add("format/mol")
    def read_mol(self):
        for c in self.get():
            formats.read_mol(c)
        logger.debug(f"done: {str(self)}")
        return self

    def write_mol(self, directory: Path = None, link: bool = True, centering: bool = True):
        if len(self.get()) != len(self.get().has_bonds()):
            self.calc_bonds()
        for c in self.get():
            formats.write_mol(c, directory, link, centering)
        Log.set_output_dir(directory)
        logger.debug(f"done: {str(self)}")
        return self

    def write_input(self, template: Path, directory=None, link: bool = True, arg: dict[str, str] = None):
        for c in self.get():
            formats.write_input(c=c, template=template, odir=directory, link=link, arg=arg)
        Log.set_output_dir(directory)
        logger.debug(f"done: {str(self)}")
        return self

    def calc_bonds(self, cov_scaling: float = 1.1, vdw_scaling: float = 1.0, aromatize=True):
        for c in self.get():
            logger.debug(f"{c.name}: calculating bonds")
            c.modeler.calc_bonds(cov_scaling, vdw_scaling)
            if aromatize:
                c.modeler.aromatize()
            else:
                logger.debug(f"{c.name}: skipped aromatic bonds detection")
        logger.debug(f"done: {str(self)}")
        return self

    def calc_symm(self, calc_all: bool = True):
        if not calc_all:
            cfs = []
            if len(self.get()) != len(self.get().has_label()):
                logger.info("labeling called automatically")
                self.labeling()
            for _confs in self.get().labels.values():
                cfs.append(_confs.get())
                if len(_confs) != 1:
                    logger.info(
                        f"symmetry information of {_confs.get().name} will be duplicated from initial confomer to all"
                    )
            cfs = Systems().bind(cfs)
        else:
            cfs = self.get()

        if len(cfs) != len(cfs.has_bonds()):
            logger.info("calc_bonds called automatically")
            self.calc_bonds()

        for c in cfs:
            mats: list[Matrix] = c.modeler.get_symmetry_matrices()
            c.data["rotamer"] = [mat for mat in mats if mat.data["type"] == "rotamer"]
            c.data["numisomer"] = [mat for mat in mats if mat.data["type"] == "numisomer"]
            c.data["has_symm"] = True

        if not calc_all:
            rot_mats_dict: dict[str, list[Matrix]] = {_c.label: _c.data["rotamer"] for _c in cfs}
            num_mats_dict: dict[str, list[Matrix]] = {_c.label: _c.data["numisomer"] for _c in cfs}
            for _label, _confs in self.get().labels.items():
                for c in _confs:
                    atoms_list = c.atoms.to_list()
                    c.data["rotamer"] = [Matrix(atoms_list).bind(_m._matrix) for _m in rot_mats_dict[_label]]
                    c.data["numisomer"] = [Matrix(atoms_list).bind(_m._matrix) for _m in num_mats_dict[_label]]
                    c.data["has_symm"] = True
        logger.debug(f"done: {str(self)}")
        return self

    def rmsd_limit(
        self,
        threshold: float = 0.01,
        all_combinations_of_confs: bool = False,
        redundant_check: int = 3,
        all_perturbation_of_rotamers: bool = False,
    ):
        if len(self.get()) != len(self.get().has_data("has_symm", True)):
            logger.info("calc_symm called automatically")
            self.calc_symm(calc_all=False)
        topology.rmsdpruning(
            self.get(),
            threshold,
            all_combinations_of_confs,
            redundant_check,
            all_perturbation_of_rotamers,
        )
        logger.debug(f"done: {str(self)}")
        return self

    def map_numbers(self, reference_box: "BoxCore" = None):
        if reference_box is None:
            reference_box = BoxCore([cs.get() for cs in self.get().labels.values()])
        rbox = reference_box.duplicate()
        rbox.calc_bonds()
        rbox.calc_symm()
        topology.map_numbers(self.get(), rbox.contents)
        logger.debug(f"done: {str(self)}")
        return self

    def calc_length(self, number_a: int, number_b: int, key: str = ""):
        number_a = int(number_a)
        number_b = int(number_b)
        for c in self.get():
            if key == "" or not isinstance(key, str):
                actual_key = "length_{}{}-{}{}".format(
                    c.atoms.get(number_a).symbol,
                    str(number_a),
                    c.atoms.get(number_b).symbol,
                    str(number_b),
                )
            else:
                actual_key = key
            c.data[actual_key] = c.atoms.get_length(number_a, number_b)
        logger.debug(f"done: {str(self)}")
        return self

    def calc_angle(self, number_a: int, number_b: int, number_c: int, key: str = "", radian: bool = False):
        number_a = int(number_a)
        number_b = int(number_b)
        number_c = int(number_c)
        for c in self.get():
            if key == "" or not isinstance(key, str):
                actual_key = "angle_{}{}-{}{}-{}{}".format(
                    c.atoms.get(number_a).symbol,
                    str(number_a),
                    c.atoms.get(number_b).symbol,
                    str(number_b),
                    c.atoms.get(number_c).symbol,
                    str(number_c),
                )
            else:
                actual_key = key
            c.data[actual_key] = c.atoms.get_angle(number_a, number_b, number_c, radian=radian)
        logger.debug(f"done: {str(self)}")
        return self

    def calc_dihedral(
        self, number_a: int, number_b: int, number_c: int, number_d: int, key: str = "", radian: bool = False
    ):
        number_a = int(number_a)
        number_b = int(number_b)
        number_c = int(number_c)
        number_d = int(number_d)
        for c in self.get():
            if key == "" or not isinstance(key, str):
                actual_key = "dihedral_{}{}-{}{}-{}{}-{}{}".format(
                    c.atoms.get(number_a).symbol,
                    str(number_a),
                    c.atoms.get(number_b).symbol,
                    str(number_b),
                    c.atoms.get(number_c).symbol,
                    str(number_c),
                    c.atoms.get(number_d).symbol,
                    str(number_d),
                )
            else:
                actual_key = key
            c.data[actual_key] = c.atoms.get_dihedral(number_a, number_b, number_c, number_d, radian=radian)
        logger.debug(f"done: {str(self)}")
        return self

    def get_average(self, keys: list[str] = [], keys_for_atoms: list[str] = []) -> Systems:
        new_ls: list[System] = []
        if len(self.get()) != len(self.get().has_distribution()):
            logger.info("all of active conformers do not have distribution data")
            logger.info("calc_distr called automatically")
            self.calc_distribution()
        for label, ref_cs in self.get().labels.items():
            c = System()
            _ref_one = ref_cs.get()
            new_ls.append(c)
            c.name = label
            c.energy = _ref_one.energy
            c.label = label
            for ref_a in _ref_one.atoms:
                c.atoms.get_new_atom().axyz = ref_a.axyz
            for ref_c in ref_cs:
                if len(ref_c.atoms) != len(c.atoms):
                    logger.error("atomic symbol does not matched")
                    raise Exception
                for ref_a in ref_c.atoms:
                    if c.atoms.get(ref_a.number).symbol != ref_a.symbol:
                        logger.error("atomic symbol does not matched")
                        raise Exception
        for key in keys:
            if len(self.get()) != len(self.get().has_data(key)):
                logger.info(f"all of active conformers do not have {key} data")
                continue
            for c in new_ls:
                for _oc in self.get().labels.get(c.name):
                    c.data[key] += _oc.data[key] * _oc.distribution
        for key in keys_for_atoms:
            for c in new_ls:
                for _a in c.atoms:
                    sum_val = 0.0
                    for _oc in self.get().labels.get(c.name):
                        tmp_val = _oc.atoms.get(_a.number).data.get(key)
                        if tmp_val is None:
                            logger.error(f"atom {_oc.atoms.get(_a.number)} of {_oc} have no {key} data")
                            break
                        sum_val += float(tmp_val) * _oc.distribution
                    else:
                        _a.data[key] = sum_val
        return Systems().bind(new_ls)

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
        for c in self.get():
            logger.info(
                f"{c.name}: moving {c.atoms.get(number_a)} and {c.atoms.get(number_b)} to make their distance {target}"
            )
            c.modeler.set_length(
                number_a, number_b, target, (bool(fix_a), bool(fix_b)), numbers_along_with_a, numbers_along_with_b
            )
        logger.debug(f"done: {str(self)}")
        return self

    def convert_to_mirror(self, centering: bool = True):
        for c in self.get():
            logger.info(f"{c.name}: converting to mirror image")
            c.modeler.mirroring(centering)
        logger.debug(f"done: {str(self)}")
        return self

    def only_minimum(self, in_label: bool = True):
        if in_label:
            if len(self.get()) != len(self.get().has_label()):
                logger.info("labeling called automatically")
                self.labeling()
            for cs in self.get().labels.values():
                sorted_cs_list: list[System] = sorted(cs, key=lambda c: c.energy)
                for c in sorted_cs_list[1:]:
                    c.deactivate("only_minimum")
        else:
            sorted_cs_list: list[System] = sorted(self.get(), key=lambda c: c.energy)
            for c in sorted_cs_list[1:]:
                c.deactivate("only_minimum")
        logger.debug(f"done: {str(self)}")
        return self
