from collections import defaultdict
from collections.abc import MutableSequence
from copy import deepcopy
from pathlib import Path
from typing import Dict, Iterator, List, Union

from accel.base.atoms import Atoms
from accel.util import FileType
from accel.util.datadict import Data
from accel.util.log import Log, logger


class Mol:
    __slots__ = [
        "path",
        "name",
        "filetype",
        "label",
        "flag",
        "defect",
        "energy",
        "atoms",
        "data",
        "cache",
        "total_charge",
        "multiplicity",
        "_initialized",
    ]

    def __init__(self, file_path=None):
        self._initialized: bool = False
        self.name: str = ""
        self.path: Path = None
        self.filetype: str = None
        if file_path is not None:
            self.path: Path = Path(file_path).absolute().resolve()
            self.name = self.path.stem
            self.filetype = FileType.analyse(self.path)
        self.label: str = ""
        self.flag: bool = True
        self.defect: str = ""
        self.energy: float = None
        self.atoms: Atoms = Atoms()
        self.data: Data = Data(self)
        self.cache: Dict = {}
        self.total_charge: int = None
        self.multiplicity: int = 1
        self._initialized = True

        if self.path is not None and not self.path.exists():
            logger.error(f"{self.path} does not exist")
            self.deactivate("file not found")
        elif self.path is not None:
            Log.input_dir = self.path.parent

    def __str__(self):
        return self.name

    def __setattr__(self, key, value):
        if key == "_initialized":
            object.__setattr__(self, key, value)
        elif not hasattr(self, "_initialized"):
            pass
        elif not self._initialized:
            pass
        elif key == "path":
            if value is not None:
                value = Path(value).absolute()
                if self.name == "":
                    self.name = value.stem
                if not value.exists():
                    logger.error(f"{value} does not exist")
                else:
                    self.filetype = FileType.analyse(value)
                logger.debug(f"{self}: set: {str(value)}")
        elif key == "name":
            value = str(value)
            logger.debug(f"{self}: name: {value}")
        elif key == "flag":
            value = bool(value)
            logger.debug(f"{self}: flag: {value}")
        elif key == "label":
            value = str(value)
            logger.debug(f"{self}: label: {value}")
        elif key == "energy":
            value = float(value)
            logger.debug(f"{self}: energy: {value}")
        elif key == "total_charge":
            value = int(value)
            logger.debug(f"{self}: total_charge: {value}")
        elif key == "multiplicity":
            value = int(value)
            logger.debug(f"{self}: multiplicity: {value}")
        object.__setattr__(self, key, value)

    @property
    def charge(self):
        if self.total_charge is not None:
            return self.total_charge
        else:
            return sum(_a.charge for _a in self.atoms)

    @charge.setter
    def charge_setter(self, value):
        self.total_charge = value

    def deactivate(self, cause: str = ""):
        logger.info(f"{self.name}: {cause}")
        if self.defect == "":
            self.defect = cause
        else:
            self.defect += f"; {cause}"
        self._initialized = False
        self.flag = False
        self._initialized = True

    def show(self):
        data_dict = {}
        data_dict["Name"] = self.name
        try:
            _path = self.path.relative_to(Path.cwd())
        except ValueError:
            _path = self.path.absolute()
        data_dict["Path"] = _path
        data_dict["FileType"] = self.filetype
        if self.flag:
            data_dict["State"] = "ACTIVE"
        else:
            data_dict["State"] = "INACTIVE"
        data_dict["Defect"] = self.defect
        data_dict["Label"] = self.label
        if self.energy is None:
            data_dict["Energy"] = "----- kcal/mol"
        else:
            data_dict["Energy"] = f"{self.energy:>5.1f} kcal/mol"
        data_dict["Charge"] = self.charge
        data_dict["Multiplicity"] = self.multiplicity
        for _key, _val in self.data._data.items():
            data_dict[f"Data[{_key}]"] = _val
        for _key, _val in data_dict.items():
            logger.info(f"{self.name}: {_key}: {str(_val)}")

    def duplicate(self) -> "Mol":
        _n = Mol()
        _n._initialized = False
        _n.path = deepcopy(self.path)
        _n.name = self.name
        _n.label = self.label
        _n.flag = self.flag
        _n.energy = self.energy
        _n.atoms = self.atoms.duplicate()
        _n.data = self.data.duplicate(_n)
        _n.cache = deepcopy(self.cache)
        _n.total_charge = self.total_charge
        _n.multiplicity = self.multiplicity
        _n._initialized = True
        return _n


class Mols(MutableSequence):
    __slots__ = ["_list"]

    def __init__(self):
        self._list: List[Mol] = []

    def _new_confomer(self, value: Union[Mol, Path, str]) -> Mol:
        if isinstance(value, Mol):
            return value
        elif isinstance(value, Path) or isinstance(value, str):
            return Mol(file_path=value)
        else:
            raise ValueError

    def __str__(self):
        return f"Box having {len(self._list)} conformers"

    def __getitem__(self, index) -> Mol:
        return self._list[index]

    def __setitem__(self, index, value: Union[Mol, Path, str]):
        new_conformer = self._new_confomer(value)
        self._list[index] = new_conformer

    def __delitem__(self, index):
        del self._list[index]

    def __len__(self):
        return len(self._list)

    def __iter__(self) -> Iterator[Mol]:
        return super().__iter__()

    def __add__(self, other: "Mols"):
        return Mols().bind(self._list + other._list)

    def show(self):
        for _index, _c in enumerate(self._list, 1):
            _c: Mol = _c
            _is_active = "ACTIVE"
            if _c.flag is False:
                _is_active = "------"
            if _c.energy is None:
                _energy = "----- kcal/mol"
            else:
                _energy = f"{_c.energy:>5.1f} kcal/mol"
            try:
                _path = _c.path.relative_to(Path.cwd())
            except ValueError:
                _path = _c.path.absolute()
            logger.info(f"{_index:>4}: {_c.name:<18}: {_is_active:<7}: {_c.label:<12}: {_energy}: {str(_path)}")

    def insert(self, index: int, value: Union[Mol, Path, str]):
        new_conformer = self._new_confomer(value)
        self._list.insert(index, new_conformer)

    def swap(self, index_a: int, index_b: int) -> "Mols":
        self._list[index_a], self._list[index_b] = (
            self._list[index_b],
            self._list[index_a],
        )
        return self

    def to_list(self) -> List[Mol]:
        return [_c for _c in self._list]

    def bind(self, data_list: List[Mol]):
        if isinstance(data_list, list):
            self._list = data_list
        else:
            self._list = list(data_list)
        return self

    def duplicate(self, parent_obj=None):
        _n = Mols()
        _n._list = [_c.duplicate() for _c in self._list]
        return _n

    @property
    def labels(self) -> Dict[str, "Mols"]:
        label_dict = defaultdict(list)
        for _c in self._list:
            label_dict[_c.label].append(_c)
        sorted_dict = {}
        for _key in sorted(label_dict.keys()):
            sorted_dict[_key] = Mols().bind(label_dict[_key])
        return sorted_dict

    @property
    def filetypes(self) -> Dict[str, "Mols"]:
        filetype_dict = defaultdict(list)
        for _c in self._list:
            if _c.filetype is None:
                filetype_dict[""].append(_c)
            else:
                filetype_dict[_c.filetype].append(_c)
        sorted_dict = {}
        for _key in sorted(filetype_dict.keys()):
            sorted_dict[_key] = Mols().bind(filetype_dict[_key])
        return sorted_dict

    def has_state(self, state: Union[bool, None] = True) -> "Mols":
        if state is not None:
            return Mols().bind([_c for _c in self._list if _c.flag is state])
        else:
            return Mols().bind([_c for _c in self._list])

    def has_bonds(self, flag: bool = True) -> "Mols":
        return Mols().bind([_c for _c in self._list if _c.atoms.has_bonds() is flag])

    def has_energy(self, min_limit: float = None, max_limit: float = None) -> "Mols":
        _cs = [_c for _c in self._list if _c.energy is not None]
        if min_limit is None and max_limit is None:
            return Mols().bind(_cs)
        if len(_cs) == 0:
            return Mols().bind(_cs)
        if min_limit is None:
            min_limit = min(_c.energy for _c in _cs)
        else:
            min_limit = float(min_limit)
        if max_limit is None:
            max_limit = max(_c.energy for _c in _cs)
        else:
            max_limit = float(max_limit)
        _cs = [_c for _c in _cs if min_limit <= _c.energy and _c.energy <= max_limit]
        return Mols().bind(_cs)

    def has_label(self, label: str = None) -> "Mols":
        if label is None:
            return Mols().bind([_c for _c in self._list if _c.label != ""])
        else:
            _cs = self.labels.get(label)
            if _cs is None:
                return Mols()
            return _cs

    def has_data(self, key: str, value=None) -> "Mols":
        if value is None:
            _cs = [_c for _c in self._list if key in _c.data.keys()]
        else:
            _cs = [_c for _c in self._list if _c.data.get(key) is value]
        return Mols().bind(_cs)

    def has_filetype(self, filetype: str) -> "Mols":
        return Mols().bind([_c for _c in self._list if _c.filetype is filetype])

    def get(self, identifier=None) -> Mol:
        if len(self._list) == 0:
            return None
        elif identifier is None:
            if len(self.has_energy()) == len(self.has_state()):
                return sorted(self.to_list(), key=lambda _c: _c.energy)[0]
            else:
                return sorted(self.to_list(), key=lambda _c: _c.name)[0]
        elif isinstance(identifier, int):
            if 1 <= identifier and identifier <= len(self):
                return self[identifier - 1]
            else:
                return None
        elif isinstance(identifier, str):
            for _c in self._list:
                if identifier in _c.name:
                    return _c
            else:
                return None
        return None

    def sorted(self, key="energy"):
        if key == "energy":
            new_list = sorted(self._list, key=lambda _c: _c.energy)
        elif key == "name":
            new_list = sorted(self._list, key=lambda _c: _c.name)
        elif key == "label":
            new_list = sorted(self._list, key=lambda _c: _c.label)
        else:
            new_list = sorted(self._list, key=lambda _c: _c.data.get(key))
        return Mols().bind(new_list)
