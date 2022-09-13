from collections import defaultdict
from collections.abc import MutableSequence
from copy import deepcopy
from pathlib import Path
from typing import Dict, Iterable, Iterator, Union

from accel.base.atoms import Atoms
from accel.base.modeler import Modeler
from accel.util import FileType
from accel.util.datadict import Data
from accel.util.log import Log, logger


class System:
    __slots__ = [
        "path",
        "name",
        "filetype",
        "label",
        "state",
        "history",
        "energy",
        "atoms",
        "data",
        "cache",
        "total_charge",
        "multiplicity",
        "distribution",
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
        self.state: bool = True
        self.history: str = ""
        self.energy: float = None
        self.atoms: Atoms = Atoms()
        self.data: Data = Data(self)
        self.cache: Dict = {}
        self.total_charge: int = None
        self.multiplicity: int = 1
        self.distribution: float = None
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
                logger.debug(f"{self}: path: {str(value)}")
            else:
                logger.debug(f"{self}: path: {value}")
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
            if value is not None:
                value = float(value)
            logger.debug(f"{self}: energy: {value}")
        elif key == "total_charge":
            if value is not None:
                value = int(value)
            logger.debug(f"{self}: total_charge: {value}")
        elif key == "multiplicity":
            value = int(value)
            logger.debug(f"{self}: multiplicity: {value}")
        elif key == "distribution":
            if value is not None:
                value = float(value)
            logger.debug(f"{self}: distribution: {value}")
        object.__setattr__(self, key, value)

    @property
    def charge(self):
        if self.total_charge is not None:
            return self.total_charge
        else:
            try:
                return sum(a.charge for a in self.atoms)
            except TypeError:
                return None

    @charge.setter
    def charge(self, value):
        self.total_charge = value

    def deactivate(self, reason: str = ""):
        logger.info(f"{self.name}: {reason}")
        if self.history == "":
            self.history = reason
        else:
            self.history += f"; {reason}"
        self._initialized = False
        self.state = False
        self._initialized = True
        return self

    def show(self):
        data_dict = {}
        data_dict["name"] = self.name
        if self.path is not None:
            try:
                _path = self.path.relative_to(Path.cwd())
            except ValueError:
                _path = self.path.absolute()
        else:
            _path = "---"
        data_dict["path"] = _path
        data_dict["fileType"] = self.filetype
        if self.state:
            data_dict["state"] = "ACTIVE"
        else:
            data_dict["state"] = "INACTIVE"
        data_dict["history"] = self.history
        data_dict["label"] = self.label
        if self.energy is None:
            data_dict["energy"] = "----- kcal/mol"
        else:
            data_dict["energy"] = f"{self.energy:>5.1f} kcal/mol"
        data_dict["charge"] = self.charge
        data_dict["multiplicity"] = self.multiplicity
        for _key, _val in self.data._data.items():
            data_dict[f"data[{_key}]"] = _val
        for _key, _val in data_dict.items():
            logger.info(f"{self.name}: {_key}: {str(_val)}")
        return self

    @property
    def modeler(self) -> Modeler:
        return Modeler(self.atoms)

    def duplicate(self) -> "System":
        n = System()
        n._initialized = False
        n.path = deepcopy(self.path)
        n.name = self.name
        n.label = self.label
        n.state = self.state
        n.energy = self.energy
        n.atoms = self.atoms.duplicate()
        n.data = self.data.duplicate(n)
        n.cache = deepcopy(self.cache)
        n.total_charge = self.total_charge
        n.multiplicity = self.multiplicity
        n._initialized = True
        return n


class Systems(MutableSequence):
    __slots__ = ["_list"]

    def __init__(self, contents: Iterable[System] = None):
        self._list: list[System] = []
        if contents is None:
            pass
        elif isinstance(contents, Iterable):
            for c in contents:
                if isinstance(c, System):
                    self._list.append(c)
                elif isinstance(c, Path) or isinstance(c, str):
                    self._list.append(System(c))
                else:
                    logger.error("Systems accepts only Iterable[System-like]")
        else:
            logger.error("Systems accepts only Iterable[System-like]")
            raise ValueError

    def _new_confomer(self, value: Union[System, Path, str]) -> System:
        if isinstance(value, System):
            return value
        elif isinstance(value, Path) or isinstance(value, str):
            return System(file_path=value)
        else:
            raise ValueError

    def __str__(self):
        return f"Box having {len(self._list)} conformers"

    def __getitem__(self, index) -> System:
        return self._list[index]

    def __setitem__(self, index, value: Union[System, Path, str]):
        new_conformer = self._new_confomer(value)
        self._list[index] = new_conformer

    def __delitem__(self, index):
        del self._list[index]

    def __len__(self):
        return len(self._list)

    def __iter__(self) -> Iterator[System]:
        return super().__iter__()

    def __add__(self, other: "Systems"):
        return Systems().bind(self._list + other._list)

    def show(self):
        for idx, c in enumerate(self._list, 1):
            _is_active = "ACTIVE"
            if c.state is False:
                _is_active = "------"
            if c.energy is None:
                _energy = "----- kcal/mol"
            else:
                _energy = f"{c.energy:>5.1f} kcal/mol"
            if c.path is not None:
                try:
                    _path = c.path.relative_to(Path.cwd())
                except ValueError:
                    _path = c.path.absolute()
            else:
                _path = "---"
            logger.info(f"{idx:>4}: {c.name:<18}: {_is_active:<7}: {c.label:<12}: {_energy}: {str(_path)}")

    def insert(self, index: int, value: Union[System, Path, str]):
        new_conformer = self._new_confomer(value)
        self._list.insert(index, new_conformer)

    def swap(self, index_a: int, index_b: int) -> "Systems":
        self._list[index_a], self._list[index_b] = (
            self._list[index_b],
            self._list[index_a],
        )
        return self

    def to_list(self) -> list[System]:
        return [c for c in self._list]

    def bind(self, data_list: list[System]):
        if isinstance(data_list, list):
            self._list = data_list
        else:
            self._list = list(data_list)
        return self

    def duplicate(self, parent_obj=None) -> "Systems":
        n = Systems()
        n._list = [c.duplicate() for c in self._list]
        return n

    @property
    def labels(self) -> dict[str, "Systems"]:
        label_dict = defaultdict(list)
        for c in self._list:
            label_dict[c.label].append(c)
        sorted_dict = {}
        for key in sorted(label_dict.keys()):
            sorted_dict[key] = Systems().bind(label_dict[key])
        return sorted_dict

    @property
    def filetypes(self) -> dict[str, "Systems"]:
        filetype_dict = defaultdict(list)
        for c in self._list:
            if c.filetype is None:
                filetype_dict[""].append(c)
            else:
                filetype_dict[c.filetype].append(c)
        sorted_dict = {}
        for key in sorted(filetype_dict.keys()):
            sorted_dict[key] = Systems().bind(filetype_dict[key])
        return sorted_dict

    def has_state(self, state: Union[bool, None] = True) -> "Systems":
        if state is not None:
            return Systems().bind([_c for _c in self._list if _c.state is state])
        else:
            return Systems().bind([_c for _c in self._list])

    def has_bonds(self, flag: bool = True) -> "Systems":
        return Systems().bind([_c for _c in self._list if _c.atoms.has_bonds() is flag])

    def has_energy(self, min_limit: float = None, max_limit: float = None) -> "Systems":
        c_list = [c for c in self._list if c.energy is not None]
        if min_limit is None and max_limit is None:
            return Systems().bind(c_list)
        if len(c_list) == 0:
            return Systems().bind(c_list)
        if min_limit is None:
            min_limit = min(_c.energy for _c in c_list)
        else:
            min_limit = float(min_limit)
        if max_limit is None:
            max_limit = max(_c.energy for _c in c_list)
        else:
            max_limit = float(max_limit)
        c_list = [_c for _c in c_list if min_limit <= _c.energy and _c.energy <= max_limit]
        return Systems().bind(c_list)

    def has_label(self, label: str = None) -> "Systems":
        if label is None:
            return Systems().bind([c for c in self._list if c.label != ""])
        else:
            cs = self.labels.get(label)
            if cs is None:
                return Systems()
            return cs

    def has_distribution(self, min_limit: float = None, max_limit: float = None) -> "Systems":
        cs = [c for c in self._list if c.distribution is not None]
        if min_limit is None and max_limit is None:
            return Systems().bind(cs)
        if len(cs) == 0:
            return Systems().bind(cs)
        if min_limit is None:
            min_limit = min(_c.distribution for _c in cs)
        else:
            min_limit = float(min_limit)
        if max_limit is None:
            max_limit = max(_c.distribution for _c in cs)
        else:
            max_limit = float(max_limit)
        cs = [c for c in cs if min_limit <= c.distribution and c.distribution <= max_limit]
        return Systems().bind(cs)

    def has_data(self, key: str, value=None) -> "Systems":
        if value is None:
            c_list = [_c for _c in self._list if key in _c.data.keys()]
        else:
            c_list = [_c for _c in self._list if _c.data.get(key) == value]
        return Systems().bind(c_list)

    def has_filetype(self, filetype: str) -> "Systems":
        return Systems().bind([_c for _c in self._list if _c.filetype == filetype])

    def get(self, identifier=None) -> System:
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
            for c in self._list:
                if identifier in c.name:
                    return c
            else:
                return None
        return None

    def sorted(self, key="energy") -> "Systems":
        if key == "energy":
            new_list = sorted(self._list, key=lambda c: c.energy)
        elif key == "name":
            new_list = sorted(self._list, key=lambda c: c.name)
        elif key == "label":
            new_list = sorted(self._list, key=lambda c: c.label)
        else:
            new_list = sorted(self._list, key=lambda c: c.data.get(key))
        return Systems().bind(new_list)
