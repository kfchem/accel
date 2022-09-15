import copy
import math
from collections.abc import MutableSequence
from typing import Iterable, Iterator, MutableMapping, Sequence, Set, Tuple

import numpy as np
from accel.util.constants import Elements
from accel.util.datadict import Data
from accel.util.log import logger


class BondType:
    undefined = -1
    none = 0
    single = 1
    double = 2
    triple = 3
    aromatic = 4
    single_or_double = 5
    single_or_aromatic = 6
    double_or_aromatic = 7
    any = 8
    contact = 9

    @classmethod
    def decode(cls, value):
        for _k, _v in cls.__dict__.items():
            if _v == value:
                return _k.upper()
        return None

    def encode(cls, value):
        for _k, _v in cls.__dict__.items():
            if _k == value:
                return _v
        return None


class Atom:
    __slots__ = [
        "_symbol",
        "x",
        "y",
        "z",
        "data",
        "cache",
        "charge",
        "parent",
    ]

    def __init__(self, symbol: str = None, x: float = 0.0, y: float = 0.0, z: float = 0.0, parent: "Atoms" = None):
        if symbol is None:
            self._symbol: str = ""
        else:
            self._symbol: str = Elements.canonicalize(symbol)
        self.x: float = float(x)
        self.y: float = float(y)
        self.z: float = float(z)
        self.data = Data(self)
        self.cache = {}
        self.charge: int = None
        self.parent: Atoms = parent

    def __str__(self):
        if self.parent is None:
            return f"{self.symbol}({self.x}, {self.y}, {self.z})"
        else:
            return f"{self.number}{self.symbol}"

    @property
    def number(self):
        if self.parent is None:
            logger.error(f"{self} has no parent atoms")
            return None
        return self.parent._list.index(self) + 1

    def _get_bonding_atom(self, bond_types: list[int]):
        if self.parent is None:
            logger.error(f"{self} has no parent atoms")
            return []
        if self.parent.bonds is None:
            logger.error(f"{self.parent} has no bonding informations")
            return []
        self_num = self.number
        bonds_list = [
            self.parent.get(_index + 1)
            for _index, _b_type in enumerate(self.parent.bonds._matrix[self_num - 1])
            if _b_type in bond_types
        ]
        return bonds_list

    @property
    def bonds(self) -> list["Atom"]:
        return self._get_bonding_atom(
            [BondType.undefined, BondType.single, BondType.double, BondType.triple, BondType.aromatic, BondType.any]
        )

    @property
    def single(self) -> list["Atom"]:
        return self._get_bonding_atom([BondType.single])

    @property
    def double(self) -> list["Atom"]:
        return self._get_bonding_atom([BondType.double])

    @property
    def triple(self) -> list["Atom"]:
        return self._get_bonding_atom([BondType.triple])

    @property
    def aromatic(self) -> list["Atom"]:
        return self._get_bonding_atom([BondType.aromatic])

    @property
    def contacts(self) -> list["Atom"]:
        return self._get_bonding_atom([BondType.contact])

    @property
    def symbol(self):
        return self._symbol

    @symbol.setter
    def symbol(self, value):
        self._symbol = Elements.canonicalize(value)

    @property
    def xyz(self):
        return [self.x, self.y, self.z]

    @xyz.setter
    def xyz(self, value: list):
        self.x = float(value[0])
        self.y = float(value[1])
        self.z = float(value[2])

    @property
    def axyz(self):
        return [self.symbol, self.x, self.y, self.z]

    @axyz.setter
    def axyz(self, value: list):
        self.symbol = value[0]
        self.x = float(value[1])
        self.y = float(value[2])
        self.z = float(value[3])

    def move(self, vector=(0.0, 0.0, 0.0)):
        self.x = self.x + float(vector[0])
        self.y = self.y + float(vector[1])
        self.z = self.z + float(vector[2])
        return self

    def duplicate(self, atoms=None) -> "Atom":
        n = Atom(parent=atoms)
        n._symbol = self._symbol
        n.x = self.x
        n.y = self.y
        n.z = self.z
        n.data = self.data.duplicate(n)
        n.cache = copy.deepcopy(self.cache)
        n.charge = self.charge
        return n

    def show(self):
        data_dict = {}
        data_dict["symbol"] = self.symbol
        data_dict["number"] = self.number
        data_dict["x"] = self.x
        data_dict["y"] = self.y
        data_dict["z"] = self.y
        data_dict["charge"] = self.charge
        data_dict["bonds"] = [str(_a) for _a in self.bonds]
        data_dict["single"] = [str(_a) for _a in self.single]
        data_dict["double"] = [str(_a) for _a in self.double]
        data_dict["triple"] = [str(_a) for _a in self.triple]
        data_dict["aromatic"] = [str(_a) for _a in self.aromatic]
        data_dict["contacts"] = [str(_a) for _a in self.contacts]
        for _key, _val in self.data._data.items():
            data_dict[f"Data[{_key}]"] = _val
        for _key, _val in data_dict.items():
            logger.info(f"{str(self)}: {_key}: {str(_val)}")
        return self


class Bonds(MutableMapping):
    __slots__ = ["_matrix"]

    def __init__(self, number_of_atoms: int):
        self._matrix: np.ndarray = np.zeros((number_of_atoms, number_of_atoms), dtype=int)

    def to_dict(self) -> dict[Tuple[int], int]:
        bonds_dict = {}
        for index_a in range(len(self._matrix)):
            for index_b in range(len(self._matrix)):
                if index_a >= index_b:
                    continue
                _type = self._matrix[index_a, index_b]
                if _type not in [BondType.none, BondType.undefined, BondType.contact]:
                    bonds_dict[(index_a + 1, index_b + 1)] = self._matrix[index_a, index_b]
        return bonds_dict

    def to_set(self) -> Set[Tuple[int]]:
        return {(_b[0], _b[1]) for _b in self.keys()}

    def __getitem__(self, key) -> int:
        if type(key) is tuple and len(key) == 2:
            number_a = int(key[0])
            number_b = int(key[1])
        else:
            raise KeyError
        return self._matrix[number_a - 1, number_b - 1]

    def __setitem__(self, key, value):
        if type(key) is tuple and len(key) == 2:
            number_a = int(key[0])
            number_b = int(key[1])
        else:
            raise KeyError
        self._matrix[number_a - 1, number_b - 1] = value
        self._matrix[number_b - 1, number_a - 1] = value

    def __delitem__(self, key):
        if type(key) is tuple and len(key) == 2:
            number_a = int(key[0])
            number_b = int(key[1])
        else:
            raise KeyError
        self._matrix[number_a - 1, number_b - 1] = BondType.none
        self._matrix[number_b - 1, number_a - 1] = BondType.none

    def __iter__(self):
        return self.to_dict().__iter__()

    def __len__(self):
        return len(self.to_dict())

    def swap(self, number_a: int, number_b: int):
        index_a = number_a - 1
        index_b = number_b - 1
        self._matrix[index_a], self._matrix[index_b] = self._matrix[index_b], self._matrix[index_a]
        for arr in self._matrix:
            arr[index_a], arr[index_b] = arr[index_b], arr[index_a]

    def insert(self, number: int, value: int = 0):
        index = number - 1
        self._matrix = np.insert(self._matrix, index, value, axis=0)
        self._matrix = np.insert(self._matrix, index, value, axis=1)

    def delete(self, number: int):
        index = number - 1
        self._matrix = np.delete(self._matrix, index, axis=0)
        self._matrix = np.delete(self._matrix, index, axis=1)

    def bind(self, ndarray: np.ndarray):
        shape = ndarray.shape
        if len(shape) == 2 and shape == self._matrix.shape:
            self._matrix = ndarray
            return self
        else:
            raise ValueError

    def show(self):
        for _b, _type in self.to_dict().items():
            logger.info(f"{_b[0]:>3} -{_b[1]:>3}: {BondType.decode(_type)}")
        return self


class Atoms(MutableSequence):
    __slots__ = ["_list", "bonds"]

    def __init__(self, contents: Iterable[Atom] = None) -> None:
        self._list: list[Atom] = []
        self.bonds: Bonds = None
        if contents is None:
            pass
        elif isinstance(contents, Iterable):
            for a in contents:
                if isinstance(a, Atom):
                    if a.parent is not None:
                        self._list.append(a.duplicate(self))
                        continue
                    self._list.append(a)
                elif isinstance(a, Sequence) and len(a) == 4:
                    self._list.append(Atom(*a, parent=self))
                else:
                    logger.error("Atoms accepts only Iterable[Atom-like]")
        else:
            logger.error("Atoms accepts only Iterable[Atom-like]")
            raise ValueError

    def _new_atom_instance(self, value) -> Atom:
        if isinstance(value, Atom):
            if value.parent is None:
                new_atom = value
            else:
                new_atom = value.duplicate()
            new_atom.parent = self
        elif isinstance(value, Sequence) and len(value) == 4:
            new_atom = Atom(*value, parent=self)
        else:
            raise ValueError
        return new_atom

    def __str__(self):
        return f"Atoms ({self.__len__} atoms)"

    def __getitem__(self, index):
        return self._list[index]

    def __setitem__(self, index, value):
        self._list[index] = self._new_atom_instance(value)

    def __delitem__(self, index):
        del self._list[index]

    def __len__(self):
        return len(self._list)

    def __iter__(self) -> Iterator[Atom]:
        return super().__iter__()

    def clear(self) -> None:
        self.bonds = None
        return self._list.clear()

    def insert(self, index: int, value) -> Atom:
        new_atom = self._new_atom_instance(value)
        self._list.insert(index, new_atom)
        if self.bonds is not None:
            self.bonds.insert(index + 1)
        return new_atom

    def swap(self, number_a: int, number_b: int, bonds: bool = True) -> "Atoms":
        self._list[number_a - 1], self._list[number_b - 1] = self._list[number_b - 1], self._list[number_a - 1]
        if bonds:
            self.bonds.swap(number_a, number_b)
        return self

    def to_list(self) -> list[Atom]:
        return [_a for _a in self._list]

    def get(self, value) -> Atom:
        if isinstance(value, Atom):
            if value in self._list:
                return value
            else:
                return None
        if isinstance(value, int):
            if 1 <= value and value <= len(self):
                return self._list[value - 1]
            else:
                return None
        return None

    def get_new_atom(self, value=None) -> Atom:
        if value is None:
            new_atom = Atom(parent=self)
        else:
            new_atom = self._new_atom_instance(value)
        self._list.append(new_atom)
        return new_atom

    def _copy_charges(self, ref_atoms: "Atoms"):
        for _a in self._list:
            _a.charge = ref_atoms.get(_a.number).charge

    def duplicate(self) -> "Atoms":
        _n = Atoms()
        _n._list = [_a.duplicate(_n) for _a in self._list]
        if self.bonds is not None:
            _n.init_bonds()
            _n.bonds.bind(copy.deepcopy(self.bonds._matrix))
        return _n

    def has_bonds(self) -> bool:
        if self.bonds is None:
            return False
        else:
            return True

    def init_bonds(self):
        self.bonds = Bonds(len(self._list))
        return self

    def add_bond(self, number_a: int, number_b: int, bond_type: int):
        if self.bonds is None:
            self.init_bonds()
        number_a = int(number_a)
        number_b = int(number_b)
        if number_a > number_b:
            number_a, number_b = number_b, number_a
        if number_a < 1 or len(self._list) < number_b:
            logger.error(f"could not add bond between {number_a} and {number_b}")
            raise ValueError
        self.bonds[number_a, number_b] = bond_type
        return self

    @property
    def mw(self) -> float:
        logger.error("not implemented yet: returns incorrect values")
        mol_weight = 0
        for a in self._list:
            if a.symbol == "H":
                mol_weight += 1
                continue
            mol_weight += Elements.get_element(a.symbol)["number"] * 2
        return mol_weight

    def show(self):
        for a in self._list:
            logger.info(f"{str(a):>4}: {a.x:>8.4f}: {a.y:>7.4f}: {a.z:>7.4f}")
        return self

    def bind(self, value: list):
        self._list = value
        return self

    def get_length(self, number_a: int, number_b: int) -> float:
        a_ = self.get(number_a).xyz
        b_ = self.get(number_b).xyz
        d_ = [float(a_[i]) - float(b_[i]) for i in range(3)]
        return math.sqrt(sum(x**2 for x in d_))

    def get_angle(self, number_a: int, number_b: int, number_c: int, radian=False) -> float:
        va = np.array(self.get(number_a).xyz)
        vb = np.array(self.get(number_b).xyz)
        vc = np.array(self.get(number_c).xyz)
        vba = vb - va
        vbc = vb - vc
        dba = np.linalg.norm(vba)
        dbc = np.linalg.norm(vbc)
        angle = np.arccos(np.sum(vba * vbc) / (dba * dbc))
        if not radian:
            angle = np.rad2deg(angle)
        return float(angle)

    def get_dihedral(self, number_a: int, number_b: int, number_c: int, number_d: int, radian=False) -> float:
        va = np.array(self.get(number_a).xyz)
        vb = np.array(self.get(number_b).xyz)
        vc = np.array(self.get(number_c).xyz)
        vd = np.array(self.get(number_d).xyz)
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
        if not radian:
            angle = np.rad2deg(angle)
        return float(angle)
