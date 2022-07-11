import copy
from collections import deque
from collections.abc import MutableSequence
from typing import Dict, Iterator, List, MutableMapping, Sequence, Set, Tuple

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
        "_atoms",
    ]

    def __init__(self, parent_atoms: "Atoms" = None):
        self._symbol: str = ""
        self.x: float = 0.0
        self.y: float = 0.0
        self.z: float = 0.0
        self.data = Data(self)
        self.cache = {}
        self.charge: int = None
        self._atoms: Atoms = parent_atoms

    def __str__(self):
        if self._atoms is None:
            return f"{self.symbol}({self.x}, {self.y}, {self.z})"
        else:
            return f"{self.number}{self.symbol}"

    @property
    def number(self):
        if self._atoms is None:
            logger.error(f"{self} has no parent atoms")
            return None
        return self._atoms._list.index(self) + 1

    def _get_bonding_atom(self, bond_types: List[int]):
        if self._atoms is None:
            logger.error(f"{self} has no parent atoms")
            return []
        if self._atoms.bonds is None:
            logger.error(f"{self._atoms} has no bonding informations")
            return []
        self_num = self.number
        bonds_list = [
            self._atoms.get(_index + 1)
            for _index, _b_type in enumerate(self._atoms.bonds._matrix[self_num - 1])
            if _b_type in bond_types
        ]
        return bonds_list

    @property
    def bonds(self) -> List["Atom"]:
        return self._get_bonding_atom(
            [BondType.undefined, BondType.single, BondType.double, BondType.triple, BondType.aromatic, BondType.any]
        )

    @property
    def single(self) -> List["Atom"]:
        return self._get_bonding_atom([BondType.single])

    @property
    def double(self) -> List["Atom"]:
        return self._get_bonding_atom([BondType.double])

    @property
    def triple(self) -> List["Atom"]:
        return self._get_bonding_atom([BondType.triple])

    @property
    def aromatic(self) -> List["Atom"]:
        return self._get_bonding_atom([BondType.aromatic])

    @property
    def contacts(self) -> List["Atom"]:
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
        self.y = self.y + float(vector[0])
        self.z = self.z + float(vector[0])
        return self

    def duplicate(self, parent_atoms=None) -> "Atom":
        _n = Atom(parent_atoms)
        _n._symbol = self._symbol
        _n.x = self.x
        _n.y = self.y
        _n.z = self.z
        _n.data = self.data.duplicate(_n)
        _n.cache = self.cache
        _n.charge = self.charge
        return _n

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

    def to_dict(self) -> Dict[Tuple[int], int]:
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


def _get_maps(
    atoms_a: "Atoms",
    atoms_b: "Atoms",
    known_pairs: List[Tuple[int]] = [],
    exclude_known: bool = False,
    terminal_first: bool = False,
) -> List[List[Tuple[Atom]]]:

    if len(known_pairs) == 0:
        known_check = False
        given_known_as = []
        given_known_bs = []
    else:
        known_check = True
        given_known_as = [atoms_a.get(_pr[0]) for _pr in known_pairs]
        given_known_bs = [atoms_b.get(_pr[1]) for _pr in known_pairs]
        if None in given_known_as or None in given_known_bs:
            logger.error(f"invalid known_pairs: {known_pairs}")
            return []

    def _should_exclude(atom_a: Atom, atom_b: Atom):
        if atom_a in given_known_as and atom_b in given_known_bs:
            if given_known_as.index(atom_a) == given_known_bs.index(atom_b):
                if exclude_known:
                    return True
            else:
                return True
        elif atom_a in given_known_as or atom_b in given_known_bs:
            return True
        return False

    initial_pairs: List[Tuple[Atom]] = []

    if terminal_first:
        for atom_a in atoms_a:
            if atom_a.symbol == "H":
                continue
            for atom_b in atoms_b:
                if atom_b.symbol == "H":
                    continue
                if atom_a.symbol != atom_b.symbol:
                    continue
                if known_check:
                    if _should_exclude(atom_a, atom_b):
                        continue
                if len([_a for _a in atom_a.bonds if _a.symbol != "H"]) > 1:
                    if len([_b for _b in atom_b.bonds if _b.symbol != "H"]) > 1:
                        continue
                initial_pairs.append((atom_a, atom_b))

    if len(initial_pairs) == 0:
        for atom_a in atoms_a:
            if atom_a.symbol == "H":
                continue
            for atom_b in atoms_b:
                if atom_b.symbol == "H":
                    continue
                if atom_a.symbol != atom_b.symbol:
                    continue
                if len([_a for _a in atom_a.bonds if _a.symbol == "H"]) != len(
                    [_b for _b in atom_b.bonds if _b.symbol == "H"]
                ):
                    continue
                if known_check:
                    if _should_exclude(atom_a, atom_b):
                        continue
                initial_pairs.append((atom_a, atom_b))

    if len(initial_pairs) == 0:
        for atom_a in atoms_a:
            if atom_a.symbol == "H":
                continue
            for atom_b in atoms_b:
                if atom_b.symbol == "H":
                    continue
                if atom_a.symbol != atom_b.symbol:
                    continue
                if known_check:
                    if _should_exclude(atom_a, atom_b):
                        continue
                initial_pairs.append((atom_a, atom_b))

    if len(initial_pairs) == 0:
        _appended_known_idxs = None
        for atom_a in atoms_a:
            for atom_b in atoms_b:
                if known_check:
                    if _should_exclude(atom_a, atom_b):
                        continue
                    _known_a_idxs = [given_known_as.index(_a) for _a in atom_a.bonds if _a in given_known_as]
                    _known_a_idxs = sorted(_known_a_idxs)
                    _known_b_idxs = [given_known_bs.index(_b) for _b in atom_b.bonds if _b in given_known_bs]
                    _known_b_idxs = sorted(_known_b_idxs)
                    if len(_known_a_idxs) == 0 or len(_known_b_idxs) == 0:
                        continue
                    if _known_a_idxs == _known_b_idxs:
                        if _appended_known_idxs is None or _appended_known_idxs == _known_a_idxs:
                            _appended_known_idxs = _known_a_idxs
                            initial_pairs.append((atom_a, atom_b))

    if len(initial_pairs) == 0:
        for atom_a in atoms_a:
            for atom_b in atoms_b:
                if known_check:
                    if _should_exclude(atom_a, atom_b):
                        continue
                if len(atom_a.bonds) == 0 and len(atom_b.bonds) == 0:
                    initial_pairs.append((atom_a, atom_b))

    if len(initial_pairs) == 0:
        for atom_a in atoms_a:
            for atom_b in atoms_b:
                if known_check:
                    if _should_exclude(atom_a, atom_b):
                        continue
                initial_pairs.append((atom_a, atom_b))

    def _is_proper_bonding(atom_a: Atom, atom_b: Atom, side_pairs: List[Tuple[Atom]]):
        atom_a_bonds = atom_a.bonds
        atom_b_bonds = atom_b.bonds
        # shoud swap here
        if len(atom_a_bonds) != len(atom_b_bonds):
            return True
        # and here
        if len([a for a in atom_a_bonds if a.symbol == "H"]) != len([a for a in atom_b_bonds if a.symbol == "H"]):
            return False
        a_side_pairs = [_ab[0] for _ab in side_pairs]
        b_side_pairs = [_ab[1] for _ab in side_pairs]
        a_index_list = sorted([a_side_pairs.index(_a) for _a in atom_a_bonds if _a in a_side_pairs])
        b_index_list = sorted([b_side_pairs.index(_a) for _a in atom_b_bonds if _a in b_side_pairs])
        if a_index_list != b_index_list:
            return False
        return True

    initial_chains: List[List[Tuple[Atom]]] = [[]]
    for initial_pair in initial_pairs:
        stack = [[initial_pair]]
        while stack:
            side_pairs: List[Tuple[Atom]] = stack.pop()
            if len(side_pairs) > len(initial_chains[-1]):
                initial_chains = [side_pairs]
                logger.debug(f"initial_chains updated: {[(_pr[0].number, _pr[1].number) for _pr in side_pairs]}")
            elif len(side_pairs) == len(initial_chains[-1]):
                initial_chains.append(side_pairs)
                logger.debug(f"initial_chains appended: {[(_pr[0].number, _pr[1].number) for _pr in side_pairs]}")
            a_list = [_pr[0] for _pr in side_pairs]
            b_list = [_pr[1] for _pr in side_pairs]
            for next_a in side_pairs[-1][0].bonds:
                if next_a.symbol == "H" or next_a in a_list:
                    continue
                for next_b in side_pairs[-1][1].bonds:
                    if next_b.symbol == "H" or next_b in b_list:
                        continue
                    if next_a.symbol != next_b.symbol:
                        continue
                    if known_check:
                        if _should_exclude(next_a, next_b):
                            continue
                    if not _is_proper_bonding(next_a, next_b, side_pairs):
                        continue
                    stack.append(side_pairs + [(next_a, next_b)])

    if len(initial_chains) == 1 and len(initial_chains[0]) == 0:
        return []

    canonical_initial_chains: Dict[List[set], Tuple[Tuple[Atom]]] = {}
    for chain in initial_chains:
        key = tuple(sorted([(_pr[0].number, _pr[1].number) for _pr in chain], key=lambda t: t[0]))
        canonical_initial_chains[key] = chain

    for chain in canonical_initial_chains.values():
        logger.debug(f"canonical_initial_chains: {[(_pr[0].number, _pr[1].number) for _pr in chain]}")

    min_invalid = None
    h_matched_chains: List[List[Tuple[Atom]]] = [[]]
    for chain in canonical_initial_chains.values():
        invalid_atoms = 0
        for pair in chain:
            hs_of_a = [_a for _a in pair[0].bonds if _a.symbol == "H"]
            hs_of_b = [_b for _b in pair[1].bonds if _b.symbol == "H"]
            if len(hs_of_a) != len(hs_of_b):
                invalid_atoms += 1
        if min_invalid is None or min_invalid > invalid_atoms:
            h_matched_chains = [chain]
            min_invalid = invalid_atoms
        elif min_invalid == invalid_atoms:
            h_matched_chains.append(chain)

    for chain in h_matched_chains:
        logger.debug(f"h_matched_chains: {[(_pr[0].number, _pr[1].number) for _pr in chain]}")

    def _det_ez(a: Atom, b: Atom, c: Atom, d: Atom) -> bool:
        _vab = np.array(a.xyz) - np.array(b.xyz)
        _vcb = np.array(c.xyz) - np.array(b.xyz)
        _vdc = np.array(d.xyz) - np.array(c.xyz)
        _pvac = np.cross(_vab, _vcb)
        _pvbd = np.cross(_vdc, _vcb)
        _angle = np.arccos(np.sum(_pvac * _pvbd) / (np.linalg.norm(_pvac) * np.linalg.norm(_pvbd)))
        if np.sum(_pvac * np.cross(_pvbd, _vcb)) < 0:
            _angle = -_angle
        _angle = float(np.rad2deg(_angle))
        if _angle > 90 or _angle < -90:
            return True
        return False

    def _det_chirality(a: Atom, b: Atom, c: Atom, d: Atom) -> bool:
        neighbors_xyz = np.array([_atom.xyz for _atom in [b, c, d]]) - np.array(a.xyz)
        if np.linalg.det(neighbors_xyz) > 0:
            return True
        return False

    extended_chains: List[List[Tuple[Atom]]] = []
    for chain in h_matched_chains:
        stack = deque([_pr for _pr in chain])
        assigned_as = given_known_as + [_pr[0] for _pr in chain]
        assigned_bs = given_known_bs + [_pr[1] for _pr in chain]

        def _new_assign(new_pair: Tuple[Atom], root_pair: Tuple[Atom]):
            assigned_as.append(new_pair[0])
            assigned_bs.append(new_pair[1])
            stack.append((new_pair[0], new_pair[1]))
            stack.appendleft(root_pair)

        def _get_large(atoms: List[Atom]) -> List[Atom]:
            trees = [[{_a}] for _a in atoms]
            for _ in range(16):
                for tree in trees:
                    tree.append(set())
                    for _a in tree[-2]:
                        tree[-1].update(_a.bonds)

                dup_atoms = trees[0][-1]
                for tree in trees:
                    dup_atoms = dup_atoms & tree[-1]
                for tree in trees:
                    tree[-1] = tree[-1] - dup_atoms
                weights = []
                for tree in trees:
                    weights.append(len(tree[-1]))
                if max(weights) == 0:
                    return [list(tree[0])[0] for tree in trees]
                alive_trees = []
                for tree_idx, weight in enumerate(weights):
                    if weight == max(weights):
                        alive_trees.append(trees[tree_idx])
                if len(alive_trees) == 1:
                    return list(alive_trees[0][0])
                trees = alive_trees
            else:
                return [list(tree[0])[0] for tree in alive_trees]

        while stack:
            root_pair: Tuple[Atom] = stack.popleft()
            next_as = [_a for _a in root_pair[0].bonds if _a not in assigned_as]
            next_bs = [_b for _b in root_pair[1].bonds if _b not in assigned_bs]
            if len(next_as) == 0 or len(next_bs) == 0:
                continue
            a_sbls = [_a.symbol for _a in next_as]
            b_sbls = [_b.symbol for _b in next_bs]

            mono_subs = [_s for _s in a_sbls if a_sbls.count(_s) == 1 and b_sbls.count(_s) == 1]
            if len(mono_subs) >= 1:
                _new_assign((next_as[a_sbls.index(mono_subs[0])], next_bs[b_sbls.index(mono_subs[0])]), root_pair)
                continue

            tri_subs = [_s for _s in a_sbls if a_sbls.count(_s) >= 3 and b_sbls.count(_s) >= 3]
            if len(tri_subs) >= 1:
                large_as = _get_large([_a for _a in next_as if _a.symbol == tri_subs[0]])
                large_bs = _get_large([_b for _b in next_bs if _b.symbol == tri_subs[0]])
                _new_assign((large_as[0], large_bs[0]), root_pair)
                continue

            di_subs = [_s for _s in a_sbls if min(a_sbls.count(_s), b_sbls.count(_s)) == 2]
            if len(di_subs) >= 1:
                large_as = _get_large([_a for _a in next_as if _a.symbol == di_subs[0]])
                large_bs = _get_large([_b for _b in next_bs if _b.symbol == di_subs[0]])
                if len(large_as) == 1 and len(large_bs) == 1:
                    _new_assign((large_as[0], large_bs[0]), root_pair)
                    continue
                if len(large_as) == 1 or len(large_bs) == 1:
                    continue
                known_as = [_a for _a in root_pair[0].bonds if _a in assigned_as]
                corre_bs = [assigned_bs[assigned_as.index(_a)] for _a in known_as]
                known_bs = [_b for _b in root_pair[1].bonds if _b in assigned_bs]
                if len(known_as) == 2 and len(known_bs) == 2:
                    flag_a = _det_chirality(known_as[0], known_as[1], large_as[0], large_as[1])
                    flag_b = _det_chirality(corre_bs[0], corre_bs[1], large_bs[0], large_bs[1])
                    if flag_a is flag_b:
                        _new_assign((large_as[0], large_bs[0]), root_pair)
                    else:
                        _new_assign((large_as[0], large_bs[1]), root_pair)
                    continue
                if len(known_as) == 1 and len(known_bs) == 1:
                    two_bond_pair = None
                    for _a in [_a for _a in known_as[0].bonds if _a in assigned_as and _a != root_pair[0]]:
                        for _b in [_b for _b in corre_bs[0].bonds if _b in assigned_bs and _b != root_pair[1]]:
                            if assigned_as.index(_a) == assigned_bs.index(_b):
                                two_bond_pair = (_a, _b)
                                break
                    if two_bond_pair is None:
                        continue
                    flag_a = _det_ez(two_bond_pair[0], known_as[0], root_pair[0], large_as[0])
                    flag_b = _det_ez(two_bond_pair[1], corre_bs[0], root_pair[1], large_bs[0])
                    if flag_a is flag_b:
                        _new_assign((large_as[0], large_bs[0]), root_pair)
                    else:
                        _new_assign((large_as[0], large_bs[1]), root_pair)
                    continue

        extended_chains.append([(_a, _b) for _a, _b in zip(assigned_as, assigned_bs)])

    canonical_extended_chains: Dict[Tuple[Tuple[int]], List[Tuple[Atom]]] = {}
    for chain in extended_chains:
        key = tuple(sorted([(_pr[0].number, _pr[1].number) for _pr in chain], key=lambda t: t[0]))
        canonical_extended_chains[key] = chain

    for chain in canonical_extended_chains.values():
        logger.debug(f"canonical_extended_chains: {[(_pr[0].number, _pr[1].number) for _pr in chain]}")

    recursive_extended_chains: List[List[Tuple[Atom]]] = []
    for chain in canonical_extended_chains.values():
        if len(chain) < min(len(atoms_a), len(atoms_b)):
            new_maps = _get_maps(atoms_a, atoms_b, known_pairs=chain, exclude_known=True, terminal_first=False)
            if len(new_maps) == 0:
                recursive_extended_chains.append(chain)
            else:
                for new_map in new_maps:
                    recursive_extended_chains.append(chain + new_map)
        else:
            recursive_extended_chains.append(chain)

    canonical_recursive_extended_chains: Dict[Tuple[Tuple[int]], List[Tuple[Atom]]] = {}
    for chain in recursive_extended_chains:
        if exclude_known:
            chain = chain[len(known_pairs) :]
        key = tuple(sorted([(_pr[0].number, _pr[1].number) for _pr in chain], key=lambda t: t[0]))
        canonical_recursive_extended_chains[key] = chain

    for chain in canonical_recursive_extended_chains.values():
        logger.debug(f"canonical_recursive_extended_chains: {[(_pr[0].number, _pr[1].number) for _pr in chain]}")

    return canonical_recursive_extended_chains.values()


def _order_maps(atoms_a: "Atoms", atoms_b: "Atoms", atom_maps: List[List[Tuple[Atom]]]) -> List[List[Tuple[Atom]]]:
    return atom_maps


class Atoms(MutableSequence):
    __slots__ = ["_list", "bonds"]

    def __init__(self) -> None:
        self._list: List[Atom] = []
        self.bonds: Bonds = None

    def _new_atom_instance(self, value) -> Atom:
        if isinstance(value, Atom):
            if value._atoms is None:
                new_atom = value
            else:
                new_atom = value.duplicate()
            new_atom._atoms = self
        elif isinstance(value, Sequence) and len(value) == 4:
            new_atom = Atom(self)
            new_atom.axyz = value
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

    def to_list(self) -> List[Atom]:
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
            new_atom = Atom(self)
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

    @property
    def mw(self) -> float:
        logger.error("not implemented yet: returns incorrect values")
        mol_weight = 0
        for _a in self._list:
            if _a.symbol == "H":
                mol_weight += 1
                continue
            mol_weight += Elements.get_element(_a.symbol)["number"] * 2
        return mol_weight

    def show(self):
        for _a in self._list:
            logger.info(f"{str(_a):>4}: {_a.x:>8.4f}: {_a.y:>7.4f}: {_a.z:>7.4f}")
        return self

    def bind(self, value: List):
        self._list = value
        return self

    def get_maps(
        self,
        target: "Atoms",
        known_pairs: List[Tuple[int]] = [],
        terminal_first: bool = False,
    ) -> List[List[int]]:
        atoms_map = _get_maps(target, self, known_pairs=known_pairs, terminal_first=terminal_first)
        atoms_map = _order_maps(target, self, atoms_map)
        for chain in atoms_map:
            logger.debug(f"maps (target, self): {[(_pr[0].number, _pr[1].number) for _pr in chain]}")
        return_chains = []
        for chain in atoms_map:
            chain_target = [_pr[0] for _pr in chain]
            chain_self = [_pr[1] for _pr in chain]
            return_nums = []
            for atom_self in self:
                try:
                    atom_target = chain_target[chain_self.index(atom_self)]
                    return_nums.append(atom_target.number)
                except ValueError:
                    return_nums.append(0)
            return_chains.append(return_nums)
        return return_chains

    def get_mapped(
        self,
        target: "Atoms",
        known_pairs: List[Tuple[int]] = [],
        terminal_first: bool = False,
    ) -> List["Atoms"]:
        atoms_map = _get_maps(target, self, known_pairs=known_pairs, terminal_first=terminal_first)
        atoms_map = _order_maps(target, self, atoms_map)
        for chain in atoms_map:
            logger.debug(f"maps (target, self): {[(_pr[0].number, _pr[1].number) for _pr in chain]}")
        return_list = []
        for chain in atoms_map:
            chain_target = [_pr[0] for _pr in chain]
            chain_self = [_pr[1] for _pr in chain]
            return_atoms = Atoms()
            mapped_numbers = []
            for atom_target in target:
                try:
                    atom_self = chain_self[chain_target.index(atom_target)]
                except ValueError:
                    logger.error("could not map atoms properly")
                    return []
                return_atoms.append(atom_self)
                mapped_numbers.append(atom_self.number)
            if self.has_bonds():
                return_atoms.init_bonds()
                for fromto, val in self.bonds.to_dict().items():
                    try:
                        _from = mapped_numbers.index(fromto[0]) + 1
                        _to = mapped_numbers.index(fromto[1]) + 1
                    except ValueError:
                        continue
                    return_atoms.bonds[_from, _to] = val
            return_list.append(return_atoms)
        return return_list
