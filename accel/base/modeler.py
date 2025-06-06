import itertools
import math
from collections import defaultdict, deque
from collections.abc import MutableSequence
from statistics import mean
from typing import Iterable, Sequence

import numpy as np

from accel.base.atoms import Atom, Atoms, BondType
from accel.base.tools import float_to_str
from accel.util.constants import Elements
from accel.util.log import logger
from accel.util.matrix import Matrix


class Modeler:
    def __init__(self, atoms: Atoms = None) -> None:
        if atoms is None:
            self.atoms = Atoms()
        elif isinstance(atoms, Atoms):
            self.atoms = atoms
        elif isinstance(atoms, Iterable):
            atoms_list = []
            for a in atoms:
                if not isinstance(a, Atom):
                    logger.error("not allowed input for Modeler")
                    raise ValueError
                else:
                    atoms_list.append(a)
            self.atoms = Atoms()
            self.atoms.bind(atoms_list)
        else:
            logger.error("not allowed input for Modeler")
            raise ValueError

    def mirroring(self, centering=True):
        if centering:
            center = [mean([a.xyz[i] for a in self.atoms]) for i in range(3)]
            prec = max(max(len(str(a.xyz[i]).split(".")[1]) for a in self.atoms) for i in range(3))
            center = [round((-1) * v, prec) for v in center]
        for a in self.atoms:
            xyz = [(-1) * a.x, (-1) * a.y, (-1) * a.z]
            if centering:
                xyz = [float_to_str(round(_v - center[i], prec)) for i, _v in enumerate(xyz)]
            xyz = [float(float_to_str(_v)) for _v in xyz]
            a.x = xyz[0]
            a.y = xyz[1]
            a.z = xyz[2]
        return self

    def set_length(
        self,
        atom_a: int,
        atom_b: int,
        target_length: float,
        fixed_atom: tuple[bool, bool] = (False, False),
        move_along_with_a: Sequence[int] = (),
        move_along_with_b: Sequence[int] = (),
    ):
        vect = [0.0, 0.0, 0.0]
        for i in range(3):
            vect[i] = self.atoms.get(atom_b).xyz[i] - self.atoms.get(atom_a).xyz[i]
        distance = math.sqrt(sum(x**2 for x in vect))

        def _move_atoms(atoms_list, vect_factor):
            for atom_no in atoms_list:
                self.atoms.get(atom_no).xyz = [
                    _val + (vect_factor * vect[i] * (distance - target_length) / distance)
                    for i, _val in enumerate(self.atoms.get(atom_no).xyz)
                ]

        if fixed_atom == (False, False):
            _move_atoms([atom_a] + list(move_along_with_a), 0.5)
            _move_atoms([atom_b] + list(move_along_with_b), -0.5)
        elif fixed_atom == (False, True):
            _move_atoms([atom_a] + list(move_along_with_a), 1.0)
        elif fixed_atom == (True, False):
            _move_atoms([atom_b] + list(move_along_with_b), -1.0)
        else:
            raise ValueError
        return self

    def aromatize(self, single_threshold: int = 1.01, double_threshold: int = 1.01, max_depth: int = 18):
        a_list = self.atoms.to_list()
        npxyz = [[a.x, a.y, a.z] for a in a_list]
        npdist_mat = np.expand_dims(npxyz, axis=1) - np.expand_dims(npxyz, axis=0)
        npdist_mat = np.sqrt(np.sum(npdist_mat**2, axis=-1))

        npcov = [Elements.get_element(_a.symbol)["single"] for _a in a_list]
        npcov = [_v if _v is not None else np.nan for _v in npcov]
        npcov_mat = np.expand_dims(npcov, axis=1) + np.expand_dims(npcov, axis=0)
        npcov_mat = np.signbit(npdist_mat - single_threshold * npcov_mat)
        # single_threshold should be like 0.9 to reduce calculation cost

        npdbl = [Elements.get_element(_a.symbol)["double"] for _a in a_list]
        npdbl = [_v if _v is not None else np.nan for _v in npdbl]
        npdbl_mat = np.expand_dims(npdbl, axis=1) + np.expand_dims(npdbl, axis=0)
        npdbl_mat = np.signbit((double_threshold * npdbl_mat) - npdist_mat)
        # double_threshold should be like 1.1 to reduce calculation cost

        aromatic_mat: list[list[bool]] = (npcov_mat & npdbl_mat).tolist()

        def _path_to_other(start_idx: int, end_idx: int, route: list[int], ring_list: list[list[int]]):
            route = route + [start_idx]
            for other_idx in [_idx for _idx, _flag in enumerate(aromatic_mat[start_idx]) if _flag is True]:
                if other_idx == end_idx and other_idx != route[-2]:
                    ring_list.append(route)
                    continue
                if other_idx in route:
                    continue
                if len(route) > max_depth:
                    continue
                _path_to_other(other_idx, end_idx, route, ring_list)

        aromatic_atom = [False for _ in range(len(a_list))]
        aromatic_rings = []
        for cheking_idx in range(len(a_list)):
            if aromatic_atom[cheking_idx]:
                continue
            ring_list = []
            _path_to_other(cheking_idx, cheking_idx, [], ring_list)
            if len(ring_list) == 0:
                continue
            for ring in ring_list:
                electrons = 0
                for atom_idx in ring:
                    _symbol = self.atoms[atom_idx].symbol
                    if _symbol == "C":
                        electrons += 1
                    elif _symbol == "O":
                        electrons += 2
                    elif _symbol == "S":
                        electrons += 2
                    elif _symbol == "N":
                        if len(self.atoms[atom_idx].bonds) == 3:
                            electrons += 2
                        else:
                            electrons += 1
                    else:
                        electrons += 1
                if ((electrons - 2) % 4) != 0:
                    continue
                # check planer
                for _idx in ring:
                    aromatic_atom[_idx] = True
                aromatic_rings.append(ring)
        aromatic_bonds_set = set()
        for ring in aromatic_rings:
            for _id in range(len(ring) - 1):
                aromatic_bonds_set.add((ring[_id], ring[_id + 1]))
            aromatic_bonds_set.add((ring[0], ring[-1]))
        aromatic_bonds_set = {(id_a, id_b) for id_a, id_b in aromatic_bonds_set if id_a < id_b}
        for id_a, id_b in aromatic_bonds_set:
            self.atoms.bonds[id_a + 1, id_b + 1] = BondType.aromatic
        return self

    def calc_bonds(self, cov_scaling=1.1, vdw_scaling=1.0, double_scaling=1.05, triple_scaling=1.05):
        a_list = self.atoms.to_list()
        npxyz = [[_a.x, _a.y, _a.z] for _a in a_list]
        npdist_mat = np.expand_dims(npxyz, axis=1) - np.expand_dims(npxyz, axis=0)
        npdist_mat = np.sqrt(np.sum(npdist_mat**2, axis=-1))

        npcov = [Elements.get_element(_a.symbol)["single"] for _a in a_list]
        npcov = [_v if _v is not None else np.nan for _v in npcov]
        npcov_mat = np.expand_dims(npcov, axis=1) + np.expand_dims(npcov, axis=0)
        npcov_mat = np.signbit(npdist_mat - cov_scaling * npcov_mat)

        npvdw = [Elements.get_element(_a.symbol)["vdw"] for _a in a_list]
        npvdw = [_v if _v is not None else np.nan for _v in npvdw]
        npvdw_mat = np.expand_dims(npvdw, axis=1) + np.expand_dims(npvdw, axis=0)
        npvdw_mat = np.signbit(npdist_mat - vdw_scaling * npvdw_mat) ^ npcov_mat

        npdbl = [Elements.get_element(_a.symbol)["double"] for _a in a_list]
        npdbl = [_v if _v is not None else np.nan for _v in npdbl]
        npdbl_mat = np.expand_dims(npdbl, axis=1) + np.expand_dims(npdbl, axis=0)
        npdbl_mat = np.signbit(npdist_mat - double_scaling * npdbl_mat)

        nptri = [Elements.get_element(_a.symbol)["triple"] for _a in a_list]
        nptri = [_v if _v is not None else np.nan for _v in nptri]
        nptri_mat = np.expand_dims(nptri, axis=1) + np.expand_dims(nptri, axis=0)
        nptri_mat = np.signbit(npdist_mat - triple_scaling * nptri_mat)

        npsgl_mat = npcov_mat ^ npdbl_mat
        npdbl_mat = npdbl_mat ^ nptri_mat

        mat_dict: dict[int, np.ndarray] = {
            BondType.single: npsgl_mat,
            BondType.double: npdbl_mat,
            BondType.triple: nptri_mat,
            BondType.contact: npvdw_mat,
        }

        self.atoms.init_bonds()

        for _tyep, _mat in mat_dict.items():
            for index_a in range(len(a_list)):
                for index_b in range(len(a_list)):
                    if index_a >= index_b:
                        continue
                    if _mat[index_a, index_b]:
                        self.atoms.bonds[index_a + 1, index_b + 1] = _tyep
        return self

    def calc_stereo(self, with_num_chirality=True):
        logger.debug("calc_stereo is under developement, might be wrong assignment")

        def set_stereo(a: Atom, stereo: str):
            logger.debug(f"{a}: {stereo}")
            a.stereo.add(stereo)

        def get_ez(a: Atom, b: Atom, c: Atom, d: Atom) -> str:
            vab = np.array(a.xyz) - np.array(b.xyz)
            vcb = np.array(c.xyz) - np.array(b.xyz)
            vdc = np.array(d.xyz) - np.array(c.xyz)
            pvac = np.cross(vab, vcb)
            pvbd = np.cross(vdc, vcb)
            angle = np.arccos(np.sum(pvac * pvbd) / (np.linalg.norm(pvac) * np.linalg.norm(pvbd)))
            if np.sum(pvac * np.cross(pvbd, vcb)) < 0:
                angle = -angle
            angle = float(np.rad2deg(angle))
            if angle > 90 or angle < -90:
                return "E"
            return "Z"

        def get_rs(a: Atom, b: Atom, c: Atom, d: Atom) -> str:
            neighbors_xyz = np.array([atom_.xyz for atom_ in [a, b, c]]) - np.array(d.xyz)
            if np.linalg.det(neighbors_xyz) > 0:
                return "S"
            return "R"

        def det_allene_chirality(a: Atom, b: Atom, c: Atom):
            if len(b.single) >= 1 or len(b.double) != 2:
                return None
            a_singles = a.single
            c_singles = c.single
            if len(a_singles) != 2 or len(c_singles) != 2:
                return None
            a_subs = Chains.get_ordered_substituents([Chains([(a, a_single)]) for a_single in a_singles])
            c_subs = Chains.get_ordered_substituents([Chains([(a, c_single)]) for c_single in c_singles])
            if None in [a_subs, c_subs]:
                return None
            chirality = get_rs(*[chain[0][1] for chain in a_subs + c_subs])
            if chirality is None:
                return None
            else:
                return chirality + "a"

        def det_allene_num_chirality(a: Atom, b: Atom, c: Atom):
            if len(b.single) >= 1 or len(b.double) != 2:
                return None
            a_singles = a.single
            c_singles = c.single
            if len(a_singles) != 2 or len(c_singles) != 2:
                return None
            ordered_atoms: list[Atom] = sorted(a_singles, key=lambda a: a.number)
            ordered_atoms.extend(sorted(c_singles, key=lambda a: a.number))
            chirality = get_rs(*ordered_atoms)
            if chirality is None:
                return None
            else:
                return "n" + chirality + "a"

        for root_atom in self.atoms:
            root_atom.stereo = set()
        for root_atom in self.atoms:
            root_atom_bonds = root_atom.bonds
            if len(root_atom_bonds) == 2:
                db_bond = root_atom.double
                if len(db_bond) == 2:
                    allene_chirality = det_allene_chirality(db_bond[0], root_atom, db_bond[1])
                    if allene_chirality is not None:
                        set_stereo(db_bond[0], allene_chirality)
                        set_stereo(root_atom, allene_chirality)
                        set_stereo(db_bond[1], allene_chirality)
                if len(db_bond) == 2 and with_num_chirality:
                    allene_chirality = det_allene_num_chirality(db_bond[0], root_atom, db_bond[1])
                    if allene_chirality is not None:
                        set_stereo(db_bond[0], allene_chirality)
                        set_stereo(root_atom, allene_chirality)
                        set_stereo(db_bond[1], allene_chirality)
                if len(db_bond) != 1:
                    continue
                if len(root_atom.single) != 1:
                    continue
                if len(db_bond[0].single) not in (2, 1):
                    continue
                if with_num_chirality:
                    db_next = sorted(db_bond[0].single, key=lambda a: a.number, reverse=True)[0]
                    chiral = "n" + get_ez(root_atom.single[0], root_atom, db_bond[0], db_next) + "i"
                    set_stereo(root_atom, chiral)
                ba_singles = Chains.get_ordered_substituents([Chains([(db_bond[0], a)]) for a in db_bond[0].single])
                if ba_singles is None:
                    continue
                elif len(ba_singles) in (2, 1):
                    orderd_atoms = [root_atom.single[0], root_atom, db_bond[0], ba_singles[0][0][1]]
                    chiral = get_ez(*orderd_atoms) + "i"
                    set_stereo(root_atom, chiral)
                else:
                    logger.error(f"error chiral handling of ez {root_atom}")
                    root_atom.stereo.add("error")
            if len(root_atom_bonds) == 3:
                # aromatic handling
                db_bond = root_atom.double
                if len(db_bond) != 1:
                    continue
                if len(db_bond[0].single) not in (2, 1):
                    continue
                if with_num_chirality:
                    root_atom_single = root_atom.single
                    if len(root_atom_single) == 0:
                        logger.error(f"unexpected 0 length single bonds of {root_atom}")
                        continue
                    root_next = sorted(root_atom_single, key=lambda a: a.number, reverse=True)[0]
                    db_next = sorted(db_bond[0].single, key=lambda a: a.number, reverse=True)[0]
                    chiral = "n" + get_ez(root_next, root_atom, db_bond[0], db_next)
                    if len(db_bond[0].single) == 1:
                        chiral += "i"
                    set_stereo(root_atom, chiral)
                ra_singles = Chains.get_ordered_substituents([Chains([(root_atom, a)]) for a in root_atom.single])
                ba_singles = Chains.get_ordered_substituents([Chains([(db_bond[0], a)]) for a in db_bond[0].single])
                if None in [ra_singles, ba_singles]:
                    continue
                if len(ra_singles) == 2 and len(ba_singles) == 2:
                    orderd_atoms = [ra_singles[0][0][1], root_atom, db_bond[0], ba_singles[0][0][1]]
                    chiral = get_ez(*orderd_atoms)
                    set_stereo(root_atom, chiral)
                elif len(ra_singles) == 2 and len(ba_singles) == 1:
                    orderd_atoms = [ra_singles[0][0][1], root_atom, db_bond[0], ba_singles[0][0][1]]
                    chiral = get_ez(*orderd_atoms) + "i"
                    set_stereo(root_atom, chiral)
                else:
                    logger.error(f"error chiral handling of ez {root_atom}")
                    root_atom.stereo.add("error")
            if len(root_atom_bonds) == 4:
                if with_num_chirality:
                    chiral = "n" + get_rs(*sorted(root_atom_bonds, key=lambda a: a.number, reverse=True))
                    set_stereo(root_atom, chiral)
                ordered_subs = Chains.get_ordered_substituents([Chains([(root_atom, a)]) for a in root_atom_bonds])
                if ordered_subs is None:
                    continue
                if len(ordered_subs) == 4:
                    orderd_atoms: list[Atom] = [sub[0][1] for sub in ordered_subs]
                    chiral = get_rs(*orderd_atoms)
                    set_stereo(root_atom, chiral)
                else:
                    logger.error(f"error chiral handling of rs {root_atom}")
                    root_atom.stereo.add("error")

    def merge(self, atoms: Atoms):
        p_len = len(self.atoms)
        for a in atoms:
            self.atoms.append(a)
        for b_ft, b_val in atoms.bonds.to_dict().items():
            self.atoms.bonds[b_ft[0] + p_len, b_ft[1] + p_len] = b_val
        logger.debug(f"merged {len(atoms)} atoms")
        return self

    def incorporate(self, atoms: Atoms, padding: float = 2.0):
        padding = float(padding)
        t_atoms = atoms.duplicate()
        s_xyzs = [a.xyz for a in self.atoms]
        max_xyz = [max([xyz[i] for xyz in s_xyzs]) for i in range(3)]
        t_xyzs = [a.xyz for a in t_atoms]
        min_xyz = [min([xyz[i] for xyz in t_xyzs]) for i in range(3)]
        vec = [max_xyz[i] - min_xyz[i] + padding for i in range(3)]
        logger.debug(f"incoporate vector: {vec}")
        for a in t_atoms:
            a.move(vec)
        self.merge(t_atoms)
        return self

    def _map_to(self, target: Atoms):
        return self

    def get_splitted(self) -> list["Atoms"]:
        atom_ll: list[list[int]] = [[_a.number] for _a in self.atoms]
        bonding_types = [
            BondType.single,
            BondType.double,
            BondType.triple,
            BondType.aromatic,
            BondType.single_or_double,
            BondType.single_or_aromatic,
            BondType.double_or_aromatic,
        ]
        for _bond, _type in self.atoms.bonds.to_dict().items():
            if _type not in bonding_types:
                continue
            new_l = []
            for atom_l in atom_ll:
                if _bond[0] in atom_l or _bond[1] in atom_l:
                    new_l.extend(atom_l)
                    atom_l.clear()
            atom_ll.append(new_l)
        atom_ll: list[list[int]] = [sorted(list(atom_l)) for atom_l in atom_ll if atom_l != []]
        ret_list: list["Atoms"] = []
        for atom_nums in atom_ll:
            return_atoms = Atoms()
            for _num in atom_nums:
                return_atoms.append(self.atoms.get(_num))
            if self.atoms.has_bonds():
                return_atoms.init_bonds()
                for fromto, val in self.atoms.bonds.to_dict().items():
                    try:
                        _from = atom_nums.index(fromto[0]) + 1
                        _to = atom_nums.index(fromto[1]) + 1
                    except ValueError:
                        continue
                    return_atoms.bonds[_from, _to] = val
            ret_list.append(return_atoms)
        return ret_list

    def get_maps(
        self,
        target: "Atoms",
        known_pairs: list[tuple[int]] = [],
        terminal_first: bool = None,
        options: dict = {},
    ) -> list[list[int]]:
        atoms_map = _get_maps(target, self.atoms, known_pairs=known_pairs, terminal_first=terminal_first)
        atoms_map = _correct_maps(target, self.atoms, atoms_map)
        atoms_map = _order_maps(target, self.atoms, atoms_map, options=options)
        for chain in atoms_map:
            logger.debug(f"maps (target, self): {[(pr[0].number, pr[1].number) for pr in chain]}")
        return_chains = []
        for chain in atoms_map:
            chain_target = [pr[0] for pr in chain]
            chain_self = [pr[1] for pr in chain]
            return_nums = []
            for atom_self in self.atoms:
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
        known_pairs: list[tuple[int]] = [],
        terminal_first: bool = None,
        options: dict = {},
    ) -> list["Atoms"]:
        atoms_map = _get_maps(target, self.atoms, known_pairs=known_pairs, terminal_first=terminal_first)
        atoms_map = _correct_maps(target, self.atoms, atoms_map)
        atoms_map = _order_maps(target, self.atoms, atoms_map, options=options)
        for chain in atoms_map:
            logger.debug(f"maps (target, self): {[(pr[0].number, pr[1].number) for pr in chain]}")
        return_list = []
        for chain in atoms_map:
            chain_target = [pr[0] for pr in chain]
            chain_self = [pr[1] for pr in chain]
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
            if self.atoms.has_bonds():
                return_atoms.init_bonds()
                for fromto, val in self.atoms.bonds.to_dict().items():
                    try:
                        _from = mapped_numbers.index(fromto[0]) + 1
                        _to = mapped_numbers.index(fromto[1]) + 1
                    except ValueError:
                        continue
                    return_atoms.bonds[_from, _to] = val
            return_list.append(return_atoms)
        return return_list

    def get_symmetry_matrices(self) -> list[Matrix]:
        def _reset_visited_flags(atoms: Atoms):
            for a in atoms:
                a.cache["unvisited"] = True
                a.cache["ref_unvisited"] = True
            return None

        def _del_visited_flags(atoms: Atoms):
            for a in atoms:
                del a.cache["unvisited"]
                del a.cache["ref_unvisited"]
            return None

        _reset_visited_flags(self.atoms)
        for a in self.atoms:
            # generate tree by BFS
            atom_trees = []
            for root in a.bonds:
                a.cache["unvisited"] = False
                prop_tree = defaultdict(list)
                stack = []
                stack.append([root, 0])
                while not (stack == []):
                    poped_stack = stack.pop(0)
                    node: Atom = poped_stack[0]
                    distance: int = poped_stack[1]
                    node.cache["unvisited"] = False
                    prop_tree[distance].append(node)
                    if node.symbol == "H":
                        continue
                    for _nb in node.bonds:
                        if _nb.cache["unvisited"]:
                            stack.append([_nb, distance + 1])
                _reset_visited_flags(self.atoms)
                prop_list = []
                for distance in sorted(prop_tree.keys()):
                    prop = defaultdict(int)
                    for node in prop_tree[distance]:
                        prop["bonds"] += len(node.bonds)
                        prop[node.symbol] += 1
                    prop_list.append(dict(prop))
                atom_trees.append({"root_atom": root, "prop_list": prop_list})
            # logger.debug('BFS on {}'.format(_a))

            # generate candidates of pair of isomorphic tree
            candidate_pair = []
            for idx_ref, tree_ref in enumerate(atom_trees):
                for idx_target, tree_target in enumerate(atom_trees[idx_ref + 1 :], idx_ref + 1):
                    for prop_ref, prop_target in zip(tree_ref["prop_list"], tree_target["prop_list"]):
                        for _key in prop_ref:
                            if prop_ref[_key] != prop_target.get(_key):
                                break
                        else:
                            continue
                        break
                    else:
                        candidate_pair.append(
                            {"ref_root": tree_ref["root_atom"], "tar_root": tree_target["root_atom"]}
                        )
            # logger.debug('generated candidates of atom mapping on {}'.format(_a))

            # check candidates by DFS
            def _recursive_dfs(node_ref: Atom, node_tar: Atom, pair_atoms: list):
                if (node_ref.symbol != node_tar.symbol) or (len(node_ref.bonds) != len(node_tar.bonds)):
                    return False
                node_ref.cache["ref_unvisited"] = False
                node_tar.cache["unvisited"] = False
                if node_ref == node_tar:
                    return True
                for _nb_ref in [_n for _n in node_ref.bonds if _n.cache["ref_unvisited"]]:
                    if not _nb_ref.cache["ref_unvisited"]:
                        break
                    for _nb_tar in [_n for _n in node_tar.bonds if _n.cache["unvisited"]]:
                        if not _nb_tar.cache["unvisited"]:
                            break
                        if _recursive_dfs(_nb_ref, _nb_tar, pair_atoms):
                            break
                    else:
                        node_ref.cache["ref_unvisited"] = True
                        node_tar.cache["unvisited"] = True
                        return False
                pair_atoms.append([node_ref, node_tar])
                return True

            for _pair in candidate_pair:
                pair_atoms_list = []
                a.cache["unvisited"] = False
                a.cache["ref_unvisited"] = False
                if _recursive_dfs(_pair["ref_root"], _pair["tar_root"], pair_atoms_list):
                    _pair["pair_list"] = pair_atoms_list
                _reset_visited_flags(self.atoms)

            # logger.debug('checked candidates of atom mapping on {}'.format(_a))

            pairs_list = [_pair for _pair in candidate_pair if "pair_list" in _pair]
            for _ps in pairs_list:
                logger.debug(
                    "{}: root {} and {}: map {}".format(
                        a,
                        _ps["ref_root"],
                        _ps["tar_root"],
                        [[str(k) for k in i] for i in _ps["pair_list"]],
                    )
                )

            # generate matrix
            for _ps in pairs_list:
                _mat = np.identity(len(self.atoms))
                for ref_atom, tar_atom in _ps["pair_list"]:
                    _ref_idx = ref_atom.number - 1
                    _tar_idx = tar_atom.number - 1
                    _mat[[_ref_idx, _tar_idx]] = _mat[[_tar_idx, _ref_idx]]
                _ps["pair_matrix"] = _mat

            if len(pairs_list) != 0:
                a.cache["isomeric_subs_list"] = [
                    {
                        "root_a": _p["ref_root"],
                        "root_b": _p["tar_root"],
                        "list": _p["pair_list"],
                        "matrix": _p["pair_matrix"],
                    }
                    for _p in pairs_list
                ]

        rotamer_mat_list = []
        numisomer_mat_list = []
        for a in self.atoms:
            if "isomeric_subs_list" not in a.cache:
                continue
            # extract rotamers
            if [len(a.bonds), len(a.cache["isomeric_subs_list"])] in [[2, 1], [3, 1]]:
                # aromatic handling should be included
                if len(a.double) != 0:
                    for double_bond_atom in a.double:
                        if double_bond_atom in [
                            a.cache["isomeric_subs_list"][0]["root_a"],
                            a.cache["isomeric_subs_list"][0]["root_a"],
                        ]:
                            break
                    else:
                        numisomer_mat_list.append(a.cache["isomeric_subs_list"][0]["matrix"])
                        continue
                rot_mat = a.cache["isomeric_subs_list"][0]["matrix"]
                rot_mat = rot_mat.astype(int)
                rotamer_mat_list.append(rot_mat)
            elif [len(a.bonds), len(a.cache["isomeric_subs_list"])] in [[4, 3], [4, 6], [3, 3]]:
                for _ps in a.cache["isomeric_subs_list"][1:]:
                    rot_mat = np.dot(_ps["matrix"], a.cache["isomeric_subs_list"][0]["matrix"])
                    rot_mat = rot_mat.astype(int)
                    rotamer_mat_list.append(rot_mat)
                for _ps in a.cache["isomeric_subs_list"]:
                    numisomer_mat_list.append(_ps["matrix"])
                if [len(a.bonds), len(a.cache["isomeric_subs_list"])] in [[3, 3]]:
                    rotamer_mat_list.append(a.cache["isomeric_subs_list"][0]["matrix"])
            else:
                for _ps in a.cache["isomeric_subs_list"]:
                    numisomer_mat_list.append(_ps["matrix"])

        # generate identical collections of rot_mat
        rotamer_mat_list = list({_m.tobytes(): _m for _m in rotamer_mat_list}.values())
        numisomer_mat_list = list({_m.tobytes(): _m for _m in numisomer_mat_list}.values())

        # check chirality here
        def _invalid_chirality(atoms: Atoms, modified_xyz: np.ndarray) -> int:
            invalid_count = 0
            original_xyz = np.array([a.xyz for a in atoms])
            for _a in atoms:
                if len(_a.bonds) != 4:
                    continue
                _b_indexs = sorted(a.number - 1 for a in _a.bonds)
                original_matirix = original_xyz[_b_indexs[1:]] - original_xyz[_b_indexs[0]]
                modified_matirix = modified_xyz[_b_indexs[1:]] - modified_xyz[_b_indexs[0]]
                if (np.linalg.det(original_matirix) > 0) != (np.linalg.det(modified_matirix) > 0):
                    invalid_count += 1
            return invalid_count

        def _get_invalid_chiral_atoms(atoms: Atoms, modified_xyz: np.ndarray) -> list[Atom]:
            original_xyz = np.array([a.xyz for a in atoms])
            ret_list = []
            for _a in atoms:
                if len(_a.bonds) != 4:
                    continue
                _b_indexs = sorted(a.number - 1 for a in _a.bonds)
                original_matirix = original_xyz[_b_indexs[1:]] - original_xyz[_b_indexs[0]]
                modified_matirix = modified_xyz[_b_indexs[1:]] - modified_xyz[_b_indexs[0]]
                if (np.linalg.det(original_matirix) > 0) != (np.linalg.det(modified_matirix) > 0):
                    ret_list.append(_a)
            return ret_list

        # for TMS, cyclic_chiral_check has error
        def _cyclic_chiral_check_legacy(atoms: Atoms, rotamer_mat_list, numisomer_mat_list):
            unchanged_flag = True
            for rm_idx, rot_mat in enumerate(rotamer_mat_list):
                original_xyz = np.array([at.xyz for at in atoms])
                number_of_invalid = _invalid_chirality(atoms, np.dot(rot_mat, original_xyz))
                logger.debug(f"number of invalid chirality for rotemer matrix {rm_idx} is {number_of_invalid}")
                if number_of_invalid == 0:
                    continue
                for num_mat in reversed(numisomer_mat_list):
                    t_mat = np.dot(num_mat, rot_mat)
                    if _invalid_chirality(atoms, np.dot(t_mat, original_xyz)) < number_of_invalid:
                        rotamer_mat_list[rm_idx] = t_mat
                        unchanged_flag = False
                        logger.debug("chirality in rotamer matrix was partially corrected")
                        break
                else:
                    logger.debug("chirality in rotamer matrix temporarily corrected")
            return unchanged_flag

        def _cyclic_chiral_check(atoms: Atoms, rotamer_mat_list, numisomer_mat_list):
            unchanged_flag = True
            for rm_idx, rot_mat in enumerate(rotamer_mat_list):
                original_xyz = np.array([at.xyz for at in atoms])
                invalid_atoms = _get_invalid_chiral_atoms(atoms, np.dot(rot_mat, original_xyz))
                if len(invalid_atoms) == 0:
                    continue
                logger.debug(
                    f"invalid chiral atom is detected: {[str(ta) for ta in invalid_atoms]} in mat_idx {rm_idx}"
                )
                for isosub in invalid_atoms[0].cache["isomeric_subs_list"]:
                    t_mat = np.dot(isosub["matrix"], rot_mat)
                    iv_atoms = _get_invalid_chiral_atoms(atoms, np.dot(t_mat, original_xyz))
                    if len(iv_atoms) < len(invalid_atoms):
                        rotamer_mat_list[rm_idx] = t_mat
                        unchanged_flag = False
                        logger.debug("chirality in rotamer matrix was partially corrected")
                        break
                else:
                    logger.debug("chirality in rotamer matrix not corrected")
            return unchanged_flag

        if len(rotamer_mat_list) != 0:
            if not _cyclic_chiral_check(self.atoms, rotamer_mat_list, numisomer_mat_list):
                for _ in range(len(rotamer_mat_list) * len(numisomer_mat_list)):
                    if _cyclic_chiral_check(self.atoms, rotamer_mat_list, numisomer_mat_list):
                        logger.debug("chirality in rotamer matrix was successfully corrected")
                        break
                else:
                    for _ in range(len(rotamer_mat_list) * len(numisomer_mat_list)):
                        if _cyclic_chiral_check_legacy(self.atoms, rotamer_mat_list, numisomer_mat_list):
                            logger.debug("chirality in rotamer matrix was successfully corrected by legacy protocol")
                            break
                    else:
                        logger.error("chirality in rotamer matrix was not completely corrected")

        # generate identical collections of rot_mat again
        rotamer_mat_list = {_m.tobytes(): _m for _m in rotamer_mat_list}.values()
        numisomer_mat_list = {_m.tobytes(): _m for _m in numisomer_mat_list}.values()

        # added matrix property
        rotamer_mat_list = [Matrix(self.atoms.to_list()).bind(_m) for _m in rotamer_mat_list]
        numisomer_mat_list = [Matrix(self.atoms.to_list()).bind(_m) for _m in numisomer_mat_list]

        for mat in rotamer_mat_list:
            mat.data["type"] = "rotamer"
        for mat in numisomer_mat_list:
            mat.data["type"] = "numisomer"

        _del_visited_flags(self.atoms)

        return rotamer_mat_list + numisomer_mat_list

    def _superimpose(self, target: Atoms):
        return self

    def _get_rmsd(self, target: Atoms) -> float:
        return 0.0


def _get_maps(
    atoms_a: "Atoms",
    atoms_b: "Atoms",
    known_pairs: list[tuple[int]] = [],
    terminal_first: bool = None,
    _recursion_counter: int = 0,
) -> list[list[tuple[Atom]]]:
    initial_pairs: list[tuple[Atom]] = []

    if _recursion_counter != 0:
        logger.debug(f"recursively called: depth {_recursion_counter}")

    if terminal_first is None and len(atoms_a) > 256:
        logger.debug("terminal_first is automatically activated")
        terminal_first = True

    given_known_as = [atoms_a.get(pr[0]) for pr in known_pairs]
    given_known_bs = [atoms_b.get(pr[1]) for pr in known_pairs]
    if None in given_known_as or None in given_known_bs:
        logger.error(f"invalid known_pairs: {known_pairs}")
        return []

    ip_settings: list[list[bool]] = [
        [True, False, True, True, False, False],
        [False, True, True, True, False, False],
        [False, False, True, True, False, False],
        [False, False, False, True, True, False],
        [False, False, False, True, False, True],
        [False, False, False, True, False, False],
        [False, False, False, False, True, False],
        [False, False, False, False, False, True],
        [False, False, False, False, False, False],
    ]

    if len(known_pairs) == 0:
        ip_settings.pop(3)

    if _recursion_counter == 0 and len(known_pairs) != 0:
        ip_settings.pop(2)
        ip_settings.pop(1)

    if not terminal_first:
        ip_settings.pop(0)

    for ip_setting in ip_settings:
        if len(initial_pairs) != 0:
            break
        ip_terminal_first_check = ip_setting[0]
        ip_check_hs = ip_setting[1]
        ip_exclude_h = ip_setting[2]
        ip_check_symbol = ip_setting[3]
        ip_known_bond_check = ip_setting[4]
        ip_isolated_check = ip_setting[5]

        appended_known_idxs = None
        for atom_a in atoms_a:
            if atom_a in given_known_as:
                continue
            if ip_exclude_h and atom_a.symbol == "H":
                continue
            atom_a_bonds = atom_a.bonds
            if ip_known_bond_check:
                known_a_idxs = sorted([given_known_as.index(_a) for _a in atom_a_bonds if _a in given_known_as])
                if len(known_a_idxs) == 0:
                    continue
            for atom_b in atoms_b:
                if atom_b in given_known_bs:
                    continue
                if ip_exclude_h and atom_b.symbol == "H":
                    continue
                if ip_check_symbol and atom_a.symbol != atom_b.symbol:
                    continue
                atom_b_bonds = atom_b.bonds
                if (
                    ip_terminal_first_check
                    and terminal_first
                    and len([_a for _a in atom_a_bonds if _a.symbol != "H"]) > 1
                    and len([_b for _b in atom_b_bonds if _b.symbol != "H"]) > 1
                ):
                    continue
                if ip_check_hs and len([_a for _a in atom_a_bonds if _a.symbol == "H"]) != len(
                    [_b for _b in atom_b_bonds if _b.symbol == "H"]
                ):
                    continue
                if ip_isolated_check and not (len(atom_a_bonds) == 0 and len(atom_b_bonds) == 0):
                    continue
                if ip_known_bond_check:
                    known_b_idxs = sorted([given_known_bs.index(_b) for _b in atom_b_bonds if _b in given_known_bs])
                    if len(known_b_idxs) == 0:
                        continue
                    if known_a_idxs != known_b_idxs:
                        continue
                    if not (appended_known_idxs is None or appended_known_idxs == known_a_idxs):
                        continue
                    appended_known_idxs = known_a_idxs
                initial_pairs.append((atom_a, atom_b))

    def _isproper_bonding(atom_a: Atom, atom_b: Atom, side_pairs: list[tuple[Atom]]):
        atom_a_bonds = atom_a.bonds
        atom_b_bonds = atom_b.bonds
        # I swaped these two if sentence, but I not sure it was needed.
        if len([a for a in atom_a_bonds if a.symbol == "H"]) != len([a for a in atom_b_bonds if a.symbol == "H"]):
            return False
        if len(atom_a_bonds) != len(atom_b_bonds):
            return True
        a_side_pairs = [_ab[0] for _ab in side_pairs]
        b_side_pairs = [_ab[1] for _ab in side_pairs]
        a_index_list = sorted([a_side_pairs.index(_a) for _a in atom_a_bonds if _a in a_side_pairs])
        b_index_list = sorted([b_side_pairs.index(_a) for _a in atom_b_bonds if _a in b_side_pairs])
        if a_index_list != b_index_list:
            # trap 1,2-migration to return True
            a_diff_idx = list(set(a_index_list) - set(b_index_list))
            b_diff_idx = list(set(b_index_list) - set(a_index_list))
            if len(a_diff_idx) == 1 and len(b_diff_idx) == 1:
                if (
                    a_side_pairs[a_diff_idx[0]] in a_side_pairs[b_diff_idx[0]].bonds
                    and b_side_pairs[b_diff_idx[0]] in b_side_pairs[a_diff_idx[0]].bonds
                ):
                    return True
            return False
        return True

    logger.debug(f"generating initial_chains from {len(initial_pairs)} initial_pairs")
    initial_chains: list[list[tuple[Atom]]] = [[]]
    for pair_num, initial_pair in enumerate(initial_pairs, 1):
        if pair_num % 50 == 0:
            logger.debug(f"processing: {pair_num}/{len(initial_pairs)}")
        stack = [[initial_pair]]
        while stack:
            side_pairs: list[tuple[Atom]] = stack.pop()
            if len(side_pairs) > len(initial_chains[-1]):
                initial_chains = [side_pairs]
            elif len(side_pairs) == len(initial_chains[-1]):
                initial_chains.append(side_pairs)
            a_list = [pr[0] for pr in side_pairs]
            b_list = [pr[1] for pr in side_pairs]
            for next_a in side_pairs[-1][0].bonds:
                if next_a.symbol == "H" or next_a in a_list:
                    continue
                if next_a in given_known_as:
                    continue
                for next_b in side_pairs[-1][1].bonds:
                    if next_b.symbol == "H" or next_b in b_list:
                        continue
                    if next_a.symbol != next_b.symbol:
                        continue
                    if next_b in given_known_bs:
                        continue
                    if not _isproper_bonding(next_a, next_b, side_pairs):
                        continue
                    stack.append(side_pairs + [(next_a, next_b)])

    if len(initial_chains) == 1 and len(initial_chains[0]) == 0:
        return []

    canonical_initial_chains: dict[list[set], tuple[tuple[Atom]]] = {}
    for chain in initial_chains:
        key = tuple(sorted([(pr[0].number, pr[1].number) for pr in chain], key=lambda t: t[0]))
        canonical_initial_chains[key] = chain

    min_invalid = None
    h_matched_canonical_initial_chains: list[list[tuple[Atom]]] = [[]]
    for chain in canonical_initial_chains.values():
        invalid_atoms = 0
        for pair in chain:
            hs_of_a = [a_ for a_ in pair[0].bonds if a_.symbol == "H"]
            hs_of_b = [b_ for b_ in pair[1].bonds if b_.symbol == "H"]
            if len(hs_of_a) != len(hs_of_b):
                invalid_atoms += 1
        if min_invalid is None or min_invalid > invalid_atoms:
            h_matched_canonical_initial_chains = [chain]
            min_invalid = invalid_atoms
        elif min_invalid == invalid_atoms:
            h_matched_canonical_initial_chains.append(chain)

    if (
        4 <= len(h_matched_canonical_initial_chains)
        and len(h_matched_canonical_initial_chains) <= 8
        and max([len(chain) for chain in h_matched_canonical_initial_chains]) == 1
    ):
        logger.debug(f"generating combination_chains from {len(h_matched_canonical_initial_chains)} chains")
        sole_chains = [chain[0] for chain in h_matched_canonical_initial_chains]
        max_combs: list[list[tuple[Atom]]] = []
        for length in reversed(list(range(1, len(sole_chains) + 1))):
            for comb in itertools.combinations(sole_chains, length):
                if len({cb[0] for cb in comb}) != length:
                    continue
                if len({cb[1] for cb in comb}) != length:
                    continue
                max_combs.append(list(comb))
            if len(max_combs) != 0:
                break
        if len(max_combs) == 0:
            logger.debug("not found appropriate combination")
        else:
            h_matched_canonical_initial_chains = max_combs

    def _det_ez(a: Atom, b: Atom, c: Atom, d: Atom) -> bool:
        vab = np.array(a.xyz) - np.array(b.xyz)
        vcb = np.array(c.xyz) - np.array(b.xyz)
        vdc = np.array(d.xyz) - np.array(c.xyz)
        pvac = np.cross(vab, vcb)
        pvbd = np.cross(vdc, vcb)
        angle = np.arccos(np.sum(pvac * pvbd) / (np.linalg.norm(pvac) * np.linalg.norm(pvbd)))
        if np.sum(pvac * np.cross(pvbd, vcb)) < 0:
            angle = -angle
        angle = float(np.rad2deg(angle))
        if angle > 90 or angle < -90:
            return True
        return False

    def _det_chirality(a: Atom, b: Atom, c: Atom, d: Atom) -> bool:
        neighbors_xyz = np.array([atom_.xyz for atom_ in [b, c, d]]) - np.array(a.xyz)
        if np.linalg.det(neighbors_xyz) > 0:
            return True
        return False

    def _get_most_dihedral(a_list: list[Atom], b: Atom, c: Atom, d: Atom) -> tuple[Atom, float]:
        max_dihedral_a = a_list[0]
        max_dihedral_angle = 0
        vb = np.array(b.xyz)
        vc = np.array(c.xyz)
        vd = np.array(d.xyz)
        vcb = vc - vb
        vdc = vd - vc
        for a in a_list:
            va = np.array(a.xyz)
            vab = va - vb
            pvac: np.ndarray = np.cross(vab, vcb)
            pvbd: np.ndarray = np.cross(vdc, vcb)
            dac = np.linalg.norm(pvac)
            dbd = np.linalg.norm(pvbd)
            angle = np.arccos(np.sum(pvac * pvbd) / (dac * dbd))
            if np.sum(pvac * np.cross(pvbd, vcb)) < 0:
                angle = -angle
            angle = np.absolute(np.rad2deg(angle))
            if max_dihedral_angle <= angle:
                max_dihedral_a = a
                max_dihedral_angle = angle
        return (max_dihedral_a, max_dihedral_angle)

    logger.debug(f"extending {len(h_matched_canonical_initial_chains)} chains")
    extended_chains: list[list[tuple[Atom]]] = []
    for chain in h_matched_canonical_initial_chains:
        stack = deque([pr for pr in chain])
        if _recursion_counter == 0 and len(known_pairs) != 0:
            stack.extendleft(zip(given_known_as, given_known_bs))
        assigned_as = given_known_as + [pr[0] for pr in chain]
        assigned_bs = given_known_bs + [pr[1] for pr in chain]

        def _new_assign(new_pair: tuple[Atom], root_pair: tuple[Atom]):
            assigned_as.append(new_pair[0])
            assigned_bs.append(new_pair[1])
            stack.append((new_pair[0], new_pair[1]))
            stack.appendleft(root_pair)

        def _get_large(atoms: list[Atom]) -> list[Atom]:
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
                if max([len(tree[-1]) for tree in trees]) == 0:
                    return [list(tree[0])[0] for tree in trees]
                symbol_idxs_dict: dict[int, list[list[int]]] = {}
                all_symbol_set: set[int] = set()
                for idx, tree in enumerate(trees):
                    symbol_idxs_dict[idx] = [Elements.symbols.index(a.symbol) for a in tree[-1]]
                    all_symbol_set.update(symbol_idxs_dict[idx])
                for s_idx in sorted(list(all_symbol_set), reverse=True):
                    symbol_idxs_counts = {
                        idx: symbol_idxs.count(s_idx) for idx, symbol_idxs in symbol_idxs_dict.items()
                    }
                    max_count = max(symbol_idxs_counts.values())
                    symbol_idxs_dict = {
                        idx: symbol_idxs
                        for idx, symbol_idxs in symbol_idxs_dict.items()
                        if symbol_idxs_counts[idx] == max_count
                    }
                    if len(symbol_idxs_dict) == 1:
                        return list(trees[list(symbol_idxs_dict.keys())[0]][0])
                trees = [trees[idx] for idx in symbol_idxs_dict.keys()]
            else:
                return [list(tree[0])[0] for tree in trees]

        while stack:
            root_pair: tuple[Atom] = stack.popleft()
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
                if len(large_as) == 1 and len(large_bs) == 1:
                    _new_assign((large_as[0], large_bs[0]), root_pair)
                    continue
                not_large_as = [a for a in next_as if a not in large_as]
                not_large_bs = [b for b in next_bs if b not in large_bs]
                if len(not_large_as) == 1 and len(not_large_bs) == 1:
                    if not_large_as[0].symbol == not_large_bs[0].symbol:
                        _new_assign((not_large_as[0], not_large_bs[0]), root_pair)
                        continue
                known_a1_idxs = [assigned_as.index(a) for a in root_pair[0].bonds if a in assigned_as]
                known_b1_idxs = [assigned_bs.index(b) for b in root_pair[1].bonds if b in assigned_bs]
                known_1idxs = list(set(known_a1_idxs) & set(known_b1_idxs))
                if len(known_1idxs) == 0:
                    _new_assign((large_as[0], large_bs[0]), root_pair)
                    continue
                pre_large_a = large_as[0]
                pre_large_b = large_bs[0]
                max_dihedral = 0.0
                for known_1idx in known_1idxs:
                    known_a2_idxs = [
                        assigned_as.index(a)
                        for a in assigned_as[known_1idx].bonds
                        if (a in assigned_as) and (a != root_pair[0])
                    ]
                    known_b2_idxs = [
                        assigned_bs.index(b)
                        for b in assigned_bs[known_1idx].bonds
                        if (b in assigned_bs) and (b != root_pair[1])
                    ]
                    known_2idxs = list(set(known_a2_idxs) & set(known_b2_idxs))
                    for known_2idx in known_2idxs:
                        pre_a, dihedral_a = _get_most_dihedral(
                            large_as, root_pair[0], assigned_as[known_1idx], assigned_as[known_2idx]
                        )
                        pre_b, dihedral_b = _get_most_dihedral(
                            large_bs, root_pair[1], assigned_bs[known_1idx], assigned_bs[known_2idx]
                        )
                        if max(dihedral_a, dihedral_b) >= max_dihedral:
                            max_dihedral = max(dihedral_a, dihedral_b)
                            pre_large_a = pre_a
                            pre_large_b = pre_b
                _new_assign((pre_large_a, pre_large_b), root_pair)
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
                if len(root_pair[0].bonds) != len(root_pair[1].bonds):
                    # I'm not sure if this block is needed.
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

    canonical_extended_chains: dict[tuple[tuple[int]], list[tuple[Atom]]] = {}
    for chain in extended_chains:
        key = tuple(sorted([(pr[0].number, pr[1].number) for pr in chain], key=lambda t: t[0]))
        canonical_extended_chains[key] = chain

    for chain in canonical_extended_chains.values():
        logger.debug(f"canonical_extended_chains: {[(pr[0].number, pr[1].number) for pr in chain]}")

    extracted_canonical_extended_chains: list[list[tuple[Atom]]] = []
    max_h_len = max(
        [len([None for pr in chain if pr[0].symbol != "H"]) for chain in canonical_extended_chains.values()]
    )
    for chain in canonical_extended_chains.values():
        if max_h_len > len([None for pr in chain if pr[0].symbol != "H"]):
            logger.debug(f"excluding canonical_extended_chains: {[(pr[0].number, pr[1].number) for pr in chain]}")
            continue
        extracted_canonical_extended_chains.append(chain)

    appended_ecec: list[list[tuple[Atom]]] = [ecec[len(known_pairs) :] for ecec in extracted_canonical_extended_chains]
    while len(known_pairs) != 0 and len(appended_ecec) >= 4 and len(set([len(ecec) for ecec in appended_ecec])) == 1:
        unit_bs = set([tuple(sorted([p[1] for p in pl], key=lambda a: a.number)) for pl in appended_ecec])
        if len(unit_bs) == 1:
            break
        unit_as = set([tuple(sorted([p[0] for p in pl], key=lambda a: a.number)) for pl in appended_ecec])
        if len(unit_as) == 1:
            break

        is_isolated = True
        all_unit: set[tuple[Atom]] = unit_as | unit_bs
        for unit in all_unit:
            unit_bonding_atoms = []
            for a in unit:
                unit_bonding_atoms.extend(a.bonds)
            if True in [a in (given_known_as + given_known_bs) for a in unit_bonding_atoms]:
                is_isolated = False
                break
        if not is_isolated:
            break

        logger.debug(f"{max((len(unit_as), len(unit_bs)))} isolated molecules detected: initiating rmsd matching")
        total_ecec: list[tuple[Atom]] = list(zip(given_known_as, given_known_bs))
        while appended_ecec:
            ecec_rmsds: list[float] = []
            for ecec in appended_ecec:
                tested_ecec = total_ecec + ecec
                a_xyzs = np.array([tpr[0].xyz for tpr in tested_ecec])
                b_xyzs = np.array([tpr[1].xyz for tpr in tested_ecec])
                ecec_rmsds.append(_kabsch(a_xyzs, b_xyzs))
            added_app_ecec = appended_ecec.pop(ecec_rmsds.index(min(ecec_rmsds)))
            while len(added_app_ecec) == 1:
                if added_app_ecec[0][0].symbol != "O":
                    break
                if added_app_ecec[0][1].symbol != "O":
                    break
                hs_in_a = added_app_ecec[0][0].bonds
                if [a.symbol for a in hs_in_a] != ["H", "H"]:
                    break
                hs_in_b = added_app_ecec[0][1].bonds
                if [a.symbol for a in hs_in_b] != ["H", "H"]:
                    break
                logger.debug(f"detected H2O without H at O{(added_app_ecec[0][0].number,added_app_ecec[0][0].number)}")
                potential_h2os = [
                    [(hs_in_a[0], hs_in_b[0]), (hs_in_a[1], hs_in_b[1])],
                    [(hs_in_a[0], hs_in_b[1]), (hs_in_a[1], hs_in_b[0])],
                ]
                h2o_rmsds = []
                for hs_h2o in potential_h2os:
                    tested_ecec = total_ecec + added_app_ecec + hs_h2o
                    a_xyzs = np.array([tpr[0].xyz for tpr in tested_ecec])
                    b_xyzs = np.array([tpr[1].xyz for tpr in tested_ecec])
                    h2o_rmsds.append(_kabsch(a_xyzs, b_xyzs))
                added_app_ecec.extend(potential_h2os[h2o_rmsds.index(min(h2o_rmsds))])
                break
            total_ecec.extend(added_app_ecec)
            added_appended_ecec_atoms_a = [p[0] for p in added_app_ecec]
            added_appended_ecec_atoms_b = [p[1] for p in added_app_ecec]
            new_appended_ecec = []
            for ecec in appended_ecec:
                if True in [p[0] in added_appended_ecec_atoms_a for p in ecec]:
                    continue
                if True in [p[1] in added_appended_ecec_atoms_b for p in ecec]:
                    continue
                new_appended_ecec.append(ecec)
            appended_ecec = new_appended_ecec
        extracted_canonical_extended_chains = [total_ecec]
        break

    recursive_extended_chains: list[list[tuple[Atom]]] = []
    for chain in extracted_canonical_extended_chains:
        if len(chain) < min(len(atoms_a), len(atoms_b)):
            new_maps = _get_maps(
                atoms_a, atoms_b, known_pairs=chain, terminal_first=False, _recursion_counter=_recursion_counter + 1
            )
            if len(new_maps) == 0:
                recursive_extended_chains.append(chain)
            else:
                for new_map in new_maps:
                    recursive_extended_chains.append(chain + new_map)
        else:
            recursive_extended_chains.append(chain)

    canonical_recursive_extended_chains: dict[tuple[tuple[int]], list[tuple[Atom]]] = {}
    for chain in recursive_extended_chains:
        if len(known_pairs) != 0 and _recursion_counter != 0:
            chain = chain[len(known_pairs) :]
        key = tuple(sorted([(pr[0].number, pr[1].number) for pr in chain], key=lambda t: t[0]))
        canonical_recursive_extended_chains[key] = chain

    for chain in canonical_recursive_extended_chains.values():
        logger.debug(f"canonical_recursive_extended_chains: {[(pr[0].number, pr[1].number) for pr in chain]}")

    return list(canonical_recursive_extended_chains.values())


def _correct_maps(
    atoms_a: "Atoms",
    atoms_b: "Atoms",
    atom_maps: list[list[tuple[Atom]]],
    max_loop: int = 4,
    max_loop_rest_que: int = 8192,
    max_diff_bonds: int = 8,
) -> list[list[tuple[Atom]]]:
    logger.debug(f"correcting {len(atom_maps)} chains")
    total_atom_maps_dict: dict[tuple[tuple[int]], list[tuple[Atom]]] = {}
    for atom_map in atom_maps:
        key = tuple(sorted([(mp[0].number, mp[1].number) for mp in atom_map], key=lambda t: t[0]))
        total_atom_maps_dict[key] = atom_map
    for initial_atom_map in atom_maps:
        que_atom_maps = [initial_atom_map]
        loop_counter = 0
        while que_atom_maps:
            atom_map = que_atom_maps.pop(0)
            potential_swaps: list[list[Atom]] = []
            a_diff_bonds: list[tuple[Atom]] = []
            for ab_map in atom_map:
                bs_in_map = [p[1] for p in atom_map]
                a_diff_bd_from_b = [atom_map[bs_in_map.index(b)][0] for b in ab_map[1].bonds if b in bs_in_map]
                for a_diff_bd in set([a for a in ab_map[0].bonds]) ^ set(a_diff_bd_from_b):
                    a_diff_bonds.append((ab_map[0], a_diff_bd))
            a_diff_bonds = list(set([tuple(sorted(adb, key=lambda a: a.number)) for adb in a_diff_bonds]))
            if len(a_diff_bonds) > max_diff_bonds:
                logger.debug(f"len(a_diff_bonds) exceeded max_diff_bonds: {len(a_diff_bonds)}: exiting loop")
                continue
            if len(a_diff_bonds) >= 2:
                all_adb = [adb[0] for adb in a_diff_bonds] + [adb[1] for adb in a_diff_bonds]
                imaginary_bonds_candidate = [a for a in all_adb if all_adb.count(a) == 1]
                while len(imaginary_bonds_candidate) == 2:
                    if len([a for a in imaginary_bonds_candidate if a in [m[0] for m in initial_atom_map]]) == 2:
                        break
                    imaginary_bonds = tuple(sorted(imaginary_bonds_candidate, key=lambda a: a.number))
                    logger.debug(f"imaginary_bond appended: {tuple(a.number for a in imaginary_bonds)}")
                    a_diff_bonds.append(imaginary_bonds)
                    break
            loop_rest_que: list[dict] = [
                {
                    "loop_atoms": (adb[0],),
                    "rest_bonds": [b for b in a_diff_bonds if b is not adb],
                    "end_atom": adb[1],
                }
                for adb in a_diff_bonds
            ]
            stored_loop_atoms = []
            while loop_rest_que:
                if len(loop_rest_que) > max_loop_rest_que:
                    logger.error(f"loop_rest_que reached {max_loop_rest_que}: exiting loop")
                    break
                q = loop_rest_que.pop(0)
                loop_atoms: tuple[Atom] = q["loop_atoms"]
                rest_bonds: list[tuple[Atom]] = q["rest_bonds"]
                if loop_atoms[-1] is q["end_atom"]:
                    stored_loop_atoms.append(loop_atoms)
                for checking_bond in rest_bonds:
                    if loop_atoms[-1] is checking_bond[0]:
                        new_terminal = checking_bond[1]
                    elif loop_atoms[-1] is checking_bond[1]:
                        new_terminal = checking_bond[0]
                    else:
                        continue
                    loop_rest_que.append(
                        {
                            "loop_atoms": loop_atoms + (new_terminal,),
                            "rest_bonds": [b for b in rest_bonds if b is not checking_bond],
                            "end_atom": q["end_atom"],
                        }
                    )
            for loop_atoms in stored_loop_atoms:
                if len(set([a.symbol for a in loop_atoms[::2]])) == 1:
                    potential_swaps.append([a for a in loop_atoms[::2]])
                if len(set([a.symbol for a in loop_atoms[1::2]])) == 1:
                    potential_swaps.append([a for a in loop_atoms[1::2]])
            potential_swaps: list[tuple[Atom]] = sorted(
                list(set([tuple(sorted(adb, key=lambda a: a.number)) for adb in potential_swaps])),
                key=lambda ap: ap[0].number,
            )
            atom_map_a = [m[0] for m in atom_map]
            temporary_atom_maps_dict = {}
            non_hydrogen_swap_counter = 0
            for pq_swap in potential_swaps:
                if len(pq_swap) != 2:
                    continue
                try:
                    p_id: int = atom_map_a.index(pq_swap[0])
                except ValueError:
                    p_id = None
                try:
                    q_id: int = atom_map_a.index(pq_swap[1])
                except ValueError:
                    q_id = None
                new_map = atom_map[:]
                if None not in (p_id, q_id):
                    new_map[p_id], new_map[q_id] = (pq_swap[0], new_map[q_id][1]), (pq_swap[1], new_map[p_id][1])
                elif p_id is not None:
                    new_map[p_id] = (pq_swap[1], new_map[p_id][1])
                elif q_id is not None:
                    new_map[q_id] = (pq_swap[0], new_map[q_id][1])
                else:
                    logger.error("error on swapping potential_swap")
                    continue
                new_map_key = tuple(sorted([(mp[0].number, mp[1].number) for mp in new_map], key=lambda t: t[0]))
                if new_map_key not in total_atom_maps_dict:
                    if pq_swap[0].symbol != "H":
                        non_hydrogen_swap_counter += 1
                    temporary_atom_maps_dict[new_map_key] = new_map
            if loop_counter == 0 or non_hydrogen_swap_counter <= 1:
                for new_map_key, new_map in temporary_atom_maps_dict.items():
                    logger.info(f"corrected_chains: {[(pr[0].number, pr[1].number) for pr in new_map]}")
                    total_atom_maps_dict[new_map_key] = new_map
                    que_atom_maps.append(new_map)
            else:
                logger.info(f"multiple non_hydrogen_swap_counter detected: {non_hydrogen_swap_counter}")
            loop_counter += 1
            if loop_counter >= max_loop:
                logger.info(f"loop_counter: {loop_counter}")
                break
    return list(total_atom_maps_dict.values())


def _order_maps(
    atoms_a: "Atoms",
    atoms_b: "Atoms",
    atom_maps: list[list[tuple[Atom]]],
    options: dict = {},
) -> list[list[tuple[Atom]]]:
    def get_normal_vector(a: Atom, b: Atom, c: Atom) -> np.ndarray:
        try:
            ord_atoms: list[Atom] = sorted([a, b, c], key=lambda x: x.cache["map_idx"])
        except TypeError:
            ord_atoms: list[Atom] = sorted(
                [t for t in [a, b, c] if t.cache["map_idx"] is not None], key=lambda x: x.cache["map_idx"]
            ) + [t for t in [a, b, c] if t.cache["map_idx"] is None]
        vop = np.array(ord_atoms[1].xyz) - np.array(ord_atoms[0].xyz)
        voq = np.array(ord_atoms[2].xyz) - np.array(ord_atoms[0].xyz)
        return np.cross(vop, voq)

    def get_reactive_face_flag(center_atom: Atom, reactive_atom: Atom) -> bool:
        center_to_reactive_vec = np.array(reactive_atom.xyz) - np.array(center_atom.xyz)
        non_rea_atoms = [a for a in center_atom.single if a is not reactive_atom]
        if len(non_rea_atoms) != 3:
            return None
        nvec = get_normal_vector(non_rea_atoms[0], non_rea_atoms[1], non_rea_atoms[2])
        st = np.linalg.norm(center_to_reactive_vec) * np.linalg.norm(nvec)
        theta = np.arccos(np.inner(center_to_reactive_vec, nvec) / st)
        if np.degrees(theta) >= 90:
            return False
        else:
            return True

    def get_relationship_flag(atom_x: Atom, atom_y: Atom) -> bool:
        x_neis = atom_x.bonds
        y_neis = atom_y.bonds
        if len(x_neis) == 3 and len(y_neis) == 3:
            xvec = get_normal_vector(x_neis[0], x_neis[1], x_neis[2])
            yvec = get_normal_vector(y_neis[0], y_neis[1], y_neis[2])
        elif len(x_neis) == 3 and len(y_neis) == 2:
            xvec = get_normal_vector(x_neis[0], x_neis[1], x_neis[2])
            yvec = get_normal_vector(y_neis[0], y_neis[1], atom_y)
        elif len(x_neis) == 2 and len(y_neis) == 3:
            xvec = get_normal_vector(x_neis[0], x_neis[1], atom_x)
            yvec = get_normal_vector(y_neis[0], y_neis[1], y_neis[2])
        else:
            return None
        st = np.linalg.norm(xvec) * np.linalg.norm(yvec)
        theta = np.arccos(np.inner(xvec, yvec) / st)
        if np.degrees(theta) >= 90:
            return False
        else:
            return True

    def get_arc_relationship_flag(atom_x: Atom, atom_y_list: list[Atom]) -> bool:
        x_neis = atom_x.bonds
        if len(x_neis) == 3 and len(atom_y_list) == 3:
            xvec = get_normal_vector(x_neis[0], x_neis[1], x_neis[2])
            vop = np.array(atom_y_list[1].xyz) - np.array(atom_y_list[0].xyz)
            voq = np.array(atom_y_list[2].xyz) - np.array(atom_y_list[0].xyz)
            yvec = np.cross(vop, voq)
        else:
            return None
        st = np.linalg.norm(xvec) * np.linalg.norm(yvec)
        theta = np.arccos(np.inner(xvec, yvec) / st)
        if np.degrees(theta) >= 90:
            return False
        else:
            return True

    def det_chirality(a: Atom, b: Atom, c: Atom, d: Atom) -> bool:
        neighbors_xyz = np.array([atom_.xyz for atom_ in [b, c, d]]) - np.array(a.xyz)
        if np.linalg.det(neighbors_xyz) > 0:
            return True
        return False

    logger.debug(f"ordering {len(atom_maps)} chains")
    evaluation_template = {
        "metathesis_point": 0,
        "reactive_sp3_carbon_penalty": 0,
        "non_hydrogen_breaking_bonds": 0,
        "heteroatoms_on_breaking_atoms": 0,
        "num_of_long_rearranged_atoms": 0,
        "inconsistent_h_count": 0,
        "pericyclic_face": 0,
        "antarafacial_penalty": 0,
        "invalid_sn2_atoms": 0,
        "boat_transition_penalty": 0,
        "syn_or_anti_addition": 0,
        "bonding_local_rmsd": 0.0,
        "local_rmsd": 0.0,
        "total_rmsd": 0.0,
    }
    evaluated_dicts = [{"map": mp} | evaluation_template for mp in atom_maps]

    for a in atoms_a:
        a.cache["original_number"] = a.number
    for b in atoms_b:
        b.cache["original_number"] = b.number
    mol_num_a = {}
    mol_num_b = {}
    for idx, atms in enumerate(Modeler(atoms_a).get_splitted()):
        for a in atms:
            mol_num_a[a.cache["original_number"]] = idx
    for idx, atms in enumerate(Modeler(atoms_b).get_splitted()):
        for b in atms:
            mol_num_b[b.cache["original_number"]] = idx

    for chain_num, evaluated_dic in enumerate(evaluated_dicts, 1):
        if chain_num % 20 == 0:
            logger.debug(f"processing: {chain_num}/{len(atom_maps)}")
        for a in atoms_a:
            a.cache["map_idx"] = None
        for a in atoms_b:
            a.cache["map_idx"] = None
        atom_map: list[tuple[Atom]] = evaluated_dic["map"]
        for idx, pr in enumerate(atom_map):
            pr[0].cache["map_idx"] = idx
            pr[1].cache["map_idx"] = idx
        bonding_pairs: list[tuple[Atom]] = []
        normal_pairs: list[tuple[Atom]] = []
        ic_map = [p[0] for p in atom_map] + [p[1] for p in atom_map]

        def ic(atom: Atom) -> Atom:
            try:
                idx = ic_map.index(atom)
            except ValueError:
                return None
            return ic_map[(idx + len(atom_map)) % len(ic_map)]

        for pair in atom_map:
            a_bonding_nums = sorted([a.number for a in pair[0].bonds])
            b_bonding_nums = sorted([ic(b).number for b in pair[1].bonds if ic(b) is not None])
            if a_bonding_nums == b_bonding_nums:
                a_to_b = True
            else:
                a_to_b = False
            b_bonding_nums = sorted([b.number for b in pair[1].bonds])
            a_bonding_nums = sorted([ic(a).number for a in pair[0].bonds if ic(a) is not None])
            if a_bonding_nums == b_bonding_nums:
                b_to_a = True
            else:
                b_to_a = False

            if a_to_b and b_to_a:
                normal_pairs.append(pair)
            else:
                bonding_pairs.append(pair)

        for pr in bonding_pairs:
            local_pairs: list[tuple[Atom]] = [pr]
            _b_bonds = pr[1].bonds
            for a_bonding in pr[0].bonds:
                _aconv = ic(a_bonding)
                if not (_aconv is not None and _aconv in _b_bonds):
                    continue
                local_pairs.append((a_bonding, _aconv))
                for a_gem in a_bonding.bonds:
                    _b_gem_bonds = _aconv.bonds
                    _a_gem_conv = ic(a_gem)
                    if _a_gem_conv is not None and _a_gem_conv in _b_gem_bonds:
                        local_pairs.append((a_gem, _a_gem_conv))
            a_xyzs = np.array([tpr[0].xyz for tpr in local_pairs])
            b_xyzs = np.array([tpr[1].xyz for tpr in local_pairs])
            evaluated_dic["bonding_local_rmsd"] += _kabsch(a_xyzs, b_xyzs)
        evaluated_dic["bonding_local_rmsd"] = round(evaluated_dic["bonding_local_rmsd"], 1)

        for pr in normal_pairs:
            local_pairs: list[tuple[Atom]] = [pr]
            _b_bonds = pr[1].bonds
            for a_bonding in pr[0].bonds:
                _aconv = ic(a_bonding)
                if not (_aconv is not None and _aconv in _b_bonds):
                    continue
                local_pairs.append((a_bonding, _aconv))
                for a_gem in a_bonding.bonds:
                    _b_gem_bonds = _aconv.bonds
                    _a_gem_conv = ic(a_gem)
                    if _a_gem_conv is not None and _a_gem_conv in _b_gem_bonds:
                        local_pairs.append((a_gem, _a_gem_conv))
            a_xyzs = np.array([tpr[0].xyz for tpr in local_pairs])
            b_xyzs = np.array([tpr[1].xyz for tpr in local_pairs])
            evaluated_dic["local_rmsd"] += _kabsch(a_xyzs, b_xyzs)
        evaluated_dic["local_rmsd"] = round(evaluated_dic["local_rmsd"], 1)

        a_xyzs = np.array([m[0].xyz for m in atom_map])
        b_xyzs = np.array([m[1].xyz for m in atom_map])
        evaluated_dic["total_rmsd"] = _kabsch(a_xyzs, b_xyzs)

        for pr in bonding_pairs:
            for i0, i1 in [(0, 1), (1, 0)]:
                for a_bonding in pr[i0].bonds:
                    if a_bonding in [pr[i0] for pr in bonding_pairs]:
                        if ic(a_bonding) is None:
                            continue
                        if ic(a_bonding) in pr[i1].bonds:
                            continue
                        for b_bonding in ic(a_bonding).bonds:
                            if ic(b_bonding) in a_bonding.bonds:
                                continue
                            if b_bonding in [pr[i1] for pr in bonding_pairs]:
                                evaluated_dic["num_of_long_rearranged_atoms"] += 1

        map_index_dict = {m[0]: idx for idx, m in enumerate(atom_map)} | {m[1]: idx for idx, m in enumerate(atom_map)}
        bonding_atoms = [p[0] for p in bonding_pairs] + [p[1] for p in bonding_pairs]
        faces_dict: dict[tuple[int], bool] = {}
        all_reactive_idx_dict: dict[int, int] = {}
        for pr in bonding_pairs + [(_pr[1], _pr[0]) for _pr in bonding_pairs]:
            if pr[0].symbol == "H" and pr[1].symbol == "H":
                pr0_single = pr[0].single
                if len(pr0_single) != 1:
                    continue
                if map_index_dict.get(pr0_single[0]) is None:
                    continue
                all_reactive_idx_dict[map_index_dict[pr[0]]] = map_index_dict[pr0_single[0]]
                faces_dict[(map_index_dict[pr[0]],)] = True
                continue
            if len(ic(pr[0]).double) != 1:
                continue
            reactive_atoms = [a for a in pr[0].single if ic(a) not in pr[1].bonds]
            if len(reactive_atoms) != 1:
                continue
            if map_index_dict.get(reactive_atoms[0]) is None:
                continue
            all_reactive_idx_dict[map_index_dict[pr[0]]] = map_index_dict[reactive_atoms[0]]
            face_flag = get_reactive_face_flag(pr[0], reactive_atoms[0])
            stack: list[Atom] = [pr[1]]
            involving_atom_number = 1
            while True:
                if involving_atom_number > 16:
                    break
                twin_atom: Atom = stack[-1]
                if twin_atom is None:
                    break
                twin_neighbors = twin_atom.double
                if len(twin_neighbors) != 1:
                    break
                if len(stack) >= 2 and twin_neighbors[0] is stack[-2]:
                    break
                rel_flag = get_relationship_flag(twin_atom, twin_neighbors[0])
                if rel_flag is None or face_flag is None:
                    face_flag = None
                elif not rel_flag:
                    face_flag = not face_flag
                my_neighbor = ic(twin_neighbors[0])
                if len(stack) >= 2 and my_neighbor is stack[-2]:
                    break
                if my_neighbor in bonding_atoms:
                    my_neighbor_reactive_atoms = [
                        a for a in my_neighbor.single if ic(a) not in twin_neighbors[0].bonds
                    ]
                    if len(my_neighbor_reactive_atoms) != 1:
                        break
                    stack.append(my_neighbor)
                    stack_as_index = tuple(map_index_dict[a] for a in stack)
                    if face_flag is None:
                        faces_dict[stack_as_index] = None
                        break
                    next_face_flag = get_reactive_face_flag(my_neighbor, my_neighbor_reactive_atoms[0])
                    if next_face_flag is None:
                        faces_dict[stack_as_index] = None
                        break
                    faces_dict[stack_as_index] = next_face_flag is face_flag
                    break
                else:
                    stack.append(my_neighbor)
                    involving_atom_number += 1

        canonical_faces_dict: dict[tuple[int], bool] = {}
        for f, syn_anti in faces_dict.items():
            face = f
            if f[0] > f[-1]:
                face = tuple(reversed(f))
            if canonical_faces_dict.get(face, syn_anti) is syn_anti:
                canonical_faces_dict[face] = syn_anti
            else:
                logger.error(f"conflicting faces detected: {face}")

        pairwise_faces: list[dict[tuple[int], bool]] = []
        for f in canonical_faces_dict:
            ridx = (all_reactive_idx_dict.get(f[0]), all_reactive_idx_dict.get(f[-1]))
            if None in ridx:
                continue
            for pair_f in canonical_faces_dict:
                if ridx == (pair_f[0], pair_f[-1]):
                    pair_f_dict = {pair_f: canonical_faces_dict[pair_f]}
                elif ridx == (pair_f[-1], pair_f[0]):
                    pair_f_dict = {tuple(reversed(pair_f)): canonical_faces_dict[pair_f]}
                else:
                    continue
                new_pair_f_dict = {f: canonical_faces_dict[f]} | pair_f_dict
                if new_pair_f_dict not in pairwise_faces:
                    pairwise_faces.append(new_pair_f_dict)

            if len(f) in (2, 3) and canonical_faces_dict[f] is not None:
                a_same = mol_num_a[atom_map[ridx[0]][0].number] == mol_num_a[atom_map[ridx[1]][0].number]
                b_same = mol_num_b[atom_map[ridx[0]][1].number] == mol_num_b[atom_map[ridx[1]][1].number]
                rearrangement_flag = False
                if not a_same and len(f) == 3:
                    mol_num = mol_num_a[atom_map[f[0]][0].number]
                    if mol_num == mol_num_a[atom_map[f[-1]][0].number]:
                        if mol_num_a[atom_map[ridx[0]][0].number] == mol_num:
                            rearrangement_flag = not rearrangement_flag
                        if mol_num_a[atom_map[ridx[1]][0].number] == mol_num:
                            rearrangement_flag = not rearrangement_flag
                if not b_same and len(f) == 3:
                    mol_num = mol_num_b[atom_map[f[0]][1].number]
                    if mol_num == mol_num_b[atom_map[f[-1]][1].number]:
                        if mol_num_b[atom_map[ridx[0]][1].number] == mol_num:
                            rearrangement_flag = not rearrangement_flag
                        if mol_num_b[atom_map[ridx[1]][1].number] == mol_num:
                            rearrangement_flag = not rearrangement_flag
                add_flags = (rearrangement_flag ^ canonical_faces_dict[f], a_same, b_same)
                if add_flags in ((True, True, True), (False, True, False), (False, False, True)):
                    evaluated_dic["syn_or_anti_addition"] -= 1

        evaluated_dic["antarafacial_penalty"] = len([sa for sa in pairwise_faces if False in sa.values()])

        for p in pairwise_faces:
            if None in p.values():
                continue
            fs = sorted(list(p.keys()), key=lambda f: len(f))
            if len(fs) != 2 or len(fs[-1]) <= 3:
                continue
            if len(fs[0]) == 1 and len(fs[-1]) % 2 != 0:
                fs_len = len(fs[-1]) + 1
            elif fs[0] == tuple(reversed(fs[-1])) and len(fs[-1]) % 2 == 0:
                fs_len = len(fs[-1])
            else:
                continue
            if (fs_len % 4 == 0) ^ p[fs[-1]]:
                evaluated_dic["pericyclic_face"] -= 1
            else:
                evaluated_dic["pericyclic_face"] += 1

        for p in pairwise_faces:
            fs = list(p.keys())
            if len(fs) != 2:
                continue
            if fs[0] != tuple(reversed(fs[1])):
                continue
            if len(fs[0]) % 2 != 0:
                continue
            if (len(fs[0]) % 4 == 0) ^ p[fs[0]]:
                evaluated_dic["pericyclic_face"] -= 1
            else:
                evaluated_dic["pericyclic_face"] += 1

        for p in pairwise_faces:
            if [len(f) for f in p.keys()] != [3, 3]:
                continue
            face_tuple = tuple(face for face in p.keys())
            for ab_idx in (0, 1):
                for face_a, face_b in (face_tuple, tuple(reversed(face_tuple))):
                    face_a_atoms: list[Atom] = [atom_map[idx][ab_idx] for idx in face_a]
                    face_a_twins: list[Atom] = [ic(a) for a in face_a_atoms]
                    face_b_atoms: list[Atom] = [atom_map[idx][ab_idx] for idx in face_b]
                    face_b_twins: list[Atom] = [ic(a) for a in face_b_atoms]
                    if face_b_atoms[0] not in face_a_atoms[0].single or face_a_atoms[0] not in face_b_atoms[0].single:
                        continue
                    face_a_flag = get_reactive_face_flag(face_a_atoms[0], face_b_atoms[0])
                    face_b_flag = get_reactive_face_flag(face_b_atoms[0], face_a_atoms[0])
                    if None in (face_a_flag, face_b_flag):
                        continue
                    arc_a_flag = get_arc_relationship_flag(face_a_twins[0], face_a_twins)
                    arc_b_flag = get_arc_relationship_flag(face_b_twins[0], face_b_twins)
                    if None in (arc_a_flag, arc_b_flag):
                        continue
                    evaluated_dic["face_flags"] = (face_a_flag, face_b_flag, arc_a_flag, arc_b_flag)
                    if not face_a_flag:
                        arc_a_flag = not arc_a_flag
                    if not face_b_flag:
                        arc_b_flag = not arc_b_flag
                    if arc_a_flag ^ arc_b_flag:
                        evaluated_dic["boat_transition_penalty"] += 1

        for pr in bonding_pairs:
            if pr[0].symbol != pr[1].symbol:
                continue
            a_singles = pr[0].single
            b_singles = pr[1].single
            if len(a_singles) != 4 or len(b_singles) != 4:
                continue
            aconv_from_a = [ic(a) for a in a_singles]
            aconv_from_b = [ic(b) for b in b_singles]
            diff_a_atoms = [a for a in a_singles if a not in aconv_from_b]
            diff_b_atoms = [b for b in b_singles if b not in aconv_from_a]
            if len(diff_a_atoms) != 1 or len(diff_b_atoms) != 1:
                continue
            if diff_a_atoms[0].symbol == "H" or diff_b_atoms[0].symbol == "H":
                continue
            if diff_a_atoms[0].symbol == "C" and diff_b_atoms[0].symbol == "C":
                continue
            a_ord_atoms = [a for a in a_singles if a not in diff_a_atoms]
            b_ord_atoms = [ic(a) for a in a_ord_atoms]
            if det_chirality(*(a_ord_atoms + diff_a_atoms)) == det_chirality(*(b_ord_atoms + diff_b_atoms)):
                evaluated_dic["invalid_sn2_atoms"] += 1

        breaking_bonds: list[tuple[Atom]] = []
        for pr in bonding_pairs:
            if pr[0].symbol != pr[1].symbol:
                continue
            a_bonds = pr[0].bonds
            b_bonds = pr[1].bonds
            aconv_from_a = [ic(a) for a in a_bonds]
            aconv_from_b = [ic(b) for b in b_bonds]
            diff_a_atoms = [a for a in a_bonds if a not in aconv_from_b]
            diff_b_atoms = [b for b in b_bonds if b not in aconv_from_a]
            for diff_atm in [a for a in a_bonds if a not in aconv_from_b]:
                breaking_bonds.append((pr[0], diff_atm))
            for diff_atm in [b for b in b_bonds if b not in aconv_from_a]:
                breaking_bonds.append((pr[1], diff_atm))
        breaking_bonds: list[tuple[Atom]] = list(
            {tuple(sorted(b, key=lambda a: a.number)): None for b in breaking_bonds}.keys()
        )
        non_hydrogen_breaking_bonds: list[tuple[Atom]] = [
            b for b in breaking_bonds if "H" != b[0].symbol and "H" != b[1].symbol
        ]
        evaluated_dic["non_hydrogen_breaking_bonds"] += len(non_hydrogen_breaking_bonds)

        for bb in non_hydrogen_breaking_bonds:
            bbb_atoms: list[Atom] = []
            bbb_atoms.extend([a for a in bb[0].bonds if a is not bb[1]] + [a for a in bb[1].bonds if a is not bb[0]])
            bbb_atoms.extend([a for a in bb[0].double if a is not bb[1]] + [a for a in bb[1].double if a is not bb[0]])
            bbb_atoms.extend(
                [a for a in bb[0].aromatic if a is not bb[1]] + [a for a in bb[1].aromatic if a is not bb[0]]
            )
            bbb_atoms.extend([a for a in bb[0].triple if a is not bb[1]] + [a for a in bb[1].triple if a is not bb[0]])
            bbb_atoms.extend([a for a in bb[0].triple if a is not bb[1]] + [a for a in bb[1].triple if a is not bb[0]])
            evaluated_dic["heteroatoms_on_breaking_atoms"] -= len([a for a in bbb_atoms if a.symbol not in ("C", "H")])

        non_hydrogen_breaking_atoms: list[Atom] = [bb[0] for bb in non_hydrogen_breaking_bonds] + [
            bb[1] for bb in non_hydrogen_breaking_bonds
        ]
        bb_dict: dict[Atom, Atom] = {
            bb[0]: bb[1] for bb in non_hydrogen_breaking_bonds if non_hydrogen_breaking_atoms.count(bb[0]) == 1
        }
        bb_dict.update(
            {bb[1]: bb[0] for bb in non_hydrogen_breaking_bonds if non_hydrogen_breaking_atoms.count(bb[1]) == 1}
        )

        metathesis_point = 0
        for start_atom in bb_dict:
            next_atom = start_atom
            for _ in range(4):
                if next_atom.symbol != "C":
                    break
                if len(next_atom.bonds) >= 4:
                    break
                next_atom = bb_dict.get(next_atom)
                if next_atom is None:
                    break
                next_atom = ic(next_atom)
                if next_atom is None:
                    break
            else:
                if next_atom is start_atom:
                    metathesis_point += 1
        if metathesis_point != 0 and (metathesis_point % 8) == 0:
            evaluated_dic["metathesis_point"] -= metathesis_point / 8

        for pr in bonding_pairs:
            if pr[0].symbol == "C" and pr[1].symbol == "C":
                if len(pr[0].bonds) == 4 and len(pr[1].bonds) == 4:
                    evaluated_dic["reactive_sp3_carbon_penalty"] += 1

        for pr in bonding_pairs:
            if pr[0].symbol != "C" or pr[1].symbol != "C":
                continue
            if len([a for a in pr[0].bonds if a.symbol == "H"]) != len([a for a in pr[1].bonds if a.symbol == "H"]):
                evaluated_dic["inconsistent_h_count"] += 1

    for opt in options:
        if opt in evaluation_template.keys():
            logger.debug(f"{opt} was modified by option value: {options[opt]}")
            for evaluated_dic in evaluated_dicts:
                evaluated_dic[opt] *= options[opt]
        else:
            logger.error(f"ordering option not found: {opt}")

    for k in reversed(evaluation_template):
        evaluated_dicts = sorted(evaluated_dicts, key=lambda evaluated_dic: evaluated_dic[k])
    return_list = []
    for idx, evaluated_dic in enumerate(evaluated_dicts):
        return_list.append(evaluated_dic["map"])
        logger.info(f"map no.{idx+1}: {[(pr[0].number, pr[1].number) for pr in evaluated_dic['map']]}")
        for k in evaluation_template:
            logger.info(f"{k}: {evaluated_dic[k]}")
    return return_list


def _kabsch(ref_xyzs: np.ndarray, tar_xyzs: np.ndarray) -> float:
    ref_xyzs = ref_xyzs - np.mean(ref_xyzs, axis=0)
    tar_xyzs = tar_xyzs - np.mean(tar_xyzs, axis=0)
    u_mat, _, v_tmat = np.linalg.svd(np.dot(tar_xyzs.T, ref_xyzs))
    if (np.linalg.det(u_mat) * np.linalg.det(v_tmat)) < 0.0:
        u_mat[:, -1] = -u_mat[:, -1]
    r_mat = np.dot(u_mat, v_tmat)
    diff_xyzs = ref_xyzs - (np.dot(tar_xyzs, r_mat))
    return np.sqrt((diff_xyzs * diff_xyzs).sum() / len(diff_xyzs))


def _rmsd_calc_loop(min_rmsd: float, max_cycle: int, rot_mats_dicts, init_rot_mat, np_ref_xyzs, np_tar_xyzs):
    prev_min_rmsd = min_rmsd
    min_rot_mat = init_rot_mat
    for _cycle in range(max_cycle):
        min_num_rotamer = -1
        for _i, _r in enumerate(rot_mats_dicts):
            _r["flag"] = not _r["flag"]
            _rot_mat = init_rot_mat
            for _true_mat in [_rot["matrix"] for _rot in rot_mats_dicts if _rot["flag"] is True]:
                _rot_mat = np.dot(_true_mat, _rot_mat)

            _rmsd = _kabsch(np_ref_xyzs, np.dot(_rot_mat, np_tar_xyzs))
            if _rmsd < prev_min_rmsd:
                prev_min_rmsd = _rmsd
                min_num_rotamer = _i
                min_rot_mat = _rot_mat
            _r["flag"] = not _r["flag"]
        if min_num_rotamer == -1:
            min_rmsd = prev_min_rmsd
            break
        else:
            logger.debug(
                "RMSD was updated using rotamer {} in cycle of {}: {:.10f}".format(
                    min_num_rotamer + 1, _cycle + 1, prev_min_rmsd
                )
            )
            rot_mats_dicts[min_num_rotamer]["flag"] = not rot_mats_dicts[min_num_rotamer]["flag"]
    else:
        logger.error("max cycle reached")
    return min_rmsd, min_rot_mat


def _cal_sym_rmsd(
    ref_atoms: Atoms,
    tar_atoms: Atoms,
    max_cycle: int,
):
    np_ref_xyzs = np.array([_a.xyz for _a in ref_atoms])
    np_tar_xyzs = np.array([_a.xyz for _a in tar_atoms])
    xyz_length = len(np_ref_xyzs)
    if xyz_length != len(np_tar_xyzs):
        raise Exception
    min_rmsd = _kabsch(np_ref_xyzs, np_tar_xyzs)
    min_rot_mat = np.identity(xyz_length)
    logger.debug(f"initial RMSD: {min_rmsd:.10f}")
    tar_conf_rotamers: list[Matrix] = tar_atoms.data["rotamer"]
    tar_conf_atoms_list = tar_atoms.to_list()
    _rotamers = [{"flag": False, "matrix": _m.ordered(tar_conf_atoms_list).to_ndarray()} for _m in tar_conf_rotamers]
    min_rmsd, min_rot_mat = _rmsd_calc_loop(
        min_rmsd,
        max_cycle,
        _rotamers,
        np.identity(xyz_length),
        np_ref_xyzs,
        np_tar_xyzs,
    )
    return min_rmsd


def _is_heavy(atom_heavy: Atom, atom_light: Atom):
    heavy_idx = Elements.symbols.index(atom_heavy.symbol)
    light_idx = Elements.symbols.index(atom_light.symbol)
    if heavy_idx > light_idx:
        return True
    elif heavy_idx < light_idx:
        return False
    return None


class Chains(MutableSequence):
    def __init__(self, chains: list[tuple[Atom]] = []) -> None:
        self._chains: list[tuple[Atom]] = chains[:]

    def __getitem__(self, idx: int) -> tuple[Atom]:
        return self._chains[idx]

    def __setitem__(self, idx: int, chain: tuple[Atom]):
        self._chains[idx] = chain

    def __delitem__(self, idx: int):
        del self._chains[idx]

    def __len__(self):
        return len(self._chains)

    def insert(self, idx: int, chain: tuple[Atom]):
        self._chains.insert(idx, chain)

    def __add__(self, other: "Chains"):
        if not isinstance(other, Chains):
            raise TypeError
        return Chains(self._chains + other._chains)

    def __iter__(self):
        return self._chains.__iter__()

    @staticmethod
    def is_same_substituent(chain_a: tuple[Atom], chain_b: tuple[Atom]):
        return (chain_a[0] is chain_b[0]) and (chain_a[1] is chain_b[1])

    @staticmethod
    def is_higher(chain_high: tuple[Atom], chain_low: tuple[Atom]):
        zipped = list(zip(chain_high, chain_low))
        for ha, la in zipped:
            is_heavy = _is_heavy(ha, la)
            if is_heavy is None:
                continue
            if is_heavy:
                return True
            else:
                return False
        if len(chain_high) > len(chain_low):
            return True
        elif len(chain_high) < len(chain_low):
            return False
        symb_dict_high: dict[str, float] = defaultdict(float)
        symb_dict_low: dict[str, float] = defaultdict(float)
        for chain, symb_dict in [(chain_high, symb_dict_high), (chain_low, symb_dict_low)]:
            for a in chain[-1].single:
                symb_dict[a.symbol] += 1
            for a in chain[-1].double:
                symb_dict[a.symbol] += 2
            for a in chain[-1].triple:
                symb_dict[a.symbol] += 3
            for a in chain[-1].aromatic:
                symb_dict[a.symbol] += 1.5
            if len(chain) >= 2:
                symb_dict[chain[-2].symbol] -= 1
        symbidx_dict_high: dict[int, float] = {Elements.symbols.index(s): v for s, v in dict(symb_dict_high).items()}
        symbidx_dict_low: dict[int, float] = {Elements.symbols.index(s): v for s, v in dict(symb_dict_low).items()}
        symbidxs: list[int] = sorted(list(set(symbidx_dict_high.keys()) | set(symbidx_dict_low.keys())), reverse=True)
        for idx in symbidxs:
            order_high = symbidx_dict_high.get(idx, 0)
            order_low = symbidx_dict_low.get(idx, 0)
            if order_high > order_low:
                return True
            elif order_high < order_low:
                return False
        return None

    @property
    def max_distance(self):
        return max([len(chain) - 1 for chain in self._chains])

    def elongated(self, accept_loop: bool = False) -> "Chains":
        new_chains = Chains()
        for chain in self._chains:
            for a in chain[-1].bonds:
                if a is chain[-2]:
                    continue
                if not accept_loop:
                    if a in chain:
                        continue
                new_chains.append(chain + (a,))
        return new_chains

    def expanded(self, accept_loop=False, max_elongation=128) -> "Chains":
        new_chains = Chains(self._chains)
        for _ in range(max_elongation):
            terminal_chains = self.elongated(accept_loop=accept_loop)
            if len(terminal_chains) == 0:
                break
            new_chains.extend(terminal_chains)
        return new_chains

    def longest(self) -> "Chains":
        max_len = self.max_distance + 1
        return Chains([chain for chain in self._chains if len(chain) == max_len])

    def heavily_terminated(self) -> "Chains":
        symbs = set([chain[-1].symbol for chain in self._chains])
        max_symbol = Elements.canonicalize(max([Elements.symbols.index(s) for s in symbs]))
        return Chains([chain for chain in self._chains if chain[-1].symbol == max_symbol])

    @property
    def terminals(self) -> list[Atom]:
        return [chain[-1] for chain in self._chains]

    @property
    def terminal_neighbors(self) -> list[list[Atom]]:
        return [chain[-1].bonds for chain in self._chains]

    @property
    def length_of_substituent(self) -> int:
        pass

    def sort(self):
        for i in range(1, len(self._chains)):
            for j in range(0, len(self._chains) - i):
                if Chains.is_higher(self._chains[j], self._chains[j + 1]) is False:
                    self._chains[j], self._chains[j + 1] = self._chains[j + 1], self._chains[j]

    @staticmethod
    def get_ordered_substituents(initial_chains_list: list["Chains"]) -> list["Chains"]:
        sorted_chains_list: list[Chains] = []
        if len(initial_chains_list) == 0:
            return sorted_chains_list
        if len(initial_chains_list) == 1:
            sorted_chains_list.append(initial_chains_list[0])
            return sorted_chains_list
        max_distance = max([chains.max_distance for chains in initial_chains_list])
        if max_distance > 64:
            logger.error("exceeded max recursive iteration")
            return None

        t_chains_list = [Chains(initial_chains._chains[:]) for initial_chains in initial_chains_list]
        for t_chains in t_chains_list:
            t_chains.sort()
        for _ in range(len(t_chains_list)):
            if len(t_chains_list) in [0, 1]:
                break
            higher_chains = [t_chains[0] for t_chains in t_chains_list]
            max_chains = [higher_chains[0]]
            for chain in higher_chains[1:]:
                for mxchain in max_chains[:]:
                    flag = Chains.is_higher(chain, mxchain)
                    if flag is True:
                        max_chains = [chain]
                        break
                    elif flag is None:
                        max_chains.append(chain)
            if len(max_chains) == 1:
                sorted_chains_list.extend(
                    [chains for chains in initial_chains_list if chains[0][1] == max_chains[0][1]]
                )
                for t_chains in t_chains_list[:]:
                    if t_chains[0][1] == max_chains[0][1]:
                        t_chains_list.remove(t_chains)
            elif len(max_chains) >= 2:
                pass
            else:
                break
            for max_chain in max_chains:
                for t_chains in t_chains_list:
                    if max_chain in t_chains._chains:
                        t_chains._chains.remove(max_chain)
            t_chains_list = [t_chains for t_chains in t_chains_list if len(t_chains) != 0]

        next_chains_list: list[Chains] = []
        elongated_flag = False
        for chains in initial_chains_list:
            for sorted_chains in sorted_chains_list:
                if chains[0][1] == sorted_chains[0][1]:
                    break
            else:
                elongated_chains = chains.elongated()
                for chain in chains:
                    if chain in elongated_chains:
                        elongated_chains.remove(chain)
                if len(elongated_chains) != 0:
                    elongated_flag = True
                    next_chains_list.append(chains + elongated_chains)
                else:
                    next_chains_list.append(chains)
        if len(next_chains_list) == 1:
            sorted_chains_list.extend(next_chains_list)
            return sorted_chains_list
        if len(next_chains_list) == 0 or elongated_flag is False:
            return None
        max_dst = max([ch.max_distance for ch in next_chains_list])
        all_term_atoms: set[Atom] = set()
        for term_atoms in [ch.longest().terminals for ch in next_chains_list if ch.max_distance == max_dst]:
            all_term_atoms = all_term_atoms | set(term_atoms)
        if len(all_term_atoms) == 1:
            return None
        ordered_chains_list = Chains.get_ordered_substituents(next_chains_list)
        if ordered_chains_list is None:
            return None
        for chains in ordered_chains_list:
            for initial_chains in initial_chains_list:
                if chains[0][1] == initial_chains[0][1]:
                    sorted_chains_list.append(initial_chains)
                    break
        return sorted_chains_list
