import math
from collections import defaultdict, deque
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
        terminal_first: bool = False,
    ) -> list[list[int]]:
        atoms_map = _get_maps(target, self.atoms, known_pairs=known_pairs, terminal_first=terminal_first)
        atoms_map = _order_maps(target, self.atoms, atoms_map)
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
        terminal_first: bool = False,
    ) -> list["Atoms"]:
        atoms_map = _get_maps(target, self.atoms, known_pairs=known_pairs, terminal_first=terminal_first)
        atoms_map = _order_maps(target, self.atoms, atoms_map)
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
            elif [len(a.bonds), len(a.cache["isomeric_subs_list"])] in [[4, 3], [4, 6]]:
                for _ps in a.cache["isomeric_subs_list"][1:]:
                    rot_mat = np.dot(_ps["matrix"], a.cache["isomeric_subs_list"][0]["matrix"])
                    rot_mat = rot_mat.astype(int)
                    rotamer_mat_list.append(rot_mat)
                for _ps in a.cache["isomeric_subs_list"]:
                    numisomer_mat_list.append(_ps["matrix"])
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

        # for TMS, cyclic_chiral_check has error
        def _cyclic_chiral_check(atoms: Atoms, rotamer_mat_list, numisomer_mat_list):
            unchanged_flag = True
            for rm_idx, rot_mat in enumerate(rotamer_mat_list):
                original_xyz = np.array([at.xyz for at in atoms])
                number_of_invalid = _invalid_chirality(atoms, np.dot(rot_mat, original_xyz))
                if number_of_invalid == 0:
                    continue
                for num_mat in numisomer_mat_list:
                    t_mat = np.dot(num_mat, rot_mat)
                    if _invalid_chirality(atoms, np.dot(t_mat, original_xyz)) < number_of_invalid:
                        rotamer_mat_list[rm_idx] = t_mat
                        unchanged_flag = False
                        logger.debug("chirality in rotamer matrix was partially corrected")
                        break
                else:
                    logger.debug("chirality in rotamer matrix temporarily corrected")
            return unchanged_flag

        if len(rotamer_mat_list) != 0:
            for _ in range(len(rotamer_mat_list)):
                if _cyclic_chiral_check(self.atoms, rotamer_mat_list, numisomer_mat_list):
                    logger.debug("chirality in rotamer matrix was successfully corrected")
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
    exclude_known: bool = False,
    terminal_first: bool = False,
) -> list[list[tuple[Atom]]]:

    if len(known_pairs) == 0:
        known_check = False
        given_known_as = []
        given_known_bs = []
    else:
        known_check = True
        given_known_as = [atoms_a.get(pr[0]) for pr in known_pairs]
        given_known_bs = [atoms_b.get(pr[1]) for pr in known_pairs]
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

    initial_pairs: list[tuple[Atom]] = []

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
        appended_known_idxs = None
        for atom_a in atoms_a:
            for atom_b in atoms_b:
                if known_check:
                    if _should_exclude(atom_a, atom_b):
                        continue
                    known_a_idxs = [given_known_as.index(_a) for _a in atom_a.bonds if _a in given_known_as]
                    known_a_idxs = sorted(known_a_idxs)
                    known_b_idxs = [given_known_bs.index(_b) for _b in atom_b.bonds if _b in given_known_bs]
                    known_b_idxs = sorted(known_b_idxs)
                    if len(known_a_idxs) == 0 or len(known_b_idxs) == 0:
                        continue
                    if known_a_idxs == known_b_idxs:
                        if appended_known_idxs is None or appended_known_idxs == known_a_idxs:
                            appended_known_idxs = known_a_idxs
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

    initial_chains: list[list[tuple[Atom]]] = [[]]
    for initial_pair in initial_pairs:
        stack = [[initial_pair]]
        while stack:
            side_pairs: list[tuple[Atom]] = stack.pop()
            if len(side_pairs) > len(initial_chains[-1]):
                initial_chains = [side_pairs]
                logger.debug(f"initial_chains updated: {[(pr[0].number, pr[1].number) for pr in side_pairs]}")
            elif len(side_pairs) == len(initial_chains[-1]):
                initial_chains.append(side_pairs)
                logger.debug(f"initial_chains appended: {[(pr[0].number, pr[1].number) for pr in side_pairs]}")
            a_list = [pr[0] for pr in side_pairs]
            b_list = [pr[1] for pr in side_pairs]
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
                    if not _isproper_bonding(next_a, next_b, side_pairs):
                        continue
                    stack.append(side_pairs + [(next_a, next_b)])

    if len(initial_chains) == 1 and len(initial_chains[0]) == 0:
        return []

    canonical_initial_chains: dict[list[set], tuple[tuple[Atom]]] = {}
    for chain in initial_chains:
        key = tuple(sorted([(pr[0].number, pr[1].number) for pr in chain], key=lambda t: t[0]))
        canonical_initial_chains[key] = chain

    for chain in canonical_initial_chains.values():
        logger.debug(f"canonical_initial_chains: {[(pr[0].number, pr[1].number) for pr in chain]}")

    min_invalid = None
    h_matched_chains: list[list[tuple[Atom]]] = [[]]
    for chain in canonical_initial_chains.values():
        invalid_atoms = 0
        for pair in chain:
            hs_of_a = [a_ for a_ in pair[0].bonds if a_.symbol == "H"]
            hs_of_b = [b_ for b_ in pair[1].bonds if b_.symbol == "H"]
            if len(hs_of_a) != len(hs_of_b):
                invalid_atoms += 1
        if min_invalid is None or min_invalid > invalid_atoms:
            h_matched_chains = [chain]
            min_invalid = invalid_atoms
        elif min_invalid == invalid_atoms:
            h_matched_chains.append(chain)

    for chain in h_matched_chains:
        logger.debug(f"h_matched_chains: {[(pr[0].number, pr[1].number) for pr in chain]}")

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

    extended_chains: list[list[tuple[Atom]]] = []
    for chain in h_matched_chains:
        stack = deque([pr for pr in chain])
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

    max_h_len = max(
        [len([None for pr in chain if pr[0].symbol != "H"]) for chain in canonical_extended_chains.values()]
    )

    recursive_extended_chains: list[list[tuple[Atom]]] = []
    for chain in canonical_extended_chains.values():
        if max_h_len > len([None for pr in chain if pr[0].symbol != "H"]):
            logger.debug(f"excluded canonical_extended_chains: {[(pr[0].number, pr[1].number) for pr in chain]}")
            continue
        if len(chain) < min(len(atoms_a), len(atoms_b)):
            new_maps = _get_maps(
                atoms_a,
                atoms_b,
                known_pairs=chain,
                exclude_known=True,
                terminal_first=False,
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
        if exclude_known:
            chain = chain[len(known_pairs) :]
        key = tuple(sorted([(pr[0].number, pr[1].number) for pr in chain], key=lambda t: t[0]))
        canonical_recursive_extended_chains[key] = chain

    for chain in canonical_recursive_extended_chains.values():
        logger.debug(f"canonical_recursive_extended_chains: {[(pr[0].number, pr[1].number) for pr in chain]}")

    return canonical_recursive_extended_chains.values()


def _order_maps(atoms_a: "Atoms", atoms_b: "Atoms", atom_maps: list[list[tuple[Atom]]]) -> list[list[tuple[Atom]]]:
    def aconv(atom: Atom, a_map: list[tuple[Atom]]) -> Atom:
        try:
            idx = [pr[0] for pr in a_map].index(atom)
        except ValueError:
            pass
        else:
            return a_map[idx][1]
        try:
            idx = [pr[1] for pr in a_map].index(atom)
        except ValueError:
            pass
        else:
            return a_map[idx][0]
        return None

    def get_normal_vector(a: Atom, b: Atom, c: Atom) -> np.ndarray:
        ord_atoms: list[Atom] = sorted([a, b, c], key=lambda x: x.cache["map_idx"])
        vop = np.array(ord_atoms[1].xyz) - np.array(ord_atoms[0].xyz)
        voq = np.array(ord_atoms[2].xyz) - np.array(ord_atoms[0].xyz)
        return np.cross(vop, voq)

    def get_reactive_face_flag(center_atom: Atom, reactive_atom: Atom) -> bool:
        center_to_reactive_vec = np.array(reactive_atom.xyz) - np.array(center_atom.xyz)
        non_rea_atoms = [a for a in center_atom.single if a is not reactive_atom]
        if len(non_rea_atoms) != 3:
            logger.error("not bearing non reactive 3 bonds")
            return True
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
        if len(x_neis) != 3 or len(y_neis) != 3:
            logger.error("not bearing 3 bonds")
            return True
        xvec = get_normal_vector(x_neis[0], x_neis[1], x_neis[2])
        yvec = get_normal_vector(y_neis[0], y_neis[1], y_neis[2])
        st = np.linalg.norm(xvec) * np.linalg.norm(yvec)
        theta = np.arccos(np.inner(xvec, yvec) / st)
        if np.degrees(theta) >= 90:
            return False
        else:
            return True

    evaluated_dicts = [{"map": mp} for mp in atom_maps]
    for evaluated_dic in evaluated_dicts:
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
        for pair in atom_map:
            a_bonding_nums = sorted([a.number for a in pair[0].bonds])
            b_bonding_nums = sorted(
                [aconv(b, atom_map).number for b in pair[1].bonds if aconv(b, atom_map) is not None]
            )
            if a_bonding_nums == b_bonding_nums:
                a_to_b = True
            else:
                a_to_b = False
            b_bonding_nums = sorted([b.number for b in pair[1].bonds])
            a_bonding_nums = sorted(
                [aconv(a, atom_map).number for a in pair[0].bonds if aconv(a, atom_map) is not None]
            )
            if a_bonding_nums == b_bonding_nums:
                b_to_a = True
            else:
                b_to_a = False

            if a_to_b and b_to_a:
                normal_pairs.append(pair)
            else:
                bonding_pairs.append(pair)

        evaluated_dic["bonding_local_rmsd"] = 0.0
        for pr in bonding_pairs:
            local_pairs: list[tuple[Atom]] = [pr]
            _b_bonds = pr[1].bonds
            for a_bonding in pr[0].bonds:
                _aconv = aconv(a_bonding, atom_map)
                if _aconv is not None and _aconv in _b_bonds:
                    local_pairs.append((a_bonding, _aconv))
            a_xyzs = np.array([tpr[0].xyz for tpr in local_pairs])
            b_xyzs = np.array([tpr[1].xyz for tpr in local_pairs])
            a_xyzs = a_xyzs - np.mean(a_xyzs, axis=0)
            b_xyzs = b_xyzs - np.mean(b_xyzs, axis=0)
            evaluated_dic["bonding_local_rmsd"] += _kabsch(a_xyzs, b_xyzs)

        evaluated_dic["bonding_tb_local_rmsd"] = 0.0
        for pr in bonding_pairs:
            local_pairs: list[tuple[Atom]] = [pr]
            _b_bonds = pr[1].bonds
            for a_bonding in pr[0].bonds:
                _aconv = aconv(a_bonding, atom_map)
                if not (_aconv is not None and _aconv in _b_bonds):
                    continue
                local_pairs.append((a_bonding, _aconv))
                for a_gem in a_bonding.bonds:
                    _b_gem_bonds = _aconv.bonds
                    _a_gem_conv = aconv(a_gem, atom_map)
                    if _a_gem_conv is not None and _a_gem_conv in _b_gem_bonds:
                        local_pairs.append((a_gem, _a_gem_conv))
            a_xyzs = np.array([tpr[0].xyz for tpr in local_pairs])
            b_xyzs = np.array([tpr[1].xyz for tpr in local_pairs])
            a_xyzs = a_xyzs - np.mean(a_xyzs, axis=0)
            b_xyzs = b_xyzs - np.mean(b_xyzs, axis=0)
            evaluated_dic["bonding_tb_local_rmsd"] += _kabsch(a_xyzs, b_xyzs)

        evaluated_dic["local_rmsd"] = 0.0
        for pr in normal_pairs:
            local_pairs: list[tuple[Atom]] = [pr]
            _b_bonds = pr[1].bonds
            for a_bonding in pr[0].bonds:
                _aconv = aconv(a_bonding, atom_map)
                if not (_aconv is not None and _aconv in _b_bonds):
                    continue
                local_pairs.append((a_bonding, _aconv))
                for a_gem in a_bonding.bonds:
                    _b_gem_bonds = _aconv.bonds
                    _a_gem_conv = aconv(a_gem, atom_map)
                    if _a_gem_conv is not None and _a_gem_conv in _b_gem_bonds:
                        local_pairs.append((a_gem, _a_gem_conv))
            a_xyzs = np.array([tpr[0].xyz for tpr in local_pairs])
            b_xyzs = np.array([tpr[1].xyz for tpr in local_pairs])
            a_xyzs = a_xyzs - np.mean(a_xyzs, axis=0)
            b_xyzs = b_xyzs - np.mean(b_xyzs, axis=0)
            evaluated_dic["local_rmsd"] += _kabsch(a_xyzs, b_xyzs)

        evaluated_dic["num_of_long_rearranged_atoms"] = 0
        for pr in bonding_pairs:
            for a_bonding in pr[0].bonds:
                if a_bonding in [pr[0] for pr in bonding_pairs]:
                    if aconv(a_bonding, atom_map) is None:
                        continue
                    if aconv(a_bonding, atom_map) in pr[1].bonds:
                        continue
                    for b_bonding in aconv(a_bonding, atom_map).bonds:
                        if aconv(b_bonding, atom_map) in a_bonding.bonds:
                            continue
                        if b_bonding in [pr[1] for pr in bonding_pairs]:
                            evaluated_dic["num_of_long_rearranged_atoms"] += 1
            for b_bonding in pr[1].bonds:
                if b_bonding in [pr[1] for pr in bonding_pairs]:
                    if aconv(b_bonding, atom_map) is None:
                        continue
                    if aconv(b_bonding, atom_map) in pr[0].bonds:
                        continue
                    for a_bonding in aconv(b_bonding, atom_map).bonds:
                        if aconv(a_bonding, atom_map) in b_bonding.bonds:
                            continue
                        if a_bonding in [pr[0] for pr in bonding_pairs]:
                            evaluated_dic["num_of_long_rearranged_atoms"] += 1
        evaluated_dic["num_of_long_rearranged_atoms"] /= 2

        evaluated_dic["sp2_face_transfer_mismatch"] = 0
        bonding_atoms = [p[0] for p in bonding_pairs] + [p[1] for p in bonding_pairs]
        for pr in bonding_pairs + [(_pr[1], _pr[0]) for _pr in bonding_pairs]:
            if pr[0].symbol != "C":
                continue
            if len(pr[0].single) != 4:
                continue
            reactive_atoms = [a for a in pr[0].single if aconv(a, atom_map) not in pr[1].bonds]
            if len(reactive_atoms) != 1:
                continue
            face_flag = get_reactive_face_flag(pr[0], reactive_atoms[0])

            stack: list[Atom] = [pr[1]]
            loop_counter = 0
            while stack:
                loop_counter += 1
                if loop_counter > 16:
                    break
                twin_atom: Atom = stack.pop()
                twin_neighbors = twin_atom.double
                if len(twin_neighbors) == 0 or twin_neighbors[0].symbol != "C":
                    continue
                if not get_relationship_flag(twin_atom, twin_neighbors[0]):
                    face_flag = not face_flag
                my_neighbor = aconv(twin_neighbors[0], atom_map)
                if my_neighbor in bonding_atoms:
                    my_neighbor_reactive_atoms = [
                        a for a in my_neighbor.single if aconv(a, atom_map) not in twin_neighbors[0].bonds
                    ]
                    if len(my_neighbor_reactive_atoms) != 1:
                        continue
                    if get_reactive_face_flag(my_neighbor, my_neighbor_reactive_atoms[0]) is not face_flag:
                        evaluated_dic["sp2_face_transfer_mismatch"] += 1
                else:
                    stack.append(my_neighbor)

    evaluated_dicts = sorted(evaluated_dicts, key=lambda d: d["local_rmsd"])
    # evaluated_dicts = sorted(evaluated_dicts, key=lambda d: round(d["bonding_local_rmsd"] * 0.5, 3))
    evaluated_dicts = sorted(evaluated_dicts, key=lambda d: round(d["bonding_tb_local_rmsd"], 1))
    evaluated_dicts = sorted(evaluated_dicts, key=lambda d: d["sp2_face_transfer_mismatch"])
    evaluated_dicts = sorted(evaluated_dicts, key=lambda d: d["num_of_long_rearranged_atoms"])
    return_list = []
    for evaluated_dic in evaluated_dicts:
        mp = [(pr[0].number, pr[1].number) for pr in evaluated_dic["map"]]
        return_list.append(evaluated_dic["map"])
        logger.info(f"mapping: {mp}")
        logger.info(f"sp2_face_transfer_mismatch: {evaluated_dic['sp2_face_transfer_mismatch']}")
        logger.info(f"num_of_long_rearranged_atoms: {evaluated_dic['num_of_long_rearranged_atoms']}")
        logger.info(f"bonding_local_rmsd: {evaluated_dic['bonding_local_rmsd']}")
        logger.info(f"bonding_tb_local_rmsd: {evaluated_dic['bonding_tb_local_rmsd']}")
        logger.info(f"local_rmsd: {evaluated_dic['local_rmsd']}")

    return return_list


def _kabsch(ref_xyzs: np.ndarray, tar_xyzs: np.ndarray) -> float:
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
    np_ref_xyzs = np_ref_xyzs - np.mean(np_ref_xyzs, axis=0)
    np_tar_xyzs = np_tar_xyzs - np.mean(np_tar_xyzs, axis=0)
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


def _order_by_cip(substitutions: list[Atom], roots: list[Atom]):
    for sub in substitutions:
        pass
    return []
