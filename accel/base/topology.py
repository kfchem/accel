import itertools
from collections import defaultdict
from typing import Dict, List

import numpy as np
from accel.base.atoms import Atom, Atoms, BondType
from accel.base.mols import Mols, Mol
from accel.util.constants import Elements
from accel.util.log import logger
from accel.util.matrix import Matrix


def aromatize(_c: Mol, single_threshold: int = 1.01, double_threshold: int = 1.01, max_depth: int = 18):
    _atoms = _c.atoms.to_list()
    npxyz = [[_a.x, _a.y, _a.z] for _a in _atoms]
    npdist_mat = np.expand_dims(npxyz, axis=1) - np.expand_dims(npxyz, axis=0)
    npdist_mat = np.sqrt(np.sum(npdist_mat**2, axis=-1))

    npcov = [Elements.get_element(_a.symbol)["single"] for _a in _atoms]
    npcov = [_v if _v is not None else np.nan for _v in npcov]
    npcov_mat = np.expand_dims(npcov, axis=1) + np.expand_dims(npcov, axis=0)
    npcov_mat = np.signbit(npdist_mat - single_threshold * npcov_mat)
    # single_threshold should be like 0.9 to reduce calculation cost

    npdbl = [Elements.get_element(_a.symbol)["double"] for _a in _atoms]
    npdbl = [_v if _v is not None else np.nan for _v in npdbl]
    npdbl_mat = np.expand_dims(npdbl, axis=1) + np.expand_dims(npdbl, axis=0)
    npdbl_mat = np.signbit((double_threshold * npdbl_mat) - npdist_mat)
    # double_threshold should be like 1.1 to reduce calculation cost

    aromatic_mat: List[List[bool]] = (npcov_mat & npdbl_mat).tolist()

    def path_to_other(start_idx: int, end_idx: int, route: List[int], ring_list: List[List[int]]):
        route = route + [start_idx]
        for other_idx in [_idx for _idx, _flag in enumerate(aromatic_mat[start_idx]) if _flag is True]:
            if other_idx == end_idx and other_idx != route[-2]:
                ring_list.append(route)
                continue
            if other_idx in route:
                continue
            if len(route) > max_depth:
                continue
            path_to_other(other_idx, end_idx, route, ring_list)

    aromatic_atom = [False for _ in range(len(_atoms))]
    aromatic_rings = []
    for cheking_idx in range(len(_atoms)):
        if aromatic_atom[cheking_idx]:
            continue
        ring_list = []
        path_to_other(cheking_idx, cheking_idx, [], ring_list)
        if len(ring_list) == 0:
            continue
        for ring in ring_list:
            electrons = 0
            for atom_idx in ring:
                _symbol = _c.atoms[atom_idx].symbol
                if _symbol == "C":
                    electrons += 1
                elif _symbol == "O":
                    electrons += 2
                elif _symbol == "S":
                    electrons += 2
                elif _symbol == "N":
                    if len(_c.atoms[atom_idx].bonds) == 3:
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
        _c.atoms.bonds[id_a + 1, id_b + 1] = BondType.aromatic


def embed_bonds(_c: Mol, cov_scaling=1.1, vdw_scaling=1.0, double_scaling=1.05, triple_scaling=1.05):
    _atoms = _c.atoms.to_list()
    npxyz = [[_a.x, _a.y, _a.z] for _a in _atoms]
    npdist_mat = np.expand_dims(npxyz, axis=1) - np.expand_dims(npxyz, axis=0)
    npdist_mat = np.sqrt(np.sum(npdist_mat**2, axis=-1))

    npcov = [Elements.get_element(_a.symbol)["single"] for _a in _atoms]
    npcov = [_v if _v is not None else np.nan for _v in npcov]
    npcov_mat = np.expand_dims(npcov, axis=1) + np.expand_dims(npcov, axis=0)
    npcov_mat = np.signbit(npdist_mat - cov_scaling * npcov_mat)

    npvdw = [Elements.get_element(_a.symbol)["vdw"] for _a in _atoms]
    npvdw = [_v if _v is not None else np.nan for _v in npvdw]
    npvdw_mat = np.expand_dims(npvdw, axis=1) + np.expand_dims(npvdw, axis=0)
    npvdw_mat = np.signbit(npdist_mat - vdw_scaling * npvdw_mat) ^ npcov_mat

    npdbl = [Elements.get_element(_a.symbol)["double"] for _a in _atoms]
    npdbl = [_v if _v is not None else np.nan for _v in npdbl]
    npdbl_mat = np.expand_dims(npdbl, axis=1) + np.expand_dims(npdbl, axis=0)
    npdbl_mat = np.signbit(npdist_mat - double_scaling * npdbl_mat)

    nptri = [Elements.get_element(_a.symbol)["triple"] for _a in _atoms]
    nptri = [_v if _v is not None else np.nan for _v in nptri]
    nptri_mat = np.expand_dims(nptri, axis=1) + np.expand_dims(nptri, axis=0)
    nptri_mat = np.signbit(npdist_mat - triple_scaling * nptri_mat)

    npsgl_mat = npcov_mat ^ npdbl_mat
    npdbl_mat = npdbl_mat ^ nptri_mat

    mat_dict: Dict[int, np.ndarray] = {
        BondType.single: npsgl_mat,
        BondType.double: npdbl_mat,
        BondType.triple: nptri_mat,
        BondType.contact: npvdw_mat,
    }

    _c.atoms.init_bonds()

    for _tyep, _mat in mat_dict.items():
        for index_a in range(len(_atoms)):
            for index_b in range(len(_atoms)):
                if index_a >= index_b:
                    continue
                if _mat[index_a, index_b]:
                    _c.atoms.bonds[index_a + 1, index_b + 1] = _tyep

    logger.debug(f"bonding information of {_c.name} was embeded")


def embed_symm(_c: Mol):
    def _reset_visited_flags(_c: Mol):
        for _a in _c.atoms:
            _a.cache["unvisited"] = True
            _a.cache["ref_unvisited"] = True
        return None

    def _del_visited_flags(_c: Mol):
        for _a in _c.atoms:
            del _a.cache["unvisited"]
            del _a.cache["ref_unvisited"]
        return None

    _reset_visited_flags(_c)
    for _a in _c.atoms:
        # generate tree by BFS
        atom_trees = []
        for root in _a.bonds:
            _a.cache["unvisited"] = False
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
            _reset_visited_flags(_c)
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
                    candidate_pair.append({"ref_root": tree_ref["root_atom"], "tar_root": tree_target["root_atom"]})
        # logger.debug('generated candidates of atom mapping on {}'.format(_a))

        # check candidates by DFS
        def recursive_dfs(node_ref: Atom, node_tar: Atom, pair_atoms: list):
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
                    if recursive_dfs(_nb_ref, _nb_tar, pair_atoms):
                        break
                else:
                    node_ref.cache["ref_unvisited"] = True
                    node_tar.cache["unvisited"] = True
                    return False
            pair_atoms.append([node_ref, node_tar])
            return True

        for _pair in candidate_pair:
            pair_atoms_list = []
            _a.cache["unvisited"] = False
            _a.cache["ref_unvisited"] = False
            if recursive_dfs(_pair["ref_root"], _pair["tar_root"], pair_atoms_list):
                _pair["pair_list"] = pair_atoms_list
            _reset_visited_flags(_c)

        # logger.debug('checked candidates of atom mapping on {}'.format(_a))

        pairs_list = [_pair for _pair in candidate_pair if "pair_list" in _pair]
        for _ps in pairs_list:
            logger.debug(
                "{}: {}: root {} and {}: map {}".format(
                    _c.name,
                    _a,
                    _ps["ref_root"],
                    _ps["tar_root"],
                    [[str(k) for k in i] for i in _ps["pair_list"]],
                )
            )

        # generate matrix
        for _ps in pairs_list:
            _mat = np.identity(len(_c.atoms))
            for ref_atom, tar_atom in _ps["pair_list"]:
                _ref_idx = ref_atom.number - 1
                _tar_idx = tar_atom.number - 1
                _mat[[_ref_idx, _tar_idx]] = _mat[[_tar_idx, _ref_idx]]
            _ps["pair_matrix"] = _mat

        if len(pairs_list) != 0:
            _a.cache["isomeric_subs_list"] = [
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
    for _a in _c.atoms:
        if "isomeric_subs_list" not in _a.cache:
            continue
        # extract rotamers
        if [len(_a.bonds), len(_a.cache["isomeric_subs_list"])] in [[2, 1], [3, 1]]:
            # aromatic handling should be included
            if len(_a.double) != 0:
                for double_bond_atom in _a.double:
                    if double_bond_atom in [
                        _a.cache["isomeric_subs_list"][0]["root_a"],
                        _a.cache["isomeric_subs_list"][0]["root_a"],
                    ]:
                        break
                else:
                    numisomer_mat_list.append(_a.cache["isomeric_subs_list"][0]["matrix"])
                    continue
            rot_mat = _a.cache["isomeric_subs_list"][0]["matrix"]
            rot_mat = rot_mat.astype(int)
            rotamer_mat_list.append(rot_mat)
        elif [len(_a.bonds), len(_a.cache["isomeric_subs_list"])] in [[4, 3]]:
            for _ps in _a.cache["isomeric_subs_list"][1:]:
                rot_mat = np.dot(_ps["matrix"], _a.cache["isomeric_subs_list"][0]["matrix"])
                rot_mat = rot_mat.astype(int)
                rotamer_mat_list.append(rot_mat)
            for _ps in _a.cache["isomeric_subs_list"]:
                numisomer_mat_list.append(_ps["matrix"])
        else:
            for _ps in _a.cache["isomeric_subs_list"]:
                numisomer_mat_list.append(_ps["matrix"])

    # generate identical collections of rot_mat
    rotamer_mat_list = {_m.tobytes(): _m for _m in rotamer_mat_list}.values()
    numisomer_mat_list = {_m.tobytes(): _m for _m in numisomer_mat_list}.values()

    # check chirality here
    def _invalid_chirality(original_conf: Mol, modified_xyz: np.ndarray) -> int:
        invalid_count = 0
        original_xyz = np.array([a.xyz for a in original_conf.atoms])
        for _a in original_conf.atoms:
            if len(_a.bonds) != 4:
                continue
            _b_indexs = sorted(a.number - 1 for a in _a.bonds)
            original_matirix = original_xyz[_b_indexs[1:]] - original_xyz[_b_indexs[0]]
            modified_matirix = modified_xyz[_b_indexs[1:]] - modified_xyz[_b_indexs[0]]
            if (np.linalg.det(original_matirix) > 0) != (np.linalg.det(modified_matirix) > 0):
                invalid_count += 1
        return invalid_count

    def _cyclic_chiral_check(_c: Mol, rotamer_mat_list, numisomer_mat_list):
        unchanged_flag = True
        for rot_mat in rotamer_mat_list:
            original_xyz = np.array([at.xyz for at in _c.atoms])
            number_of_invalid = _invalid_chirality(_c, np.dot(rot_mat, original_xyz))
            if number_of_invalid == 0:
                continue
            for num_mat in numisomer_mat_list:
                t_mat = np.dot(num_mat, rot_mat)
                if _invalid_chirality(_c, np.dot(t_mat, original_xyz)) < number_of_invalid:
                    rot_mat = t_mat
                    unchanged_flag = False
                    logger.debug(f"chirality in rotamer matrix was partially corrected: {rot_mat}")
                    break
            else:
                logger.debug(f"chirality in rotamer matrix temporarily corrected: {rot_mat}")
        return unchanged_flag

    if len(rotamer_mat_list) != 0:
        for _ in range(len(rotamer_mat_list)):
            if _cyclic_chiral_check(_c, rotamer_mat_list, numisomer_mat_list):
                logger.debug(f"{_c.name}: chirality in rotamer matrix was successfully corrected")
                break
        else:
            logger.error(f"{_c.name}: chirality in rotamer matrix was not completely corrected")

    # generate identical collections of rot_mat again
    rotamer_mat_list = {_m.tobytes(): _m for _m in rotamer_mat_list}.values()
    numisomer_mat_list = {_m.tobytes(): _m for _m in numisomer_mat_list}.values()

    # added matrix property
    _c.data["rotamer"] = [Matrix(_c.atoms.to_list()).bind(_m) for _m in rotamer_mat_list]
    _c.data["numisomer"] = [Matrix(_c.atoms.to_list()).bind(_m) for _m in numisomer_mat_list]

    _del_visited_flags(_c)


def rmsd_pruning(
    confs: Mols,
    rmsd_threshold: float = 0.01,
    for_all: bool = False,
    redundant_check: int = 3,
    all_perturbation: bool = True,
    include_numeric_isomers: bool = False,
):
    def kabsch_rmsd(ref_xyzs: np.ndarray, tar_xyzs: np.ndarray):
        u_mat, _, v_tmat = np.linalg.svd(np.dot(tar_xyzs.T, ref_xyzs))
        if (np.linalg.det(u_mat) * np.linalg.det(v_tmat)) < 0.0:
            u_mat[:, -1] = -u_mat[:, -1]
        r_mat = np.dot(u_mat, v_tmat)
        diff_xyzs = ref_xyzs - (np.dot(tar_xyzs, r_mat))
        return np.sqrt((diff_xyzs * diff_xyzs).sum() / len(diff_xyzs))

    def cyclic_rmsd(min_rmsd, max_cycle, rot_mats_dicts, init_rot_mat, np_ref_xyzs, np_tar_xyzs):
        prev_min_rmsd = min_rmsd
        min_rot_mat = init_rot_mat
        for _cycle in range(max_cycle):
            min_num_rotamer = -1
            for _i, _r in enumerate(rot_mats_dicts):
                _r["flag"] = not _r["flag"]
                _rot_mat = init_rot_mat
                for _true_mat in [_rot["matrix"] for _rot in rot_mats_dicts if _rot["flag"] is True]:
                    _rot_mat = np.dot(_true_mat, _rot_mat)

                _rmsd = kabsch_rmsd(np_ref_xyzs, np.dot(_rot_mat, np_tar_xyzs))
                if _rmsd < prev_min_rmsd:
                    prev_min_rmsd = _rmsd
                    min_num_rotamer = _i
                    min_rot_mat = _rot_mat
                _r["flag"] = not _r["flag"]
            if min_num_rotamer == -1:
                min_rmsd = prev_min_rmsd
                # logger.debug(f"end at cycle {_cycle + 1}")
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

    def cal_sym_rmsd(
        ref_conf: Mol,
        tar_conf: Mol,
        all_perturbation: bool,
        max_cycle: int,
        numisomer_swap: bool,
    ):
        np_ref_xyzs = np.array([_a.xyz for _a in ref_conf.atoms])
        np_tar_xyzs = np.array([_a.xyz for _a in tar_conf.atoms])
        xyz_length = len(np_ref_xyzs)
        if xyz_length != len(np_tar_xyzs):
            raise Exception
        np_ref_xyzs = np_ref_xyzs - np.mean(np_ref_xyzs, axis=0)
        np_tar_xyzs = np_tar_xyzs - np.mean(np_tar_xyzs, axis=0)
        min_rmsd = kabsch_rmsd(np_ref_xyzs, np_tar_xyzs)
        min_rot_mat = np.identity(xyz_length)
        logger.debug(f"first RMSD of {ref_conf.name} & {tar_conf.name}: {min_rmsd:.10f}")
        tar_conf_rotamers: List[Matrix] = tar_conf.data["rotamer"]
        tar_conf_atoms_list = tar_conf.atoms.to_list()
        _rotamers = [
            {"flag": False, "matrix": _m.ordered(tar_conf_atoms_list).to_ndarray()} for _m in tar_conf_rotamers
        ]
        if all_perturbation:
            for _i in range(len(_rotamers)):
                for _comb in itertools.combinations([_r["matrix"] for _r in _rotamers], _i + 1):
                    _rot_mat = np.identity(xyz_length)
                    for _true_mat in _comb:
                        _rot_mat = np.dot(_true_mat, _rot_mat)
                    _rmsd = kabsch_rmsd(np_ref_xyzs, np.dot(_rot_mat, np_tar_xyzs))
                    if _rmsd < min_rmsd:
                        min_rmsd = _rmsd
                        min_rot_mat = _rot_mat
                        logger.debug(f"minimum RMSD was updated: {min_rmsd:.10f}")
                        # logger.debug('rot_mat was {}'.format(_rot_mat.astype(int).tolist()))
        else:
            min_rmsd, min_rot_mat = cyclic_rmsd(
                min_rmsd,
                max_cycle,
                _rotamers,
                np.identity(xyz_length),
                np_ref_xyzs,
                np_tar_xyzs,
            )
        if numisomer_swap:
            tar_conf_numisomers: List[Matrix] = tar_conf.data["numisomer"]
            tar_conf_atoms_list = tar_conf.atoms.to_list()
            _numisomers = [
                {"flag": False, "matrix": _m.ordered(tar_conf_atoms_list).to_ndarray()} for _m in tar_conf_numisomers
            ]
            min_rmsd, _ = cyclic_rmsd(min_rmsd, max_cycle, _numisomers, min_rot_mat, np_ref_xyzs, np_tar_xyzs)
        return min_rmsd

    for label in confs.labels.keys():
        try:
            _molcfs = sorted(confs.labels[label], key=lambda t: t.energy)
            _for_all = for_all
        except TypeError:
            _molcfs = confs.labels[label]
            logger.error(f"sorting failed in {label}: for_all is flagged")
            _for_all = True

        for i in range(len(_molcfs)):
            if _molcfs[i].flag:
                redundant_counter = 0
                for q in range(i + 1, len(_molcfs)):
                    if _molcfs[q].flag:
                        _rms = cal_sym_rmsd(
                            _molcfs[i],
                            _molcfs[q],
                            all_perturbation=all_perturbation,
                            max_cycle=1024,
                            numisomer_swap=include_numeric_isomers,
                        )
                        logger.debug(f"Final RMSD of {_molcfs[i].name} & {_molcfs[q].name}: {_rms:.10f}")
                        if _rms <= rmsd_threshold:
                            _molcfs[q].deactivate("rmsd_limit")
                        elif not _for_all:
                            if redundant_counter >= redundant_check:
                                logger.debug(f"counter reached {redundant_counter}: moving to next conf")
                                break
                            else:
                                redundant_counter += 1


# not coded yet
def parse_angles(_c: Mol):
    bonds = _c.bonds_list
    ang_list = []
    for i, bond in enumerate(bonds):
        for atom in bond:
            bonds_poped = bonds.copy()
            bonds_poped.pop(i)
            next_atoms = set(sum((b for b in bonds_poped if atom in b), []))
            next_atoms.discard(atom)
            next_atoms = list(next_atoms)
            bond_i = bonds[i].copy()
            bond_i.remove(atom)
            ang_list.extend([[bond_i[0], atom, next_atom] for next_atom in next_atoms])
    for ang in ang_list:
        inv_ang = [ang[i] for i in reversed(range(3))]
        ang_list.remove(inv_ang)
    ang_list.sort()
    _c.data["angles_list"] = ang_list


# not coded yet
def parse_dihedrals(_c: Mol):
    bonds = _c.bonds_list
    dih_list = []
    for i, bond in enumerate(bonds):
        for atom in bond:
            bonds_poped = bonds.copy()
            bonds_poped.pop(i)
            next_atoms = set(sum((b for b in bonds_poped if atom in b), []))
            next_atoms.discard(atom)
            next_atoms = list(next_atoms)
            bond_i = bonds[i].copy()
            bond_i.remove(atom)
            for next_atom in next_atoms:
                next_bonds = [b for b in bonds if (next_atom in b) and (atom not in b)]
                for next_bond in next_bonds:
                    next_bond_copied = next_bond.copy()
                    next_bond_copied.remove(next_atom)
                    dih_list.append([bond_i[0], atom, next_atom, next_bond_copied[0]])
    for dih in dih_list:
        inv_dih = [dih[i] for i in reversed(range(4))]
        dih_list.remove(inv_dih)
    dih_list.sort()
    _c.data["dihedrals_list"] = dih_list


def map_numbers(confs: Mols, reference_confs: Mols):
    def _get_dihedral(atom_a: Atom, atom_b: Atom, atom_c: Atom, atom_d: Atom):
        _va = np.array(atom_a.xyz)
        _vb = np.array(atom_b.xyz)
        _vc = np.array(atom_c.xyz)
        _vd = np.array(atom_d.xyz)
        _vab = _va - _vb
        _vcb = _vc - _vb
        _vdc = _vd - _vc
        _pvac = np.cross(_vab, _vcb)
        _pvbd = np.cross(_vdc, _vcb)
        _dac = np.linalg.norm(_pvac)
        _dbd = np.linalg.norm(_pvbd)
        _angle = np.arccos(np.sum(_pvac * _pvbd) / (_dac * _dbd))
        if np.sum(_pvac * np.cross(_pvbd, _vcb)) < 0:
            _angle = -_angle
        _angle = float(np.rad2deg(_angle))
        # logger.debug('dihedral angle {}{}-{}{}-{}{}-{}{}: {} degree'.format(
        #     atom_a.symbol, atom_a.number, atom_b.symbol, atom_b.number, atom_c.symbol, atom_c.number,
        #     atom_d.symbol, atom_d.number, _angle))
        return _angle

    def _det_rs(subs: List[int], test_atoms: Atoms):
        neighbors_xyz = [test_atoms.get(_num).xyz for _num in sorted(subs)]
        neighbors_xyz_np = np.array(neighbors_xyz[1:])
        neighbors_xyz_np = neighbors_xyz_np - np.array(neighbors_xyz[0])
        if np.linalg.det(neighbors_xyz_np) > 0:
            flag = True
        else:
            flag = False
        return flag

    def _get_dihedral_nums(atom: Atom, isomeric_atom_x: Atom, isomeric_atom_y: Atom):
        atom_a: Atom = sorted([isomeric_atom_x, isomeric_atom_y], key=lambda t: t.number)[0]
        atom_b: Atom = atom
        atom_c: Atom = [_a for _a in atom_b.bonds if _a not in [isomeric_atom_x, isomeric_atom_y]][0]
        if len(atom_c.bonds) == 1:
            raise ValueError
        atom_d: Atom = sorted((_a for _a in atom_c.bonds if _a != atom_b), key=lambda t: t.number)[0]
        return [
            atom_a.number,
            atom_b.number,
            atom_c.number,
            atom_d.number,
        ]

    def _det_ez(subs: List[int], test_atoms: Atoms):
        _dihed = _get_dihedral(
            test_atoms.get(subs[0]),
            test_atoms.get(subs[1]),
            test_atoms.get(subs[2]),
            test_atoms.get(subs[3]),
        )
        if _dihed > 90 or _dihed < -90:
            _ret = True
        else:
            _ret = False
        return _ret

    for label in confs.labels:
        reference_c = reference_confs.labels[label].get()
        problematic_rs: List[Atom] = []
        problematic_ez: List[Atom] = []
        problematic_allene: List[Atom] = []
        for _a in reference_c.atoms:
            if "isomeric_subs_list" in _a.cache:
                if len(_a.bonds) == 4 and len(_a.cache["isomeric_subs_list"]) in (1, 2, 3):
                    _a.data["num_configuration"] = _det_rs([a.number for a in _a.bonds], reference_c.atoms)
                    problematic_rs.append(_a)
                if len(_a.bonds) == 3 and len(_a.double) == 1 and len(_a.cache["isomeric_subs_list"]) == 1:
                    if len(_a.double[0].double) == 2 and len(_a.double[0].bonds) == 2:
                        _end_atom = [_ea for _ea in _a.double[0].double if _ea is not _a][0]
                        _num_list = [a.number for a in _a.single] + [a.number for a in _end_atom.single]
                        _num_list.extend([a.number for a in _a.aromatic] + [a.number for a in _end_atom.aromatic])
                        if len(_num_list) != 4:
                            continue
                        _a.data["num_configuration"] = _det_rs(_num_list, reference_c.atoms)
                        problematic_allene.append(_a)
                        continue
                    try:
                        _dhed_nums = _get_dihedral_nums(
                            _a,
                            _a.cache["isomeric_subs_list"][0]["root_a"],
                            _a.cache["isomeric_subs_list"][0]["root_b"],
                        )
                        _ez = _det_ez(_dhed_nums, reference_c.atoms)
                    except ValueError:
                        continue
                    _a.data["num_configuration"] = _ez
                    problematic_ez.append(_a)
        logger.info(
            "{} mappings were detected at the sp3 centers {}".format(
                len(problematic_rs), [str(_a) for _a in problematic_rs]
            )
        )
        logger.info(
            "{} mappings were detected at the sp2 centers {}".format(
                len(problematic_ez), [str(_a) for _a in problematic_ez]
            )
        )
        if len(problematic_allene) != 0:
            logger.info(
                "{} mappings were detected at the allene chiral centers {}".format(
                    len(problematic_allene), [str(_a) for _a in problematic_allene]
                )
            )

        def resolve_mapping(
            _c: Mol, problematic_rs: List[Atom], problematic_ez: List[Atom], problematic_allene: List[Atom]
        ):
            for ref_atom in problematic_rs:
                if "valid_numeric_chirality" in _c.atoms.get(ref_atom.number).cache:
                    if _c.atoms.get(ref_atom.number).cache["valid_numeric_chirality"] in [True, False]:
                        continue
                if ref_atom.data["num_configuration"] != _det_rs([a.number for a in ref_atom.bonds], _c.atoms):
                    # [TASK] reverse flag of isomeric_subs_list will be inserted
                    for iso in ref_atom.cache["isomeric_subs_list"]:
                        _c.atoms.swap(iso["root_a"].number, iso["root_b"].number, bonds=False)
                        return_flag = ref_atom.data["num_configuration"] == _det_rs(
                            [a.number for a in ref_atom.bonds], _c.atoms
                        )
                        _c.atoms.swap(iso["root_a"].number, iso["root_b"].number, bonds=False)
                        if return_flag:
                            for pair in iso["list"]:
                                _c.atoms.swap(pair[0].number, pair[1].number)
                                logger.info(f"{_c.name}: atoms {pair[0].number} and {pair[1].number} were swaped")
                            _c.atoms.get(ref_atom.number).cache["valid_numeric_chirality"] = True
                            return False
                    logger.error(f"atom mapping of {ref_atom.symbol}{ref_atom.number} was not resolved in {_c.name}")
                    _c.atoms.get(ref_atom.number).cache["valid_numeric_chirality"] = False

            for ref_atom in problematic_ez:
                if "valid_numeric_chirality" in _c.atoms.get(ref_atom.number).cache:
                    if _c.atoms.get(ref_atom.number).cache["valid_numeric_chirality"] in [True, False]:
                        continue
                _dhed_nums = _get_dihedral_nums(
                    ref_atom,
                    ref_atom.cache["isomeric_subs_list"][0]["root_a"],
                    ref_atom.cache["isomeric_subs_list"][0]["root_b"],
                )
                if ref_atom.data["num_configuration"] != _det_ez(_dhed_nums, _c.atoms):
                    iso = ref_atom.cache["isomeric_subs_list"][0]
                    for pair in iso["list"]:
                        _c.atoms.swap(pair[0].number, pair[1].number, bonds=False)
                    return_flag = ref_atom.data["num_configuration"] == _det_ez(_dhed_nums, _c.atoms)
                    for pair in iso["list"]:
                        _c.atoms.swap(pair[0].number, pair[1].number, bonds=False)
                    if return_flag:
                        for pair in iso["list"]:
                            _c.atoms.swap(pair[0].number, pair[1].number)
                            logger.info(f"{_c.name}: atoms {pair[0].number} and {pair[1].number} were swaped")
                        _c.atoms.get(ref_atom.number).cache["valid_numeric_chirality"] = True
                        return False
                    logger.error(f"atom mapping of {ref_atom.symbol}{ref_atom.number} was not resolved in {_c.name}")
                    _c.atoms.get(ref_atom.number).cache["valid_numeric_chirality"] = False

            for ref_atom in problematic_allene:
                if "valid_numeric_chirality" in _c.atoms.get(ref_atom.number).cache:
                    if _c.atoms.get(ref_atom.number).cache["valid_numeric_chirality"] in [True, False]:
                        continue
                _end_atom = [_ea for _ea in ref_atom.double[0].double if _ea is not ref_atom][0]
                _num_list = [a.number for a in ref_atom.single] + [a.number for a in _end_atom.single]
                _num_list.extend([a.number for a in ref_atom.aromatic] + [a.number for a in _end_atom.aromatic])
                if len(_num_list) != 4:
                    logger.error("unexpected error handling allene chirality")
                    continue
                if ref_atom.data["num_configuration"] != _det_rs(_num_list, _c.atoms):
                    # [TASK] reverse flag of isomeric_subs_list will be inserted
                    for iso in ref_atom.cache["isomeric_subs_list"]:
                        _c.atoms.swap(iso["root_a"].number, iso["root_b"].number, bonds=False)
                        return_flag = ref_atom.data["num_configuration"] == _det_rs(_num_list, _c.atoms)
                        _c.atoms.swap(iso["root_a"].number, iso["root_b"].number, bonds=False)
                        if return_flag:
                            for pair in iso["list"]:
                                _c.atoms.swap(pair[0].number, pair[1].number)
                                logger.info(f"{_c.name}: atoms {pair[0].number} and {pair[1].number} were swaped")
                            _c.atoms.get(ref_atom.number).cache["valid_numeric_chirality"] = True
                            return False
                    logger.error(f"atom mapping of {ref_atom.symbol}{ref_atom.number} was not resolved in {_c.name}")
                    _c.atoms.get(ref_atom.number).cache["valid_numeric_chirality"] = False
            return True

        for _c in confs.labels[label]:
            for _ in range(len(problematic_rs) + len(problematic_ez) + len(problematic_allene) + 1):
                if resolve_mapping(_c, problematic_rs, problematic_ez, problematic_allene):
                    # [TASK] reset "valid_numeric_chirality" then check all of mapping again
                    logger.info(f"{_c.name}: successfully mapped")
                    break
            else:
                logger.error(f"{_c.name}: reached max iteration")


def order_by_cip(substitutions: List[Atom], roots: List[Atom]):
    for sub in substitutions:
        pass
    return []
