import itertools

import numpy as np
from accel.base.atoms import Atom, Atoms
from accel.base.systems import System, Systems
from accel.util.log import logger
from accel.util.matrix import Matrix


def rmsdpruning(
    confs: Systems,
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
        ref_conf: System,
        tar_conf: System,
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
        tar_conf_rotamers: list[Matrix] = tar_conf.data["rotamer"]
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
            tar_conf_numisomers: list[Matrix] = tar_conf.data["numisomer"]
            tar_conf_atoms_list = tar_conf.atoms.to_list()
            _numisomers = [
                {"flag": False, "matrix": _m.ordered(tar_conf_atoms_list).to_ndarray()} for _m in tar_conf_numisomers
            ]
            min_rmsd, _ = cyclic_rmsd(min_rmsd, max_cycle, _numisomers, min_rot_mat, np_ref_xyzs, np_tar_xyzs)
        return min_rmsd

    for label in confs.labels.keys():
        try:
            _Systemcfs = sorted(confs.labels[label], key=lambda t: t.energy)
            _for_all = for_all
        except TypeError:
            _Systemcfs = confs.labels[label]
            logger.error(f"sorting failed in {label}: for_all is flagged")
            _for_all = True

        for i in range(len(_Systemcfs)):
            if _Systemcfs[i].state:
                redundant_counter = 0
                for q in range(i + 1, len(_Systemcfs)):
                    if _Systemcfs[q].state:
                        _rms = cal_sym_rmsd(
                            _Systemcfs[i],
                            _Systemcfs[q],
                            all_perturbation=all_perturbation,
                            max_cycle=1024,
                            numisomer_swap=include_numeric_isomers,
                        )
                        logger.debug(f"Final RMSD of {_Systemcfs[i].name} & {_Systemcfs[q].name}: {_rms:.10f}")
                        if _rms <= rmsd_threshold:
                            _Systemcfs[q].deactivate("rmsd_limit")
                        elif not _for_all:
                            if redundant_counter >= redundant_check:
                                logger.debug(f"counter reached {redundant_counter}: moving to next conf")
                                break
                            else:
                                redundant_counter += 1


def map_numbers(confs: Systems, reference_confs: Systems):
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

    def _det_rs(subs: list[int], test_atoms: Atoms):
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

    def _det_ez(subs: list[int], test_atoms: Atoms):
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
        problematic_rs: list[Atom] = []
        problematic_ez: list[Atom] = []
        problematic_allene: list[Atom] = []
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
            _c: System, problematic_rs: list[Atom], problematic_ez: list[Atom], problematic_allene: list[Atom]
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


def set_chirality(c: System, center_index: int, sub_index: list[int]):
    if len(sub_index) != 4:
        raise ValueError
    else:
        sorted_index = sorted(sub_index)

    _sub_xyzs = np.array([c.atoms.get(i).xyz for i in sorted_index[1:]]) - np.array(
        [c.atoms.get(sorted_index[0]).xyz for _ in range(3)]
    )
    _ret = np.linalg.det(_sub_xyzs)
    if _ret > 0:
        _ret = 1
    elif _ret < 0:
        _ret = -1
    else:
        _ret = 0
    c.data[f"chiral_{center_index}_to_{sub_index}"] = _ret
