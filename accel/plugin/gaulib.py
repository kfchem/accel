import statistics
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import List

import numpy as np
from accel.base.boxcore import BoxCore
from accel.base.mols import Mol
from accel.base.selector import Selectors
from accel.util import Execmd, FileType
from accel.util.constants import Elements, Unit, Units
from accel.util.log import logger


# not coding yet
def read_multiple_xyz(mulcos: BoxCore):
    pass


def calc_sum_energy(mulcos: BoxCore):
    for _c in mulcos.pack():
        _c.data["g16_t_zero"] = _c.data["g16_zpc"] + _c.data["g16_scf"]
        _c.data["g16_t_energy"] = _c.data["g16_corr_to_energy"] + _c.data["g16_scf"]
        _c.data["g16_t_enthalpy"] = _c.data["g16_corr_to_enthalpy"] + _c.data["g16_scf"]
        _c.data["g16_t_gibbs"] = _c.data["g16_corr_to_gibbs"] + _c.data["g16_scf"]
        _c.energy = Units.hartree(_c.data["g16_t_gibbs"]).to_kcal_mol


def read_low_freqs(mulcos: BoxCore):
    for _c in mulcos.pack():
        with _c.path.open() as f:
            _ls = f.readlines()
        for i, _l in enumerate(_ls):
            if "1                      2                      3" in _l:
                _freq = [float(_ls[i + 2].split()[k]) for k in range(2, 5)]
                _c.data["g16_freq_1"] = _freq[0]
                _c.data["g16_freq_2"] = _freq[1]
                _c.data["g16_freq_3"] = _freq[2]


def read_vibration(mulcos: BoxCore, freq_number: int = 1):
    for _c in mulcos.pack():
        _freq = []
        with _c.path.open() as f:
            _ls = f.readlines()
        for i, _l in enumerate(_ls):
            if "Frequencies --" in _l:
                no_row = [int(s) for s in _ls[i - 2].split()]
                if freq_number in no_row:
                    column_no = no_row.index(freq_number)
                    for k, _l2 in enumerate(_ls[i:]):
                        if "Atom" in _l2:
                            row_no = i + k + 1
                            break
        for _l in _ls[row_no:]:
            if len(_l.split()) in [3, 0, 1]:
                break
            _freq.append(_l.split()[2 + column_no * 3 : 2 + (column_no + 1) * 3])
        _c.data["g16_vibration"] = _freq


def check_bonding(
    mulcos: BoxCore,
    bonding_atoms: List[int],
    additional_valid_range: int = 0,
    acceptable_invalid_atoms: int = 0,
):
    for _c in mulcos.pack():
        _dists = []
        for i, _v in enumerate(_c.data["g16_vibration"]):
            _dists.append([i + 1, np.linalg.norm(np.array(_v))])
        _rank = [_dist[0] for _dist in sorted(_dists, reverse=True, key=lambda x: x[1])]
        logger.debug(f"vib_rank: {_rank}")
        invalid_atoms = []
        for _atom in bonding_atoms:
            valid_range = len(bonding_atoms) + additional_valid_range
            _atom_rank = _rank.index(_atom) + 1
            logger.debug(f"the rank of atom {_atom} was {_atom_rank}: {_c.path.name}")
            if _atom_rank <= valid_range:
                logger.debug(f"the rank of atom {_atom} was in valid range")
            else:
                logger.info(f"the rank of atom {_atom} was in invalid range")
                invalid_atoms.append(_atom)
        if len(invalid_atoms) > acceptable_invalid_atoms:
            logger.info(f"atoms {invalid_atoms} are out of range")
            _c.flag = False


def read_nmr(mulcos: BoxCore):
    for _c in mulcos.pack():
        with _c.path.open() as f:
            _ls = f.readlines()
        for i, _l in enumerate(_ls):
            if "Anisotropy =" in _l:
                _lsp = _l.split()
                _a = _c.atoms.get(int(_lsp[0]))
                if _a.symbol != Elements.canonicalize(_lsp[1]):
                    logger.error("atomic symbol does not matched")
                    break
                else:
                    _a.data["isotropic"] = float(_l.split()[4])
                    _a.data["anisotropy"] = float(_l.split()[7])


def read_ecd(mulcos: BoxCore):
    for _c in mulcos.pack():
        with _c.path.open() as f:
            _ls = f.readlines()
        _stat_dict = defaultdict(dict)
        for i, _l in enumerate(_ls):
            if "Excited State" in _l:
                _lsp = _l.split()
                _d = _stat_dict[int(_lsp[2].replace(":", ""))]
                _d["shape"] = _lsp[3]
                _d["energy"] = float(_lsp[4])
                _d["length"] = float(_lsp[6])
                _d["f"] = float(_lsp[8][2:])
                _d["S2"] = float(_lsp[9][7:])
            if "R(velocity)" in _l:
                for _li in _ls[i + 1 :]:
                    if "rxdel" in _li:
                        break
                    _stat_dict[int(_li.split()[0])]["R_velocity"] = float(_li.split()[4])
            if "R(length)" in _l:
                for _li in _ls[i + 1 :]:
                    if "del" in _li:
                        break
                    _stat_dict[int(_li.split()[0])]["R_length"] = float(_li.split()[4])
        _c.data["ecd"] = dict(_stat_dict)


def read_spinspin(mulcos: BoxCore, key_value: str = "Total_J"):
    """
    key_value
    'FC_to_K': Fermi Contact (FC) contribution to K (Hz)
    'FC_to_J': Fermi Contact (FC) contribution to J (Hz)
    'SC_to_K': Spin-dipolar (SD) contribution to K (Hz)
    'SC_to_J': Spin-dipolar (SD) contribution to J (Hz)
    'PSO_to_K': Paramagnetic spin-orbit (PSO) contribution to K (Hz)
    'PSO_to_J': Paramagnetic spin-orbit (PSO) contribution to J (Hz)
    'DSO_to_K': Diamagnetic spin-orbit (DSO) contribution to K (Hz)
    'DSO_to_J': Diamagnetic spin-orbit (DSO) contribution to J (Hz)
    'Total_K': Total nuclear spin-spin coupling K (Hz)
    'Total_J': Total nuclear spin-spin coupling J (Hz)
    """
    _key_list = {
        "FC_to_K": "Fermi Contact (FC) contribution to K (Hz):",
        "FC_to_J": "Fermi Contact (FC) contribution to J (Hz):",
        "SC_to_K": "Spin-dipolar (SD) contribution to K (Hz):",
        "SC_to_J": "Spin-dipolar (SD) contribution to J (Hz):",
        "PSO_to_K": "Paramagnetic spin-orbit (PSO) contribution to K (Hz):",
        "PSO_to_J": "Paramagnetic spin-orbit (PSO) contribution to J (Hz):",
        "DSO_to_K": "Diamagnetic spin-orbit (DSO) contribution to K (Hz):",
        "DSO_to_J": "Diamagnetic spin-orbit (DSO) contribution to J (Hz):",
        "Total_K": "Total nuclear spin-spin coupling K (Hz):",
        "Total_J": "Total nuclear spin-spin coupling J (Hz):",
    }
    _key_str = _key_list[key_value]
    logger.info(f"{_key_str[:-1]} is used for extracting coupling data")
    for _c in mulcos.pack():
        with _c.path.open() as f:
            _ls = f.readlines()
        for i, _l in enumerate(_ls):
            if _key_str in _l:
                for _a in _c.atoms:
                    _a.data["coupling"] = {}
                col_idx = []
                for _ll in _ls[i + 1 :]:
                    if not _ll.split()[0].isdecimal():
                        break
                    elif _ll.startswith("            "):
                        col_idx = [int(k) for k in _ll.split()]
                    else:
                        for k, _jvalue in enumerate(_ll.split()[1:]):
                            _jvalue = float(_jvalue.replace("D", "E"))
                            _c.atoms.get(int(_ll.split()[0]) + 1).data["coupling"][col_idx[k]] = _jvalue
                            _c.atoms.get(col_idx[k] + 1).data["coupling"][int(_ll.split()[0])] = _jvalue


def check_resub(mulcos: BoxCore, input_suffix=".gjf", log_suffix=".log"):
    def _gen_gjf_from_final_stdorient(_p: Path):
        with _p.open() as f:
            _ls = f.readlines()
        geom_no = 0
        for i, _l in enumerate(_ls):
            if "Standard orientation:" in _l:
                geom_no = i
        if geom_no == 0:
            logger.error("could not find Standard orientation entry: " + _p.name)
            return False
        geom = []
        for _l in _ls[(geom_no + 5) :]:
            if "---------------------------------------------------------------------" in _l:
                break
            geom.append(_l.split())
        geom = [
            "{1:<2} {0[3]:>11} {0[4]:>11} {0[5]:>11}\n".format(_g, Elements.canonicalize(int(_g[1]))) for _g in geom
        ]

        with _p.with_suffix(input_suffix).open() as f:
            _gjf_ls = f.readlines()
        blank_line_index = [i for i, _gjf_l in enumerate(_gjf_ls) if _gjf_l == "\n"]
        _out_lines = _gjf_ls[0 : blank_line_index[1] + 2]
        _out_lines.extend(geom)
        _out_lines.extend(_gjf_ls[blank_line_index[2] :])

        _suffix = "_"
        for _ext in [input_suffix, log_suffix, ".chk"]:
            while _p.with_suffix(_ext + _suffix).exists():
                _suffix = _suffix + "_"

        for _ext in [input_suffix, log_suffix, ".chk"]:
            if _p.with_suffix(_ext).exists():
                _p.with_suffix(_ext).rename(_p.with_suffix(_ext + _suffix))
                logger.debug(f"{_p.with_suffix(_ext).name} was renamed as {_p.with_suffix(_ext + _suffix).name}")

        with _p.with_suffix(input_suffix).open("w") as f:
            f.writelines(_out_lines)
        logger.debug(f"{_p.with_suffix(input_suffix).name} was generated")
        return True

    class Gerr:
        def __init__(self, name="", statement="", resub=True, use_final_geom=True) -> None:
            self.name = name
            self.statement = statement
            self.resub = resub
            self.use_final_geom = use_final_geom

        def resolve(self, conf: Mol):
            logger.info(f"{_c.name}: {self.name} was detected")
            if self.resub:
                if self.use_final_geom:
                    if _gen_gjf_from_final_stdorient(conf.path.with_suffix(log_suffix)):
                        conf.data["resubmission"] = True
                    else:
                        return False
                else:
                    conf.data["resubmission"] = True
            return True

    no_log_err = Gerr("no log file error", use_final_geom=False)
    unknown_err = Gerr("unknown error", use_final_geom=False)

    link9999_err = Gerr("link9999 error", statement="Error termination request processed by link 9999")
    link9999_wo_vib = Gerr("link9999 (not vibrational) error", resub=False)
    link9999_w_vib = Gerr("link9999 (vibrational ending) error")

    normal_errs = [
        Gerr(
            "floating point exception",
            statement="Error: floating point exception, integer divide by zero",
            use_final_geom=False,
        ),
        Gerr(
            "memory allocation error",
            statement="could not allocate memory.: Cannot allocate memory",
            use_final_geom=False,
        ),
        Gerr("PCMMkU failure", statement="Inv3 failed in PCMMkU", resub=False),
        Gerr("SCF convergence failure", statement="SCF has not converged", resub=False),
        Gerr("FormBX problem", statement="FormBX had a problem"),
    ]

    for _c in mulcos.pack():
        _c.data["resubmission"] = False
        with _c.path.with_suffix(input_suffix).open("r") as f:
            gjf_ls = f.readlines()
        tot_jobs = 1
        for i, gjf_l in enumerate(gjf_ls):
            if "--link1--" in gjf_l.lower():
                tot_jobs += 1
            if "freq" in gjf_l.lower() and "opt" in gjf_l.lower():
                tot_jobs += 1

        log_p = _c.path.with_suffix(log_suffix)
        if not log_p.exists():
            no_log_err.resolve(_c)
            continue

        norm_term = 0
        err_flag = False

        with log_p.open("r") as f:
            log_ls = f.readlines()

        for log_l in reversed(log_ls):
            if "Normal termination" in log_l:
                norm_term += 1
                continue
            if norm_term > 0:
                continue
            if link9999_err.statement in log_l:
                rms_force_vs = []
                for tmp_log_l in log_ls:
                    if "RMS     Force" in tmp_log_l:
                        try:
                            rms_force_vs.append(float(tmp_log_l.split()[2]))
                        except ValueError:
                            logger.error(f"ValueError during collecting RMS forces:{tmp_log_l}")
                            rms_force_vs = []
                            break
                if len(rms_force_vs) < 12:
                    logger.error("not enough force entries (less than 12) to detect vibration")
                    link9999_err.resolve(_c)
                elif (statistics.stdev(rms_force_vs[-12::2]) < 0.000001) or (
                    statistics.stdev(rms_force_vs[-12::3]) < 0.000001
                ):
                    link9999_w_vib.resolve(_c)
                else:
                    link9999_wo_vib.resolve(_c)
                err_flag = True
                break
            for _e in normal_errs:
                if _e.statement in log_l:
                    _e.resolve(_c)
                    err_flag = True
                    break
            else:
                continue
            break

        if not err_flag:
            if norm_term == tot_jobs:
                logger.debug(f"{_c.name} terminated normally")
            else:
                unknown_err.resolve(_c)
                logger.error(f"{_c.name}: unknown termination was detected")

    if len(mulcos.pack().has_data("resubmission", True)) == 0:
        logger.info("There is no resubmission")


@FileType.add("app/g16/input", 30)
def is_g16_input(_p: Path) -> bool:
    if _p.suffix not in (".gjf", ".com", ".inp"):
        return False

    route_pass = False
    with _p.open() as _f:
        for _i, _l in enumerate(_f):
            if _l.startswith("%") and not route_pass:
                continue
            elif _l.startswith("#"):
                route_pass = True
                continue
            elif _l == "\n" and route_pass:
                return True
            if not route_pass:
                return False
            if _i > 20:
                break
    return False


@FileType.add("app/g16/output", 50)
def is_g16_output(_p: Path) -> bool:
    if _p.suffix not in (".log", ".out", ".g16", ".g09", ".g98"):
        return False
    with _p.open() as _f:
        for _i, _l in enumerate(_f):
            if "Entering Gaussian System" in _l:
                return True
            if _i > 100:
                break
    return False


def run(_c: Mol):
    _cmd = [Execmd.get("g16"), str(_c.path.resolve().absolute())]
    try:
        logger.info(f"running: {_c.name}: {''.join(_cmd)}")
        _proc = subprocess.run(_cmd, cwd=str(_c.path.parent), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        _out = _proc.stdout.decode("utf-8").split("\n")
        logger.info(f"finished: {_c.name}: {_out}")
    except subprocess.CalledProcessError:
        _c.flag = False
        logger.error(f"failed: {_c.name}: {''.join(_cmd)}")


def submit(_c: Mol):
    _cmd = [Execmd.get("g16"), str(_c.path.resolve().absolute())]
    try:
        subprocess.Popen(_cmd, cwd=str(_c.path.parent), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        logger.info(f"submited: {_c.name}: {''.join(_cmd)}")
    except subprocess.CalledProcessError:
        _c.flag = False
        logger.error(f"failed: {_c.name}: {''.join(_cmd)}")


class GauBox(BoxCore):
    @Selectors.read_atoms.add("app/g16/output")
    def read_atoms(self, source: str = "archive"):
        """
        source = 'archive' | 'standard_orientation' | 'input_orientation'
        """
        for _c in self.pack():
            with _c.path.open() as f:
                _ls = f.readlines()

            if source == "archive":
                arc_sec = ""
                for i, _l in enumerate(_ls):
                    if "1\\1\\" in _l:
                        arc_sec = ""
                        for _l in _ls[i:]:
                            arc_sec += "".join(_l[1:].splitlines())
                            if "\\\\@" in _l:
                                break
                if arc_sec == "":
                    logger.error(f"could not find archive section in {_c.path.name}")
                    _c.deactivate("read_atoms: g16")
                    continue
                arc_sec = arc_sec.split("\\\\")
                _axyz = [i.split(",") for i in arc_sec[3].split("\\")[1:]]
                if len(_axyz[0]) == 5:
                    for _l in _axyz:
                        _l.pop(1)
                elif len(_axyz[0]) == 1:
                    logger.error(f"could not resolve archive section in {_c.path.name}")
                    _c.deactivate("read_atoms: g16")
                    continue

            elif source in ("standard_orientation", "input_orientation"):
                key_words = {
                    "standard_orientation": "Standard orientation:",
                    "input_orientation": "Input orientation:",
                }
                line_number = 0
                for i, _l in enumerate(_ls):
                    if key_words[source] in _l:
                        line_number = i
                if line_number == 0:
                    logger.error(f"could not find orientation section in {_c.path.name}")
                    _c.deactivate("read_atoms: g16")
                    continue
                _axyz = []
                for _l in _ls[(line_number + 5) :]:
                    if "---------------------------------------------------------------------" in _l:
                        break
                    _axyz.append(_l.split())
                _axyz = [[Elements.canonicalize(int(_l[1])), _l[3], _l[4], _l[5]] for _l in _axyz]

            _c.atoms.clear()
            for _l in _axyz:
                _c.atoms.append(_l)
            logger.debug(f"read xyz of {_c.path.name} from {source}")
        logger.debug(f"done: {len(self.pack())}/{len(self.mols)} confomers")
        return self

    @Selectors.read_atoms.add("app/g16/input")
    def read_atoms_from_input(self):
        for _c in self.pack():
            with _c.path.open() as f:
                _ls = f.readlines()
            blank_line_count = 0
            atomic_coord = []
            for i, _l in enumerate(_ls):
                if blank_line_count == 2:
                    atomic_coord.append(_l.replace(",", " ").replace("/", " ").replace("\t", " ").split())
                elif blank_line_count == 3:
                    break
                if _l == "\n":
                    blank_line_count += 1
            chg_mult = atomic_coord.pop(0)
            atomic_coord.pop()
            try:
                if len(chg_mult) != 2:
                    raise ValueError
                _c.total_charge = int(chg_mult[0])
                _c.multiplicity = int(chg_mult[1])
                uni_len = list({len(_l) for _l in atomic_coord})
                if len(uni_len) != 1 or uni_len[0] != 4:
                    logger.error("currently cannot read Z-matrix")
                    raise ValueError
                _c.atoms.clear()
                _c.atoms.extend(atomic_coord)
            except ValueError:
                _c.deactivate("read_atoms: g16")
                logger.error(f"could not read coordinate of {_c.path.name}")
                continue
            else:
                logger.debug(f"read atoms of {_c.path.name}")
        logger.debug(f"done: {len(self.pack())}/{len(self.mols)} confomers")
        return self

    @Selectors.read_energy.add("app/g16/output")
    def read_scf(self, calculation_method: str = None):
        for _c in self.pack():
            if calculation_method is None:
                scf_key = "SCF Done:"
            else:
                scf_key = "SCF Done:  E(" + calculation_method + ")"
                logger.info(f"{scf_key} is used for search SCF")
            with _c.path.open() as f:
                _ls = f.readlines()
            _n = 0
            for i, _l in enumerate(_ls):
                if scf_key in _l:
                    _n = i
            if _n == 0:
                logger.error("could not find the scf entry")
                _c.deactivate("read_scf")
                continue
            _c.data["g16_scf"] = float(_ls[_n].split()[4])
            _c.energy = Units.hartree(_c.data["g16_scf"]).to_kcal_mol
        logger.debug(f"done: {len(self.pack())}/{len(self.mols)} confomers")
        return self

    @Selectors.read_correction.add("app/g16/output")
    def read_correction(self):
        for _c in self.pack():
            with _c.path.open() as f:
                _ls = f.readlines()
            for i, _l in enumerate(_ls):
                if "Zero-point correction" in _l:
                    _c.data["g16_zpc"] = float(_ls[i].split()[2])
                if "Thermal correction to Energy" in _l:
                    _c.data["g16_corr_to_energy"] = float(_ls[i].split()[4])
                if "Thermal correction to Enthalpy" in _l:
                    _c.data["g16_corr_to_enthalpy"] = float(_ls[i].split()[4])
                if "Thermal correction to Gibbs Free Energy" in _l:
                    _c.data["g16_corr_to_gibbs"] = float(_ls[i].split()[6])
        logger.debug(f"done: {len(self.pack())}/{len(self.mols)} confomers")
        return self

    @Selectors.check_end.add("app/g16/output")
    def check_end(self, input_suffix: str = ".gjf"):
        for _c in self.pack():
            if _c.path.with_suffix(input_suffix).exists():
                with _c.path.with_suffix(input_suffix).open("r") as f:
                    _gjf_ls = f.readlines()
                _tot_jobs = 1
                for i, _gjf_l in enumerate(_gjf_ls):
                    if "--link1--" in _gjf_l.lower():
                        _tot_jobs += 1
                    if "freq" in _gjf_l.lower() and "opt" in _gjf_l.lower():
                        _tot_jobs += 1
                with _c.path.open() as f:
                    _ls = f.readlines()
                for _l in _ls:
                    if "Normal termination" in _l:
                        _tot_jobs -= 1
                if _tot_jobs == 0:
                    logger.debug(f"{_c.path.name} was terminated normally")
                else:
                    logger.info(f"{_c.path.name} was NOT terminated normally")
                    _c.deactivate("check_end")
            else:
                logger.debug(f"{_c.path.name}: could not find input: check if Normal termination is in last 10 lines")
                with _c.path.open() as f:
                    _ls = f.readlines()
                for _l in _ls[-10:]:
                    if "Normal termination" in _l:
                        logger.debug(f"{_c.path.name} was terminated normally")
                        break
                else:
                    logger.info(f"{_c.path.name} was NOT terminated normally")
                    _c.deactivate("check_end")
        logger.debug(f"done: {len(self.pack())}/{len(self.mols)} confomers")
        return self

    @Selectors.check_freq.add("app/g16/output")
    def check_freq(self, im: int = 0):
        for _c in self.pack():
            _freq = []
            with _c.path.open() as f:
                _ls = f.readlines()
            for i, _l in enumerate(_ls):
                if "1                      2                      3" in _l:
                    _freq = [float(_ls[i + 2].split()[k]) for k in range(2, 5)]
                    _c.data["g16_img_freq"] = len([_fr for _fr in _freq if _fr < 0])
                    if im != _c.data["g16_img_freq"]:
                        logger.info(
                            "number of imaginary frequency of {} was {}".format(_c.name, _c.data["g16_img_freq"])
                        )
                        _c.deactivate("check_freq")
            if _freq == []:
                logger.error(f"could not find frequency in {_c.path.name}")
        logger.debug(f"done: {len(self.pack())}/{len(self.mols)} confomers")
        return self

    def read_vibration(self, freq_number: int = 1):
        read_vibration(self, freq_number)
        logger.debug(f"done: {len(self.pack())}/{len(self.mols)} confomers")
        return self

    def check_bonding(
        self, bonding_atoms: List[int], additional_valid_range: int = 0, acceptable_invalid_atoms: int = 0
    ):
        if len(self.pack()) != len(self.pack().has_data("g16_vibration")):
            logger.info("read_vibration called automatically")
            self.read_vibration()
        check_bonding(self, bonding_atoms, additional_valid_range, acceptable_invalid_atoms)
        logger.debug(f"done: {len(self.pack())}/{len(self.mols)} confomers")
        return self

    def read_nmr(self):
        read_nmr(self)
        logger.debug(f"done: {len(self.pack())}/{len(self.mols)} confomers")
        return self

    def read_ecd(self):
        read_ecd(self)
        logger.debug(f"done: {len(self.pack())}/{len(self.mols)} confomers")
        return self

    def read_coupling(self, key_value: str = "Total_J"):
        read_spinspin(self, key_value)
        logger.debug(f"done: {len(self.pack())}/{len(self.mols)} confomers")
        return self

    def check_resub(self, input_suffix=".gjf", output_suffix=".log"):
        check_resub(self, input_suffix, output_suffix)
        logger.debug(f"done: {len(self.pack())}/{len(self.mols)} confomers")
        return self

    @Selectors.calc_free_energy.add("app/g16/output")
    def calc_free_energy(self, keys: List[str] = ["g16_scf", "g16_corr_to_gibbs"], unit: Unit = Units.hartree):
        return self.calc_energy(keys=keys, unit=unit)

    def is_input(self):
        for _c in self.pack():
            if not is_g16_input(_c.path):
                _c.flag = False
        logger.debug(f"done: {len(self.pack())}/{len(self.mols)} confomers")
        return self

    def is_output(self):
        for _c in self.pack():
            if not is_g16_output(_c.path):
                _c.flag = False
        logger.debug(f"done: {len(self.pack())}/{len(self.mols)} confomers")
        return self

    @Selectors.run.add("app/g16/input")
    def run(self):
        for _c in self.pack():
            run(_c)
        logger.debug(f"done: {len(self.pack())}/{len(self.mols)} confomers")
        return self

    @Selectors.submit.add("app/g16/input")
    def submit(self):
        for _c in self.pack():
            submit(_c)
        logger.debug(f"done: {len(self.pack())}/{len(self.mols)} confomers")
        return self
