import statistics
import subprocess
from collections import defaultdict
from pathlib import Path

import numpy as np

from accel.base.boxcore import BoxCore
from accel.base.selector import Selectors
from accel.base.systems import System, Systems
from accel.util import Execmd, FileType
from accel.util.constants import Elements, Unit, Units
from accel.util.log import logger


# not coding yet
def read_multiple_xyz(box: BoxCore):
    pass


def calc_sum_energy(box: BoxCore):
    for c in box.get():
        c.data["g16_t_zero"] = c.data["g16_zpc"] + c.data["g16_scf"]
        c.data["g16_t_energy"] = c.data["g16_corr_to_energy"] + c.data["g16_scf"]
        c.data["g16_t_enthalpy"] = c.data["g16_corr_to_enthalpy"] + c.data["g16_scf"]
        c.data["g16_t_gibbs"] = c.data["g16_corr_to_gibbs"] + c.data["g16_scf"]
        c.energy = Units.hartree(c.data["g16_t_gibbs"]).to_kcal_mol


def read_low_freqs(box: BoxCore):
    for c in box.get():
        with c.path.open() as f:
            ls = f.readlines()
        for i, l_ in enumerate(ls):
            if "1                      2                      3" in l_:
                freq = [float(ls[i + 2].split()[k]) for k in range(2, 5)]
                c.data["g16_freq_1"] = freq[0]
                c.data["g16_freq_2"] = freq[1]
                c.data["g16_freq_3"] = freq[2]


def read_vibration(box: BoxCore, freq_number: int = 1):
    for c in box.get():
        freq = []
        with c.path.open() as f:
            ls = f.readlines()
        for i, line in enumerate(ls):
            if "Frequencies --" in line:
                no_row = [int(s) for s in ls[i - 2].split()]
                if freq_number in no_row:
                    column_no = no_row.index(freq_number)
                    for k, l_following in enumerate(ls[i:]):
                        if "Atom" in l_following:
                            row_no = i + k + 1
                            break
        for line in ls[row_no:]:
            if len(line.split()) in [3, 0, 1]:
                break
            freq.append(line.split()[2 + column_no * 3 : 2 + (column_no + 1) * 3])
        c.data["g16_vibration"] = freq


def check_bonding(
    box: BoxCore,
    bonding_atoms: list[int],
    additional_valid_range: int = 0,
    acceptable_invalid_atoms: int = 0,
):
    for c in box.get():
        distances = []
        for i, vib in enumerate(c.data["g16_vibration"]):
            distances.append([i + 1, np.linalg.norm(np.array(vib))])
        rank = [distance[0] for distance in sorted(distances, reverse=True, key=lambda x: x[1])]
        logger.debug(f"vib_rank: {rank}")
        invalid_atoms = []
        for a in bonding_atoms:
            valid_range = len(bonding_atoms) + additional_valid_range
            atom_rank = rank.index(a) + 1
            logger.debug(f"the rank of atom {a} was {atom_rank}: {c.path.name}")
            if atom_rank <= valid_range:
                logger.debug(f"the rank of atom {a} was in valid range")
            else:
                logger.info(f"the rank of atom {a} was in invalid range")
                invalid_atoms.append(a)
        if len(invalid_atoms) > acceptable_invalid_atoms:
            logger.info(f"atoms {invalid_atoms} are out of range")
            c.state = False


def read_nmr(box: BoxCore):
    for c in box.get():
        with c.path.open() as f:
            ls = f.readlines()
        for line in ls:
            if "Anisotropy =" in line:
                splitted_line = line.split()
                a = c.atoms.get(int(splitted_line[0]))
                if a.symbol != Elements.canonicalize(splitted_line[1]):
                    logger.error("atomic symbol does not matched")
                    break
                else:
                    a.data["isotropic"] = float(line.split()[4])
                    a.data["anisotropy"] = float(line.split()[7])


def read_ecd(box: BoxCore):
    for c in box.get():
        with c.path.open() as f:
            ls = f.readlines()
        state_dicts = defaultdict(dict)
        for i, line in enumerate(ls):
            if "Excited State" in line:
                splitted_line = line.split()
                d = state_dicts[int(splitted_line[2].replace(":", ""))]
                d["shape"] = splitted_line[3]
                d["energy"] = float(splitted_line[4])
                d["length"] = float(splitted_line[6])
                d["f"] = float(splitted_line[8][2:])
                d["S2"] = float(splitted_line[9][7:])
            if "R(velocity)" in line:
                for l_following in ls[i + 1 :]:
                    if "rxdel" in l_following:
                        break
                    state_dicts[int(l_following.split()[0])]["R_velocity"] = float(l_following.split()[4])
            if "R(length)" in line:
                for l_following in ls[i + 1 :]:
                    if "del" in l_following:
                        break
                    state_dicts[int(l_following.split()[0])]["R_length"] = float(l_following.split()[4])
        c.data["ecd"] = dict(state_dicts)


def read_spinspin(box: BoxCore, key_value: str = "Total_J"):
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
    key_dict = {
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
    key_string = key_dict[key_value]
    logger.info(f"{key_string[:-1]} is used for extracting coupling data")
    for c in box.get():
        with c.path.open() as f:
            ls = f.readlines()
        for i, line in enumerate(ls):
            if key_string in line:
                for a in c.atoms:
                    a.data["coupling"] = {}
                col_idx = []
                for l_following in ls[i + 1 :]:
                    if not l_following.split()[0].isdecimal():
                        break
                    elif l_following.startswith("            "):
                        col_idx = [int(k) for k in l_following.split()]
                    else:
                        for k, j_value in enumerate(l_following.split()[1:]):
                            j_value = float(j_value.replace("D", "E"))
                            c.atoms.get(int(l_following.split()[0]) + 1).data["coupling"][col_idx[k]] = j_value
                            c.atoms.get(col_idx[k] + 1).data["coupling"][int(l_following.split()[0])] = j_value


def check_resub(box: BoxCore, input_suffix=".gjf", log_suffix=".log"):
    def _gen_gjf_from_final_stdorient(p: Path):
        with p.open() as f:
            ls = f.readlines()
        geometory_position = 0
        for i, line in enumerate(ls):
            if "Standard orientation:" in line:
                geometory_position = i
        if geometory_position == 0:
            logger.error("could not find Standard orientation entry: " + p.name)
            return False
        geometories = []
        for line in ls[(geometory_position + 5) :]:
            if "---------------------------------------------------------------------" in line:
                break
            geometories.append(line.split())
        geometories = [
            "{1:<2} {0[3]:>11} {0[4]:>11} {0[5]:>11}\n".format(g, Elements.canonicalize(int(g[1])))
            for g in geometories
        ]

        with p.with_suffix(input_suffix).open() as f:
            gjf_ls = f.readlines()
        blank_line_indexes = [i for i, gjf_line in enumerate(gjf_ls) if gjf_line == "\n"]
        out_lines = gjf_ls[0 : blank_line_indexes[1] + 2]
        out_lines.extend(geometories)
        out_lines.extend(gjf_ls[blank_line_indexes[2] :])

        suffix = "_"
        for ext_ in [input_suffix, log_suffix, ".chk"]:
            while p.with_suffix(ext_ + suffix).exists():
                suffix = suffix + "_"

        for ext_ in [input_suffix, log_suffix, ".chk"]:
            if p.with_suffix(ext_).exists():
                p.with_suffix(ext_).rename(p.with_suffix(ext_ + suffix))
                logger.debug(f"{p.with_suffix(ext_).name} was renamed as {p.with_suffix(ext_ + suffix).name}")

        with p.with_suffix(input_suffix).open("w") as f:
            f.writelines(out_lines)
        logger.debug(f"{p.with_suffix(input_suffix).name} was generated")
        return True

    class Gerr:
        def __init__(self, name="", statement="", resub=True, use_final_geom=True) -> None:
            self.name = name
            self.statement = statement
            self.resub = resub
            self.use_final_geom = use_final_geom

        def resolve(self, c: System):
            logger.info(f"{c.name}: {self.name} was detected")
            if self.resub:
                if self.use_final_geom:
                    if _gen_gjf_from_final_stdorient(c.path.with_suffix(log_suffix)):
                        c.data["resubmission"] = True
                    else:
                        return False
                else:
                    c.data["resubmission"] = True
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

    for c in box.get():
        c.data["resubmission"] = False
        with c.path.with_suffix(input_suffix).open("r") as f:
            gjf_ls = f.readlines()
        tot_jobs = 1
        for i, gjf_l in enumerate(gjf_ls):
            if "--link1--" in gjf_l.lower():
                tot_jobs += 1
            if "freq" in gjf_l.lower() and "opt" in gjf_l.lower():
                tot_jobs += 1

        log_p = c.path.with_suffix(log_suffix)
        if not log_p.exists():
            no_log_err.resolve(c)
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
                    link9999_err.resolve(c)
                elif (statistics.stdev(rms_force_vs[-12::2]) < 0.000001) or (
                    statistics.stdev(rms_force_vs[-12::3]) < 0.000001
                ):
                    link9999_w_vib.resolve(c)
                else:
                    link9999_wo_vib.resolve(c)
                err_flag = True
                break
            for _e in normal_errs:
                if _e.statement in log_l:
                    _e.resolve(c)
                    err_flag = True
                    break
            else:
                continue
            break

        if not err_flag:
            if norm_term == tot_jobs:
                logger.debug(f"{c.name} terminated normally")
            else:
                unknown_err.resolve(c)
                logger.error(f"{c.name}: unknown termination was detected")

    if len(box.get().has_data("resubmission", True)) == 0:
        logger.info("There is no resubmission")


@FileType.add("app/g16/input", 30)
def is_g16_input(p: Path) -> bool:
    if p.suffix not in (".gjf", ".com", ".inp"):
        return False

    route_pass = False
    with p.open() as _f:
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
def is_g16_output(p: Path) -> bool:
    if p.suffix not in (".log", ".out", ".g16", ".g09", ".g98"):
        return False
    with p.open() as f:
        for i, line in enumerate(f):
            if "Entering Gaussian System" in line:
                return True
            if i > 100:
                break
    return False


def run(c: System):
    cmd = [Execmd.get("g16"), str(c.path.resolve().absolute())]
    try:
        logger.info(f"running: {c.name}: {''.join(cmd)}")
        proc = subprocess.run(cmd, cwd=str(c.path.parent), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        _out = proc.stdout.decode("utf-8").split("\n")
        logger.info(f"finished: {c.name}: {_out}")
    except subprocess.CalledProcessError:
        c.state = False
        logger.error(f"failed: {c.name}: {''.join(cmd)}")


def submit(c: System):
    cmd = [Execmd.get("g16"), str(c.path.resolve().absolute())]
    try:
        subprocess.Popen(cmd, cwd=str(c.path.parent), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        logger.info(f"submited: {c.name}: {''.join(cmd)}")
    except subprocess.CalledProcessError:
        c.state = False
        logger.error(f"failed: {c.name}: {''.join(cmd)}")


class GauBox(BoxCore):
    @Selectors.read_atoms.add("app/g16/output")
    def read_atoms(self, source: str = "archive"):
        """
        source = 'archive' | 'standard_orientation' | 'input_orientation'
        """
        for c in self.get():
            with c.path.open() as f:
                ls = f.readlines()
            charge = None
            multiplicity = None
            if source == "archive":
                arc_sec = ""
                for i, line in enumerate(ls):
                    if "1\\1\\" in line:
                        arc_sec = ""
                        for line in ls[i:]:
                            arc_sec += "".join(line[1:].splitlines())
                            if "\\\\@" in line:
                                break
                if arc_sec == "":
                    logger.error(f"could not find archive section in {c.path.name}")
                    c.deactivate("read_atoms: g16")
                    continue
                arc_sec = arc_sec.split("\\\\")
                splitted_arc_sec = arc_sec[3].split("\\")
                charge = int(splitted_arc_sec[0].split(",")[0])
                multiplicity = int(splitted_arc_sec[0].split(",")[1])
                axyzs = [i.split(",") for i in splitted_arc_sec[1:]]
                if len(axyzs[0]) == 5:
                    for line in axyzs:
                        line.pop(1)
                elif len(axyzs[0]) == 1:
                    logger.error(f"could not resolve archive section in {c.path.name}")
                    c.deactivate("read_atoms: g16")
                    continue

            elif source in ("standard_orientation", "input_orientation"):
                key_words = {
                    "standard_orientation": "Standard orientation:",
                    "input_orientation": "Input orientation:",
                }
                line_idx = None
                chg_mult_idx = None
                for i, line in enumerate(ls):
                    if key_words[source] in line:
                        line_idx = i
                    if "Charge =" in line:
                        chg_mult_idx = i
                if line_idx is None:
                    logger.error(f"could not find orientation section in {c.path.name}")
                    c.deactivate("read_atoms: g16")
                    continue
                if chg_mult_idx is not None:
                    charge = int(ls[chg_mult_idx].split()[2])
                    multiplicity = int(ls[chg_mult_idx].split()[5])
                axyzs = []
                for line in ls[(line_idx + 5) :]:
                    if "---------------------------------------------------------------------" in line:
                        break
                    axyzs.append(line.split())
                axyzs = [[Elements.canonicalize(int(axyz[1])), axyz[3], axyz[4], axyz[5]] for axyz in axyzs]

            c.atoms.clear()
            for line in axyzs:
                c.atoms.append(line)
            if charge is not None:
                c.charge = charge
            if multiplicity is not None:
                c.multiplicity = multiplicity
            logger.debug(f"read {source} from {c.path.name}")
        logger.debug(f"done: {str(self)}")
        return self

    @Selectors.read_atoms.add("app/g16/input")
    def read_atoms_from_input(self):
        for c in self.get():
            with c.path.open() as f:
                ls = f.readlines()
            blank_line_count = 0
            atomic_coord = []
            for line in ls:
                if blank_line_count == 2:
                    atomic_coord.append(line.replace(",", " ").replace("/", " ").replace("\t", " ").split())
                elif blank_line_count == 3:
                    break
                if line == "\n":
                    blank_line_count += 1
            chg_mult = atomic_coord.pop(0)
            atomic_coord.pop()
            try:
                if len(chg_mult) != 2:
                    raise ValueError
                c.charge = int(chg_mult[0])
                c.multiplicity = int(chg_mult[1])
                uni_len = list({len(_l) for _l in atomic_coord})
                if len(uni_len) != 1 or uni_len[0] != 4:
                    logger.error("currently cannot read Z-matrix")
                    raise ValueError
                c.atoms.clear()
                c.atoms.extend(atomic_coord)
            except ValueError:
                c.deactivate("read_atoms: g16")
                logger.error(f"could not read coordinate of {c.path.name}")
                continue
            else:
                logger.debug(f"read atoms of {c.path.name}")
        logger.debug(f"done: {str(self)}")
        return self

    @Selectors.read_energy.add("app/g16/output")
    def read_scf(self, calculation_method: str = None):
        for c in self.get():
            if calculation_method is None:
                scf_key = "SCF Done:"
            else:
                scf_key = "SCF Done:  E(" + calculation_method + ")"
                logger.info(f"{scf_key} is used for search SCF")
            with c.path.open() as f:
                ls = f.readlines()
            position_num = 0
            for i, line in enumerate(ls):
                if scf_key in line:
                    position_num = i
            if position_num == 0:
                logger.error("could not find the scf entry")
                c.deactivate("read_scf")
                continue
            c.data["g16_scf"] = float(ls[position_num].split()[4])
            c.energy = Units.hartree(c.data["g16_scf"]).to_kcal_mol
        logger.debug(f"done: {str(self)}")
        return self

    @Selectors.read_thermal.add("app/g16/output")
    def read_thermal(self):
        for c in self.get():
            with c.path.open() as f:
                ls = f.readlines()
            for i, line in enumerate(ls):
                if "Zero-point correction" in line:
                    c.data["g16_zpc"] = float(ls[i].split()[2])
                if "Thermal correction to Energy" in line:
                    c.data["g16_corr_to_energy"] = float(ls[i].split()[4])
                if "Thermal correction to Enthalpy" in line:
                    c.data["g16_corr_to_enthalpy"] = float(ls[i].split()[4])
                if "Thermal correction to Gibbs Free Energy" in line:
                    c.data["g16_corr_to_gibbs"] = float(ls[i].split()[6])
        logger.debug(f"done: {str(self)}")
        return self

    @Selectors.check_end.add("app/g16/output")
    def check_end(self, input_suffix: str = ".gjf"):
        for c in self.get():
            if c.path.with_suffix(input_suffix).exists():
                with c.path.with_suffix(input_suffix).open("r") as f:
                    gjf_ls = f.readlines()
                total_jobs_num = 1
                for gjf_l in gjf_ls:
                    if "--link1--" in gjf_l.lower():
                        total_jobs_num += 1
                    if "freq" in gjf_l.lower() and "opt" in gjf_l.lower():
                        total_jobs_num += 1
                with c.path.open() as f:
                    ls = f.readlines()
                for line in ls:
                    if "Normal termination" in line:
                        total_jobs_num -= 1
                if total_jobs_num == 0:
                    logger.debug(f"{c.path.name} was terminated normally")
                else:
                    logger.info(f"{c.path.name} was NOT terminated normally")
                    c.deactivate("check_end")
            else:
                logger.debug(f"{c.path.name}: could not find input: check if Normal termination is in last 10 lines")
                with c.path.open() as f:
                    ls = f.readlines()
                for line in ls[-10:]:
                    if "Normal termination" in line:
                        logger.debug(f"{c.path.name} was terminated normally")
                        break
                else:
                    logger.info(f"{c.path.name} was NOT terminated normally")
                    c.deactivate("check_end")
        logger.debug(f"done: {str(self)}")
        return self

    @Selectors.check_freq.add("app/g16/output")
    def check_freq(self, im: int = 0):
        for c in self.get():
            freqs = []
            with c.path.open() as f:
                ls = f.readlines()
            for i, line in enumerate(ls):
                if "1                      2                      3" in line:
                    freqs = [float(ls[i + 2].split()[k]) for k in range(2, 5)]
                    c.data["g16_img_freq"] = len([freq for freq in freqs if freq < 0])
                    if im != c.data["g16_img_freq"]:
                        logger.info(
                            "number of imaginary frequency of {} was {}".format(c.name, c.data["g16_img_freq"])
                        )
                        c.deactivate("check_freq")
            if freqs == []:
                logger.error(f"could not find frequency in {c.path.name}")
        logger.debug(f"done: {str(self)}")
        return self

    def read_vibration(self, freq_number: int = 1):
        read_vibration(self, freq_number)
        logger.debug(f"done: {str(self)}")
        return self

    def check_bonding(
        self, bonding_atoms: list[int], additional_valid_range: int = 0, acceptable_invalid_atoms: int = 0
    ):
        if len(self.get()) != len(self.get().has_data("g16_vibration")):
            logger.info("read_vibration called automatically")
            self.read_vibration()
        check_bonding(self, bonding_atoms, additional_valid_range, acceptable_invalid_atoms)
        logger.debug(f"done: {str(self)}")
        return self

    def read_nmr(self):
        read_nmr(self)
        logger.debug(f"done: {str(self)}")
        return self

    def read_ecd(self):
        read_ecd(self)
        logger.debug(f"done: {str(self)}")
        return self

    def read_coupling(self, key_value: str = "Total_J"):
        read_spinspin(self, key_value)
        logger.debug(f"done: {str(self)}")
        return self

    def check_resub(self, input_suffix=".gjf", output_suffix=".log"):
        check_resub(self, input_suffix, output_suffix)
        logger.debug(f"done: {str(self)}")
        return self

    def get_resub(self, input_suffix=".gjf", output_suffix=".log") -> Systems:
        check_resub(self, input_suffix, output_suffix)
        logger.debug(f"done: {str(self)}")
        return self.get().has_data("resubmission", True)

    @Selectors.calc_free_energy.add("app/g16/output")
    def calc_free_energy(self, keys: list[str] = ["g16_scf", "g16_corr_to_gibbs"], unit: Unit = Units.hartree):
        return self.calc_energy(keys=keys, unit=unit)

    def is_input(self):
        for c in self.get():
            if not is_g16_input(c.path):
                c.state = False
        logger.debug(f"done: {str(self)}")
        return self

    def is_output(self):
        for c in self.get():
            if not is_g16_output(c.path):
                c.state = False
        logger.debug(f"done: {str(self)}")
        return self

    @Selectors.run.add("app/g16/input")
    def run(self):
        for c in self.get():
            run(c)
        logger.debug(f"done: {str(self)}")
        return self

    @Selectors.submit.add("app/g16/input")
    def submit(self):
        for c in self.get():
            submit(c)
        logger.debug(f"done: {str(self)}")
        return self

    def get_irc(self) -> Systems:
        ret_cs = []
        for c in self.get():
            logger.debug(f"reading irc of {c.path.name}")
            with c.path.open() as f:
                ls = f.readlines()
            rxn_direction = "NAN"
            ts_appended = False
            t_charge = c.total_charge
            for i, line in enumerate(ls):
                if "Rxn path following direction =" in line:
                    rxn_direction = line.split()[5]
                if "Charge =" in line:
                    t_charge = int(line.split()[2])
                if "Cartesian Coordinates (Ang):" in line or "Input orientation:" in line:
                    geom = []
                    for _l in ls[(i + 5) :]:
                        if _l.find("---------------------------------------------------------------------") > -1:
                            break
                        axyz = []
                        axyz.append(int(_l.split()[1]))
                        axyz.append(float(_l.split()[-3]))
                        axyz.append(float(_l.split()[-2]))
                        axyz.append(float(_l.split()[-1]))
                        geom.append(axyz)
                if "NET REACTION COORDINATE UP TO THIS POINT =" in line:
                    rxn_coord = float(line.split()[8])
                    if rxn_direction == "Reverse":
                        rxn_coord = (-1) * rxn_coord
                if "Corrected End Point Energy =" in line:
                    scf_energy = float(line.split()[5])
                if "# OF POINTS ALONG THE PATH" in line:
                    num_of_points = int(line.split()[7])
                    new_c = c.duplicate()
                    new_c.label = new_c.name
                    new_c.name = f"{new_c.name}_{rxn_direction[0]}{num_of_points:03}"
                    new_c.energy = Units.hartree(scf_energy).to_kcal_mol
                    new_c.data["g16_irc_scf"] = scf_energy
                    new_c.data["rxn_coordinate"] = rxn_coord
                    new_c.total_charge = t_charge
                    new_c.atoms.clear()
                    new_c.atoms.extend(geom)
                    ret_cs.append(new_c)
                if "Current Structure is TS" in line and not ts_appended:
                    new_c = c.duplicate()
                    new_c.label = new_c.name
                    new_c.name = f"{new_c.name}_TS"
                    bf_ls = ls[:i]
                    while bf_ls:
                        _l = bf_ls.pop()
                        if "SCF Done:" in _l:
                            new_c.data["g16_irc_scf"] = float(_l.split()[4])
                            new_c.energy = Units.hartree(new_c.data["g16_irc_scf"]).to_kcal_mol
                            break
                    new_c.data["rxn_coordinate"] = 0
                    new_c.total_charge = t_charge
                    new_c.atoms.clear()
                    new_c.atoms.extend(geom)
                    ret_cs.append(new_c)
                    ts_appended = True
        logger.debug(f"done: {str(self)}")
        return Systems(ret_cs)

    def get_trj(self) -> Systems:
        ret_cs = []
        for c in self.get():
            logger.debug(f"reading trajectory of {c.path.name}")
            with c.path.open() as f:
                ls = f.readlines()
            t_charge = c.total_charge
            for i, line in enumerate(ls):
                if "Charge =" in line:
                    t_charge = int(line.split()[2])
                if "Input orientation:" in line:
                    geom = []
                    for _l in ls[(i + 5) :]:
                        if _l.find("---------------------------------------------------------------------") > -1:
                            break
                        axyz = []
                        axyz.append(int(_l.split()[1]))
                        axyz.append(float(_l.split()[-3]))
                        axyz.append(float(_l.split()[-2]))
                        axyz.append(float(_l.split()[-1]))
                        geom.append(axyz)
                if "SCF Done:" in line:
                    scf_energy = float(line.split()[4])
                    new_c = c.duplicate()
                    new_c.label = new_c.name
                    new_c.name = f"{new_c.name}_{len(ret_cs)+1:03}"
                    new_c.energy = Units.hartree(scf_energy).to_kcal_mol
                    new_c.data["g16_scf"] = scf_energy
                    new_c.data["trj_num"] = len(ret_cs) + 1
                    new_c.total_charge = t_charge
                    new_c.atoms.clear()
                    new_c.atoms.extend(geom)
                    ret_cs.append(new_c)
        logger.debug(f"done: {str(self)}")
        return Systems(ret_cs)
