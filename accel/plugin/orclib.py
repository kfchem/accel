import subprocess
from pathlib import Path

from accel.base.boxcore import BoxCore
from accel.base.selector import Selectors
from accel.base.systems import System
from accel.util import Execmd, FileType, Units
from accel.util.log import logger


def check_optimized(box: BoxCore):
    for c in box.get():
        with c.path.open() as f:
            _ls = f.readlines()
        optimized = False
        for _l in _ls:
            if "*** OPTIMIZATION RUN DONE ***" in _l:
                optimized = True
        if optimized:
            logger.debug(f"ORCA: {c.path.name} was optimized successfully")
        else:
            logger.info(f"ORCA: {c.path.name} was NOT optimized successfully")
            c.state = False


def read_energy(c: System):
    with c.path.open() as f:
        ls = f.readlines()
    position_idx = 0
    for i, line in enumerate(ls):
        if "FINAL SINGLE POINT ENERGY" in line:
            position_idx = i
    if position_idx != 0:
        logger.debug("Orca: {}: {}".format(c.path.name, ls[position_idx].replace("\n", "")))
        c.energy = Units.hartree(float(ls[position_idx].split()[4])).to_kcal_mol
    else:
        logger.error(f"Orca: {c.path.name}: the energy entry was not found")
        c.deactivate("read_energy: orca")


@FileType.add("app/orca/input", 40)
def is_orca_input(p: Path) -> bool:
    if p.suffix not in (".com", ".inp", ".inp"):
        return False
    with p.open() as f:
        for line in f:
            if line.startswith("#"):
                continue
            elif line.startswith("!"):
                return True
            else:
                return False
    return False


@FileType.add("app/orca/output", 60)
def is_orca_output(p: Path) -> bool:
    if p.suffix not in (".log", ".out"):
        return False
    with p.open() as f:
        for i, line in enumerate(f):
            if "* O   R   C   A *" in line:
                return True
            if i > 100:
                break
    return False


def read_atoms_from_xyz(c: System):
    xyz_path = c.path.with_suffix(".xyz")
    if not xyz_path.exists():
        c.deactivate("read_atoms: orca from xyz file: not exists")
        return None
    with xyz_path.open() as f:
        ls = f.readlines()
    axyz = [line.split() for line in ls[2:] if len(line.split()) == 4]
    if len(axyz) != int(ls[0]):
        c.deactivate("read_atoms: orca from xyz file")
        return None
    c.atoms.clear()
    for line in axyz:
        c.atoms.append(line)


def run(c: System):
    cmd_txts = [
        Execmd.get("orca"),
        str(c.path.resolve().absolute()),
        str(c.path.resolve().absolute().with_suffix(".out")),
    ]
    try:
        logger.info(f"running: {c.name}: {''.join(cmd_txts)}")
        proc = subprocess.run(cmd_txts, cwd=str(c.path.parent), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        _out = proc.stdout.decode("utf-8").split("\n")
        logger.info(f"finished: {c.name}: {_out}")
    except subprocess.CalledProcessError:
        c.state = False
        logger.error(f"failed: {c.name}: {''.join(cmd_txts)}")


def submit(c: System):
    cmd_txts = [
        Execmd.get("orca"),
        str(c.path.resolve().absolute()),
        str(c.path.resolve().absolute().with_suffix(".out")),
    ]
    try:
        subprocess.Popen(cmd_txts, cwd=str(c.path.parent), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        logger.info(f"submited: {c.name}: {''.join(cmd_txts)}")
    except subprocess.CalledProcessError:
        c.state = False
        logger.error(f"failed: {c.name}: {''.join(cmd_txts)}")


class OrcBox(BoxCore):
    @Selectors.check_end.add("app/orca/output")
    def check_end(self):
        for c in self.get():
            with c.path.open() as f:
                ls = f.readlines()
            terminated_normally = False
            for line in ls:
                if "****ORCA TERMINATED NORMALLY****" in line:
                    terminated_normally = True
            if terminated_normally:
                logger.debug(f"ORCA: {c.path.name} was terminated normally")
            else:
                logger.info(f"ORCA: {c.path.name} was NOT terminated normally")
                c.deactivate("check_end")
        logger.debug(f"done: {str(self)}")
        return self

    def check_optimized(self):
        check_optimized(self)
        logger.debug(f"done: {str(self)}")
        return self

    @Selectors.read_energy.add("app/orca/output")
    def read_energy(self):
        for c in self.get():
            read_energy(c)
        logger.debug(f"done: {str(self)}")
        return self

    @Selectors.read_atoms.add("app/orca/output")
    def read_atoms_from_xyz(self):
        for c in self.get():
            read_atoms_from_xyz(c)
        logger.debug(f"done: {str(self)}")
        return self

    def is_input(self):
        for c in self.get():
            if not is_orca_input(c.path):
                c.state = False
        logger.debug(f"done: {str(self)}")
        return self

    def is_output(self):
        for c in self.get():
            if not is_orca_output(c.path):
                c.state = False
        logger.debug(f"done: {str(self)}")
        return self

    @Selectors.run.add("app/orca/input")
    def run(self):
        for c in self.get():
            run(c)
        logger.debug(f"done: {str(self)}")
        return self

    @Selectors.submit.add("app/orca/input")
    def submit(self):
        for c in self.get():
            submit(c)
        logger.debug(f"done: {str(self)}")
        return self
