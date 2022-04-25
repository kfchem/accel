import subprocess
from pathlib import Path

from accel.base.boxcore import BoxCore
from accel.base.mols import Mol
from accel.base.selector import Selectors
from accel.util import Execmd, FileType, Units
from accel.util.log import logger


def check_optimized(mulcos: BoxCore):
    for _c in mulcos.pack():
        with _c.path.open() as f:
            _ls = f.readlines()
        _flag = False
        for _l in _ls:
            if "*** OPTIMIZATION RUN DONE ***" in _l:
                _flag = True
        if _flag:
            logger.debug(f"ORCA: {_c.path.name} was optimized successfully")
        else:
            logger.info(f"ORCA: {_c.path.name} was NOT optimized successfully")
            _c.flag = False


def read_energy(_c: Mol):
    with _c.path.open() as f:
        _ls = f.readlines()
    _n = 0
    for i, _l in enumerate(_ls):
        if "FINAL SINGLE POINT ENERGY" in _l:
            _n = i
    if _n != 0:
        logger.debug("Orca: {}: {}".format(_c.path.name, _ls[_n].replace("\n", "")))
        _c.energy = Units.hartree(float(_ls[_n].split()[4])).to_kcal_mol
    else:
        logger.error(f"Orca: {_c.path.name}: the energy entry was not found")
        _c.deactivate("read_energy: orca")


@FileType.add("app/orca/input", 40)
def is_orca_input(_p: Path) -> bool:
    if _p.suffix not in (".com", ".inp", ".inp"):
        return False

    with _p.open() as _f:
        for _i, _l in enumerate(_f):
            if _l.startswith("#"):
                continue
            elif _l.startswith("!"):
                return True
            else:
                return False
    return False


@FileType.add("app/orca/output", 60)
def is_orca_output(_p: Path) -> bool:
    if _p.suffix not in (".log", ".out"):
        return False
    with _p.open() as _f:
        for _i, _l in enumerate(_f):
            if "* O   R   C   A *" in _l:
                return True
            if _i > 100:
                break
    return False


def read_atoms_from_xyz(_c: Mol):
    _xyz_path = _c.path.with_suffix(".xyz")
    if not _xyz_path.exists():
        _c.deactivate("read_atoms: orca from xyz file: not exists")
        return None
    with _xyz_path.open() as f:
        _ls = f.readlines()
    _axyz = [_l.split() for _l in _ls[2:] if len(_l.split()) == 4]
    if len(_axyz) != int(_ls[0]):
        _c.deactivate("read_atoms: orca from xyz file")
        return None
    _c.atoms.clear()
    for _l in _axyz:
        _c.atoms.append(_l)


def run(_c: Mol):
    _cmd = [
        Execmd.get("orca"),
        str(_c.path.resolve().absolute()),
        str(_c.path.resolve().absolute().with_suffix(".out")),
    ]
    try:
        logger.info(f"running: {_c.name}: {''.join(_cmd)}")
        _proc = subprocess.run(_cmd, cwd=str(_c.path.parent), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        _out = _proc.stdout.decode("utf-8").split("\n")
        logger.info(f"finished: {_c.name}: {_out}")
    except subprocess.CalledProcessError:
        _c.flag = False
        logger.error(f"failed: {_c.name}: {''.join(_cmd)}")


def submit(_c: Mol):
    _cmd = [
        Execmd.get("orca"),
        str(_c.path.resolve().absolute()),
        str(_c.path.resolve().absolute().with_suffix(".out")),
    ]
    try:
        subprocess.Popen(_cmd, cwd=str(_c.path.parent), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        logger.info(f"submited: {_c.name}: {''.join(_cmd)}")
    except subprocess.CalledProcessError:
        _c.flag = False
        logger.error(f"failed: {_c.name}: {''.join(_cmd)}")


class OrcBox(BoxCore):
    @Selectors.check_end.add("app/orca/output")
    def check_end(self):
        for _c in self.pack():
            with _c.path.open() as f:
                _ls = f.readlines()
            _flag = False
            for _l in _ls:
                if "****ORCA TERMINATED NORMALLY****" in _l:
                    _flag = True
            if _flag:
                logger.debug(f"ORCA: {_c.path.name} was terminated normally")
            else:
                logger.info(f"ORCA: {_c.path.name} was NOT terminated normally")
                _c.deactivate("check_end")
        logger.debug(f"done: {len(self.pack())}/{len(self.mols)} confomers")
        return self

    def check_optimized(self):
        check_optimized(self)
        logger.debug(f"done: {len(self.pack())}/{len(self.mols)} confomers")
        return self

    @Selectors.read_energy.add("app/orca/output")
    def read_energy(self):
        for _c in self.pack():
            read_energy(_c)
        logger.debug(f"done: {len(self.pack())}/{len(self.mols)} confomers")
        return self

    @Selectors.read_atoms.add("app/orca/output")
    def read_atoms_from_xyz(self):
        for _c in self.pack():
            read_atoms_from_xyz(_c)
        logger.debug(f"done: {len(self.pack())}/{len(self.mols)} confomers")
        return self

    def is_input(self):
        for _c in self.pack():
            if not is_orca_input(_c.path):
                _c.flag = False
        logger.debug(f"done: {len(self.pack())}/{len(self.mols)} confomers")
        return self

    def is_output(self):
        for _c in self.pack():
            if not is_orca_output(_c.path):
                _c.flag = False
        logger.debug(f"done: {len(self.pack())}/{len(self.mols)} confomers")
        return self

    @Selectors.run.add("app/orca/input")
    def run(self):
        for _c in self.pack():
            run(_c)
        logger.debug(f"done: {len(self.pack())}/{len(self.mols)} confomers")
        return self

    @Selectors.submit.add("app/orca/input")
    def submit(self):
        for _c in self.pack():
            submit(_c)
        logger.debug(f"done: {len(self.pack())}/{len(self.mols)} confomers")
        return self
