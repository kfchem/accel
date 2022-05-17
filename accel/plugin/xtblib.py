from pathlib import Path

from accel.base.boxcore import BoxCore
from accel.base.selector import Selectors
from accel.util import FileType, Units
from accel.util.log import logger


def read_total_free_energy(mulcos: BoxCore):
    for _c in mulcos.mols:
        with _c.path.open(encoding="utf-8") as f:
            _ls = f.readlines()
        _n = 0
        for i, _l in enumerate(_ls):
            if "TOTAL FREE ENERGY" in _l:
                _n = i
        if _n != 0:
            logger.debug("xTB: {}: {}".format(_c.path.name, _ls[_n].replace("\n", "")))
            _c.energy = Units.hartree(float(_ls[_n].split()[4])).to_kcal_mol
        else:
            logger.error(f"xTB: {_c.path.name}: the energy entry was not found")
            _c.flag = False


@FileType.add("app/xtb/output", 50)
def is_xtb_output(_p: Path) -> bool:
    if _p.suffix not in (".log", ".out", ".xtb"):
        return False
    with _p.open() as _f:
        for _i, _l in enumerate(_f):
            if "|                           x T B                           |" in _l:
                return True
            if _i > 20:
                break
    return False


class XtbBox(BoxCore):
    @Selectors.check_end.add("app/xtb/output")
    def check_end(self):
        for _c in self.mols:
            with _c.path.open(encoding="utf-8") as f:
                _ls = f.readlines()
            _flag = False
            for _l in _ls:
                if "normal termination of xtb" in _l:
                    _flag = True
            if _flag:
                logger.debug(f"xTB: {_c.path.name} was terminated normally")
            else:
                logger.info(f"xTB: {_c.path.name} was NOT terminated normally")
                _c.deactivate("check_end")
        logger.debug(f"done: {str(self)}")
        return self

    @Selectors.read_energy.add("app/xtb/output")
    def read_energy(self):
        for _c in self.mols:
            with _c.path.open(encoding="utf-8") as f:
                _ls = f.readlines()
            _n = 0
            for i, _l in enumerate(_ls):
                if "TOTAL ENERGY" in _l:
                    _n = i
            if _n != 0:
                logger.debug("xTB: {}: {}".format(_c.path.name, _ls[_n].replace("\n", "")))
                _c.energy = Units.hartree(float(_ls[_n].split()[3])).to_kcal_mol
            else:
                logger.error(f"xTB: {_c.path.name}: the energy entry was not found")
                _c.deactivate("read_energy: xtb")
        logger.debug(f"done: {str(self)}")
        return self

    def read_free_energy(self):
        read_total_free_energy(self)
        logger.debug(f"done: {str(self)}")
        return self
