from pathlib import Path

from accel.base.boxcore import BoxCore
from accel.base.selector import Selectors
from accel.util import FileType, Units
from accel.util.log import logger


def read_total_free_energy(box: BoxCore):
    for c in box.get():
        with c.path.open(encoding="utf-8") as f:
            ls = f.readlines()
        position_idx = 0
        for i, line in enumerate(ls):
            if "TOTAL FREE ENERGY" in line:
                position_idx = i
        if position_idx != 0:
            logger.debug("xTB: {}: {}".format(c.path.name, ls[position_idx].replace("\n", "")))
            c.energy = Units.hartree(float(ls[position_idx].split()[4])).to_kcal_mol
        else:
            logger.error(f"xTB: {c.path.name}: the energy entry was not found")
            c.state = False


@FileType.add("app/xtb/output", 50)
def is_xtb_output(p: Path) -> bool:
    if p.suffix not in (".log", ".out", ".xtb"):
        return False
    with p.open() as f:
        for i, line in enumerate(f):
            if "|                           x T B                           |" in line:
                return True
            if i > 20:
                break
    return False


class XtbBox(BoxCore):
    @Selectors.check_end.add("app/xtb/output")
    def check_end(self):
        for c in self.get():
            with c.path.open(encoding="utf-8") as f:
                ls = f.readlines()
            terminated_nomally = False
            for line in ls:
                if "normal termination of xtb" in line:
                    terminated_nomally = True
            if terminated_nomally:
                logger.debug(f"xTB: {c.path.name} was terminated normally")
            else:
                logger.info(f"xTB: {c.path.name} was NOT terminated normally")
                c.deactivate("check_end")
        logger.debug(f"done: {str(self)}")
        return self

    @Selectors.read_energy.add("app/xtb/output")
    def read_energy(self):
        for c in self.get():
            with c.path.open(encoding="utf-8") as f:
                ls = f.readlines()
            position_idx = 0
            for i, line in enumerate(ls):
                if "TOTAL ENERGY" in line:
                    position_idx = i
            if position_idx != 0:
                logger.debug("xTB: {}: {}".format(c.path.name, ls[position_idx].replace("\n", "")))
                c.energy = Units.hartree(float(ls[position_idx].split()[3])).to_kcal_mol
            else:
                logger.error(f"xTB: {c.path.name}: the energy entry was not found")
                c.deactivate("read_energy: xtb")
        logger.debug(f"done: {str(self)}")
        return self

    def read_free_energy(self):
        read_total_free_energy(self)
        logger.debug(f"done: {str(self)}")
        return self
