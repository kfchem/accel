from accel.base.boxcore import BoxCore
from accel.plugin.gaulib import GauBox
from accel.plugin.maelib import MaeBox
from accel.plugin.orclib import OrcBox
from accel.plugin.pbslib import PbsBox
from accel.plugin.txtlib import TxtBox
from accel.plugin.xtblib import XtbBox


class Plugins:
    def __init__(self, box: BoxCore) -> None:
        self.box = box

    @property
    def gau(self):
        return GauBox(self.box)

    @property
    def orc(self):
        return OrcBox(self.box)

    @property
    def pbs(self):
        return PbsBox(self.box)

    @property
    def txt(self):
        return TxtBox(self.box)

    @property
    def mae(self):
        return MaeBox(self.box)

    @property
    def xtb(self):
        return XtbBox(self.box)
