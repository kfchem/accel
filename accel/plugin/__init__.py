from accel.base.boxcore import BoxCore
from accel.plugin.gaulib import GauPlugin
from accel.plugin.orclib import OrcPlugin
from accel.plugin.pbslib import PbsPlugin
from accel.plugin.txtlib import TxtPlugin


class Plugins:
    def __init__(self, box: BoxCore) -> None:
        self.box = box

    @property
    def gau(self):
        return GauPlugin(self.box)

    @property
    def orc(self):
        return OrcPlugin(self.box)

    @property
    def pbs(self):
        return PbsPlugin(self.box)

    @property
    def txt(self):
        return TxtPlugin(self.box)
