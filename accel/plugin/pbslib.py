import subprocess
import time
from pathlib import Path

from accel.base.boxcore import BoxCore
from accel.base.mols import Mol
from accel.base.selector import Selectors
from accel.util import Execmd, FileType
from accel.util.log import logger


def que_submit(_c: Mol):
    try:
        _jid = subprocess.run([Execmd.get("qsub"), str(_c.path)], cwd=str(_c.path.parent), stdout=subprocess.PIPE)
        _jid = _jid.stdout.decode("utf-8").replace("\n", "")
        logger.info(f"qsub: {_c.path.name} was submitted")
    except subprocess.CalledProcessError:
        _c.flag = False
        _jid = "Submission Error"
        logger.error(f"qsub: failed submission of {_c.path.name}")
    _c.data["jobid"] = _jid


def que_wait(mulcos: BoxCore, interval_time=10):
    logger.info("qsub: waiting completion of tasks")
    _jid_list = []
    for _c in mulcos.mols:
        _jid_list.append(_c.data["jobid"])
    while len(_jid_list) != 0:
        time.sleep(interval_time)
        _proc = subprocess.run([Execmd.get("qstat")], stdout=subprocess.PIPE)
        _plist = _proc.stdout.decode("utf-8").split("\n")
        _plist = [_l.split()[0] for _l in _plist[2:] if len(_l) != 0]
        for _jid in _jid_list:
            if not (_jid in _plist):
                _jid_list.remove(_jid)
                logger.info(f"Job ID {_jid} was completed")


@FileType.add("app/pbs/jobscript", 50)
def is_pbs_jobscript(_p: Path) -> bool:
    if _p.suffix not in (".sh", ".qsh", ".qsub"):
        return False
    with _p.open() as _f:
        for _i, _l in enumerate(_f):
            if "#PBS" in _l:
                return True
            if _i > 10:
                break
    return False


class PbsBox(BoxCore):
    @Selectors.submit.add("app/pbs/jobscript")
    def submit(self):
        for _c in self.mols:
            que_submit(_c)
        logger.debug(f"done: {str(self)}")
        return self

    def wait(self, interval_time=10):
        que_wait(self, interval_time)
        logger.debug(f"done: {str(self)}")
        return self

    @Selectors.run.add("app/pbs/jobscript")
    def run(self):
        for _c in self.mols:
            que_submit(_c)
        que_wait(self)
        logger.debug(f"done: {str(self)}")
        return self
