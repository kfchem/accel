import subprocess
import time
from pathlib import Path

from accel.base.boxcore import BoxCore
from accel.base.selector import Selectors
from accel.base.systems import System
from accel.util import Execmd, FileType
from accel.util.log import logger


def que_submit(c: System):
    try:
        joc_id = subprocess.run([Execmd.get("qsub"), str(c.path)], cwd=str(c.path.parent), stdout=subprocess.PIPE)
        joc_id = joc_id.stdout.decode("utf-8").replace("\n", "")
        logger.info(f"qsub: {c.path.name} was submitted")
    except subprocess.CalledProcessError:
        c.state = False
        joc_id = "Submission Error"
        logger.error(f"qsub: failed submission of {c.path.name}")
    c.data["jobid"] = joc_id


def que_wait(box: BoxCore, interval_time=10):
    logger.info("qsub: waiting completion of tasks")
    job_ids = []
    for c in box.get():
        job_ids.append(c.data["jobid"])
    while len(job_ids) != 0:
        time.sleep(interval_time)
        proc = subprocess.run([Execmd.get("qstat")], stdout=subprocess.PIPE)
        proc_list = proc.stdout.decode("utf-8").split("\n")
        proc_list = [_l.split()[0] for _l in proc_list[2:] if len(_l) != 0]
        for job_id in job_ids:
            if not (job_id in proc_list):
                job_ids.remove(job_id)
                logger.info(f"Job ID {job_id} was completed")


@FileType.add("app/pbs/jobscript", 50)
def is_pbs_jobscript(p: Path) -> bool:
    if p.suffix not in (".sh", ".qsh", ".qsub"):
        return False
    with p.open() as f:
        for i, line in enumerate(f):
            if "#PBS" in line:
                return True
            if i > 10:
                break
    return False


class PbsBox(BoxCore):
    @Selectors.submit.add("app/pbs/jobscript")
    def submit(self):
        for c in self.get():
            que_submit(c)
        logger.debug(f"done: {str(self)}")
        return self

    def wait(self, interval_time=10):
        que_wait(self, interval_time)
        logger.debug(f"done: {str(self)}")
        return self

    @Selectors.run.add("app/pbs/jobscript")
    def run(self):
        for c in self.get():
            que_submit(c)
        que_wait(self)
        logger.debug(f"done: {str(self)}")
        return self
