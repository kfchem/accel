import datetime
import subprocess
import sys
from pathlib import Path

init_path = Path(__file__).parent / "accel" / "__init__.py"

with init_path.open() as f:
    _ls = f.readlines()

for _l in _ls:
    if "__version__ = " in _l:
        old_ver = _l.split()[2]
        break

new_ver = str(datetime.date.today().strftime("%Y%m%d")) + "01"


def rec_ver(_ver: str):
    if old_ver[:8] == _ver[:8] and int(old_ver[2:]) >= int(_ver[2:]):
        _ver = str(int(_ver) + 1)
        return rec_ver(_ver)
    else:
        return _ver


new_ver = rec_ver(new_ver)

for i, _l in enumerate(_ls):
    if "__version__ = " in _l:
        _ls[i] = f"__version__ = {new_ver}\n"
        break

with init_path.open("w", newline="\n") as f:
    f.writelines(_ls)

try:
    _gvout = subprocess.run(["python", "setup.py", "bdist_wheel"])
except subprocess.CalledProcessError:
    print("Error")
    sys.exit()

try:
    _gvout = subprocess.run(["python", "setup.py", "clean", "--all"])
except subprocess.CalledProcessError:
    print("Error")
    sys.exit()
