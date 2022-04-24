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

old_ver = old_ver.replace("'", "").replace('"', "")
old_ver_list = [int(s) for s in old_ver.split(".")]
old_ver_list[-1] += 1
new_ver = ".".join([str(s) for s in old_ver_list])

for i, _l in enumerate(_ls):
    if "__version__ = " in _l:
        _ls[i] = f'__version__ = "{new_ver}"\n'
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
