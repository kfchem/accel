from pathlib import Path

from accel.util.log import logger


def make_dir(dir_name: Path):
    _p = Path(dir_name)
    if not (_p.exists()):
        _p.mkdir(parents=True)
        logger.info(f"{_p.name} was created")
    return _p


def change_dir(current_path: Path, new_dir: Path, new_name: str = None) -> Path:
    if new_dir is None or new_dir == "":
        # logger.info("directory not specified: used the same directory")
        _p = current_path
    elif Path(new_dir).is_dir():
        _p = Path(new_dir).joinpath(current_path.name)
    elif Path(new_dir).is_file():
        logger.debug("new_dir: file path was refered as directory")
        _p = Path(new_dir).parent.joinpath(current_path.name)
    else:
        _p = make_dir(new_dir).joinpath(current_path.name)
    if not (new_name is None or new_name == ""):
        _p = _p.with_name(new_name)
    return _p


def float_to_str(float_num):
    _fstr = str(float_num)
    if "e" in _fstr:
        _dig, _exp = _fstr.split("e")
        _dig = _dig.replace(".", "").replace("-", "")
        _exp = int(_exp)
        _zeros = "0" * (abs(int(_exp)) - 1)
        _sign = "-" if float_num < 0 else ""
        if _exp > 0:
            _fstr = f"{_sign}{_dig}{_zeros}.0"
        else:
            _fstr = f"{_sign}0.{_zeros}{_dig}"
    return _fstr
