from pathlib import Path

from accel.util.log import logger


def make_dir(dir_name: Path):
    p = Path(dir_name)
    if not (p.exists()):
        p.mkdir(parents=True)
        logger.info(f"{p.name} was created")
    return p


def change_dir(current_path: Path, new_dir: Path, new_name: str = None) -> Path:
    if new_dir is None or new_dir == "":
        # logger.info("directory not specified: used the same directory")
        p = current_path
    elif Path(new_dir).is_dir():
        p = Path(new_dir).joinpath(current_path.name)
    elif Path(new_dir).is_file():
        logger.debug("new_dir: file path was refered as directory")
        p = Path(new_dir).parent.joinpath(current_path.name)
    else:
        p = make_dir(new_dir).joinpath(current_path.name)
    if not (new_name is None or new_name == ""):
        p = p.with_name(new_name)
    return p


def float_to_str(float_num):
    fstr = str(float_num)
    if "e" in fstr:
        dig_, exp_ = fstr.split("e")
        dig_ = dig_.replace(".", "").replace("-", "")
        exp_ = int(exp_)
        zeros_ = "0" * (abs(int(exp_)) - 1)
        sign_ = "-" if float_num < 0 else ""
        if exp_ > 0:
            fstr = f"{sign_}{dig_}{zeros_}.0"
        else:
            fstr = f"{sign_}0.{zeros_}{dig_}"
    return fstr
