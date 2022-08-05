from pathlib import Path

from accel.base.systems import System
from accel.base.tools import change_dir, float_to_str
from accel.util.log import logger


def replace_key(c: System, lines: list[str]):
    rls = "".join(lines)
    while True:
        if "#NAME#" in rls:
            rls = rls.replace("#NAME#", c.name)
            continue
        if "#DATA[" in rls:
            key = rls.split("#DATA[")[1].split("]#")[0]
            rls = rls.replace("#DATA[" + key + "]#", c.data[key])
            continue
        if "#AXYZ#" in rls:
            _xyz = [
                "{:<2} {:>15} {:>15} {:>15}".format(
                    _a.symbol,
                    float_to_str(_a.x),
                    float_to_str(_a.y),
                    float_to_str(_a.z),
                )
                for _a in c.atoms
            ]
            rls = rls.replace("#AXYZ#", "\n".join(_xyz))
            continue
        if "#ATOMS#" in rls:
            rls = rls.replace("#ATOMS#", str(len(c.atoms)))
            continue
        if "#CHG#" in rls:
            rls = rls.replace("#CHG#", str(c.charge))
            continue
        if "#MULT#" in rls:
            rls = rls.replace("#MULT#", str(c.multiplicity))
            continue
        if "#ENRGY#" in rls:
            rls = rls.replace("#ENRGY#", str(c.energy))
            continue
        if "#PATH#" in rls:
            rls = rls.replace("#PATH#", str(c.path))
            continue
        if "#LABEL#" in rls:
            rls = rls.replace("#LABEL#", str(c.label))
            continue
        break
    rls = [_l + "\n" for _l in rls.split("\n")]
    return rls[:-1]


def replace_arg(c: System, lines: list[str], arg: dict[str, str]):
    rls = "".join(lines)
    while True:
        for key, val in arg.items():
            if f"#{key}#" in rls:
                rls = rls.replace(f"#{key}#", val)
                break
        else:
            break
    rls = [_l + "\n" for _l in rls.split("\n")]
    return rls[:-1]


def write_input(c: System, template: Path, odir=None, link=False, arg: dict[str, str] = None):
    tp = Path(template).resolve()
    if not tp.exists():
        logger.error(f"{str(tp)} not exist")
        raise ValueError
    with tp.open("r") as f:
        ls = f.readlines()
    ls = replace_key(c, ls)
    if arg is not None:
        arg = {str(_k): str(_v) for _k, _v in arg.items()}
        ls = replace_arg(c, ls, arg)
    p = change_dir(c.path, odir, c.name).with_suffix(tp.suffix)
    with p.open("w", newline="\n") as f:
        f.writelines(ls)
    logger.info(f"{p.name} created")
    if link:
        c.path = p
