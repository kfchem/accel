from pathlib import Path

from accel.base.systems import System
from accel.base.tools import change_dir, float_to_str
from accel.util.log import logger


def replace_key(system: System, lines: list[str]):
    rls = "".join(lines)
    while True:
        if "#NAME#" in rls:
            rls = rls.replace("#NAME#", system.name)
            continue
        if "#DATA[" in rls:
            key = rls.split("#DATA[")[1].split("]#")[0]
            rls = rls.replace("#DATA[" + key + "]#", system.data[key])
            continue
        if "#AXYZ#" in rls:
            _xyz = [
                "{:<2} {:>15} {:>15} {:>15}".format(
                    _a.symbol,
                    float_to_str(_a.x),
                    float_to_str(_a.y),
                    float_to_str(_a.z),
                )
                for _a in system.atoms
            ]
            rls = rls.replace("#AXYZ#", "\n".join(_xyz))
            continue
        if "#ATOMS#" in rls:
            rls = rls.replace("#ATOMS#", str(len(system.atoms)))
            continue
        if "#CHG#" in rls:
            rls = rls.replace("#CHG#", str(system.charge))
            continue
        if "#MULT#" in rls:
            rls = rls.replace("#MULT#", str(system.multiplicity))
            continue
        if "#ENRGY#" in rls:
            rls = rls.replace("#ENRGY#", str(system.energy))
            continue
        if "#PATH#" in rls:
            rls = rls.replace("#PATH#", str(system.path))
            continue
        if "#LABEL#" in rls:
            rls = rls.replace("#LABEL#", str(system.label))
            continue
        break
    rls = [_l + "\n" for _l in rls.split("\n")]
    return rls[:-1]


def replace_arg(system: System, lines: list[str], arg: dict[str, str]):
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


def write_input(system: System, template: Path, odir=None, link=False, arg: dict[str, str] = None):
    tp = Path(template).resolve()
    if not tp.exists():
        logger.error(f"{str(tp)} not exist")
        raise ValueError
    with tp.open("r") as f:
        ls = f.readlines()
    ls = replace_key(system, ls)
    if arg is not None:
        arg = {str(_k): str(_v) for _k, _v in arg.items()}
        ls = replace_arg(system, ls, arg)
    p = change_dir(system.path, odir, system.name).with_suffix(tp.suffix)
    with p.open("w", newline="\n") as f:
        f.writelines(ls)
    logger.info(f"{p.name} created")
    if link:
        system.path = p
