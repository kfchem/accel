from typing import Any

from accel.util.log import logger


class FuncSelector:
    def __init__(self) -> None:
        self.funcs: dict[str, Any] = {}

    def add(self, filetype: str):
        def _add_dec(func):
            self.funcs[filetype] = func
            return func

        return _add_dec

    def __call__(self, filetype, *args, **kwargs):
        func = self.funcs.get(filetype)
        if func is None:
            logger.error(f"could not detect the file type {filetype}")
            raise ValueError
        return func(*args, **kwargs)


class Selectors:
    read_atoms = FuncSelector()
    read_energy = FuncSelector()
    read_thermal = FuncSelector()
    check_end = FuncSelector()
    check_freq = FuncSelector()
    calc_free_energy = FuncSelector()
    submit = FuncSelector()
    run = FuncSelector()
