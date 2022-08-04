from pathlib import Path
from typing import List

from accel.base.boxcore import BoxCore
from accel.base.text import replace_key
from accel.base.tools import change_dir
from accel.util.log import logger


def _line_formatter(new_line):
    if isinstance(new_line, str):
        ret_lines = [str(_l) + "\n" for _l in new_line.split("\n")]
    else:
        raise TypeError
    return ret_lines


class TxtBox(BoxCore):
    def read_text(self):
        for _c in self.get():
            with _c.path.open() as f:
                _ls = f.readlines()
            _c.data["txt"] = _ls
        return self

    def write_text(self, directory: Path = None, change_path: bool = False, suffix: str = None):
        for _c in self.get():
            _p = change_dir(_c.path, directory, _c.name)
            if suffix is not None:
                _p = _p.with_suffix(suffix)
            with _p.open("w", newline="\n") as f:
                f.writelines(_c.data["txt"])
            logger.info(f"{_p.name} created")
            if change_path:
                _c.path = _p
        return self

    def delete_lines_by_key(self, keyword: str):
        for _c in self.get():
            _ls: list[str] = _c.data["txt"]
            for i, _l in enumerate(_ls):
                if keyword in _l:
                    _ls.pop(i)
        return self

    def delete_lines(self, numbers: list[int]):
        for _c in self.get():
            _ls: list[str] = _c.data["txt"]
            for _i in numbers:
                _ls.pop(int(_i - 1))
        return self

    def insert_lines(self, number: int, new_lines=""):
        new_lines = _line_formatter(new_lines)
        for _c in self.get():
            _ls: list[str] = _c.data["txt"]
            _ls = _ls[: number - 1] + new_lines + _ls[number - 2 :]
        return self

    def append_lines(self, new_lines=""):
        new_lines = _line_formatter(new_lines)
        for _c in self.get():
            _ls: list[str] = _c.data["txt"]
            _ls.extend(new_lines)
        return self

    def replace_lines(self, number: int, new_lines=""):
        new_lines = _line_formatter(new_lines)
        for _c in self.get():
            _ls: list[str] = _c.data["txt"]
            _ls[number - 1] = "".join(new_lines)
        return self

    def parse_keys(self):
        for _c in self.get():
            _ls: list[str] = _c.data["txt"]
            _ls = replace_key(_c, _ls)
        return self
