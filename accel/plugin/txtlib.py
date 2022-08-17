from pathlib import Path

from accel.base.boxcore import BoxCore
from accel.base.formats import replace_key
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
        for c in self.get():
            with c.path.open() as f:
                ls = f.readlines()
            c.data["txt"] = ls
        return self

    def write_text(self, directory: Path = None, change_path: bool = False, suffix: str = None):
        for c in self.get():
            p = change_dir(c.path, directory, c.name)
            if suffix is not None:
                p = p.with_suffix(suffix)
            with p.open("w", newline="\n") as f:
                f.writelines(c.data["txt"])
            logger.info(f"{p.name} created")
            if change_path:
                c.path = p
        return self

    def delete_lines_by_key(self, keyword: str):
        for c in self.get():
            ls: list[str] = c.data["txt"]
            for i, line in enumerate(ls):
                if keyword in line:
                    ls.pop(i)
        return self

    def delete_lines(self, numbers: list[int]):
        for c in self.get():
            ls: list[str] = c.data["txt"]
            for num in numbers:
                ls.pop(int(num - 1))
        return self

    def insert_lines(self, number: int, new_lines=""):
        new_lines = _line_formatter(new_lines)
        for c in self.get():
            ls: list[str] = c.data["txt"]
            ls = ls[: number - 1] + new_lines + ls[number - 2 :]
        return self

    def append_lines(self, new_lines=""):
        new_lines = _line_formatter(new_lines)
        for c in self.get():
            ls: list[str] = c.data["txt"]
            ls.extend(new_lines)
        return self

    def replace_lines(self, number: int, new_lines=""):
        new_lines = _line_formatter(new_lines)
        for c in self.get():
            ls: list[str] = c.data["txt"]
            ls[number - 1] = "".join(new_lines)
        return self

    def parse_keys(self):
        for c in self.get():
            ls: list[str] = c.data["txt"]
            ls = replace_key(c, ls)
        return self
