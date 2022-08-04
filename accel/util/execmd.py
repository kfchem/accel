from pathlib import Path
from typing import Union

from accel.util.log import logger


class Execmd:
    _path_dict: dict[str, str] = {}

    @classmethod
    def get(cls, key: str) -> str:
        exe_string = cls._path_dict.get(str(key))
        if exe_string is None:
            return key
        else:
            return exe_string

    @classmethod
    def add(cls, key: str, exe_cmd: Union[Path, str]) -> None:
        cls._path_dict[str(key)] = str(exe_cmd)
        logger.info(f"set executable file: {key}: {cls._path_dict[str(key)]}")
