from pathlib import Path
from typing import Dict, Union

from accel.util.log import logger


class Execmd:
    _path_dict: dict[str, str] = {}

    @classmethod
    def get(cls, key: str) -> str:
        _exe = cls._path_dict.get(str(key))
        if _exe is None:
            return key
        else:
            return _exe

    @classmethod
    def add(cls, key: str, exe_cmd: Union[Path, str]) -> None:
        cls._path_dict[str(key)] = str(exe_cmd)
        logger.info(f"set executable file: {key}: {cls._path_dict[str(key)]}")
