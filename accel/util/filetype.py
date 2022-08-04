from pathlib import Path
from typing import Any

from accel.util.log import logger


class FileType:
    funcs: dict[str, Any] = {}
    metric: dict[str, int] = {}

    @classmethod
    def add(cls, filetype: str, metric: int = 50):
        def _add_dec(func):
            cls.funcs[filetype] = func
            cls.metric[filetype] = min([int(metric), 100])
            return func

        return _add_dec

    @classmethod
    def analyse(cls, filepath: Path):
        for filetype, _ in sorted(cls.metric.items(), key=lambda x: x[1]):
            if cls.funcs[filetype](filepath):
                return filetype
        logger.error(f"could not detect filetype of {filepath.name}")
