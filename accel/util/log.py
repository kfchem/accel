import atexit
import datetime
import logging
import tempfile
from pathlib import Path

import numpy as np
from accel import __version__


class Log:
    logger = logging.getLogger("accel")
    logger.setLevel(logging.DEBUG)

    input_dir: Path = None
    output_dir: Path = None

    @classmethod
    def set_input_dir(cls, directory):
        if directory is None:
            return None
        directory = Path(directory)
        if directory.exists() and directory.is_dir():
            cls.input_dir = directory

    @classmethod
    def set_output_dir(cls, directory):
        if directory is None:
            return None
        directory = Path(directory)
        if directory.exists() and directory.is_dir():
            cls.output_dir = directory

    temp_file = tempfile.TemporaryFile("a+")
    temp_handler = logging.StreamHandler(temp_file)
    temp_handler.setLevel(logging.DEBUG)
    temp_handler.setFormatter(logging.Formatter("%(asctime)s: %(levelname)s: %(funcName)s: %(message)s"))
    logger.addHandler(temp_handler)
    logger.info(f"ACCeL version: {__version__}")
    logger.debug(f"Numpy version:  {np.__version__}")

    stderr_handler = logging.StreamHandler()
    stderr_handler.setLevel("DEBUG")
    stderr_handler.setFormatter(logging.Formatter("%(funcName)s: %(message).119s"))
    logger.addHandler(stderr_handler)

    file_handler: logging.FileHandler = None

    @classmethod
    def console(cls, show=True, level="DEBUG"):
        if show is True:
            cls.stderr_handler.setLevel(level)
            if cls.stderr_handler not in cls.logger.handlers:
                cls.logger.addHandler(cls.stderr_handler)
                cls.logger.debug("stderr activated")
        elif show is False:
            if cls.stderr_handler in cls.logger.handlers:
                cls.logger.removeHandler(cls.stderr_handler)
                cls.logger.debug("stderr deactivated")

    @classmethod
    def flush(cls, filepath: Path = None):
        if cls.temp_handler is None:
            cls.logger.error("logs were already flushed")
            return None
        cls.logger.removeHandler(cls.temp_handler)
        cls.temp_handler = None
        cls.temp_file.seek(0)
        _ls = cls.temp_file.readlines()
        cls.temp_file.close()
        with filepath.open("a") as f:
            f.writelines(_ls)

    @classmethod
    def file(cls, filepath=None, level="DEBUG"):
        if filepath is None:
            _name = datetime.date.today().strftime("%Y%m%d")
            _p = Path().cwd().joinpath(f"{_name}.mcl")
        else:
            _p = Path(filepath)
            if not (_p.parent.exists()):
                _p.parent.mkdir(parents=True)
        if cls.file_handler is None:
            cls.flush(Path(_p))
        elif cls.file_handler in cls.logger.handlers:
            cls.logger.removeHandler(cls.file_handler)
        cls.file_handler = logging.FileHandler(str(_p))
        cls.file_handler.setLevel(logging.DEBUG)
        cls.file_handler.setFormatter(logging.Formatter("%(asctime)s: %(levelname)s: %(funcName)s: %(message)s"))
        cls.logger.addHandler(cls.file_handler)
        cls.logger.debug(f"set: {str(_p.absolute())}")

    @classmethod
    def write(cls, text=""):
        logger.info(f"log: {text}")

    @classmethod
    def force_flush(cls):
        if cls.file_handler is not None:
            return None
        if cls.input_dir is None and cls.output_dir is None:
            cls.logger.removeHandler(cls.temp_handler)
            cls.temp_handler = None
            cls.temp_file.close()
            return None
        _dir = Path.cwd()
        if cls.input_dir is not None:
            _dir = cls.input_dir
        if cls.output_dir is not None:
            _dir = cls.output_dir
        _dir = Path(_dir)
        _name = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
        _p = _dir / f"{_name}.mcl"
        _p = _p.absolute()
        cls.logger.info(f"log: {str(_p)}")
        cls.flush(_p)


atexit.register(Log.force_flush)

logger = Log.logger
