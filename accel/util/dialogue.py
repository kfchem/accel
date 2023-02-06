import sys
from pathlib import Path

from accel.util.log import logger


class Dialog:
    @classmethod
    def directory(cls, title: str = "Input Directory") -> Path:
        print(title)
        p = input()
        if p is None or p == "":
            logger.error("Dialog returned None")
            sys.exit()
        logger.debug(f"Dialog returned {p}")
        return Path(p)

    @classmethod
    def file(cls, title="Input File", filetypes=[("any", "*.*")]) -> Path:
        print(f"{title}: filetypes {filetypes}")
        p = input()
        if p is None or p == "":
            logger.error("Dialog returned None")
            sys.exit()
        logger.debug(f"Dialog returned {p}")
        return Path(p)

    @classmethod
    def yesno(cls, title="Question", message="") -> bool:
        print(f"{title}: {message}: yes or no")
        answer = input()
        if answer.lower() in ("yes", "y"):
            _ret = True
        elif answer.lower() in ("no", "n"):
            _ret = False
        else:
            logger.error("Dialog returned None")
            sys.exit()
        logger.debug(f"Dialog returned {_ret}")
        return _ret

    @classmethod
    def value(cls, title="Question", message="") -> str:
        print(f"{title}: {message}: input value")
        answer = input()
        if answer != "" or answer is not None:
            logger.debug(f"Dialog returned {answer}")
            return answer
        else:
            logger.error("Dialog returned None")
            sys.exit()
