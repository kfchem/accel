import sys
from pathlib import Path
from tkinter import PhotoImage, Tk, filedialog, messagebox, simpledialog

from accel.util.log import logger


class Dialog:
    root = None

    @classmethod
    def initialize(cls):
        if not cls.root:
            cls.root = Tk()
            cls.root.withdraw()
            cls.root.iconphoto(False, PhotoImage(file=Path(__file__).parent.joinpath("logo32.png")))

    @classmethod
    def directory(cls, title: str = "Select Directory") -> Path:
        cls.initialize()
        p = filedialog.askdirectory(title=title)
        if p is None or p == "":
            logger.error("Dialog returned None")
            sys.exit()
        logger.debug(f"Dialog returned {p}")
        return Path(p)

    @classmethod
    def file(cls, title="Select File", filetypes=[("any", "*.*")]) -> Path:
        cls.initialize()
        p = filedialog.askopenfilename(title=title, filetypes=filetypes)
        if p is None or p == "":
            logger.error("Dialog returned None")
            sys.exit()
        logger.debug(f"Dialog returned {p}")
        return Path(p)

    @classmethod
    def yesno(cls, title="Question", message="") -> bool:
        cls.initialize()
        answer = messagebox.askquestion(title=title, message=message)
        if answer == "yes":
            _ret = True
        elif answer == "no":
            _ret = False
        else:
            logger.error("Dialog returned None")
            sys.exit()
        logger.debug(f"Dialog returned {_ret}")
        return _ret

    @classmethod
    def value(cls, title="Question", message="") -> str:
        cls.initialize()
        answer = simpledialog.askstring(title=title, prompt=message)
        if answer != "" or answer is not None:
            logger.debug(f"Dialog returned {answer}")
            return answer
        else:
            logger.error("Dialog returned None")
            sys.exit()
