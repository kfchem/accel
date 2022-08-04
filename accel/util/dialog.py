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
            icon_data = "iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAAACXBIWXMAAACZAAAAmQFwj2B8AAAAGXRFWHRTb2Z0d2FyZQB3d3cuaW\
                5rc2NhcGUub3Jnm+48GgAAAUBJREFUOI1jYBgFFANGGOM/AwPja0/tbFE2ltsQkb/XGDdeeYyu4X+oJSfDnx9mDH//cn5l+BcKN+Cdt3b577+M2\
                8V3XLnEwMDA8N9XX5uBiUEG087/3xgEP53+/pYn9cf//w8gQoFGmt999VbgcuZ/BwcOFH6AgcIvP73j/xkYGFn+Oziw/Pr7bhkLy/8sBgYGBoaE\
                Bg4+HsZMdh5uS1Ye3nN//jOsZqxNv4vs1d///s39859hEhsDw3/G7wEGbUf/cTu6/FXdwMD4z42BgUmdmZX5Ep+C0iUOPj71/6wsjoz//j/4///\
                /BhYWxnWPT82y+cfwv4SJTU2VcfXqvyySv7QkPzCxMjD8+8vLwMxWycD56ezf1av/vodZ2dDAIsEmZf3/54+Av3+ZV6rrxb6d9XBfssOi1X8JRN\
                AoIBYAAK0Tda2iIfAiAAAAAElFTkSuQmCC"
            cls.root.iconphoto(False, PhotoImage(data=icon_data))

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
