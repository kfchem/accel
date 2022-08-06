from pathlib import Path
from typing import Iterable, Union

from accel.base.boxcore import BoxCore
from accel.base.selector import FuncSelector, Selectors
from accel.base.systems import System
from accel.plugin import Plugins
from accel.util.log import logger


def adaptive_function_caller(box: "Box", selector: FuncSelector, filetype, **options):
    if filetype is not None:
        try:
            selector(filetype, Box().bind(box.get()), **options)
        except ValueError:
            for c in box.get():
                c.deactivate("could not find an appropriate function")
    else:
        for ft, confs in box.get().filetypes.items():
            try:
                selector(ft, Box().bind(confs), **options)
            except ValueError:
                for c in confs:
                    c.deactivate("could not find an appropriate function")
    logger.debug(f"done: {str(box)}")


class Box(BoxCore):
    def __init__(self, contents: Iterable[Union[System, Path, str]] = None):
        if isinstance(contents, BoxCore):
            self.contents = contents.contents
            self.data = contents.data
        else:
            super().__init__(contents=contents)

    @property
    def plugin(self):
        return Plugins(self)

    # function selector
    def check_end(self, filetype=None, **options):
        adaptive_function_caller(self, Selectors.check_end, filetype, **options)
        return self

    def read_atoms(self, filetype=None, **options):
        adaptive_function_caller(self, Selectors.read_atoms, filetype, **options)
        return self

    def read_energy(self, filetype=None, **options):
        adaptive_function_caller(self, Selectors.read_energy, filetype, **options)
        return self

    def read_thermal(self, filetype=None, **options):
        adaptive_function_caller(self, Selectors.read_thermal, filetype, **options)
        return self

    def check_freq(self, filetype=None, **options):
        adaptive_function_caller(self, Selectors.check_freq, filetype, **options)
        return self

    def run(self, filetype=None, **options):
        adaptive_function_caller(self, Selectors.run, filetype, **options)
        return self

    def submit(self, filetype=None, **options):
        adaptive_function_caller(self, Selectors.submit, filetype, **options)
        return self

    def calc_free_energy(self, filetype=None, **options):
        adaptive_function_caller(self, Selectors.calc_free_energy, filetype, **options)
        return self
