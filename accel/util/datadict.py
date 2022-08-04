from collections.abc import MutableMapping
from copy import deepcopy
from typing import Any

from accel.util.log import logger


class Data(MutableMapping):
    __slots__ = ["_data", "parent"]

    def __init__(self, parent=None):
        self._data: dict[str, Any] = {}
        self.parent = parent

    def set_data(self, key, value):
        logger.debug(f"{str(self.parent)}: {key}: {value}")
        return self._data.__setitem__(key, value)

    def __setitem__(self, key, value):
        return self.set_data(key, value)

    def __getitem__(self, key) -> Any:
        return self._data[key]

    def __delitem__(self, key):
        del self._data[key]

    def __iter__(self):
        return self._data.__iter__()

    def __len__(self):
        return len(self._data)

    def duplicate(self, parent=None):
        n = Data()
        n._data = deepcopy(self._data)
        n.parent = parent
        return n
