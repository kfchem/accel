from collections.abc import MutableMapping
from copy import deepcopy
from typing import Any, Dict

from accel.util.log import logger


class Data(MutableMapping):
    __slots__ = ["_data", "_parent"]

    def __init__(self, parent_obj=None):
        self._data: dict[str, Any] = {}
        self._parent = parent_obj

    def set_data(self, key, value):
        logger.debug(f"{str(self._parent)}: {key}: {value}")
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

    def duplicate(self, parent_obj=None):
        _n = Data()
        _n._data = deepcopy(self._data)
        _n._parent = parent_obj
        return _n
