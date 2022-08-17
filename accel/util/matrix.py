from copy import deepcopy
from typing import Any

import numpy as np
from accel.util.datadict import Data
from accel.util.log import logger


class Matrix:
    __slots__ = ["indexes", "_matrix", "symmetric", "data"]

    def __init__(self, indexes: list[Any]) -> None:
        self.indexes: list[Any] = indexes
        self._matrix: np.ndarray = None
        self.symmetric: bool = False
        self.data: Data = Data(self)

    def bind(self, ndarray: np.ndarray) -> "Matrix":
        if ndarray.shape == (len(self.indexes), len(self.indexes)):
            self._matrix = ndarray
            return self
        else:
            raise ValueError

    def __getitem__(self, key):
        if type(key) is tuple and len(key) == 2:
            number_a = int(key[0])
            number_b = int(key[1])
        else:
            raise KeyError
        return self._matrix[number_a - 1, number_b - 1]

    def __setitem__(self, key, value):
        if type(key) is tuple and len(key) == 2:
            number_a = int(key[0])
            number_b = int(key[1])
        else:
            raise KeyError
        self._matrix[number_a - 1, number_b - 1] = value
        self._matrix[number_b - 1, number_a - 1] = value

    def __len__(self):
        return len(self.indexes)

    def index(self, index):
        try:
            return self.indexes.index(index)
        except ValueError:
            logger.error(f"could not find {index} in the matrix")
            raise ValueError

    @classmethod
    def zeros(cls, indexes: list[Any]):
        return Matrix(indexes).bind(np.zeros((len(indexes), len(indexes))))

    @classmethod
    def ones(cls, indexes: list[Any]):
        return Matrix(indexes).bind(np.ones((len(indexes), len(indexes))))

    @classmethod
    def identity(cls, indexes: list[Any]):
        return Matrix(indexes).bind(np.identity(len(indexes)))

    def ordered(self, indexes: list) -> "Matrix":
        try:
            idx_list = [self.indexes.index(idx) for idx in indexes]
        except ValueError:
            logger.error("could not resolve indexes in the matrix")
            raise ValueError
        return Matrix(indexes).bind(self._matrix[idx_list][:, idx_list])

    def to_ndarray(self) -> np.ndarray:
        return deepcopy(self._matrix)

    def swap(self, index_a, index_b):
        idx_a = self.index(index_a)
        idx_b = self.index(index_b)
        self._matrix[idx_a], self._matrix[idx_b] = self._matrix[idx_b], self._matrix[idx_a]
        for arr in self._matrix:
            arr[idx_a], arr[idx_b] = arr[idx_b], arr[idx_a]

    def insert(self, index, value=None):
        idx = self.indexes.index(index)
        self._matrix = np.insert(self._matrix, idx, value, axis=0)
        self._matrix = np.insert(self._matrix, idx, value, axis=1)

    def delete(self, index):
        idx = self.indexes.index(index)
        self._matrix = np.delete(self._matrix, idx, axis=0)
        self._matrix = np.delete(self._matrix, idx, axis=1)
