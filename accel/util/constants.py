import pickle
from pathlib import Path
from typing import Any, Dict, List, Union

from accel.util.log import logger


class Elements:
    with Path(__file__).parent.joinpath("elements.pkl").open("rb") as f:
        # the following data were taken from
        # https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page) and references therein
        data: List[Dict[str, Union[str, float]]] = pickle.load(f)
    symbols: List[str] = [_e["symbol"] for _e in data]
    lower_symbols = [sym.lower() for sym in symbols]

    @classmethod
    def canonicalize(cls, input):
        if isinstance(input, int) and input < len(cls.symbols):
            canonical_symbol = cls.symbols[input]
            # logger.debug('resolved interger {} as atomic number'.format(input))
        elif isinstance(input, str):
            try:
                input = str(input).strip()
                canonical_symbol = cls.symbols[cls.lower_symbols.index(input.lower())]
            except ValueError:
                logger.error(f"could not find atom symbol {input}")
                canonical_symbol = ""
        else:
            logger.error(f"could not find atom symbol {input}")
            canonical_symbol = ""
        return canonical_symbol

    @classmethod
    def get_element(cls, input: Union[int, str]) -> Dict[str, Any]:
        if isinstance(input, int) and input < len(cls.data):
            return cls.data[input]
        elif isinstance(input, str):
            return cls.data[cls.lower_symbols.index(input.lower())]


# part of the following constants values were taken from https://physics.nist.gov/cuu/Constants/
class CONSTANTS:
    hartree_to_J = 4.3597447222071 * (10 ** -18)
    hartree_to_eV = 27.211386245988
    cal_to_J = 4.184
    boltzman_constant_in_J_over_K = 1.380649 * (10 ** -23)
    avogadros_constant_in_over_mol = 6.02214076 * (10 ** 23)
    gas_constant_in_J_over_mol_K = boltzman_constant_in_J_over_K * avogadros_constant_in_over_mol


class Converter:
    def __init__(self, value_in_J_mol):
        self.value_in_J_mol = float(value_in_J_mol)

    @property
    def to_hartree(self):
        return self.value_in_J_mol / Units.hartree.scale_to_J_mol

    @property
    def to_kJ_mol(self):
        return self.value_in_J_mol / Units.kJ_mol.scale_to_J_mol

    @property
    def to_J_mol(self):
        return self.value_in_J_mol / Units.J_mol.scale_to_J_mol

    @property
    def to_kcal_mol(self):
        return self.value_in_J_mol / Units.kcal_mol.scale_to_J_mol

    @property
    def to_eV(self):
        return self.value_in_J_mol / Units.eV.scale_to_J_mol


class Unit:
    def __init__(self, scale_to_J_mol):
        self.scale_to_J_mol = float(scale_to_J_mol)

    def __call__(self, value):
        return Converter(value * self.scale_to_J_mol)


class Units:
    hartree = Unit(CONSTANTS.hartree_to_J * CONSTANTS.avogadros_constant_in_over_mol)
    kJ_mol = Unit(10 ** 3)
    kcal_mol = Unit(CONSTANTS.cal_to_J * 10 ** 3)
    J_mol = Unit(1)
    eV = Unit(CONSTANTS.hartree_to_J * CONSTANTS.avogadros_constant_in_over_mol / CONSTANTS.hartree_to_eV)
