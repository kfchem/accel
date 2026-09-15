"""Energy processing: relative energies, Boltzmann populations, energy sums, units."""

import math

import pytest
from builders import METHANE, box_of

from accel.util.constants import Units

# CODATA 2018 reference values, written out independently of accel.util.constants.
HARTREE_IN_KCAL_MOL = 627.5094740631
EV_IN_KCAL_MOL = 23.060547830619
GAS_CONSTANT = 8.314462618  # J / (mol K)


def boltzmann(energies_kcal, temperature=298.15):
    rel = [e - min(energies_kcal) for e in energies_kcal]
    w = [math.exp(-r * 4184.0 / (GAS_CONSTANT * temperature)) for r in rel]
    return [x / sum(w) for x in w]


def test_unit_conversions():
    assert Units.hartree(1.0).to_kcal_mol == pytest.approx(HARTREE_IN_KCAL_MOL, rel=1e-10)
    assert Units.eV(1.0).to_kcal_mol == pytest.approx(EV_IN_KCAL_MOL, rel=1e-9)
    assert Units.hartree(1.0).to_eV == pytest.approx(27.211386245988, rel=1e-12)
    assert Units.kJ_mol(4.184).to_kcal_mol == pytest.approx(1.0, rel=1e-12)
    assert Units.kcal_mol(1.0).to_J_mol == pytest.approx(4184.0, rel=1e-12)


def test_calc_rel_energy_per_label_overwrites_energy():
    box = box_of([10.0, 12.5, 3.0, 4.0], labels=["a", "a", "b", "b"])
    box.calc_rel_energy()
    assert [c.energy for c in box.contents] == [0.0, 2.5, 0.0, 1.0]


def test_calc_rel_energy_global():
    box = box_of([10.0, 12.5, 3.0, 4.0], labels=["a", "a", "b", "b"])
    box.calc_rel_energy(in_label=False)
    assert [c.energy for c in box.contents] == [7.0, 9.5, 0.0, 1.0]


def test_calc_rel_energy_ignores_inactive_systems():
    box = box_of([1.0, 2.0, 5.0], labels=["a", "a", "a"])
    box.contents[0].deactivate("x")
    box.calc_rel_energy()
    assert [c.energy for c in box.contents] == [1.0, 0.0, 3.0]


def test_calc_distribution_matches_boltzmann_formula():
    energies = [0.0, 0.5, 1.2, 3.0]
    box = box_of(energies, labels=["m"] * 4)
    box.calc_distribution()
    got = [c.distribution for c in box.contents]
    assert got == pytest.approx(boltzmann(energies), rel=1e-9)
    assert sum(got) == pytest.approx(1.0, rel=1e-12)
    # relative energies are stashed in the (transient) cache
    assert [c.cache["relative_energy"] for c in box.contents] == [0.0, 0.5, 1.2, 3.0]


def test_calc_distribution_temperature():
    energies = [0.0, 1.0]
    box = box_of(energies, labels=["m", "m"])
    box.calc_distribution(temperature=500.0)
    assert [c.distribution for c in box.contents] == pytest.approx(boltzmann(energies, 500.0), rel=1e-9)


def test_calc_distribution_normalises_per_label():
    box = box_of([0.0, 1.0, 100.0, 101.0], labels=["a", "a", "b", "b"])
    box.calc_distribution()
    got = [c.distribution for c in box.contents]
    assert got[0] + got[1] == pytest.approx(1.0)
    assert got[2] + got[3] == pytest.approx(1.0)
    assert got[0] == pytest.approx(got[2])


def test_calc_distribution_global_and_inactive():
    box = box_of([0.0, 1.0, 2.0], labels=["a", "b", "b"])
    box.contents[2].deactivate("x")
    box.calc_distribution(in_label=False)
    got = [c.distribution for c in box.contents]
    assert got[:2] == pytest.approx(boltzmann([0.0, 1.0]), rel=1e-9)
    assert got[2] is None


def test_calc_energy_sums_data_keys_and_converts_units():
    box = box_of([None, None])
    for c, (scf, corr) in zip(box.contents, [(-40.0, 0.027), (-39.9, 0.030)]):
        c.data["scf"] = scf
        c.data["corr"] = corr
    box.calc_energy(keys=["scf", "corr"], unit=Units.hartree)
    assert box.contents[0].energy == pytest.approx((-40.0 + 0.027) * HARTREE_IN_KCAL_MOL, rel=1e-12)
    assert box.contents[1].energy == pytest.approx((-39.9 + 0.030) * HARTREE_IN_KCAL_MOL, rel=1e-12)


def test_get_average_is_boltzmann_weighted_per_label():
    energies = [0.0, 1.0, 5.0]
    box = box_of(energies, labels=["m", "m", "n"], atoms=METHANE)
    for c, v in zip(box.contents, [1.0, 3.0, 7.0]):
        c.data["x"] = v
    averaged = box.get_average(keys=["x"])
    w = boltzmann([0.0, 1.0])
    assert [c.name for c in averaged] == ["m", "n"]
    assert [c.label for c in averaged] == ["m", "n"]
    assert averaged[0].data["x"] == pytest.approx(1.0 * w[0] + 3.0 * w[1], rel=1e-9)
    assert averaged[1].data["x"] == pytest.approx(7.0)
    assert averaged[0].energy == 0.0
    assert len(averaged[0].atoms) == len(METHANE)
