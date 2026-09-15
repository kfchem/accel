"""System / Systems data-model baseline."""

import pytest
from builders import METHANE, system_with_atoms, write_xyz

from accel import System, Systems


def test_system_defaults():
    c = System()
    assert c.name == ""
    assert c.path is None
    assert c.filetype is None
    assert c.label == ""
    assert c.state is True
    assert c.history == ""
    assert c.energy is None
    assert c.total_charge is None
    assert c.charge == 0
    assert c.multiplicity == 1
    assert c.distribution is None
    assert len(c.atoms) == 0
    assert len(c.data) == 0


def test_system_setattr_coerces_types():
    c = System()
    c.energy = "1.5"
    c.state = 0
    c.label = 3
    c.name = 12
    c.total_charge = "-1"
    c.multiplicity = "2"
    c.distribution = "0.25"
    assert c.energy == 1.5 and isinstance(c.energy, float)
    assert c.state is False
    assert c.label == "3"
    assert c.name == "12"
    assert c.total_charge == -1
    assert c.multiplicity == 2
    assert c.distribution == 0.25


def test_charge_falls_back_to_sum_of_atomic_formal_charges():
    c = system_with_atoms("ion_pair", [("N", 0, 0, 0), ("O", 2, 0, 0), ("C", 4, 0, 0)])
    c.atoms[0].charge = 1
    c.atoms[1].charge = -1
    assert c.total_charge is None
    assert c.charge == 0
    c.atoms[1].charge = None
    assert c.charge == 1


def test_explicit_charge_overrides_atomic_charges():
    c = system_with_atoms("x", [("N", 0, 0, 0)])
    c.atoms[0].charge = 1
    c.charge = -2
    assert c.total_charge == -2
    assert c.charge == -2


def test_deactivate_records_reason_history():
    c = System()
    c.deactivate("first").deactivate("second")
    assert c.state is False
    assert c.history == "first; second"


def test_system_from_file_sets_name_path_filetype(tmp_path):
    p = write_xyz(tmp_path / "mol_a_1.xyz", METHANE)
    c = System(p)
    assert c.name == "mol_a_1"
    assert c.path == p.resolve()
    assert c.filetype == "format/xyz"
    assert c.state is True


def test_system_from_missing_file_is_inactive(tmp_path):
    c = System(tmp_path / "missing.xyz")
    assert c.state is False
    assert c.history == "file not found"


def test_setting_path_redetects_filetype(tmp_path):
    p = write_xyz(tmp_path / "a.xyz", METHANE)
    q = tmp_path / "a.sdf"
    q.write_text("")
    c = System(p)
    c.path = q
    assert c.filetype == "format/mol"
    assert c.name == "a"


def test_duplicate_is_deep_copy():
    c = system_with_atoms("x", METHANE, energy=1.0, label="m")
    c.data["k"] = [1, 2]
    c.charge = 1
    c.multiplicity = 2
    d = c.duplicate()
    assert (d.name, d.label, d.energy, d.charge, d.multiplicity, d.state) == ("x", "m", 1.0, 1, 2, True)
    d.atoms[0].x = 9.0
    d.data["k"].append(3)
    assert c.atoms[0].x == 0.0
    assert c.data["k"] == [1, 2]


@pytest.mark.quirk
def test_duplicate_drops_filetype_history_and_distribution(tmp_path):
    c = System(write_xyz(tmp_path / "a.xyz", METHANE))
    c.distribution = 0.5
    c.deactivate("reason")
    d = c.duplicate()
    assert d.path == c.path
    assert d.state is False
    assert d.filetype is None
    assert d.history == ""
    assert d.distribution is None


def _systems():
    a = system_with_atoms("a_1", energy=2.0, label="a")
    b = system_with_atoms("b_1", energy=1.0, label="b")
    c = system_with_atoms("a_2", energy=5.0, label="a")
    u = system_with_atoms("u", energy=0.5)
    return Systems([a, b, c, u])


def test_labels_groups_and_sorts_keys():
    labels = _systems().labels
    assert list(labels.keys()) == ["", "a", "b"]
    assert [c.name for c in labels["a"]] == ["a_1", "a_2"]
    assert isinstance(labels["a"], Systems)


def test_has_state():
    ss = _systems()
    ss[1].deactivate("x")
    assert [c.name for c in ss.has_state(True)] == ["a_1", "a_2", "u"]
    assert [c.name for c in ss.has_state(False)] == ["b_1"]
    assert len(ss.has_state(None)) == 4


def test_has_label():
    ss = _systems()
    assert [c.name for c in ss.has_label()] == ["a_1", "b_1", "a_2"]
    assert [c.name for c in ss.has_label("a")] == ["a_1", "a_2"]
    assert len(ss.has_label("missing")) == 0


def test_has_energy_bounds_are_inclusive():
    ss = _systems()
    ss.append(System())
    assert len(ss.has_energy()) == 4
    assert [c.name for c in ss.has_energy(1.0, 2.0)] == ["a_1", "b_1"]
    assert [c.name for c in ss.has_energy(min_limit=2.0)] == ["a_1", "a_2"]


def test_has_data():
    ss = _systems()
    ss[0].data["k"] = 1
    ss[1].data["k"] = 2
    assert [c.name for c in ss.has_data("k")] == ["a_1", "b_1"]
    assert [c.name for c in ss.has_data("k", 2)] == ["b_1"]


def test_get_by_index_name_and_default():
    ss = _systems()
    assert ss.get(1).name == "a_1"
    assert ss.get(4).name == "u"
    assert ss.get(0) is None
    assert ss.get(5) is None
    assert ss.get("b_").name == "b_1"
    assert ss.get("zzz") is None
    # with energies on every system, get() returns the lowest-energy one
    assert ss.get().name == "u"
    assert Systems().get() is None


def test_get_default_falls_back_to_name_order_without_energies():
    ss = Systems([system_with_atoms("b"), system_with_atoms("a")])
    assert ss.get().name == "a"


def test_sorted():
    ss = _systems()
    assert [c.name for c in ss.sorted()] == ["u", "b_1", "a_1", "a_2"]
    assert [c.name for c in ss.sorted("name")] == ["a_1", "a_2", "b_1", "u"]
    assert [c.label for c in ss.sorted("label")] == ["", "a", "a", "b"]


def test_add_concatenates():
    ss = _systems()
    total = ss + Systems([System()])
    assert len(total) == 5
    assert total[0] is ss[0]


def test_views_share_system_objects():
    ss = _systems()
    view = ss.has_label("a")
    view[0].energy = 42.0
    assert ss[0].energy == 42.0
