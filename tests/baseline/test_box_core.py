"""Box collection behaviour: construction, chaining, state, labels, metadata."""

import csv

import pytest
from builders import METHANE, box_of, system_with_atoms, write_xyz

from accel import Box, Systems
from accel.plugin.gaulib import GauBox


def _xyz_dir(tmp_path, names):
    for n in names:
        write_xyz(tmp_path / f"{n}.xyz", METHANE)
    return tmp_path


def test_box_from_glob_is_sorted(tmp_path):
    _xyz_dir(tmp_path, ["mol_b_1", "mol_a_2", "mol_a_1"])
    box = Box(tmp_path / "*.xyz")
    assert [c.name for c in box.contents] == ["mol_a_1", "mol_a_2", "mol_b_1"]
    assert all(c.filetype == "format/xyz" for c in box.contents)


def test_box_from_directory_and_path_list(tmp_path):
    _xyz_dir(tmp_path, ["a", "b"])
    assert len(Box(tmp_path)) == 2
    assert len(Box([tmp_path / "a.xyz", str(tmp_path / "b.xyz")])) == 2


def test_box_accepts_system_objects():
    systems = [system_with_atoms("x"), system_with_atoms("y")]
    box = Box(systems)
    assert box.contents[0] is systems[0]
    assert len(Box(Systems(systems))) == 2


@pytest.mark.quirk
def test_box_of_box_shares_contents():
    box = box_of([0.0, 1.0])
    alias = Box(box)
    assert alias.contents is box.contents
    assert box.plugin.gau.contents is box.contents
    assert isinstance(box.plugin.gau, GauBox)


def test_add_box_adds_only_active_systems_by_reference():
    src = box_of([0.0, 1.0])
    src.contents[1].deactivate("x")
    dst = Box().add(src)
    assert len(dst) == 1
    assert dst.contents[0] is src.contents[0]


CHAIN_METHODS = [
    ("labeling", {}),
    ("set_data", {"key": "k", "value": 1}),
    ("set_state", {"state": True}),
    ("set_label", {"label": "x"}),
    ("zero_fill", {}),
    ("count", {}),
    ("show", {}),
    ("energy_limit", {}),
    ("calc_rel_energy", {}),
    ("calc_distribution", {}),
    ("calc_energy", {"keys": []}),
    ("only_minimum", {}),
    ("calc_bonds", {}),
    ("calc_symm", {}),
    ("rmsd_limit", {}),
    ("calc_length", {"number_a": 1, "number_b": 2}),
    ("calc_angle", {"number_a": 2, "number_b": 1, "number_c": 3}),
    ("calc_dihedral", {"number_a": 2, "number_b": 1, "number_c": 3, "number_d": 4}),
    ("convert_to_mirror", {"centering": False}),
    ("modify_length", {"number_a": 1, "number_b": 2, "target": 1.2}),
]


@pytest.mark.parametrize("name,kwargs", CHAIN_METHODS, ids=[m for m, _ in CHAIN_METHODS])
def test_chain_methods_return_the_same_box(name, kwargs):
    box = box_of([0.0, 1.0], labels=["m", "m"], atoms=METHANE)
    assert getattr(box, name)(**kwargs) is box


def test_get_variants():
    box = box_of([0.0, 1.0, 2.0], labels=["a", "b", "a"])
    box.contents[2].deactivate("x")
    assert [c.name for c in box.get()] == ["c1", "c2"]
    assert [c.name for c in box.get(True)] == ["c1", "c2"]
    assert [c.name for c in box.get(False)] == ["c3"]
    assert [c.name for c in box.get(None)] == ["c1", "c2", "c3"]
    # a string selects a label and includes inactive systems
    assert [c.name for c in box.get("a")] == ["c1", "c3"]
    assert str(box) == "Box: active/all = 2/3"
    assert len(box) == 3


def test_labeling_default_and_custom():
    box = Box([system_with_atoms("KEF20958_a_256"), system_with_atoms("x")])
    box.labeling()
    assert [c.label for c in box.contents] == ["KEF20958_a", "x"]
    box.labeling(separator="_", index_list=[0])
    assert [c.label for c in box.contents] == ["KEF20958", "x"]


def test_bulk_setters_apply_to_inactive_systems_too():
    box = box_of([0.0, 1.0])
    box.contents[1].deactivate("x")
    box.set_label("L").set_data("k", "v")
    assert [c.label for c in box.contents] == ["L", "L"]
    assert [c.data["k"] for c in box.contents] == ["v", "v"]
    box.set_state(True)
    assert [c.state for c in box.contents] == [True, True]
    # history is kept when re-activated
    assert box.contents[1].history == "x"


def test_zero_fill():
    box = Box([system_with_atoms("KEF20958_a_2"), system_with_atoms("KEF20958_a_12")])
    box.zero_fill(digit=3)
    assert [c.name for c in box.contents] == ["KEF20958_a_002", "KEF20958_a_012"]


def test_zero_fill_skips_inactive():
    box = Box([system_with_atoms("m_a_2"), system_with_atoms("m_a_3")])
    box.contents[1].deactivate("x")
    box.zero_fill()
    assert [c.name for c in box.contents] == ["m_a_002", "m_a_3"]


def test_duplicate_is_independent():
    box = box_of([0.0, 1.0], atoms=METHANE)
    dup = box.duplicate()
    assert isinstance(dup, Box)
    dup.contents[0].energy = 9.0
    dup.contents[1].deactivate("x")
    assert box.contents[0].energy == 0.0
    assert box.contents[1].state is True


def test_export_data_csv(tmp_path):
    box = box_of([1.5, 2.0], labels=["m", "m"])
    box.contents[0].data["note"] = "hello"
    box.contents[1].data["long"] = "x" * 300
    box.contents[1].deactivate("why")
    box.export_data(tmp_path / "result.txt")
    with (tmp_path / "result.csv").open(newline="") as f:
        rows = list(csv.reader(f))
    assert rows[0] == [
        "name",
        "path",
        "filetype",
        "state",
        "history",
        "label",
        "energy",
        "distribution",
        "charge",
        "multiplicity",
        "note",
        "long",
    ]
    assert rows[1] == ["c1", "", "", "True", "", "m", "1.5", "", "0", "1", "hello", ""]
    assert rows[2] == ["c2", "", "", "False", "why", "m", "2.0", "", "0", "1", "", "######"]


@pytest.mark.quirk
def test_search_without_directory_or_suffix_raises(tmp_path):
    box = Box(write_xyz(tmp_path / "a.xyz", METHANE))
    with pytest.raises(UnboundLocalError):
        box.search()


def test_search_by_suffix(tmp_path):
    box = Box(write_xyz(tmp_path / "a.xyz", METHANE))
    (tmp_path / "a.log").write_text("")
    (tmp_path / "b.xyz").write_text("")
    box.add(tmp_path / "b.xyz")
    box.search(suffix=".log")
    assert box.contents[0].path.name == "a.log"
    assert box.contents[1].state is False
    assert box.contents[1].history == "search_dir"
