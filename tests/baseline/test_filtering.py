"""energy_limit filtering and active-state semantics."""

import pytest
from builders import box_of


def _states(box):
    return [c.state for c in box.contents]


def test_energy_limit_threshold_is_exclusive():
    # a system exactly at the threshold is removed (relative energy >= threshold)
    box = box_of([0.0, 2.999, 3.0, 4.0], labels=["m"] * 4)
    box.energy_limit(3.0)
    assert _states(box) == [True, True, False, False]
    assert [c.history for c in box.contents] == ["", "", "energy_limit", "energy_limit"]


def test_energy_limit_default_threshold_is_3_kcal_mol():
    box = box_of([0.0, 2.9, 3.1], labels=["m"] * 3)
    box.energy_limit()
    assert _states(box) == [True, True, False]


def test_energy_limit_reference_is_per_label_minimum():
    box = box_of([10.0, 12.0, 14.0, 0.0, 1.0, 5.0], labels=["a", "a", "a", "b", "b", "b"])
    box.energy_limit(3.0)
    assert _states(box) == [True, True, False, True, True, False]


def test_energy_limit_global():
    box = box_of([10.0, 12.0, 14.0, 0.0, 1.0, 5.0], labels=["a", "a", "a", "b", "b", "b"])
    box.energy_limit(3.0, in_label=False)
    assert _states(box) == [False, False, False, True, True, False]


def test_unlabelled_systems_form_one_group():
    box = box_of([10.0, 0.0, 2.0])
    box.energy_limit(3.0)
    assert _states(box) == [False, True, True]


def test_energy_limit_max_limit_keeps_n_lowest_per_label():
    box = box_of([0.3, 0.0, 0.2, 0.1, 5.0, 0.0], labels=["a", "a", "a", "a", "b", "b"])
    box.energy_limit(10.0, max_limit=2)
    assert _states(box) == [False, True, False, True, True, True]
    assert box.contents[0].history == "max limit in energy_limit"


def test_energy_limit_ignores_inactive_systems_for_reference():
    box = box_of([0.0, 2.0, 4.0, 6.0], labels=["m"] * 4)
    box.contents[0].deactivate("manual")
    box.energy_limit(3.0)
    assert _states(box) == [False, True, True, False]
    assert box.contents[0].history == "manual"


def test_energy_limit_accepts_string_threshold():
    box = box_of([0.0, 2.0, 4.0], labels=["m"] * 3)
    box.energy_limit("3")
    assert _states(box) == [True, True, False]


def test_filtering_is_cumulative_and_never_removes_systems():
    box = box_of([0.0, 1.0, 2.0, 4.0], labels=["m"] * 4)
    box.energy_limit(3.0).energy_limit(1.5)
    assert len(box) == 4
    assert _states(box) == [True, True, False, False]


@pytest.mark.quirk
def test_energy_limit_with_missing_energy_raises():
    box = box_of([0.0, None], labels=["m", "m"])
    with pytest.raises(TypeError):
        box.energy_limit()
