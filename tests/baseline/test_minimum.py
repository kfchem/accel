"""only_minimum selection."""

from builders import box_of


def test_only_minimum_per_label():
    box = box_of([3.0, 1.0, 2.0, 7.0, 5.0], labels=["a", "a", "a", "b", "b"])
    box.only_minimum()
    assert [c.name for c in box.get()] == ["c2", "c5"]
    assert [c.history for c in box.contents] == ["only_minimum", "", "only_minimum", "only_minimum", ""]


def test_only_minimum_global():
    box = box_of([3.0, 1.0, 2.0, 0.5], labels=["a", "a", "b", "b"])
    box.only_minimum(in_label=False)
    assert [c.name for c in box.get()] == ["c4"]


def test_only_minimum_ignores_inactive_systems():
    box = box_of([0.0, 1.0, 2.0], labels=["m"] * 3)
    box.contents[0].deactivate("manual")
    box.only_minimum()
    assert [c.name for c in box.get()] == ["c2"]


def test_only_minimum_tie_keeps_first_in_order():
    # sorting is stable, so with equal energies the earlier system wins
    box = box_of([1.0, 1.0, 1.0], labels=["m"] * 3)
    box.only_minimum()
    assert [c.name for c in box.get()] == ["c1"]
