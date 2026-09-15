"""Typical end-to-end Box method chains (highest-priority baseline).

These chains mirror the README quick start. Program outputs are minimal
synthetic excerpts (see builders.py); no external program is executed.
"""

from builders import (
    METHANE,
    g16_log,
    orca_out,
    rotated_permuted_methane,
    stretched_methane,
    write_xyz,
)

from accel import Box

# Hartree energies chosen so that, relative to conformer 1:
#   2 -> +0.31 kcal/mol (same geometry as 1: removed by rmsd_limit)
#   3 -> +0.63 kcal/mol (different geometry: kept)
#   4 -> +6.28 kcal/mol (removed by energy_limit at 3 kcal/mol)
CONFORMERS = [
    ("mol_a_1", METHANE, -40.000000),
    ("mol_a_2", rotated_permuted_methane(), -39.999500),
    ("mol_a_3", stretched_methane(), -39.999000),
    ("mol_a_4", METHANE, -39.990000),
]


def _states(box):
    return [(c.name, c.state, c.history) for c in box.contents]


EXPECTED_STATES = [
    ("mol_a_1", True, ""),
    ("mol_a_2", False, "rmsd_limit"),
    ("mol_a_3", True, ""),
    ("mol_a_4", False, "energy_limit"),
]


def test_readme_chain_with_orca_outputs(tmp_path):
    for name, atoms, energy in CONFORMERS:
        (tmp_path / f"{name}.out").write_text(orca_out(energy))
        write_xyz(tmp_path / f"{name}.xyz", atoms)
    template = tmp_path / "template.inp"
    template.write_text("! B3LYP def2-SVP Opt\n* xyz #CHG# #MULT#\n#AXYZ#\n*\n")

    box = (
        Box(str(tmp_path / "*.out"))
        .check_end()
        .read_atoms()
        .read_energy()
        .energy_limit()
        .rmsd_limit()
        .write_input(str(template))
    )

    assert _states(box) == EXPECTED_STATES
    assert (tmp_path / "mol_a_1.inp").exists()
    assert (tmp_path / "mol_a_3.inp").exists()
    assert not (tmp_path / "mol_a_2.inp").exists()
    assert not (tmp_path / "mol_a_4.inp").exists()
    # write_input re-links active systems to the generated inputs
    active = box.get()
    assert [c.path.name for c in active] == ["mol_a_1.inp", "mol_a_3.inp"]
    assert [c.filetype for c in active] == ["app/orca/input", "app/orca/input"]
    text = (tmp_path / "mol_a_1.inp").read_text()
    assert text.startswith("! B3LYP def2-SVP Opt\n* xyz 0 1\nC ")
    assert text.endswith("\n*\n")


def test_readme_chain_with_gaussian_outputs(tmp_path):
    for name, atoms, energy in CONFORMERS:
        (tmp_path / f"{name}.log").write_text(g16_log(atoms, scf=energy, charge=1, mult=2))
    template = tmp_path / "template.gjf"
    template.write_text("%chk=#NAME#.chk\n# B3LYP/6-31G(d) opt freq\n\n#NAME#\n\n#CHG# #MULT#\n#AXYZ#\n\n")

    box = Box(tmp_path / "*.log").check_end().read_atoms().read_energy().energy_limit().rmsd_limit()
    box.write_input(template)

    assert _states(box) == EXPECTED_STATES
    active = box.get()
    assert [c.path.name for c in active] == ["mol_a_1.gjf", "mol_a_3.gjf"]
    assert [c.filetype for c in active] == ["app/g16/input", "app/g16/input"]
    lines = (tmp_path / "mol_a_1.gjf").read_text().split("\n")
    assert lines[:6] == ["%chk=mol_a_1.chk", "# B3LYP/6-31G(d) opt freq", "", "mol_a_1", "", "1 2"]


def test_labelled_chain_filters_per_label(tmp_path):
    energies = {"mol_a_1": -40.0, "mol_a_2": -39.990, "mol_b_1": -80.0, "mol_b_2": -79.999}
    for name, energy in energies.items():
        (tmp_path / f"{name}.out").write_text(orca_out(energy))
        write_xyz(tmp_path / f"{name}.xyz", stretched_methane(1.2 if name.endswith("2") else 1.5))

    box = Box(tmp_path / "*.out").read_atoms().read_energy().labeling().energy_limit(3.0)
    assert [c.label for c in box.contents] == ["mol_a", "mol_a", "mol_b", "mol_b"]
    assert [c.state for c in box.contents] == [True, False, True, True]

    box.only_minimum()
    assert [c.name for c in box.get()] == ["mol_a_1", "mol_b_1"]


def test_analysis_only_chain_on_existing_gaussian_outputs(tmp_path):
    thermal = {"zpc": 0.044, "energy": 0.047, "enthalpy": 0.048, "gibbs": 0.027}
    for name, energy in [("m_x_1", -40.0), ("m_x_2", -39.9995)]:
        (tmp_path / f"{name}.log").write_text(
            g16_log(METHANE, scf=energy, thermal=thermal, freqs=(100.0, 200.0, 300.0))
        )

    box = (
        Box(tmp_path / "*.log")
        .check_end()
        .check_freq()
        .read_atoms()
        .read_energy()
        .read_thermal()
        .calc_free_energy()
        .labeling()
        .calc_distribution()
    )
    assert [c.state for c in box.contents] == [True, True]
    dist = [c.distribution for c in box.contents]
    assert abs(sum(dist) - 1.0) < 1e-12
    assert dist[0] > dist[1]
