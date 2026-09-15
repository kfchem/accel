"""Filetype detection and program-output parsers on minimal synthetic excerpts.

These are format-characterisation tests: they pin which lines and columns the
current parsers read. They do not validate against real program outputs.
"""

import pytest
from builders import METHANE, g16_gjf, g16_log, orca_out, stretched_methane, write_xyz, xtb_out

from accel import Box

HARTREE_IN_KCAL_MOL = 627.5094740631
THERMAL = {"zpc": 0.044793, "energy": 0.047650, "enthalpy": 0.048594, "gibbs": 0.027228}


def _one(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return Box(p)


def _xyz(c):
    return [a.xyz for a in c.atoms]


def _flat(rows):
    return [v for row in rows for v in row]


@pytest.mark.parametrize(
    "name,text,filetype",
    [
        ("a.xyz", "1\n\nH 0 0 0\n", "format/xyz"),
        ("a.sdf", "", "format/mol"),
        ("a.mol", "", "format/mol"),
        ("a.gjf", g16_gjf(METHANE), "app/g16/input"),
        ("a.com", g16_gjf(METHANE), "app/g16/input"),
        ("a.log", g16_log(METHANE), "app/g16/output"),
        ("a.inp", "! B3LYP def2-SVP\n* xyz 0 1\nH 0 0 0\n*\n", "app/orca/input"),
        ("a.out", orca_out(-40.0), "app/orca/output"),
        ("b.out", xtb_out(-4.0), "app/xtb/output"),
        ("a.sh", "#!/bin/bash\n#PBS -l nodes=1\n", "app/pbs/jobscript"),
        ("a.mae", "", "format/mae"),
        ("a.txt", "hello\n", None),
    ],
)
def test_filetype_detection(tmp_path, name, text, filetype):
    assert _one(tmp_path, name, text).contents[0].filetype == filetype


def test_unknown_filetype_deactivates_on_dispatch(tmp_path):
    box = _one(tmp_path, "a.txt", "hello\n").read_energy()
    assert box.contents[0].state is False
    assert box.contents[0].history == "could not find an appropriate function"


def test_mixed_filetypes_dispatch_per_system(tmp_path):
    (tmp_path / "a.log").write_text(g16_log(METHANE, scf=-40.0))
    (tmp_path / "b.out").write_text(orca_out(-41.0))
    (tmp_path / "c.out").write_text(xtb_out(-4.0))
    box = Box([tmp_path / "a.log", tmp_path / "b.out", tmp_path / "c.out"]).read_energy()
    got = [c.energy for c in box.contents]
    assert got == pytest.approx([-40.0 * HARTREE_IN_KCAL_MOL, -41.0 * HARTREE_IN_KCAL_MOL, -4.0 * HARTREE_IN_KCAL_MOL])


# --- Gaussian ---------------------------------------------------------------


def test_g16_read_atoms_from_archive(tmp_path):
    box = _one(tmp_path, "m.log", g16_log(METHANE, charge=1, mult=2, std_shift=1.0)).read_atoms()
    c = box.contents[0]
    assert [a.symbol for a in c.atoms] == ["C", "H", "H", "H", "H"]
    assert _xyz(c) == [[x, y, z] for _, x, y, z in METHANE]
    assert (c.charge, c.multiplicity) == (1, 2)


def test_g16_read_atoms_standard_and_input_orientation(tmp_path):
    box = _one(tmp_path, "m.log", g16_log(METHANE, charge=-1, mult=3, std_shift=1.0))
    box.read_atoms(source="standard_orientation")
    c = box.contents[0]
    assert _flat(_xyz(c)) == pytest.approx(_flat([[x + 1.0, y, z] for _, x, y, z in METHANE]))
    assert (c.charge, c.multiplicity) == (-1, 3)
    box.read_atoms(source="input_orientation")
    assert _flat(_xyz(c)) == pytest.approx(_flat([[x, y, z] for _, x, y, z in METHANE]))


def test_g16_read_energy(tmp_path):
    box = _one(tmp_path, "m.log", g16_log(METHANE, scf=-40.5183)).read_energy()
    c = box.contents[0]
    assert c.data["g16_scf"] == -40.5183
    assert c.energy == pytest.approx(-40.5183 * HARTREE_IN_KCAL_MOL, rel=1e-12)


def test_g16_read_energy_with_calculation_method(tmp_path):
    ok = _one(tmp_path, "a.log", g16_log(METHANE, scf=-40.1, method="RB3LYP"))
    ok.read_energy(calculation_method="RB3LYP")
    assert ok.contents[0].data["g16_scf"] == -40.1
    ng = _one(tmp_path, "b.log", g16_log(METHANE, scf=-40.1, method="RB3LYP"))
    ng.read_energy(calculation_method="UB3LYP")
    assert ng.contents[0].state is False
    assert ng.contents[0].history == "read_scf"


def test_g16_read_thermal_and_free_energy(tmp_path):
    box = _one(tmp_path, "m.log", g16_log(METHANE, scf=-40.5, thermal=THERMAL))
    box.read_energy().read_thermal().calc_free_energy()
    d = box.contents[0].data
    assert (d["g16_zpc"], d["g16_corr_to_energy"], d["g16_corr_to_enthalpy"], d["g16_corr_to_gibbs"]) == (
        0.044793,
        0.047650,
        0.048594,
        0.027228,
    )
    assert box.contents[0].energy == pytest.approx((-40.5 + 0.027228) * HARTREE_IN_KCAL_MOL, rel=1e-12)


def test_g16_check_end(tmp_path):
    (tmp_path / "ok.log").write_text(g16_log(METHANE))
    (tmp_path / "ng.log").write_text(g16_log(METHANE, normal=False))
    box = Box(tmp_path / "*.log").check_end()
    assert [(c.name, c.state, c.history) for c in box.contents] == [("ng", False, "check_end"), ("ok", True, "")]


def test_g16_check_end_counts_jobs_from_sibling_input(tmp_path):
    # "opt" and "freq" in one route are counted as two jobs; one Normal termination is not enough
    (tmp_path / "m.gjf").write_text(g16_gjf(METHANE, route="# B3LYP/STO-3G opt freq"))
    box = _one(tmp_path, "m.log", g16_log(METHANE)).check_end()
    assert box.contents[0].state is False


def test_g16_check_freq(tmp_path):
    log = g16_log(METHANE, freqs=(-120.5, 200.0, 300.0))
    box0 = _one(tmp_path, "a.log", log).check_freq()
    assert box0.contents[0].data["g16_img_freq"] == 1
    assert box0.contents[0].state is False
    box1 = _one(tmp_path, "b.log", log).check_freq(im=1)
    assert box1.contents[0].state is True


def test_g16_missing_scf_deactivates(tmp_path):
    text = "\n".join(line for line in g16_log(METHANE).split("\n") if "SCF Done" not in line)
    box = _one(tmp_path, "m.log", text).read_energy()
    assert box.contents[0].state is False


def test_g16_read_atoms_from_gjf(tmp_path):
    box = _one(tmp_path, "m.gjf", g16_gjf(METHANE, charge=1, mult=2)).read_atoms()
    c = box.contents[0]
    assert _xyz(c) == [[x, y, z] for _, x, y, z in METHANE]
    assert (c.charge, c.multiplicity) == (1, 2)


# --- ORCA -------------------------------------------------------------------


def test_orca_read_energy_and_check_end(tmp_path):
    (tmp_path / "ok.out").write_text(orca_out(-40.25))
    (tmp_path / "ng.out").write_text(orca_out(-40.25, normal=False))
    box = Box(tmp_path / "*.out").check_end().read_energy()
    assert [(c.name, c.state) for c in box.contents] == [("ng", False), ("ok", True)]
    assert box.contents[1].energy == pytest.approx(-40.25 * HARTREE_IN_KCAL_MOL, rel=1e-12)


def test_orca_read_atoms_uses_sibling_xyz(tmp_path):
    (tmp_path / "a.out").write_text(orca_out(-40.0))
    write_xyz(tmp_path / "a.xyz", stretched_methane())
    (tmp_path / "b.out").write_text(orca_out(-40.0))
    box = Box(tmp_path / "*.out").read_atoms()
    assert _xyz(box.contents[0]) == [[x, y, z] for _, x, y, z in stretched_methane()]
    assert box.contents[1].state is False
    assert box.contents[1].history == "read_atoms: orca from xyz file: not exists"


def test_orca_check_optimized_plugin(tmp_path):
    (tmp_path / "a.out").write_text(orca_out(-40.0, optimized=True))
    (tmp_path / "b.out").write_text(orca_out(-40.0))
    box = Box(tmp_path / "*.out")
    box.plugin.orc.check_optimized()
    assert [c.state for c in box.contents] == [True, False]


# --- xTB --------------------------------------------------------------------


def test_xtb_energy_thermal_and_free_energy(tmp_path):
    thermal = {"total_free_energy": -4.15, "zpe": 0.044, "g_rrho_wo_zpve": 0.001, "g_rrho_contrib": 0.02}
    box = _one(tmp_path, "m.out", xtb_out(-4.175218109986, free_energy=-4.155, thermal=thermal))
    box.check_end().read_energy().read_thermal()
    c = box.contents[0]
    assert c.state is True
    assert c.energy == pytest.approx(-4.175218109986 * HARTREE_IN_KCAL_MOL, rel=1e-12)
    assert (c.data["xtb_total_free_energy"], c.data["xtb_zpe"]) == (-4.15, 0.044)
    assert (c.data["xtb_g_rrho_wo_zpve"], c.data["xtb_g_rrho_contrib"]) == (0.001, 0.02)
    box.plugin.xtb.read_free_energy()
    assert c.energy == pytest.approx(-4.155 * HARTREE_IN_KCAL_MOL, rel=1e-12)


def test_xtb_check_end_abnormal(tmp_path):
    box = _one(tmp_path, "m.out", xtb_out(-4.0, normal=False)).check_end()
    assert box.contents[0].history == "check_end"
