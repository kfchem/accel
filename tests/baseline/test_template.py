"""Template expansion in write_input (and TxtBox.parse_keys, which shares the engine)."""

import pytest
from builders import RIGHT_ANGLE, write_xyz

from accel import Box


def _axyz_line(sym, x, y, z):
    return f"{sym:<2} {x:>15} {y:>15} {z:>15}"


def _box(tmp_path, name="mol_a_1"):
    return Box(write_xyz(tmp_path / f"{name}.xyz", RIGHT_ANGLE)).read_atoms().labeling()


def test_all_builtin_keywords(tmp_path):
    box = _box(tmp_path)
    c = box.contents[0]
    c.energy = 1.25
    c.data["solvent"] = "water"
    c.charge = -1
    c.multiplicity = 2
    source_path = str(c.path)
    tpl = tmp_path / "tpl.gjf"
    tpl.write_text("#NAME#|#LABEL#|#ATOMS#|#CHG#|#MULT#|#ENRGY#|#PATH#|#DATA[solvent]#\n#AXYZ#\n")

    box.write_input(tpl, directory=tmp_path / "out")

    out = tmp_path / "out" / "mol_a_1.gjf"
    expected_axyz = "\n".join(
        [
            _axyz_line("O", "0.0", "0.0", "0.0"),
            _axyz_line("H", "1.0", "0.0", "0.0"),
            _axyz_line("H", "0.0", "1.0", "0.0"),
        ]
    )
    assert out.read_text() == f"mol_a_1|mol_a|3|-1|2|1.25|{source_path}|water\n{expected_axyz}\n"
    assert b"\r\n" not in out.read_bytes()


def test_arg_substitution_and_unknown_placeholders(tmp_path):
    box = _box(tmp_path)
    tpl = tmp_path / "tpl.inp"
    tpl.write_text("! #FUNCTIONAL# #BASIS#\n%pal nprocs #NPROC# end\n#UNKNOWN#\n")
    box.write_input(tpl, directory=tmp_path / "out", arg={"FUNCTIONAL": "B3LYP", "BASIS": "def2-SVP", "NPROC": 8})
    assert (tmp_path / "out" / "mol_a_1.inp").read_text() == "! B3LYP def2-SVP\n%pal nprocs 8 end\n#UNKNOWN#\n"


def test_output_name_suffix_and_relinking(tmp_path):
    box = _box(tmp_path)
    tpl = tmp_path / "anything.gjf"
    tpl.write_text("%chk=#NAME#.chk\n# hf/sto-3g\n\n#NAME#\n\n#CHG# #MULT#\n#AXYZ#\n\n")
    box.write_input(tpl)
    c = box.contents[0]
    assert c.path == (tmp_path / "mol_a_1.gjf").absolute()
    assert c.filetype == "app/g16/input"
    assert c.name == "mol_a_1"


def test_link_false_keeps_original_path(tmp_path):
    box = _box(tmp_path)
    before = box.contents[0].path
    tpl = tmp_path / "t.inp"
    tpl.write_text("#NAME#\n")
    box.write_input(tpl, directory=tmp_path / "out", link=False)
    assert box.contents[0].path == before
    assert (tmp_path / "out" / "mol_a_1.inp").exists()


def test_only_active_systems_are_written(tmp_path):
    box = Box([write_xyz(tmp_path / "a.xyz", RIGHT_ANGLE), write_xyz(tmp_path / "b.xyz", RIGHT_ANGLE)]).read_atoms()
    box.contents[1].deactivate("x")
    tpl = tmp_path / "t.inp"
    tpl.write_text("#NAME#\n")
    box.write_input(tpl, directory=tmp_path / "out")
    assert (tmp_path / "out" / "a.inp").exists()
    assert not (tmp_path / "out" / "b.inp").exists()


def test_trailing_newline_is_preserved_or_not_added(tmp_path):
    box = _box(tmp_path)
    tpl = tmp_path / "t.inp"
    tpl.write_text("#NAME#")
    box.write_input(tpl, directory=tmp_path / "out")
    assert (tmp_path / "out" / "mol_a_1.inp").read_text() == "mol_a_1"


def test_crlf_template_is_written_with_lf(tmp_path):
    box = _box(tmp_path)
    tpl = tmp_path / "t.inp"
    tpl.write_bytes(b"#NAME#\r\nline\r\n")
    box.write_input(tpl, directory=tmp_path / "out")
    assert (tmp_path / "out" / "mol_a_1.inp").read_bytes() == b"mol_a_1\nline\n"


def test_missing_template_raises_value_error(tmp_path):
    box = _box(tmp_path)
    with pytest.raises(ValueError):
        box.write_input(tmp_path / "missing.inp")


@pytest.mark.quirk
def test_enrgy_without_energy_renders_none(tmp_path):
    box = _box(tmp_path)
    tpl = tmp_path / "t.inp"
    tpl.write_text("#ENRGY#\n")
    box.write_input(tpl, directory=tmp_path / "out")
    assert (tmp_path / "out" / "mol_a_1.inp").read_text() == "None\n"


@pytest.mark.quirk
def test_data_placeholder_requires_string_value(tmp_path):
    box = _box(tmp_path).set_data("n", 3)
    tpl = tmp_path / "t.inp"
    tpl.write_text("#DATA[n]#\n")
    with pytest.raises(TypeError):
        box.write_input(tpl, directory=tmp_path / "out")


@pytest.mark.quirk
def test_txtbox_parse_keys_joins_lines(tmp_path):
    p = tmp_path / "mol_a_1.txt"
    p.write_text("a #NAME#\nb")
    box = Box(p)
    txt = box.plugin.txt.read_text().parse_keys()
    assert txt.contents[0].data["txt"] == "a mol_a_1b"
