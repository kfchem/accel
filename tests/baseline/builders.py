"""Builders for the small synthetic inputs used by the baseline tests.

The program-output builders (Gaussian, ORCA, xTB) emit *minimal excerpts* that
contain only the lines ACCeL's current parsers look for, laid out the way those
programs print them. They are format fixtures, not real calculation outputs:
the numbers in them carry no scientific meaning, and passing tests do not prove
compatibility with every real output file.
"""

import math
from pathlib import Path

from accel import Box, System

ATOMIC_NUMBERS = {"H": 1, "C": 6, "N": 7, "O": 8}

# Idealised tetrahedral methane (C-H = 1.0896 A).
METHANE = [
    ("C", 0.0, 0.0, 0.0),
    ("H", 0.6291, 0.6291, 0.6291),
    ("H", -0.6291, -0.6291, 0.6291),
    ("H", -0.6291, 0.6291, -0.6291),
    ("H", 0.6291, -0.6291, -0.6291),
]

# Right-angle triatomic used for exact geometry checks (not a real water geometry).
RIGHT_ANGLE = [
    ("O", 0.0, 0.0, 0.0),
    ("H", 1.0, 0.0, 0.0),
    ("H", 0.0, 1.0, 0.0),
]


def rotated_permuted_methane(angle: float = 0.7, shift=(1.0, 2.0, 3.0)):
    """METHANE rotated about z, translated, with H2<->H3 and H4<->H5 swapped.

    It is the same structure as METHANE, so a symmetry-aware RMSD must be ~0.
    """
    c, s = math.cos(angle), math.sin(angle)
    atoms = []
    for i in (0, 2, 1, 4, 3):
        sym, x, y, z = METHANE[i]
        atoms.append((sym, c * x - s * y + shift[0], s * x + c * y + shift[1], z + shift[2]))
    return atoms


def stretched_methane(factor: float = 1.5):
    """METHANE with the first hydrogen pushed outward; clearly a different geometry."""
    atoms = list(METHANE)
    sym, x, y, z = atoms[1]
    atoms[1] = (sym, x * factor, y * factor, z * factor)
    return atoms


def xyz_text(atoms, comment: str = "accel test") -> str:
    lines = [str(len(atoms)), comment] + [f"{s} {x} {y} {z}" for s, x, y, z in atoms]
    return "\n".join(lines) + "\n"


def write_xyz(path: Path, atoms, comment: str = "accel test") -> Path:
    path.write_text(xyz_text(atoms, comment))
    return path


def system_with_atoms(name: str = "", atoms=(), energy=None, label: str = "") -> System:
    c = System()
    c.name = name
    c.label = label
    c.energy = energy
    for a in atoms:
        c.atoms.append(list(a))
    return c


def box_of(energies, labels=None, atoms=()) -> Box:
    """In-memory Box of Systems named c1, c2, ... with the given energies (kcal/mol)."""
    systems = []
    for i, e in enumerate(energies, start=1):
        label = "" if labels is None else labels[i - 1]
        systems.append(system_with_atoms(f"c{i}", atoms, energy=e, label=label))
    return Box(systems)


def _orientation_block(title: str, atoms, shift: float = 0.0):
    lines = [
        f"                         {title}                         ",
        " ---------------------------------------------------------------------",
        " Center     Atomic      Atomic             Coordinates (Angstroms)",
        " Number     Number       Type             X           Y           Z",
        " ---------------------------------------------------------------------",
    ]
    for i, (s, x, y, z) in enumerate(atoms, start=1):
        lines.append(f" {i:>6} {ATOMIC_NUMBERS[s]:>10} {0:>11} {x + shift:>15.6f}{y:>12.6f}{z:>12.6f}")
    lines.append(" ---------------------------------------------------------------------")
    return lines


def _g16_archive(atoms, charge: int, mult: int, scf: float):
    coords = "\\".join(f"{s},{x},{y},{z}" for s, x, y, z in atoms)
    arc = (
        "1\\1\\GINC-LOCALHOST\\SP\\RB3LYP\\STO-3G\\C1H4\\ACCEL\\01-Jan-2026\\0\\\\"
        "#P B3LYP/STO-3G\\\\accel test\\\\"
        f"{charge},{mult}\\{coords}\\\\Version=ES64L-G16RevC.01\\State=1-A1\\HF={scf}\\\\@"
    )
    return [" " + arc[i : i + 70] for i in range(0, len(arc), 70)]


def g16_log(
    atoms,
    scf: float = -40.0,
    charge: int = 0,
    mult: int = 1,
    method: str = "RB3LYP",
    thermal: dict = None,
    freqs=None,
    normal: bool = True,
    std_shift: float = 0.0,
) -> str:
    """Minimal Gaussian 16 log excerpt.

    ``std_shift`` offsets x in the Standard orientation block so tests can tell
    which geometry block a parser used.
    """
    lines = [
        " Entering Gaussian System, Link 0=g16",
        " Input=accel.gjf",
        " ----------------------",
        " #P B3LYP/STO-3G",
        " ----------------------",
        " Symbolic Z-matrix:",
        f" Charge = {charge} Multiplicity = {mult}",
    ]
    lines += _orientation_block("Input orientation:", atoms)
    lines += _orientation_block("Standard orientation:", atoms, shift=std_shift)
    lines.append(f" SCF Done:  E({method}) =  {scf:.10f}     A.U. after   10 cycles")
    if freqs is not None:
        lines += [
            "                      1                      2                      3",
            "                      A                      A                      A",
            f" Frequencies --  {freqs[0]:>10.4f}             {freqs[1]:>10.4f}             {freqs[2]:>10.4f}",
        ]
    if thermal is not None:
        lines += [
            f" Zero-point correction=                           {thermal['zpc']:.6f} (Hartree/Particle)",
            f" Thermal correction to Energy=                    {thermal['energy']:.6f}",
            f" Thermal correction to Enthalpy=                  {thermal['enthalpy']:.6f}",
            f" Thermal correction to Gibbs Free Energy=         {thermal['gibbs']:.6f}",
        ]
    lines += _g16_archive(atoms, charge, mult, scf)
    if normal:
        lines.append(" Normal termination of Gaussian 16 at Thu Jan  1 00:00:00 2026.")
    else:
        lines.append(" Error termination via Lnk1e in /opt/g16/l502.exe at Thu Jan  1 00:00:00 2026.")
    return "\n".join(lines) + "\n"


def g16_gjf(atoms, charge: int = 0, mult: int = 1, route: str = "# B3LYP/STO-3G opt") -> str:
    lines = ["%chk=accel.chk", route, "", "accel test", "", f"{charge} {mult}"]
    lines += [f"{s} {x} {y} {z}" for s, x, y, z in atoms]
    lines += ["", ""]
    return "\n".join(lines)


def orca_out(energy: float, normal: bool = True, optimized: bool = False) -> str:
    lines = [
        "",
        "                                 *****************",
        "                                 * O   R   C   A *",
        "                                 *****************",
        "",
    ]
    if optimized:
        lines.append("                    *** OPTIMIZATION RUN DONE ***")
    lines += [
        "-------------------------   --------------------",
        f"FINAL SINGLE POINT ENERGY       {energy:.12f}",
        "-------------------------   --------------------",
    ]
    if normal:
        lines.append("                             ****ORCA TERMINATED NORMALLY****")
    return "\n".join(lines) + "\n"


def xtb_out(total_energy: float, free_energy: float = None, thermal: dict = None, normal: bool = True) -> str:
    lines = [
        "      -----------------------------------------------------------      ",
        "     |                   =====================                   |     ",
        "     |                           x T B                           |     ",
        "     |                   =====================                   |     ",
        "      -----------------------------------------------------------      ",
        "",
    ]
    if thermal is not None:
        lines += [
            f"         :: total free energy        {thermal['total_free_energy']:.12f} Eh   ::",
            f"         :: zero point energy        {thermal['zpe']:.12f} Eh   ::",
            f"         :: G(RRHO) w/o ZPVE         {thermal['g_rrho_wo_zpve']:.12f} Eh   ::",
            f"         :: G(RRHO) contrib.         {thermal['g_rrho_contrib']:.12f} Eh   ::",
        ]
    lines += [
        "           -------------------------------------------------",
        f"          | TOTAL ENERGY             {total_energy:.12f} Eh   |",
    ]
    if free_energy is not None:
        lines.append(f"          | TOTAL FREE ENERGY        {free_energy:.12f} Eh   |")
    lines.append("           -------------------------------------------------")
    if normal:
        lines.append(" * finished run on 2026/01/01 at 00:00:00.000")
        lines.append(" normal termination of xtb")
    return "\n".join(lines) + "\n"
