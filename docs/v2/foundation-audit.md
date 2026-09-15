# ACCeL v2.0 foundation audit (as of v0.3.1)

This document records an audit of the current code base (v0.3.1, commit `d746482`) and the
design questions that must be answered before the v2.0 work (`Structure` / `Flow` /
incremental recomputation) starts.

**This document does not settle any specification.** Anything that touches scientific
meaning, the behaviour existing Box users depend on, Flow UX, or structure identity is
recorded as an open question (see section 14) and left for an explicit decision.

File references use `path:line` (line numbers as of v0.3.1).

---

## 1. Current architecture

```text
accel/__init__.py            exports Box, System, Systems, Atom, Atoms
accel/base/
  atoms.py      Atom / Atoms / Bonds / BondType     atoms, coordinates, bonds, geometry
  systems.py    System / Systems                    one structure + state + metadata, and lists of them
  boxcore.py    BoxCore                             the chainable collection API (~40 methods)
  box.py        Box(BoxCore)                        filetype-dispatched methods + .plugin
  selector.py   FuncSelector / Selectors            registry: filetype -> parser / execution function
  formats.py    xyz / mol IO, template engine
  topology.py   rmsdpruning, map_numbers
  modeler.py    Modeler (bonds, aromaticity, stereo, symmetry matrices, mapping), Chains
accel/plugin/
  gaulib.py (Gaussian)  orclib.py (ORCA)  xtblib.py (xTB)  maelib.py (Maestro)
  pbslib.py (PBS/qsub)  txtlib.py (text editing)
accel/util/
  constants.py (Elements[pickle], CONSTANTS, Units)  datadict.py (Data)  filetype.py (FileType)
  log.py (global logger, writes .mcl at exit)  execmd.py (executable name map)  matrix.py  dialog*.py
```

Central mechanisms:

- **Registration at import time.** `@FileType.add(...)` and `@Selectors.<op>.add(filetype)` in each
  plugin module populate global registries when the module is imported
  (`box.py` -> `accel.plugin` -> every plugin).
- **Filetype detection.** Assigning `System.path` calls `FileType.analyse`, which opens the file and
  returns the first type whose detector returns True (detectors are tried in ascending metric order).
  This file I/O happens inside `__setattr__` (`systems.py:62-105`).
- **Dispatch.** `Box.read_energy()` and friends go through `adaptive_function_caller` (`box.py:11`),
  which groups active Systems by filetype and calls the plugin's **unbound** method with a plain Box
  as `self`: `selector(ft, Box().bind(confs))`. With no matching function the systems are deactivated
  with the reason `could not find an appropriate function`.
- **State.** Every operation mutates `System` objects in place and returns `self`.

## 2. Responsibilities of System / Systems

`System` (`systems.py:14`) holds, via `__slots__`:

| Attribute | Nature | Notes |
|---|---|---|
| `atoms: Atoms` | scientific data | elements, coordinates, atomic formal charges, bonds (`atoms.bonds`), stereo |
| `total_charge` / `charge` | scientific data | `charge` falls back to the sum of atomic formal charges when `total_charge` is None (`systems.py:107`) |
| `multiplicity` | scientific data | **defaults to 1**, even when the source carries no such information |
| `energy` | result | a single scalar in kcal/mol; **whether it is an SCF energy, a free energy or a relative value is not recorded** |
| `data: Data` | results + arbitrary metadata | `g16_scf`, thermal corrections, vibrations, `rotamer` matrices, raw `sdf` blocks, `jobid`, ... all mixed together |
| `path`, `name`, `filetype` | provenance / file | `name` defaults to the file stem; assigning `path` re-detects the filetype |
| `label` | grouping | the key for per-molecule (conformer-set) operations |
| `state`, `history` | workflow state | active/inactive plus a string log of the reasons |
| `distribution` | ensemble-dependent result | Boltzmann population, depends on the whole label group |
| `cache` | scratch | `relative_energy`, `distr_factor`, ... |

`System.__setattr__` additionally coerces types, logs at DEBUG level, re-detects the filetype and
updates the global `Log.input_dir`. There is no notion of identity beyond Python object identity
and `name`.

`Systems` (`systems.py:185`) is a `MutableSequence` providing **views** that share the same System
objects (`has_state/has_label/has_energy/has_data/has_filetype/has_bonds/has_distribution`),
grouping (`labels` / `filetypes`, returned as sorted dicts), `get()` (None = lowest energy or first
by name, int = 1-based index, str = substring match), `sorted()` and `duplicate()`.

## 3. Responsibilities of Box

`BoxCore` (`boxcore.py:37`) is a collection wrapper around `contents: Systems` plus a `data` object.
Responsibilities are mixed:

| Category | Methods |
|---|---|
| Ingestion | `__init__` (glob/dir/path/System/Systems/Box), `add`, `bind` |
| State, labels, metadata | `get`, `labeling`, `set_label`, `set_state`, `set_data`, `zero_fill`, `count`, `show`, `duplicate` |
| Energy analysis | `energy_limit`, `calc_rel_energy`, `calc_distribution`, `calc_energy`, `only_minimum`, `get_average` |
| Structure algorithms | `calc_bonds`, `calc_stereo`, `calc_symm`, `rmsd_limit`, `map_numbers`, `calc_length/angle/dihedral`, `modify_length`, `convert_to_mirror` |
| File IO | `read_xyz`, `write_xyz`, `read_mol`, `write_mol`, `write_input`, `export_data`, `copy_files`, `search` |
| Filetype dispatch (Box) | `read_atoms`, `read_energy`, `read_thermal`, `check_end`, `check_freq`, `calc_free_energy`, `run`, `submit` |
| Plugins | `.plugin.gau/.orc/.xtb/.mae/.pbs/.txt` (a plugin Box sharing the same contents) |

Conventions that make chaining work:

- **Every method returns `self`** (except `get`, `get_average`, `get_resub`, `get_irc`, `get_trj`, `get_unzip`).
- **Operations apply to `self.get()` (active systems only)** - except `labeling`, `set_label`,
  `set_data` and `set_state`, which apply to all contents including inactive ones (`boxcore.py:83-106`).
- **Filtering deactivates instead of removing**, appending a reason to `history`. The full population
  is kept, so `show()` / `export_data()` can show why each structure was dropped.
- **`in_label=True` is the default**: energy operations and `only_minimum` work per label (per molecule).
  Systems without a label form a single `""` group.
- **Prerequisites run automatically**: `rmsd_limit` -> `calc_symm(calc_all=False)` -> `calc_bonds`,
  `write_mol` -> `calc_bonds`, `get_average` -> `calc_distribution`, `GauBox.check_bonding` -> `read_vibration`.
- **`Box(box)` and `.plugin.*` alias the same contents** (`boxcore.py:39`), and `add(box)` appends the
  active Systems **by reference**.

## 4. The template mechanism

Implementation: `formats.write_input` (`formats.py:205`) -> `replace_key` (`:224`) -> `replace_arg` (`:275`).

1. The template path is resolved; a missing file raises `ValueError`. The file is re-read for every System.
2. `replace_key` loops over the whole text until no keyword remains (replacement order is fixed):

| Keyword | Replaced with | Notes |
|---|---|---|
| `#NAME#` | `c.name` | |
| `#DATA[key]#` | `c.data[key]` | **raises TypeError unless the value is a str** (pinned by a quirk test) |
| `#AXYZ#` | `"{:<2} {:>15} {:>15} {:>15}"` (symbol, `float_to_str(x/y/z)`) | precision is Python's repr, no rounding |
| `#ATOMS#` | number of atoms | |
| `#CHG#` | `c.charge` (including the fallback to atomic formal charges) | |
| `#MULT#` | `c.multiplicity` (defaults to 1) | |
| `#ENRGY#` | `str(c.energy)` (kcal/mol; `"None"` when unset) | unit and kind are implicit |
| `#PATH#` | `c.path` **before** the write | |
| `#LABEL#` | `c.label` | |

3. If `arg: dict` is given, `#KEY#` is replaced with `str(value)` (user-defined keywords).
4. Output path: `change_dir(c.path, directory, c.name)` plus the template's suffix; directories are
   created as needed and newlines are always LF.
5. With `link=True` (the default), `c.path` is re-pointed at the generated input file, so the filetype
   is re-detected as an input format and a following `.run()` / `.submit()` dispatches on it.

Other properties: unknown `#XXX#` placeholders are silently left in place; a replacement value that
itself contains a keyword is expanded again (a self-referencing value can loop forever); and the
template content, its parameters and the expansion result are not recorded anywhere.
`TxtBox.parse_keys` (`txtlib.py:77`) reuses `replace_key` but joins the lines together (pinned by a
quirk test). PBS job scripts can also be produced by `write_input` (a `.sh` file containing `#PBS` is
detected as `app/pbs/jobscript`) and `Box.submit()` then dispatches to `qsub` - that is, calculation
inputs and execution scripts go through exactly the same mechanism.

## 5. Structure of calculation execution

| Target | Implementation | Behaviour |
|---|---|---|
| Gaussian | `gaulib.run/submit` (`gaulib.py:370,382`) | `run`: blocking `subprocess.run([g16, path], cwd)`. `submit`: fire-and-forget `Popen` |
| ORCA | `orclib.run/submit` (`orclib.py:85,101`) | same pattern; stdout is captured via PIPE and logged (ACCeL itself does not write the `.out` file; the meaning of the second argument is unverified) |
| PBS / OpenPBS | `pbslib.que_submit/que_wait` (`pbslib.py:12,24`) | stores the `qsub` stdout in `data["jobid"]`, then polls `qstat` every 10 s and treats a job as finished once its id disappears from the listing. No timeout and no failure detection |
| Slurm | **not implemented** | |
| Local parallelism | none (the `Popen` in `submit` is the de facto parallelism) | |
| Script generation | no dedicated mechanism (templates are used instead) | |
| Retry / error handling | `GauBox.check_resub` (`gaulib.py:181`) | classifies error signatures in the log (link9999 with oscillation detection, SCF convergence failure, memory, ...), regenerates a `.gjf` from the last Standard orientation, renames the old files with a `_` suffix, sets `data["resubmission"]` |
| Command resolution | `Execmd` (`execmd.py`) | global key -> executable mapping |

Problems: without `check=True`, `subprocess.run` never raises `CalledProcessError`, so **failures go
undetected**; `submit` keeps no process handle, so completion cannot be tracked; and the link between
an execution and its output is left to the user's workflow (run the next script over the output files).
Conversely, this means that typical v0.3 usage is already centred on *analysing existing outputs*, with
execution as a thin layer on top.

## 6. Parser structure

- All parsers scan the file line by line, re-reading the whole file on every call.
- Results are **written straight into the System** (`atoms`, `energy`, `charge`, `multiplicity`,
  `data[...]`), and failures call `deactivate`. **Parsing, validation and selection are fused together.**
- Energies are converted to kcal/mol at read time and overwrite `energy`. Gaussian keeps the raw
  hartree value in `data["g16_scf"]`, but ORCA and xTB discard the raw value.

| filetype | Detection (metric) | Registered operations |
|---|---|---|
| `format/xyz` | suffix .xyz (10) | read_atoms |
| `format/mol` | .sdf/.mol/.sd (11) | read_atoms (raw block kept in `data["sdf"]`) |
| `format/mae` / `maegz` | .mae (15) / .maegz (20) | read_atoms, read_energy (`r_mmod_Potential_Energy*`, kJ/mol) |
| `app/g16/input` | .gjf/.com/.inp + route line (30) | read_atoms (no Z-matrix support), run, submit |
| `app/orca/input` | .com/.inp + leading `!` (40) | run, submit |
| `app/g16/output` | .log/.out/... + "Entering Gaussian System" (50) | read_atoms (archive / standard / input orientation), read_energy (SCF, optional method filter), read_thermal, check_end, check_freq, calc_free_energy |
| `app/xtb/output` | .log/.out/.xtb + banner (50) | read_energy, read_thermal, check_end |
| `app/pbs/jobscript` | .sh/.qsh/.qsub + `#PBS` (50) | submit, run (submit + wait) |
| `app/orca/output` | .log/.out + banner (60) | read_atoms (sibling .xyz), read_energy, check_end |

Gaussian-only (reachable through `.plugin.gau`): `read_vibration`, `check_bonding`, `read_nmr`,
`read_ecd`, `read_coupling`, `check_resub`/`get_resub`, `get_irc`, `get_trj`, `is_input`/`is_output`.
ORCA: `check_optimized`. xTB: `read_free_energy`. Maestro: `read_maegz`, `get_unzip`, `write_mae`
(via the external `sdconvert`).

## 7. Main algorithms and where they should live

| Algorithm | Current location | Nature | Proposed v2 placement |
|---|---|---|---|
| Relative energy | `BoxCore.calc_rel_energy` | ensemble-dependent, **overwrites `energy`** | standalone algorithm (pure function) + Flow annotation keeping absolute values |
| Energy filtering | `BoxCore.energy_limit` | ensemble-dependent selection (drops dE >= threshold, `max_limit`) | standalone algorithm -> Flow selection operation |
| Boltzmann populations | `BoxCore.calc_distribution` | ensemble-dependent | standalone algorithm -> Flow analysis operation |
| Minimum selection | `BoxCore.only_minimum` | ensemble-dependent selection (ties: first wins) | standalone algorithm -> Flow selection operation |
| Boltzmann averaging | `BoxCore.get_average` | aggregation (creates new Systems) | standalone algorithm -> Flow reduction operation |
| RMSD pruning | `topology.rmsdpruning` | ensemble-dependent, **order-dependent, approximate search** (`redundant_check`) | algorithm returning decisions -> Flow selection |
| Symmetry matrices | `Modeler.get_symmetry_matrices` | derived from a single structure | standalone algorithm (memoisable derived property of a Structure) |
| Bond perception, aromaticity | `Modeler.calc_bonds/aromatize` | derived from a single structure (threshold-dependent) | standalone algorithm (result becomes a perception attribute) |
| Stereo assignment | `Modeler.calc_stereo` | single structure (logs "under development") | standalone algorithm |
| Atom-number mapping | `topology.map_numbers`, `Modeler.get_maps` | reference-dependent, **reorders atoms in place** | standalone algorithm returning a new Structure |
| Geometry measurement | `Atoms.get_length/angle/dihedral` | single structure | read-only Structure methods |
| Geometry manipulation | `Modeler.set_length/mirroring/merge/incorporate` | single structure | transforms returning a new Structure |
| State management | `System.deactivate`, `Systems.has_*` | workflow state | Flow (selections as sets) / Box compatibility layer |

## 8. Technical debt and confirmed issues (not fixed in this change)

Items marked with a star are pinned by `quirk`-marked tests; the rest were confirmed by reading the code.

- (*) `System.duplicate()` does not copy `filetype`, `history` or `distribution` (`systems.py:168`);
  `Box.duplicate()` inherits this.
- (*) `Box(box)` and `.plugin.*` alias the contents.
- (*) `search()` with no arguments raises `UnboundLocalError` (`boxcore.py:288`).
- (*) `energy_limit` and friends raise `TypeError` when an active System has no energy.
- (*) `#DATA[key]#` raises `TypeError` for non-str values; `#ENRGY#` renders `"None"` when unset.
- (*) `TxtBox.parse_keys` joins lines together (`txtlib.py:77`).
- Centring in `write_xyz` / `convert_to_mirror` rounds the centroid to the number of decimals in the
  input coordinates, so it is only approximate (pinned by tests).
- `run`/`submit` do not detect failures (no `check=True`); `submit` keeps no process handle.
- `Atom.show()` prints y in place of z (`atoms.py:174`).
- `is_mae_format`/`is_maegz_format` use `p.suffix not in (".mae")`, a substring test on a string rather
  than a tuple (`maelib.py:139,146`).
- `read_xyz` reports the failure reason as "orca from xyz file" (`formats.py:18`, copy-paste).
- `TxtBox.delete_lines_by_key` **keeps** the lines containing the keyword (the opposite of its name);
  the slicing in `insert_lines` needs review (`txtlib.py:52`).
- `check_freq` only looks at the first three frequencies of each block and repeats its verdict per block.
- `calc_symm(calc_all=False)` (the default path of `rmsd_limit`) copies the symmetry matrices of the
  label representative to every conformer of that label, **implicitly assuming identical atom ordering**
  (`boxcore.py:387-395`).
- `rmsdpruning` is an approximate search with an early break (`redundant_check=3`) whose result depends
  on the energy ordering. `include_numeric_isomers` cannot be passed from Box, and the default of
  `all_perturbation` differs between Box (False) and topology (True).
- The kind of `energy` (SCF / G / relative) and its unit are not recorded, and `calc_rel_energy`
  destroys the absolute values.
- `multiplicity` defaults to 1 and `charge` falls back to atomic formal charges even when the source
  provides no such information.
- Unused or stub code: `modeler._cal_sym_rmsd` (references `Atoms.data`, which does not exist, so it
  cannot run), `Modeler.incorporate`, `Modeler.get_mapped`, `_superimpose`, `_map_to`,
  `topology.set_chirality`, `gaulib.read_multiple_xyz`, `Atoms.mw` (logs that it returns wrong values).
- `Log` installs a DEBUG-level stderr handler at import time and writes a `.mcl` file into the last
  input directory at interpreter exit (a side effect on the user's data directories). `Log.input_dir`
  is global state rewritten by every `System()` creation.
- Global registries (`FileType`, `Selectors`) rely on import side effects; plugin import order can
  affect detection order.
- `Elements` is a pickle (`elements.pkl`): not reviewable in diffs, and the provenance of the radii is
  only stated in a comment.
- README drifts from the implementation (`write_xyz(change_path=...)` is really `link`; `write_input`'s
  `arg` is undocumented).
- `.gitignore` ignored `tests/*` (now with an exception for `tests/baseline/`).
- The pinned versions in `.pre-commit-config.yaml` do not run on current environments (verified on
  Python 3.14): black 21.11b1 fails with `ImportError: cannot import name '_unicodefun'` on click >= 8.1,
  and pyupgrade v2.29.1 fails with `AttributeError` because `ast.Str` was removed. check-yaml,
  end-of-file-fixer, trailing-whitespace and mixed-line-ending still work. Updating the versions would
  produce formatting diffs, so nothing was changed here.

---

## 9. What to carry over from Box into v2.0

### 9.1 What actually makes Box pleasant to write

1. **A single entry point**: `Box("*.log")` accepts globs, directories, path lists, Systems and Boxes alike.
2. **Verb + domain vocabulary method names**: `read_*` / `check_*` / `calc_*` / `write_*` / `*_limit` /
   `only_minimum`. No workflow jargon (task, node, connect, graph) appears anywhere.
3. **Every method returns self**, so the written order is the execution order and the script reads top to bottom.
4. **Deactivation instead of removal**: excluded structures stay with their reason and remain auditable
   through `show()` / `export_data()`.
5. **Implicit grouping by label** (`in_label=True`): per-molecule conformer sets without writing loops.
6. **Automatic filetype dispatch**: `read_energy()` works across Gaussian, ORCA, xTB and Maestro.
7. **Scientifically sensible defaults**: 3 kcal/mol, 0.01 A, 298.15 K.
8. **Prerequisites run themselves**: `rmsd_limit()` obtains bonds and symmetry on its own.
9. **Templates as the calculation definition**: the calculation lives in the user's template file and
   ACCeL only fills in the `#KEY#` placeholders.
10. **Side outputs fit inside the chain**: `write_xyz()`, `export_data()` and `count()` do not break the flow.
11. **Every decision is logged**: exclusion reasons, minimum energies and so on.

Flow should preserve points 1-10 and turn point 11 into structured provenance.

### 9.2 APIs that can move to Flow unchanged (same name, same meaning)

Selection, analysis and measurement that can be expressed purely as "input set -> output set / annotation":

`energy_limit(threshold, max_limit, in_label)`, `only_minimum(in_label)`, `rmsd_limit(...)`,
`calc_distribution(in_label, temperature)`, `calc_energy(keys, unit)`, `calc_free_energy()`,
`labeling(separator, index_list)`, `set_label`, `set_data`, `calc_length/angle/dihedral`,
`calc_bonds`, `calc_symm`, `count`, `show`, `export_data`, and `read_atoms`, `read_energy`,
`read_thermal`, `check_end`, `check_freq` as analyses of existing outputs.

### 9.3 APIs that keep their name but change meaning

| API | Meaning in Box | Proposed meaning in Flow |
|---|---|---|
| `write_input(template)` | writes a file and re-points `path` | **declares a Task** built from a Template (lazily); input generation is part of the Task |
| `read_*` | overwrites the System | parses an output artifact into a **new version of a Structure / result** |
| `calc_rel_energy` | overwrites `energy` | annotates relative values while keeping the absolute ones |
| `energy_limit` / `only_minimum` / `rmsd_limit` | sets `state` to False | returns the selected output set; exclusions are kept as provenance |
| `run` / `submit` | immediate execution / fire-and-forget | either unnecessary, or a declaration that a Task should be executed; execution happens from the plan at `start()` |
| `labeling` | rewrites labels | defines the grouping key |
| `convert_to_mirror` / `modify_length` / `map_numbers` | mutates coordinates or atom order in place | transforms producing a new Structure |
| `get()` | the active Systems | access to evaluated results |

### 9.4 APIs that should exist only in the Box compatibility layer

`bind`, direct access to `contents` / `data`, `set_state` (manual reactivation), `duplicate`
(meaningless for an immutable Flow), the aliasing behaviour of `Box(box)`, `search`, `copy_files`,
`zero_fill` (file-name tidying), the `.plugin.*` accessor, `PbsBox.wait`, `submit` (fire-and-forget),
the file renaming in `check_resub`/`get_resub`, the `TxtBox` line editing operations, and `Dialog`.

---

## 10. Structure

### 10.1 System features that can move into Structure

Element list and coordinates (`Atoms`), atomic formal charges, total charge (keeping an explicit value
distinct from a derived one), multiplicity (**with a way to express "unspecified"**), bonds (marking
whether they came from the input or from perception), geometry measurement (length/angle/dihedral),
and coordinate transforms (mirror / set_length / merge) expressed as functions returning a new Structure.

### 10.2 Box-internal data that should move into Structure (or an attached result record)

- `energy` -> several energies carrying kind, unit and provenance (which calculation, which parser);
  e.g. SCF, ZPE correction, Gibbs correction.
- The parts of `data` that are calculation results (`g16_scf`, thermal corrections, frequencies,
  vibrational modes, NMR, ECD, couplings) -> results with provenance.
- `path` -> a reference to input/output artifacts (content hash plus storage location, not a path).
- `name` -> a display name, decoupled from identity.
- `label` -> undecided: an attribute of the Structure, or metadata of the Flow grouping (see Q3).

### 10.3 What must not go into Structure

`state` / `history` (Flow selection state), `distribution` and relative energies (ensemble-dependent),
label-wise operations, file IO, template expansion, calculation execution, filetype detection and
dispatch, logging, algorithm scratch data (`atom.cache["isomeric_subs_list"]` and similar), symmetry
matrices (memoised derived data), and the `Modeler` algorithms themselves (functions taking a Structure).

---

## 11. Flow

### 11.1 Existing code reusable as Flow operations

- The cores of `energy_limit`, `only_minimum`, `calc_distribution`, `calc_rel_energy`, `calc_energy`
  and `get_average` (per-label minimum, threshold test, Boltzmann weights) can be extracted as pure
  functions taking `(id, energy, group)` and returning the kept set or the values.
- `topology.rmsdpruning` -> a function returning which structures to exclude (today it calls
  `deactivate` directly).
- `Modeler.calc_bonds/aromatize/get_symmetry_matrices/calc_stereo` and geometry measurement/manipulation.
- `formats.replace_key/replace_arg` -> a pure function Structure + Template + params -> text.
- Each parser -> a function artifact -> result record (dropping the System writes and `deactivate`).
- The error classification in `check_end` / `check_freq` / `check_resub` -> validators.

### 11.2 Changes needed to use Templates as Flow tasks

1. **Separate rendering, writing and re-linking**: today `write_input` does all three at once.
2. **A Template object** (name undecided): source path, content bytes and their hash, the list of
   placeholders used (extractable statically from the text), and the `arg` parameters.
3. **Declared dependencies**: the placeholders a template uses determine which Structure attributes the
   Task depends on (a template using `#ENRGY#` must include the energy in the input fingerprint).
4. **A strict mode**: optionally treat unknown `#X#`, non-str `#DATA[]#` and `None` values as errors
   (without changing existing Box behaviour).
5. **Output binding**: which program and which parser read the result, and which output files are expected.
6. **Distinguish execution templates** (PBS job scripts) from calculation-input templates.
7. **Cache invalidation**: by default an exact byte match of the template content plus the parameters
   plus the fingerprint of the attributes used. Whether comment/whitespace-only edits should invalidate
   is an open question (Q6).
8. **Provenance**: record the template content (or its hash plus a stored copy), the parameters, the hash
   of the rendered input, and the program version readable from the output (e.g. `Version=` in a Gaussian
   archive).

### 11.3 Coupling that currently blocks incremental recomputation

1. **In-place mutation of System** (overwriting `energy`, reordering atoms, re-pointing `path`): there
   are no versions, so previous and current results cannot be compared.
2. **Selection results live as a flag on the System** (`state`), so a selection's output cannot be kept
   or compared as a set.
3. **`write_input` re-points `path`**: the same object turns from "output structure" into "input file",
   tying identity to a path.
4. **Aliasing** (`Box(box)`, views, plugins): the blast radius of a mutation cannot be traced.
5. **Parsing, validation and selection are fused**: parsers call `deactivate` directly.
6. **Executions are not tied to their results**: `submit` keeps no handle, `run` detects no failure, and
   reading the outputs is left to the user.
7. **Global state**: `Log.input_dir`, `Execmd`, the registries.
8. **Selections are not monotone**: the reference of `energy_limit` (the minimum) and the comparison
   order of `rmsd_limit` depend on the whole ensemble, so adding a new low-energy conformer D can change
   the verdict for the existing A/B/C (including turning a kept structure into an excluded one).

### 11.4 Information needed to track differences per Structure

- **Content fingerprint**: element list (including atom order), coordinates (rounding policy to be
  decided), total charge, multiplicity, user-specified bonds, and later periodic cells, fragment
  definitions and constraints.
- **Lineage id**: the content hash of the source artifact plus the index within the file (a single
  maegz / IRC / trajectory file holds many structures), the parent Structure, and the producing Task.
- **Task identity**: template content hash, parameters, the Structure attributes it depends on, the
  program name/version, and any environment settings that affect the result.
- **Result identity**: the content hash of the output artifact plus the parser kind/version.
- **Grouping key** (label) and the kind, unit and provenance of the energy used for selection.

### 11.5 Design questions for provenance / cache / invalidation

See Q1-Q10 in section 14; the most important ones are:

- Identity defined by content or by lineage (should identical coordinates from different sources be merged?).
- Numerical comparison of coordinates (decimal text from a file vs floating-point results of computation).
- Distinguishing "the parser changed, so re-parse" from "the calculation must be re-run".
- Caching of failed calculations (deterministic failures vs transient ones).
- Cache location and portability across machines (avoiding absolute paths).
- Non-monotone selections: results for structures that leave the selection must be retained rather than
  deleted, and merely marked as "not part of the current output".

### 11.6 Internal boundaries that let a parameter change propagate as a difference

```text
User-facing Flow chain (Box vocabulary)
    |  each method returns a new immutable Flow definition with one generic operation appended
Operation definitions
    - set-level operations (selection / grouping / analysis / reduction): cheap, may be recomputed in full
    - per-structure maps (Template Tasks, parsing, single-structure derived data): expensive,
      memoised by (task key, structure key)
    v
Dependency graph (the operation sequence with its inputs and outputs)
    v
Incremental evaluator
    - compares each operation's output fingerprint (id set + annotation values) and stops when unchanged (early cutoff)
    - for map operations, the new work is "current input ids minus cached keys"
    - because selections depend on upstream results, evaluation is staged (dynamic) rather than planned once
    v
Execution plan (only the Tasks that genuinely need to run)
    v
Executor interface (submit -> future -> artifact) - Local / Parsl (Local, Slurm, PBS)
    v
Result store (content-addressed artifacts plus a task-key index; nothing is deleted)
```

The key split is: **recompute the cheap set-level operations in full every time, and reuse the expensive
per-structure work at structure granularity.** Then widening `energy_limit` from 3 to 5 kcal/mol so that
A,B,C becomes A,B,C,D,E schedules only D and E, and if `only_minimum` still selects A the output
fingerprint is unchanged, so A's expensive downstream calculation is not re-run. Users never see this
distinction: they write the same chain they write with Box.

---

## 12. Execution / Parser

### 12.1 What can be delegated to Parsl

Process launching (`run`/`submit`), submission to PBS and polling with `qstat` (handled by Parsl's
providers/executors), Slurm support (new), parallelism and block/node management, retries for transient
failures, and waiting for completion (futures).

### 12.2 Workflow logic that must stay in ACCeL

Deciding what to run (the planner), input rendering, choosing the program and parser, output validation
(`check_end`, `check_freq`), **scientific error classification and resubmission policy** (`check_resub`
regenerates a modified geometry - that is domain knowledge), cache keys and provenance, selection
semantics, grouping, working-directory and artifact naming, and executable resolution (the `Execmd` role).

### 12.3 Parser work that can move to cclib (parity must be verified)

Gaussian: coordinates (the final geometry of `atomcoords`, corresponding to the archive/orientation
blocks), atomic numbers, charge and multiplicity, SCF energy, normal-termination detection, frequencies
(`vibfreqs`), vibrational modes (`vibdisps`), thermochemistry, excited states (`etenergies`, `etoscs`,
`etrotats`), and optimisation trajectories (the `get_trj` equivalent). ORCA: final energy, coordinates
(removing the dependency on a sibling `.xyz`), normal termination.

Caveats: cclib's **unit conventions** (energy units have differed between versions) and the fact that its
thermochemistry is provided as totals (electronic energy plus thermal correction) rather than as the
corrections ACCeL stores in `g16_corr_to_gibbs`. **Parity tests comparing the existing parsers and cclib
on the same output files are a precondition for migrating.** cclib support for xTB, NMR tensors and
spin-spin couplings has not been verified.

### 12.4 Parser work that should stay ACCeL-specific

`check_resub` (error classification plus input regeneration), `check_bonding` (the vibration-rank
heuristic), SCF extraction with a method filter (`SCF Done:  E(method)`), job counting from the input
file (`--link1--`, opt+freq), Maestro `.mae`/`.maegz`, sdf/mol (including keeping the raw block),
reading Gaussian `.gjf` inputs, xyz, filetype detection (input files and job scripts are outside cclib's
scope), IRC direction and reaction coordinate, and - where cclib does not provide them - the coupling
component selection and the xTB thermochemistry breakdown.

### 12.5 Where ASE looks useful

File-format interoperability as a fallback (extxyz, cif, POSCAR, pdb, ...), especially **crystals and
periodic systems** (the lattice and periodic boundaries that a generalised Structure will need),
conversion to and from `ase.Atoms` for external tooling, and neighbour lists for large systems.
Conversion rules for units, atom order and periodicity must be covered by tests.

---

## 13. Proposed migration roadmap

The two things to design first are **(a) the identity/fingerprint specification for Structure** and
**(b) the operation contract** (separating cheap set-level operations from expensive per-structure maps,
with early cutoff on output fingerprints). They are what makes "a Box-like chain" and "incremental
recomputation" compatible, and both the `Structure` class design and the Flow method design follow from them.

| Phase | Content | Done when |
|---|---|---|
| 0 | Baseline tests / CI / CLAUDE.md / this audit (**this change**) | existing behaviour is pinned in CI |
| 1 | Design document for identity, fingerprints and the operation contract; decisions on section 14 | the decisions are recorded (ADR); no implementation |
| 2 | `Structure` (immutable value object) and `System` <-> `Structure` adapters | Box unchanged; round-trip tests |
| 3 | Turn the algorithms into pure functions (energy / Boltzmann / only_minimum / RMSD verdicts); Box just calls them | the baseline suite passes unchanged |
| 4 | Parser layer: split out "artifact -> result record" functions, add a cclib backend and parity tests | values agree on the same outputs (differences documented) |
| 5 | Template object (pure rendering, placeholder extraction, fingerprint); `write_input` uses it internally | byte-identical output |
| 6 | Flow definition model (immutable chain, Box vocabulary) and an in-memory evaluator for **analysis-only flows** | same results as Box from existing outputs alone |
| 7 | Incremental engine: result store, per-structure memoisation, early cutoff, staged evaluation, tested with a **deterministic fake executor** | the 3->5 kcal/mol and minimum-unchanged/changed scenarios are covered by automated tests |
| 8 | Executor interface and a local implementation, then the Parsl backend (Local/Slurm/PBS) | the same scenarios pass locally |
| 9 | Box compatibility integration (Box using Structure and the algorithms internally) and Flow UX polish | the baseline suite passes unchanged |

Main differences from the originally sketched order: the identity specification comes before implementing
Structure (Phase 1); the parser split comes before Template/Flow, because analysis-only Flows depend on
pure parser functions; and Parsl is connected only after the engine is complete behind a fake executor,
so that the engine's correctness is verified independently of any scheduler.

---

## 14. Open questions (deliberately not settled here)

- **Q1 Structure identity**: defined by content (elements, coordinates, charge, multiplicity), by lineage
  (source file plus index plus parent), or both? May structurally identical conformers from different
  sources be treated as one?
- **Q2 Coordinate comparison precision**: round coordinates for the fingerprint (e.g. 1e-6 A), or treat the
  decimal text from the file as authoritative? How should near-boundary cases be handled?
- **Q3 The role of `label`**: an attribute of Structure, or a Flow grouping definition? In v0.3 it is
  derived from the file name.
- **Q4 Atom order**: are two otherwise identical molecules with permuted atoms the same structure?
  (They differ as calculation inputs but are equivalent for comparison.)
- **Q5 Kind of energy**: how should a Flow specify and record which energy (SCF / G / relative) a selection
  such as `energy_limit` uses? v0.3 simply uses "whatever was written last".
- **Q6 Template change detection**: exact byte match, or ignore whitespace and comment-only differences?
- **Q7 Parser updates**: should there be a mechanism to re-parse without re-running calculations when the
  parser version changes?
- **Q8 Caching failures**: should deterministic failures (e.g. SCF non-convergence) be cached or retried?
  Should a `check_resub`-style automatic resubmission be the Flow default?
- **Q9 UX for parameter changes**: (a) edit the value in the script and re-run, relying on a persistent
  cache (closest to how Box is used today), or (b) named steps with a parameter-override API? This audit
  recommends (a) as the primary route with (b) as a supplement.
- **Q10 Explicit execution in Flow**: should `run()` be required after `write_input(template)`, or should
  `write_input` -> `read_*` implicitly establish the Task?
- **Q11 Existing quirks**: should the starred items in section 8 (the `duplicate` omissions, `search()`,
  and so on) be fixed in v0.x, or preserved in the v2 Box compatibility layer?
- **Q12 Defaults**: should "multiplicity defaults to 1" and "total charge derived from atomic formal
  charges" remain the defaults in Structure, or should "unspecified" become mandatory?
