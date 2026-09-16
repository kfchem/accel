# CLAUDE.md

Guidance for AI agents and contributors working on ACCeL.

## What ACCeL is

ACCeL is a research library for computational chemistry: it manages molecular
structures and conformer ensembles, generates quantum-chemistry inputs from
templates, and parses/analyses calculation results (energies, filtering, RMSD,
Boltzmann populations, ...). Users are researchers whose results end up in
papers. **Scientific reproducibility comes before everything else.**

## Commands

```bash
python -m pip install -e .                 # numpy is the only runtime dependency
python -m pip install -r requirements-dev.txt
python -m pytest                           # baseline suite in tests/baseline (see pytest.ini)
ruff check .                               # minimal static checks (ruff.toml)
```

CI (`.github/workflows/ci.yml`) runs ruff, byte-compilation, an import smoke test
of the installed package, and the baseline tests on Python 3.9-3.14 (Linux) and
Windows. External programs (Gaussian, ORCA, xTB, qsub) are never required in tests.

## Layout (v0.3)

- `accel/base/systems.py` - `System` (one structure + its state/metadata) and `Systems` (list with filter views)
- `accel/base/atoms.py` - `Atom`, `Atoms`, `Bonds` (NxN bond-type matrix)
- `accel/base/boxcore.py` - `BoxCore`: the chainable collection API (energy, filtering, IO, structure ops)
- `accel/base/box.py` - `Box`: adds filetype-dispatched methods (`read_atoms`, `read_energy`, `run`, ...) and `.plugin`
- `accel/base/selector.py`, `accel/util/filetype.py` - import-time registries: filetype detection and filetype -> parser dispatch
- `accel/base/formats.py` - xyz/mol IO and the template engine (`write_input`, `#NAME#` etc.)
- `accel/base/topology.py`, `accel/base/modeler.py` - RMSD pruning, atom mapping, bonds, symmetry, stereo
- `accel/plugin/*lib.py` - program plugins (Gaussian, ORCA, xTB, Maestro, PBS, text)
- `tests/baseline/` - characterisation tests that pin current behaviour

## Design principles (v2.0 and beyond)

1. **Reproducibility first.** Never change a numerical result, threshold, default,
   unit, rounding rule or selection order silently. If existing behaviour looks
   scientifically wrong but its intent is unclear, do not "fix" it by guessing:
   pin it with a `quirk`-marked test, document it, and ask.
2. **`Structure` is coming in v2.0.** It represents a general collection of atoms:
   a single molecule, several molecules, complexes, transition states, ion pairs,
   clusters, crystals, ... It is *not* a rename of `System` and *not* a rename of `Box`.
3. **Do not break the existing `Box` API.** Box stays as the compatibility layer for
   existing users; no removal, no deprecation, no signature or behaviour changes
   without an explicit decision. The baseline suite must keep passing.
4. **Flow inherits Box's writing style.** Short, top-to-bottom method chains using
   Box vocabulary (`read_energy()`, `energy_limit(3.0)`, `rmsd_limit()`,
   `only_minimum()`, `write_input(template)`). Never require users to wire graphs
   (`task.connect(...)`, `graph.add_node(...)`).
5. **Flow's internals are not bound to Box's mutable implementation.** Reuse Box's
   *names and UX*, not its in-place mutation model.
6. **Flow definitions are immutable.** Each chained call returns a new Flow definition.
7. **Flow is lazily evaluated.** Chaining only records operations, parameters,
   templates and dependencies; heavy work starts on an explicit call.
8. **Incremental recomputation is the core requirement.** After a parameter change
   (e.g. `energy_limit` 3 -> 5 kcal/mol) only the genuinely new work is done.
9. **Caching works at the level of structures, not just nodes.** A filter going from
   [A, B, C] to [A, B, C, D, E] reuses A, B, C downstream and schedules only D, E.
10. **No downstream recomputation unless an output actually changed.** Re-evaluating
    an upstream node is not a reason to recompute its dependants (early cutoff):
    if `only_minimum` still picks A, A's expensive follow-up is reused.
11. **Reuse existing results as much as possible**, including results produced
    outside ACCeL (analysis-only use must never require running a calculation).
12. **Templates are the primary way to define calculations.** Keep and extend the
    `#KEYWORD#` template approach. Do not hard-code calculation kinds
    (opt/freq/SP/IRC/NMR/TD-DFT, ...) as separate ACCeL methods.
13. **Parsl is an execution backend**, used to run tasks ACCeL has already decided
    to run (Local / Slurm / PBS). Structure diffs, selections, provenance and
    cache validity stay in ACCeL.
14. **Parsers: cclib first**, with ACCeL-specific extensions only where cclib lacks
    the data; keep existing parsers until parity is demonstrated by tests.
15. **Separate calculation from analysis.** Input generation, execution, result
    parsing and scientific analysis must be usable independently.
16. **Structure identity / cache validity have scientific meaning.** Do not decide
    them implicitly (e.g. by path or filename alone); agree them with the user first.

## Working rules

- Work on a branch named `agent/...`; never commit, push or merge to `master` directly; never force-push.
- Keep diffs focused: no mass reformatting (the code base uses black-style, line
  length 119), no unrelated refactors, no drive-by fixes of scientific code.
- Tests must not need Gaussian/ORCA/xTB/schedulers. Use minimal synthetic excerpts
  (see `tests/baseline/builders.py`) and do not assert scientific values that cannot
  be derived independently.
- `accel.util.log` writes a `.mcl` log file at interpreter exit into the last input
  directory; keep test data in `tmp_path`.

## Communication and GitHub workflow

Permanent rules; they apply to every session.

### Language

- **User communication in Claude is in Japanese.** Design discussion, questions,
  confirmations, presented options, progress updates, explanations of implemented
  work and final reports are written in Japanese in the Claude conversation.
- **All GitHub-facing content must be in English**: branch names, commit messages,
  Pull Request titles, descriptions and comments, Issue titles, descriptions and
  comments, release notes, and any other text published on GitHub. Never write
  Japanese on GitHub.
- **Documentation committed to the repository is English too** (README, `CLAUDE.md`,
  everything under `docs/`, code comments), since it is published on GitHub.

### Where discussion happens

- **Design discussions and clarification requests happen in the Claude conversation,
  not in Pull Requests or Issues.** Never open a PR or Issue just to ask a question,
  and never use a PR body as a design-discussion log.
- **Only finished work is published.** Audits, design notes, interim reports and
  anything that is not part of the deliverable stay local in `.local/` (gitignored);
  do not commit them. Keep everything that does reach GitHub minimal.

### Decisions

- **Routine implementation decisions are made autonomously**, based on the existing
  code, tests, documentation and the v2.0 design principles above.
- **Discuss with the user before finalising important decisions**, in Japanese in the
  Claude conversation: anything affecting scientific meaning or results, existing Box
  compatibility, Flow UX, Structure identity, cache / invalidation / incremental
  recomputation semantics, Template semantics, the public API, or the overall v2.0
  architecture.

### Pull Requests

- **A Pull Request is a review-ready implementation unit** for final human review,
  not a place for design discussion or interim confirmation.
- One task = one `agent/` working branch = one final Pull Request.
- A finished PR description is concise English and covers only **Purpose** and
  **Summary of changes**.
- Never merge your own Pull Request. Never close or delete existing Pull Requests or
  branches without an explicit instruction from the user.
