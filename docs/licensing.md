# Licensing and copyright headers

potpourri is distributed under the **MIT License**. The canonical text is
[`LICENSE`](https://github.com/RWTH-IAEW/opf-potpourri/blob/main/LICENSE) at
the repository root; `LICENSES/MIT.txt` is a byte-identical copy in the layout
the [REUSE specification](https://reuse.software/) expects, and a unit test
keeps the two from drifting apart.

This page is the contributor-facing policy. It is also the audit record: what
was checked in September 2026, what changed, and what is still open.

---

## The first-party header

Every Python file in the repository starts with this, as **real comments**
before anything else:

```python
# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT
```

Points that matter:

- **It must be a comment.** The same text inside a module docstring, a string
  literal or a documentation example does **not** count, and the checker
  rejects it. Tooling reads the physical header, so that is where the notice
  has to be.
- **The holder wording is copied verbatim from `LICENSE`.** Do not shorten it
  to "IAEW" or "RWTH" — the long form is the rights holder as declared.
- **One line, however long.** Wrapping the notice across two comment lines
  truncates it in generated SPDX documents, so it stays on one line. Ruff
  exempts `SPDX-` pragma comments from `E501`, so this does not fight the
  79-character limit.
- Order: copyright line(s) first, then a bare `#`, then the licence
  identifier. A blank line separates the header from the module docstring.

### Copyright years

The year range is **`2023-2026`** on every first-party file, and `LICENSE`
carries the same span: 2023 is the repository's earliest revision, 2026 the
most recent. It is a single project-wide range rather than a per-file one,
because a 2026 reorganisation rewrote the path history of most files —
`git log --follow` reports 2026 for files whose own notices say 2023, so
per-file git dates would *understate* the real authorship.

Extend the end year when the project moves into a new year, in `LICENSE`,
`LICENSES/MIT.txt` and `FIRST_PARTY_COPYRIGHT` together. Do not stamp files
individually.

### Authorship is not copyright

Copyright in this repository is **institutional**. Files do not carry a
personal copyright line, and neither should new ones.

Several example scripts used to end their docstring with a personal claim:

```
(c) 2023, Steffen Kortmann
```

The maintainer confirmed in September 2026 that copyright is institutional,
so those claims were retired in favour of the holder named in `LICENSE`, and
the name kept as plain attribution:

```
Author: Steffen Kortmann (2023)
```

`Author:` is attribution and carries no rights claim, so it is welcome
anywhere. A `(c) YEAR, Name` line in prose *is* a copyright claim, and the
checker reports it wherever it appears — a copyright notice belongs in the
SPDX header, not in text that `help()` prints.

An **outside** contributor keeps their own copyright; see
[Contributions from outside IAEW](#contributions-from-outside-iaew). The
checker therefore tolerates extra `SPDX-FileCopyrightText` lines, while the
fixer refuses to touch a file that has one — whose claim it is, is a
question for a person.

---

## Copied or adapted third-party code

There is currently **no third-party Python source in this repository**, so the
rules below are for the next time someone adds some.

Importing a package is not copying it. `import pandapower` or calling
`from_mpc()` creates no obligation beyond the dependency's own licence.
Pasting or porting *source code* does.

When you bring in copied or adapted code:

1. **Keep every upstream notice.** Copyright lines, licence blocks,
   disclaimers, `@author` tags and attribution requirements stay in the file.
   Adding an SPDX identifier never justifies deleting a full licence block.
2. **Do not relabel it MIT.** The upstream licence governs that file. Record
   the real SPDX expression, even where it complicates the build.
3. **Add the licence text** to `LICENSES/<SPDX-ID>.txt` (for example
   `LICENSES/BSD-3-Clause.txt`) so the text ships with the project, and add the
   file to `license-files` in `pyproject.toml` if it is bundled in the
   distribution.
4. **Record the provenance**: upstream project, source path, version or commit,
   licence at *that* version, and the required notices. Today's default branch
   is not evidence for what a file was licensed under when it was copied, and
   not every file in a project uses that project's top-level licence.
5. **Register it** in `LICENSE_EXCEPTIONS` in `tools/license_headers.py`, with
   the SPDX expression and a one-line evidence note. The checker then requires
   *that* expression for the file instead of MIT, and the fixer refuses to
   touch it.

For a file that is genuinely part first-party and part upstream, keep the
attribution for the upstream portion and work out the correct SPDX expression
from the actual terms. Do not write `MIT AND BSD-3-Clause` mechanically, and do
not use `OR` — that implies a choice nobody granted.

---

## Running the checks

```bash
# Read-only gate. Covers every tracked and newly added Python file.
python tools/check_license_headers.py

# Preview the header the policy would write (writes nothing).
python tools/fix_license_headers.py

# Apply it: set APPLY = True in that file, then re-run.
```

The gate also runs:

- on every commit, via the `license-headers` pre-commit hook;
- in CI, in the `licensing` job of `.github/workflows/ci.yml`, which
  additionally cross-checks the files with the
  [`reuse`](https://reuse.software/) reference implementation.

Both use the same code. **CI never repairs a file** — it only reports. The
fixer is the only thing that writes, it previews by default, and it refuses any
file whose header declares a licence or a holder the policy cannot verify.

Exit codes: `0` compliant, `1` policy violations, `2` discovery or tool error.
A broken scan is an error, never a pass.

To regenerate the per-file audit record:

```bash
python tools/licensing_inventory.py   # -> docs/licensing-inventory.csv
```

### What the checker does and does not prove

It verifies **placement and consistency**: that a machine-readable notice is in
the Python header, that the SPDX expression parses (validated with
`packaging.licenses`, PEP 639), that it matches the recorded policy for that
file, and that nothing is missing, empty, duplicated or contradictory.

It cannot verify **ownership**. A syntactically valid `MIT` proves only that
someone typed a valid identifier. Provenance is decided by a human reading the
file and its history, and recorded here. Keep the two apart: a green check is
not a legal clearance.

---

## When provenance or ownership is unclear

**Leave the declaration alone and escalate.** Do not guess a holder, invent a
year, write "and contributors" to paper over a gap, or pick a licence because
it makes the build pass. An author list, a `git` commit author or an
institutional affiliation is not by itself proof of copyright ownership.

Open a pull request that changes nothing but adds the file and the open
question to this page, and ask the maintainers to route it to RWTH's licensing
contact. An unresolved question recorded in the open is the correct outcome;
a passing scan that hides it is not.

---

## Audit record — September 2026

Scope: commit `f222137` on `main`, all **88** tracked Python files plus the 5
added by this work (93 total). No `.pyi` or `.pyw` files exist. Covered
`src/potpourri/` (39), `scripts/` (22), `tests/` (27, including the three empty
`__init__.py` files), `.github/scripts/` (1) and `tools/` (4).

| Outcome | Files |
| --- | --- |
| First-party, header added | 93 |
| Third-party or mixed-origin | 0 |
| Already compliant before the audit | 0 |
| Outside editable scope | submodule `benchmarks/pglib-opf` (no Python) |
| Unresolved | 0 files blocked; see [Still open](#still-open) |

Provenance was checked rather than assumed. References to MATPOWER, PYPOWER and
pandapower were each traced:

- Mentions in docstrings and comments describe a **convention** (for example
  the `-1/x` DC susceptance, or `ANGMIN`/`ANGMAX` column semantics). Prose is
  not copied code.
- `from_mpc` and `makeYbus` appear as **imports** of the installed dependency.
- `src/potpourri/models/basemodel.py` defines `_BR_B = 4` and `_BR_G = 23`,
  two column offsets into pandapower's runtime `ppc` array. These are
  interface constants for reading a dependency's data structure.
- `.github/scripts/zenodo_release.py` is first-party and uses only the
  standard library.

No file in the repository carries an upstream copyright line, licence block,
`@author` tag or "derived from" note. The submodule `benchmarks/pglib-opf`
(pinned at `v23.07`, `dc6be4b`) is dual-licensed — CC-BY-4.0 for the case data,
MIT for the software, © 2017 *A Library of IEEE PES Power Grid Benchmarks* — and
contains **no Python files**. It was inspected read-only and left untouched, as
were its submodule reference and contents.

Nothing third-party is bundled in the distribution: the wheel and sdist contain
first-party Python modules, `LICENSE` and `LICENSES/MIT.txt`, and no data
files. Dependencies are installed separately and their licences are not
aggregated into potpourri's own metadata.

The build previously declared `setuptools>=69` while using PEP 639's
`license` and `license-files` fields, which setuptools only supports from
77.0.0 on; below that it silently emits legacy `License:` metadata instead of
`License-Expression:`. The floor was raised, and a real build was inspected to
confirm the wheel carries `License-Expression: MIT` and both licence files.

### Decisions taken, September 2026

Two of the three questions the audit raised were settled by the maintainer
on 2026-09-21 and are now implemented:

1. **Copyright is institutional.** The 15 personal `(c) YEAR, Steffen
   Kortmann` claims were retired in favour of the holder named in `LICENSE`.
   The names stay as `Author:` attribution, which asserts no rights. The
   checker now reports a personal copyright claim in prose rather than
   mirroring it.
2. **`LICENSE` carries the project's real span.** It read
   "Copyright (c) 2024" while the repository's revisions run 2023-2026;
   it and `LICENSES/MIT.txt` now read `2023-2026`, matching the headers.

These record a maintainer decision about this institute's own work. They are
not a legal opinion, and neither is the passing check.

### Still open

3. **Contributor rights.** `pyproject.toml` and `CITATION.cff` list seven
   authors, all IAEW. Under the decision above their work is institutional,
   so the listed authors are covered — but nothing in the repository records
   the **inbound** terms for contributions from outside IAEW, and the public
   mirror accepts pull requests. See below.

### Contributions from outside IAEW

Until the maintainers choose otherwise, contributions are accepted on the
ordinary open-source basis: **inbound = outbound**. Opening a pull request
means you offer your contribution under the MIT License in `LICENSE`, and
you keep your own copyright in what you wrote. No assignment to RWTH is
asked for, and none is implied.

A substantial outside contribution may therefore add its own
`SPDX-FileCopyrightText` line alongside the institutional one. The checker
accepts that; the fixer refuses to rewrite such a file. Do not remove
someone else's copyright line to make a file look uniform.

Three mechanisms exist if the maintainers want something stronger than the
statement above, in increasing order of weight:

| Mechanism | What it gives | Cost |
| --- | --- | --- |
| Inbound = outbound (current) | A clear licence for contributions | A sentence in `CONTRIBUTING.md` |
| DCO sign-off | Each contributor certifies they may submit the work; a per-commit audit trail | A `Signed-off-by` line, plus a CI check |
| CLA | Explicit grant, possibly assignment, to RWTH | Legal drafting, a signature process, friction for drive-by fixes |

For a research tool of this size, inbound = outbound is the norm and is what
JOSS reviewers expect to see stated. A DCO is a reasonable middle step if
provenance of outside patches becomes a concern. A CLA needs RWTH's legal
office and would need to come from them, not from this file.
