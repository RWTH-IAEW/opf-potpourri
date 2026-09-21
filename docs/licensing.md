# Licensing

potpourri is distributed under the **MIT License**. The canonical text
is [`LICENSE`](https://github.com/RWTH-IAEW/opf-potpourri/blob/main/LICENSE)
at the repository root; `LICENSES/MIT.txt` is a byte-identical copy in
the layout the [REUSE specification](https://reuse.software/) expects,
and a unit test keeps the two from drifting apart.

Copyright is held by the Institute for High Voltage Equipment and
Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen
University.

## The header every Python file carries

```python
# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT
```

- **It must be a comment.** The same text inside a docstring or a
  string literal does not count, and the checker rejects it.
- **Copy the holder wording verbatim** from `LICENSE`; do not shorten
  it to "IAEW" or "RWTH".
- **One line, however long.** Wrapping truncates the holder in
  generated SPDX documents. Ruff exempts `SPDX-` pragmas from `E501`.
- **Years are `2023-2026`**, the span of the repository's own history,
  the same in `LICENSE` and in every file. Extend the end year in
  `LICENSE`, `LICENSES/MIT.txt` and `FIRST_PARTY_COPYRIGHT` together;
  do not stamp files individually.

### Authorship is not copyright

Copyright here is **institutional**, so no file carries a personal
copyright line. Attribution is welcome, and is written as

```
Author: Steffen Kortmann (2023)
```

`Author:` claims no rights. A `(c) YEAR, Name` line in prose *is* a
copyright claim, and the checker reports it wherever it appears — a
copyright notice belongs in the SPDX header.

An **outside contributor keeps their own copyright**, so a file may
legitimately carry a second `SPDX-FileCopyrightText` line. The checker
accepts that; the fixer refuses to touch such a file.

## Contributing code

Contributions are accepted **inbound = outbound**: opening a pull
request offers your contribution under the MIT License, and you keep
your own copyright in what you wrote. No assignment to RWTH is asked
for, and none is implied.

### Code you did not write

Importing a package is not copying it — `import pandapower` creates no
obligation beyond that package's own licence. Pasting or porting
*source code* does. When you bring some in:

1. **Keep every upstream notice**: copyright lines, licence blocks,
   disclaimers, `@author` tags. Adding an SPDX identifier never
   justifies deleting a licence block.
2. **Do not relabel it MIT.** Record the real SPDX expression.
3. **Add the licence text** to `LICENSES/<SPDX-ID>.txt`, and to
   `license-files` in `pyproject.toml` if it ships in the distribution.
4. **Record the provenance**: upstream project, path, version or
   commit, and the licence *at that version* — today's default branch
   is not evidence for what a file was licensed under when it was
   copied.
5. **Register it** in `LICENSE_EXCEPTIONS` in
   `tools/license_headers.py`, with its evidence. The checker then
   requires that expression instead of MIT, and the fixer refuses to
   touch the file.

For a part-first-party, part-upstream file, work the SPDX expression
out from the actual terms. Do not write `MIT AND BSD-3-Clause`
mechanically, and do not use `OR` — that implies a choice nobody
granted.

## Running the check

```bash
python tools/check_license_headers.py   # read-only gate
python tools/fix_license_headers.py     # preview; APPLY = True to write
```

It also runs in pre-commit and in CI, which additionally cross-checks
with the [`reuse`](https://reuse.software/) reference tool. **CI never
repairs a file.** Exit codes: `0` compliant, `1` violations, `2` a
discovery or tool error — a broken scan is never reported as a pass.

The check verifies **placement and consistency**: that a notice sits in
the physical header, that the SPDX expression parses, and that nothing
is missing, duplicated or contradictory. It cannot verify *ownership*.
A green check is not a legal clearance.

## If provenance or ownership is unclear

**Leave the declaration alone and escalate.** Do not guess a holder,
invent a year, write "and contributors" to paper over a gap, or pick a
licence because it makes the build pass. An author list or a `git`
commit author is not, by itself, proof of copyright ownership.

## Audit record

The Python sources were audited in September 2026. Every tracked file
was inspected: all are first-party, and **no third-party or
mixed-origin Python source was found**, so `LICENSE_EXCEPTIONS` is
empty. References to MATPOWER, PYPOWER and pandapower were each traced
and turned out to be prose describing a convention, imports of an
installed dependency, or two `ppc` column offsets — not copied code.
No file carries an upstream copyright line, licence block or `@author`
tag.

The `benchmarks/pglib-opf` submodule (v23.07) is dual-licensed —
CC-BY-4.0 for the case data, MIT for the software, © 2017 *A Library of
IEEE PES Power Grid Benchmarks* — and contains no Python. It was
inspected read-only. Nothing third-party is bundled in the
distribution.

Per-file results are in `docs/licensing-inventory.csv`, regenerated by
`python tools/licensing_inventory.py`.

Two questions the audit raised were settled by the maintainer and are
reflected above: copyright is institutional, so 15 personal claims were
retired to `Author:` attribution; and `LICENSE` carries the 2023-2026
span rather than a single year.
