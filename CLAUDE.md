# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

`potpourri` is a Python tool for **multi-period Optimal Power Flow (OPF)** in distribution grids. It wraps [Pyomo](https://pyomo.readthedocs.io/) for optimization modeling over [pandapower](https://pandapower.readthedocs.io/) network objects, supporting AC/DC formulations and flexible resources (batteries, heat pumps, PV, wind).

## Setup

```bash
conda env create -f environment.yaml   # creates `potpourri_env`
conda activate potpourri_env
pip install -e .
```

Solvers (IPOPT, GLPK, CBC, Gurobi) must be installed separately. The Dockerfile shows how to compile IPOPT 3.14.20.

## Development Commands

```bash
ruff check .              # lint
ruff format .             # format
pytest                    # run all tests
pytest -m "not integration"   # skip solver-dependent tests
```

## Architecture

### Data flow

```
pandapower Network
  → Basemodel.__init__()       # extracts buses, lines, loads into Pyomo sets/params
  → [AC|DC] power flow mixin   # adds power flow equations
  → OPF mixin                  # adds operational constraints + objectives
  → .solve(solver='ipopt')     # calls Pyomo SolverFactory
  → pyo_to_net()               # writes Pyomo vars back to net.res_* DataFrames
```

### Class hierarchy (single-period, `src/potpourri/models/`)

- `Basemodel` — creates the `ConcreteModel`, maps pandapower DataFrames to Pyomo sets/parameters, provides `solve()`.
- `AC` / `DC` — extend Basemodel with power-flow equations (complex vs. linearised).
- `OPF` — adds generator/load limits, line loading limits, and objective functions.
- `ACOPF_base` — multiple-inherits `AC + OPF` for a full AC OPF.
- `HC_ACOPF` — hosting-capacity variant.

### Multi-period models (`src/potpourri/models_multi_period/`)

`Basemodel_multi_period` adds a time dimension and simbench profile integration. Device modules are instantiated as **mix-in objects** that attach their own Pyomo constraints/variables to an existing multi-period model:

```python
mpopf = ACOPF_multi_period(net, toT=24)
# Devices take the net (not the model) and attach in a separate get_all() call
battery = Battery_multi_period(mpopf.net, T=24, scenario=1)
battery.get_all(mpopf.model)
mpopf.add_OPF()
mpopf.add_voltage_deviation_objective()
mpopf.solve(solver='ipopt')
```

Device modules (all suffixed `_multi_period`, in `src/potpourri/technologies/`):
`Battery`, `Heatpump` (lower-case `p`), `PV`, `Windpower`, `Demand`, `Sgens`,
`Shunts`, `Generator`, and the `Flexibility` base class. There is no EV module.

### Supporting modules

- `net_augmentation/prepare_net.py` — adds missing pandapower columns, scales profiles before model construction.
- `plotting/plot_functions.py` — visualises network state and results.
- `pyo_to_net[_multi_period].py` — post-processing: reads Pyomo solution and writes to `net.res_*`.
- `init_pyo_from_pp_res[_multi_period].py` — warm-starts Pyomo variables from a prior pandapower power-flow result.

## Key conventions

- **Deep copy** the pandapower network before passing it to a model to avoid mutation.
- Pyomo components (Sets, Params, Vars, Constraints) are added to `self.model` inside each class.
- `solve(to_net=True)` (the default) populates `net.res_*` itself. Multi-period
  writes **one** time step — the last of the horizon — because `net.res_*` has
  no time dimension; use `map_to_net(t)` for any other step.
- Example scripts in `scripts/` are the primary usage examples (see `scripts/README.md`).
- **Pyomo style** follows the
  [MO-book style guide](https://mobook.github.io/MO-book/notebooks/appendix/pyomo-style-guide-update.html):
  `import pyomo.environ as pyo` (never a star import — ruff's F403/F405
  enforce it), component **decorators** rather than `rule=`, `domain=`
  rather than `within=` on a `Var`, ALL_CAPS set names. The decorated
  function's name *becomes* the component name, and those names are
  public API — tests and scripts index into `model.line_lim_from`, so
  renaming one is a breaking change, not a style edit. Decorate on the
  expression that already holds the model (`@self.model.Constraint(...)`,
  or `@model....` where the method takes it as a parameter); do not bind
  `model = self.model` at the top of a method, because `create_model()`
  overrides replace `self.model` partway through. Maths-derived names
  (`qsG`, `pTlv`, `SLmax`, `l`) stay as they are. Full policy in
  `docs/contributing-pyomo.md`.
- **Docstrings are Google-style Markdown**, rendered by mkdocstrings —
  `$v_b$` for maths, `` [`X`][potpourri.a.b.X] `` for cross-references,
  never reStructuredText roles. Document units, sign conventions, ppc-vs-
  pandapower indexing, side effects on `self.model`/`self.net`, and what a
  Pyomo rule returns (an expression or a `(lo, expr, hi)` tuple, never a
  bool). Enforced by `interrogate src/potpourri` and `ruff check`; the full
  policy is in `docs/contributing-docs.md`.
- **Every Python file needs an SPDX header** — an `SPDX-FileCopyrightText`
  line plus an `SPDX-License-Identifier` line naming MIT, as real comments at
  the very top, new files included. Enforced in pre-commit and CI by
  `python tools/check_license_headers.py`. Never relabel copied third-party
  code as MIT or invent a copyright holder to make the check pass; the policy
  and the open ownership questions live in `docs/licensing.md`.
