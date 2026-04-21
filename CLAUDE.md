# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

`potpourri` is a Python tool for **multi-period Optimal Power Flow (OPF)** in distribution grids. It wraps [Pyomo](https://pyomo.readthedocs.io/) for optimization modeling over [pandapower](https://pandapower.readthedocs.io/) network objects, supporting AC/DC formulations and flexible resources (batteries, EVs, heat pumps, PV, wind).

## Setup

```bash
conda env create -f environment.yaml   # creates `potpourri_env`
conda activate potpourri_env
pip install -e .
```

Solvers (IPOPT, GLPK, CBC, Gurobi) must be installed separately. The Dockerfile shows how to compile IPOPT 3.14.16.

## Development Commands

```bash
black .       # format
flake8 .      # lint
pytest        # run tests (no formal test suite yet; see scripts/ for examples)
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
battery = Battery_multi_period(mpopf)   # attaches battery constraints to mpopf
mpopf.solve(solver='ipopt')
```

Device modules: `Battery`, `EVs`, `HeatPump`, `PV`, `Windpower`, `Demand`, `Sgens`, `Flexibility`.

### Supporting modules

- `net_augmentation/prepare_net.py` — adds missing pandapower columns, scales profiles before model construction.
- `plotting/plot_functions.py` — visualises network state and results.
- `pyo_to_net[_multi_period].py` — post-processing: reads Pyomo solution and writes to `net.res_*`.
- `init_pyo_from_pp_res[_multi_period].py` — warm-starts Pyomo variables from a prior pandapower power-flow result.

## Key conventions

- **Deep copy** the pandapower network before passing it to a model to avoid mutation.
- Pyomo components (Sets, Params, Vars, Constraints) are added to `self.model` inside each class.
- `pyo_to_net` must be called after `solve()` to populate `net.res_*` DataFrames.
- Tutorial scripts in `tutorials/` and analysis scripts in `scripts/` are the primary usage examples.
