# potpourri

Multi-**P**eriod **O**p**t**imal **Po**wer Flow for Distrib**u**tion G**r**ids with Sto**r**age Appl**i**cation

> *Potpourri — piece of music composed from various popular smaller works or melodies*

## Overview

`potpourri` is a Python library for AC/DC Optimal Power Flow (OPF) in distribution grids, with support for multi-period planning and flexible resources (batteries, EVs, heat pumps, PV, wind). It wraps [Pyomo](https://pyomo.readthedocs.io/) for optimization modelling over [pandapower](https://pandapower.readthedocs.io/) network objects.

---

## Documentation

The project documentation is built with [MkDocs](https://www.mkdocs.org/).

To serve the documentation locally:

```bash
mkdocs serve
```

This starts a local development server, usually at:

```
(opf-potpourri) (base) kortmann@kortmann:~/opf-potpourri$ mkdocs serve
INFO    -  Building documentation...
INFO    -  Cleaning site directory
INFO    -  Documentation built in 0.17 seconds
INFO    -  [11:01:44] Serving on http://127.0.0.1:8000/
```

If MkDocs is not installed yet, install the documentation dependencies first:

```bash
pip install -e .[docs]
```

then run:

```bash
mkdocs serve
```

## Installation

Requires Python ≤ 3.12 and Conda.

```bash
conda env create -f environment.yaml
conda activate potpourri_env
pip install -e .
```

To update an existing environment:

```bash
conda env update -f environment.yaml
```

A Dockerfile is also provided for a fully containerised setup including compiled IPOPT 3.14.16, SHOT, and CBC solvers.

---

## Usage

### Single-period AC OPF

```python
import simbench as sb
from potpourri.models.ACOPF_base import ACOPF

net = sb.get_simbench_net("1-LV-rural1--0-sw")
opf = ACOPF(net)
opf.add_OPF()
opf.add_voltage_deviation_objective()
opf.solve(solver='ipopt', print_solver_output=False)
# results available in net.res_bus, net.res_line, etc.
```

### Multi-period AC OPF with flexible resources

```python
from potpourri.models_multi_period.ACOPF_multi_period import ACOPF_multi_period
from potpourri.models_multi_period.battery_multi_period import Battery_multi_period

opf = ACOPF_multi_period(net, toT=96, fromT=0, num_vehicles=5)
battery = Battery_multi_period(opf)   # attaches battery constraints
opf.add_OPF()
opf.solve(solver='ipopt')
```

See `tutorials/opf.py` and `tutorials/mpopf.py` for full examples.

---

## Architecture

### Single-period models (`src/potpourri/models/`)

```
Basemodel     creates ConcreteModel, maps pandapower → Pyomo sets/params, solve()
  ├── AC      full AC power flow (KCL/KVL, voltage magnitudes, reactive power)
  ├── DC      linearised DC power flow (no reactive power)
  └── OPF     operational constraints (P/Q limits, line loading, voltage bounds)

ACOPF    = AC + OPF   (multiple inheritance)
DCOPF    = DC + OPF
HC_ACOPF = ACOPF + binary variables for hosting-capacity analysis
```

### Multi-period models (`src/potpourri/models_multi_period/`)

```
Basemodel_multi_period    adds time index T, integrates simbench load/gen profiles
  └── ACOPF_multi_period  (AC_multi_period + OPF_multi_period)

Flexibility_multi_period  abstract base for all flexible devices
  ├── Battery_multi_period
  ├── EV_multi_period / Electric_vehicle_multi_period
  ├── HeatPump_multi_period
  ├── PV_multi_period
  ├── Windpower_multi_period
  ├── Demand_multi_period
  ├── Sgens_multi_period
  └── Generator_multi_period
```

Flexible devices are **composed, not inherited** — each is instantiated separately and attaches its own Pyomo Sets/Params/Vars/Constraints to the parent model.

### Data flow

```
pandapower net
  → Basemodel.__init__()    pp.runpp(), extract admittance data
  → create_model()          Pyomo ConcreteModel + sets/params/vars
  → add_OPF()               unfix controllable vars, add limits/objectives
  → .solve(solver)          SolverFactory → NLP/MIP
  → pyo_to_net()            write solution back to net.res_*
```

---

## External Dependencies

| Package | Role |
|---|---|
| `pandapower >= 3.1.2` | Network data model, initial power flow |
| `pyomo >= 6.9.4` | Optimization modelling |
| `simbench >= 1.6.1` | Benchmark networks and time-series profiles |
| `numpy`, `pandas`, `numba` | Numerical/data processing |
| `matplotlib` | Plotting |

### Solvers

| Solver | Type | Notes |
|---|---|---|
| IPOPT | NLP (primary) | installed via conda-forge; Docker compiles 3.14.16 from source |
| GLPK | LP/MILP | conda-forge |
| CBC | MILP | Docker only |
| SHOT | MINLP | Docker only |
| Gurobi | MIP/QP | optional; requires licence (`gurobipy`) |
| NEOS | Remote | `SolverManagerFactory('neos')` |

---

## Development

```bash
black .     # format
flake8 .    # lint
pytest      # run tests
```

Validation and analysis scripts are in `scripts/`. See the [Scripts documentation](docs/scripts/validation.md) for details.

---

## Authors

- Steffen Kortmann — Institute for High Voltage Equipment and Grids, Digitalisation and Energy Economics (IAEW), RWTH Aachen University
- Andreas Bong — IAEW, RWTH Aachen University
- Simon Braun — IAEW, RWTH Aachen University
- Alexander Och — IAEW, RWTH Aachen University
- Farah Nasr — IAEW, RWTH Aachen University
- Philip Kvesic — IAEW, RWTH Aachen University

