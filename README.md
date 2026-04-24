# potpourri

**Multi-Period Optimal Power Flow for Distribution Grids with Storage Application**

> *Potpourri — piece of music composed from various popular smaller works or melodies*

[![CI](https://github.com/RWTH-IAEW/opf-potpourri/actions/workflows/ci.yml/badge.svg)](https://github.com/RWTH-IAEW/opf-potpourri/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

`potpourri` is a Python library for AC/DC Optimal Power Flow (OPF) in distribution grids, with support for multi-period planning and flexible resources (batteries, EVs, heat pumps, PV, wind). It wraps [Pyomo](https://pyomo.readthedocs.io/) for optimisation modelling over [pandapower](https://pandapower.readthedocs.io/) network objects.

**Repository:** <https://github.com/RWTH-IAEW/opf-potpourri>

---

## Documentation

The project documentation is built with [MkDocs](https://www.mkdocs.org/).

To serve the documentation locally:

```bash
pip install -e .[docs]
mkdocs serve          # usually available at http://127.0.0.1:8000/
```

---

## Installation

Requires Python 3.9–3.12 and Conda (or Mamba).

```bash
conda env create -f environment.yaml
conda activate potpourri_env
pip install -e .
```

To update an existing environment:

```bash
conda env update -f environment.yaml --prune
```

A Dockerfile is provided for a fully containerised setup with IPOPT 3.14.16
compiled from source, CBC, and SHOT solvers.

---

## Quick start

### Single-period AC OPF

```python
import simbench as sb
from potpourri.models.ACOPF_base import ACOPF

net = sb.get_simbench_net("1-LV-rural1--0-sw")
opf = ACOPF(net)
opf.add_OPF()
opf.add_voltage_deviation_objective()
opf.solve(solver="ipopt", print_solver_output=False)

# results available in net.res_bus, net.res_line, net.res_sgen, ...
print(opf.net.res_bus[["vm_pu", "va_degree"]])
```

### Multi-period AC OPF with battery storage

```python
import simbench as sb
from potpourri.models_multi_period.ACOPF_multi_period import ACOPF_multi_period
from potpourri.technologies.battery import Battery_multi_period

net = sb.get_simbench_net("1-LV-urban6--0-sw")
opf = ACOPF_multi_period(net, toT=96, fromT=0)   # 96 × 15 min = 1 day

battery = Battery_multi_period(opf.net, T=96, scenario=1)
battery.get_all(opf.model)

opf.add_OPF()
opf.add_voltage_deviation_objective()
opf.solve(solver="ipopt")
```

See `scripts/` for runnable examples covering each feature area.

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
Basemodel_multi_period    adds time index T, integrates SimBench profiles
  └── ACOPF_multi_period  (AC_multi_period + OPF_multi_period)

Flexibility_multi_period  abstract base for all flexible devices
  ├── Battery_multi_period
  ├── HeatPump_multi_period
  ├── PV_multi_period
  ├── Windpower_multi_period
  ├── Demand_multi_period
  ├── Sgens_multi_period
  └── Generator_multi_period
```

Flexible devices are **composed, not inherited** — each is instantiated
separately and attaches its own Pyomo Sets/Params/Vars/Constraints to the
parent model.

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

## External dependencies

| Package | Role |
|---|---|
| `pandapower >= 2.13` | Network data model, initial power flow |
| `pyomo >= 6.7` | Optimisation modelling |
| `simbench >= 1.4` | Benchmark networks and time-series profiles |
| `numpy`, `pandas` | Numerical / data processing |
| `matplotlib` | Plotting |

### Solvers

`potpourri` does not bundle any solvers. Install at least one before
calling `solve()`.

| Solver | Type | Install |
|---|---|---|
| **IPOPT** | NLP — AC OPF | `conda install -c conda-forge ipopt` |
| **GLPK** | LP / MIP — DC OPF | `conda install -c conda-forge glpk` |
| **CBC** | LP / MIP | `conda install -c conda-forge coincbc` |
| **Gurobi** | LP / MIP / NLP | `pip install gurobipy` (licence required) |
| **NEOS** | Remote (free) | `opf.solve(solver='neos', neos_opt='ipopt')` |

IPOPT and GLPK are included automatically when you create the environment
from `environment.yaml`.

---

## Development

```bash
pip install -e ".[dev]"   # installs ruff, pytest, pytest-cov, pre-commit
ruff check .              # lint
ruff format .             # format
pytest                    # run tests
pytest -m "not integration"   # skip solver-dependent tests
```

Analysis and example scripts are in `scripts/`. See `scripts/README.md`
for an overview of what each example demonstrates.

---

## Authors

- Steffen Kortmann — IAEW, RWTH Aachen University
- Andreas Bong — IAEW, RWTH Aachen University
- Simon Braun — IAEW, RWTH Aachen University
- Alexander Och — IAEW, RWTH Aachen University
- Farah Nasr — IAEW, RWTH Aachen University
- Philip Kvesic — IAEW, RWTH Aachen University
- Nina Stumberger — IAEW, RWTH Aachen University

---

## Citation

If you use `potpourri` in your research, please cite it using the metadata in
[`CITATION.cff`](CITATION.cff). A BibTeX entry will be available once a
Zenodo DOI is registered for the release.

---

## License

`potpourri` is released under the [MIT License](LICENSE).
