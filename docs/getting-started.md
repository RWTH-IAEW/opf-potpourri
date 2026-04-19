# Getting Started

## Requirements

- Python ≤ 3.12
- Conda (Miniconda or Anaconda)
- At least one supported solver (IPOPT recommended)

## Installation

### 1. Create the Conda environment

The `environment.yaml` file installs all Python dependencies and the GLPK solver:

```bash
conda env create -f environment.yaml
conda activate potpourri_env
```

To update an existing environment:

```bash
conda env update -f environment.yaml
```

### 2. Install the package

```bash
pip install -e .
```

### 3. Install a solver

POTPOURRI requires an external solver. IPOPT is recommended for AC OPF:

```bash
conda install -c conda-forge ipopt   # free, nonlinear
conda install -c conda-forge glpk    # free, linear/MIP
```

For commercial mixed-integer problems, Gurobi can be installed via:

```bash
pip install gurobipy   # requires a valid Gurobi licence
```

A Dockerfile is provided for a fully containerised environment that includes IPOPT 3.14.16 compiled from source, CBC, and the SHOT solver:

```bash
docker build -t potpourri:latest .
```

## Verify installation

```python
import pandapower as pp
import simbench as sb
from potpourri.models.ACOPF_base import ACOPF

net = sb.get_simbench_net("1-LV-rural1--0-sw")
opf = ACOPF(net)
opf.add_OPF()
opf.add_voltage_deviation_objective()
opf.solve(solver="ipopt", print_solver_output=False)
print(net.res_bus.head())
```

If the solver returns an optimal solution and `net.res_bus` is populated, the installation is working correctly.

## Development setup

```bash
black .     # format code
flake8 .    # lint
pytest      # run tests
```

Validation scripts are available in `scripts/`. `scripts/test_classbased.py` compares Pyomo solutions against pandapower power flow results.
