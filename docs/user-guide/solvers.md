# Solving Models

Every POTPOURRI model (single-period and multi-period) exposes a `solve()` method that wraps Pyomo's `SolverFactory` and handles result mapping back to the pandapower network.

## Basic usage

```python
opf.solve(solver='ipopt', print_solver_output=False)
```

After a successful solve, `net.res_bus`, `net.res_line`, etc. are populated automatically (unless `to_net=False`).

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `solver` | `str` | `'ipopt'` | Solver name: `'ipopt'`, `'glpk'`, `'cbc'`, `'gurobi'`, `'mindtpy'`, or `'neos'` |
| `print_solver_output` | `bool` | `False` | Stream solver log to stdout |
| `to_net` | `bool` | `True` | Write solution back to `net.res_*` DataFrames |
| `load_solutions` | `bool` | `True` | Load variable values from solver into Pyomo model |
| `time_limit` | `int` | `600` | Wall-clock time limit in seconds |
| `max_iter` | `int` | `None` | Maximum solver iterations (passed as `max_iter` option) |
| `mip_solver` | `str` | `'gurobi'` | MIP sub-solver for `mindtpy` |
| `init_strategy` | `str` | `'rNLP'` | Initialisation strategy for `mindtpy` |
| `neos_opt` | `str` | `'ipopt'` | Solver requested from the NEOS server |

## Local solvers

Install one of the following before calling `solve()`:

| Solver | Problem type | Install |
|--------|-------------|---------|
| **IPOPT** | NLP (AC OPF) | `conda install -c conda-forge ipopt` |
| **GLPK** | LP / MIP (DC OPF) | `conda install -c conda-forge glpk` |
| **CBC** | LP / MIP | `conda install -c conda-forge coincbc` |
| **Gurobi** | LP / MIP / NLP | `pip install gurobipy` (licence required) |

### Continuous NLP (AC OPF)

```python
opf.solve(solver='ipopt', time_limit=300)
```

### Linear / MIP (DC OPF, hosting capacity binary relaxation)

```python
dcopf.solve(solver='glpk')
```

### Mixed-integer nonlinear (HC_ACOPF with binary variables)

Use the `mindtpy` decomposition algorithm, which alternates between an NLP and a MIP sub-solver:

```python
hc.solve(
    solver='mindtpy',
    mip_solver='gurobi',   # or 'glpk', 'cbc'
    max_iter=50,
    init_strategy='rNLP',
)
```

## NEOS — remote solver

[NEOS](https://neos-server.org/) is a free public optimisation server. It accepts Pyomo models over the network and returns results without requiring a local solver installation.

### Requirements

NEOS requires a registered e-mail address, passed via the `NEOS_EMAIL` environment variable:

```bash
export NEOS_EMAIL="your@email.address"
```

Or set it in Python before calling `solve()`:

```python
import os
os.environ["NEOS_EMAIL"] = "your@email.address"
```

### Usage

```python
opf.solve(solver='neos', neos_opt='ipopt')
```

`neos_opt` selects the solver on the NEOS side. Common choices are `'ipopt'` (NLP), `'knitro'` (NLP), `'cplex'` (MIP), and `'couenne'` (MINLP).

!!! note
    NEOS solves run on shared public infrastructure. For large multi-period problems or time-sensitive work, a local solver installation is strongly preferred.

## Checking the result

The raw Pyomo result object is stored on the model instance:

```python
import pyomo.environ as pe

opf.solve(solver='ipopt')

print(opf.results.solver.termination_condition)   # TerminationCondition.optimal
print(pe.value(opf.model.obj))                     # objective value
```

If the solve is **infeasible** or the solver **times out**, `to_net` mapping is skipped and a warning is printed. Increase `time_limit` or relax constraints if this happens.
