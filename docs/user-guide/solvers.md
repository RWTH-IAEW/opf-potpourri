# Solving Models

Every `potpourri` model (single-period and multi-period) exposes a `solve()` method that wraps Pyomo's `SolverFactory` and handles result mapping back to the pandapower network.

## Basic usage

```python
opf.solve(solver='ipopt', print_solver_output=False)
```

After a successful solve, `net.res_bus`, `net.res_line`, etc. are populated automatically (unless `to_net=False`).

### Multi-period: the solve is warm-started

`solve()` seeds every state variable from a pandapower power flow at each time
step before handing the model to the solver, controlled by `warm_start` (default
`True`).

This is a correctness matter, not a speed optimisation. A cold-started
multi-period AC OPF begins with `v = 1`, every angle at zero and no branch-flow
values at all, so Kirchhoff's laws are violated at every bus by the full nodal
injection — and IPOPT does not always recover. On a 12-step midday window of
`1-LV-rural1--0-sw` it reported a locally infeasible point even though curtailing
the PV to zero is available and feasible. With the seed it converges. (It is also
markedly faster: the unit suite went from ~8 to ~3 minutes.)

The seed does not need to be near the optimum — an uncurtailed and a fully
curtailed seed converge to the same solution — it only has to satisfy the power
flow.

```python
mpopf.solve(solver="ipopt")                    # seeded (default)
mpopf.solve(solver="ipopt", warm_start=False)  # keep existing initial values
mpopf.warm_start_from_pf()                     # seed explicitly
```

Pass `warm_start=False` when the variables already hold something better, such
as a previous solution you want to re-solve from.

!!! note "Compare the model's objective, not a hand-rolled metric"

    `add_voltage_deviation_objective()` is **not** `Σ(v − 1)²` over every bus.
    It leaves the slack out of that term and instead penalises the slack's
    distance from its base-case magnitude `v_b0`:

    ```
    Σ_{b ∉ b0} Σ_t (v[b,t] − 1)²  +  Σ_{b ∈ b0} Σ_t (v[b,t] − v_b0[b])²
    ```

    With `free_slack_vm=True` the slack magnitude moves, so a re-derived
    `Σ(v − 1)²` over *all* buses measures a different quantity and can move in
    the **opposite direction** to the value actually being minimised. Read
    `pyo.value(model.obj_v_deviation)` when comparing two runs.

    The AC OPF is also nonconvex, so IPOPT returns a local solution — worth
    keeping in mind, though on the cases tested here adding a device does
    improve the objective as expected.

### Multi-period: one time step at a time

The pandapower result tables have no time dimension, so a multi-period solve
cannot write a whole horizon into them. `to_net=True` therefore maps the
**last** time step of the horizon; pass an integer to pick a different one:

```python
mpopf.solve(solver='ipopt', to_net=2)   # write time step 2
```

To read several steps, call `map_to_net()` for each and consume `net.res_*`
in between:

```python
for t in mpopf.model.T:
    mpopf.map_to_net(t)
    print(t, mpopf.net.res_bus.vm_pu.max())
```

`map_to_net()` returns the step it wrote, and raises `ValueError` for a step
outside `model.T`.

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `solver` | `str` | `'ipopt'` | Solver name: `'ipopt'`, `'glpk'`, `'cbc'`, `'gurobi'`, `'gurobi_direct_minlp'`, `'mindtpy'`, or `'neos'` |
| `print_solver_output` | `bool` | `False` | Stream solver log to stdout |
| `to_net` | `bool` \| `int` | `True` | Write solution back to `net.res_*` DataFrames. Multi-period: `True` writes the last time step, an `int` writes that step |
| `warm_start` | `bool` | `True` | Multi-period only: seed every state variable from a per-step power flow before solving |
| `load_solutions` | `bool` | `True` | Load variable values from solver into Pyomo model |
| `time_limit` | `int` | `600` | Wall-clock time limit in seconds. Honoured by `mindtpy` and by `gurobi*` (sent as `TimeLimit`); **ignored by IPOPT, GLPK and CBC** |
| `max_iter` | `int` | `None` | Maximum solver iterations (sent as `max_iter`, or as `IterationLimit` for `gurobi*`) |
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
| **Gurobi** | LP / MIP / MINLP | `pip install "gurobipy>=12"` (licence required) |

### Continuous NLP (AC OPF)

```python
opf.solve(solver='ipopt', time_limit=300)
```

### Linear / MIP (DC OPF, hosting capacity binary relaxation)

```python
dcopf.solve(solver='glpk')
```

### Mixed-integer nonlinear (HC_ACOPF, discrete tap changers)

Two routes are available.

**Gurobi global MINLP (recommended).** Pyomo's `gurobi_direct_minlp` interface maps Pyomo's `sin`, `cos`, `exp`, `log` and `sqrt` onto `gurobipy.nlfunc`, so the polar-form AC power flow goes to Gurobi unchanged and is solved by global spatial branch-and-bound:

```python
hc.solve(solver='gurobi_direct_minlp', time_limit=300)
```

Requires **Pyomo >= 6.10** and **gurobipy >= 12**. The older `gurobi`, `gurobi_direct` and `gurobi_persistent` interfaces are limited to expressions of degree 2 and reject the AC power flow with `DegreeError` — use them only for `DCOPF`.

Unsupported functions: `asin`, `acos`, `atan`, `sinh`, `cosh` and their inverses are not in Pyomo's dispatcher and will fail. The AC models use only `sin`/`cos`, so this does not affect them.

**MindtPy decomposition.** Alternates between an NLP and a MIP sub-solver:

```python
hc.solve(
    solver='mindtpy',
    mip_solver='gurobi',   # or 'glpk', 'cbc'
    max_iter=50,
    init_strategy='rNLP',
)
```

Outer approximation is only globally valid for *convex* MINLP. The AC OPF is nonconvex, so MindtPy returns a local solution and may terminate `feasible` rather than `optimal`.

### Choosing between them

Measured on SimBench `1-LV-rural1--0-sw` (15 buses), 120 s limit:

| Model | `gurobi_direct_minlp` | `mindtpy` |
|-------|----------------------|-----------|
| ACOPF (continuous) | same optimum as IPOPT, but hits the time limit proving global optimality | — (use IPOPT: ~10 s) |
| ACOPF + discrete tap | `maxTimeLimit`, obj `7.7e-5` | `feasible`, obj `8.2e-3` (~100× worse) |
| HC_ACOPF (14 binaries) | `optimal` in ~8 s | fails on uninitialised `pLfrom` |

Rules of thumb: use **IPOPT** for continuous AC OPF — Gurobi finds the same optimum but burns the rest of the budget closing the global gap. Use **`gurobi_direct_minlp`** once integer variables are present. Always pass `time_limit`, and check `termination_condition`: `maxTimeLimit` means the incumbent is feasible but not proven globally optimal.

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
