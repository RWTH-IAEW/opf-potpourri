# OPF diagnostics

An OPF that does not solve tells you almost nothing by itself. `infeasible`
is a verdict, not an explanation, and a solver has no idea that constraint
`line_lim_from[23]` is the from-side current limit of the cable feeding your
worst-supplied street.

`diagnose()` closes that gap. It reads the network, the Pyomo model, the
solver's verdict and the solution, and reports what it finds in terms of the
pandapower objects you built — while keeping the Pyomo identifier alongside,
so you can always get back to the equation.

```python
from potpourri.models.ACOPF_base import ACOPF

opf = ACOPF(net)
opf.add_OPF()
opf.add_voltage_deviation_objective()
opf.solve(solver="ipopt")

report = opf.diagnose()
print(report)
```

It is safe to call at any point: before `add_OPF()`, after a failed solve,
or after a successful one. Checks that cannot run in the current state say
so rather than raising, and nothing it does modifies your model.

---

## What you get back

```
potpourri OPF diagnostics
==============================================================

Model
  formulation ......... AC (ACOPF)
  free variables ...... 21
  equalities .......... 20
  inequalities ........ 12

Network
  buses ............... 4
  lines ............... 2

Solver
  termination ......... optimal
  status .............. ok

Balance
  external grid ....... 0.025 MW
  load ................ 0.06 MW
  losses .............. 0 MW
  residual ............ 2.1e-15 MW

Findings: 9 info, 0 warning, 0 error
```

A finding looks like this — the network object first, the Pyomo component
kept for the trail back:

```
[ERROR] THERMAL net.line[12] "Feeder 3-4"
  loading 108.3 %, above the element's rating
  value: 108.3 %
  limits: max 100
  violation: 8.3 % (8.30 %)
  pyomo: line_lim_from[12]
```

---

## Levels

Not every check is worth running every time.

| level | adds | cost |
|---|---|---|
| `basic` | network data, bounds, islands, adequacy, model size, solver verdict, solution violations | no solver, no power flow |
| `standard` *(default)* | binding limits, plausibility, scaling, power-flow cross-check | one power flow |
| `deep` | structural singularity, Jacobian conditioning, feasibility relaxation | extra solver calls |

```python
report = opf.diagnose(level="basic")     # safe inside a loop
report = opf.diagnose(level="deep")      # when you are stuck
```

`basic` never calls a solver, so it is cheap enough to run after every solve.
`deep` can take longer than the original solve.

---

## Working with the findings

The printed report is for reading. For anything else, the findings are typed:

```python
report.errors                  # only the ERROR-level findings
report.warnings
report.by_code("BUS_VOLTAGE_HIGH")
report.by_category(DiagnosticCategory.THERMAL)

report.to_dict()               # JSON-safe, no Pyomo objects
report.to_dataframe()          # one row per finding
```

Filter on `code`, never on `message`. Codes are stable; the prose is not.

Each finding carries, where it applies: the pandapower element and its name,
the Pyomo component and index, the measured value, the bounds, the violation
in absolute and relative terms, the unit, and the time step for multi-period
models.

---

## The five situations, and what to look for

### The OPF does not solve

Start at `basic`. The checks that need no solver find the contradictions
that need no solver:

- `BOUND_MIN_ABOVE_MAX` — a limit pair that nothing can satisfy.
- `NET_ISLAND_WITHOUT_SOURCE` — load with no supply and no reference.
- `ADEQUACY_INSUFFICIENT_GENERATION` — the most the sources can produce is
  below the least the loads can take.
- `VOLTAGE_SETPOINT_OUTSIDE_BAND` — a slack held outside its own band.

Any of these is enough on its own to make the problem infeasible.

### The solver says "infeasible"

Read the termination condition carefully, because **"locally infeasible" is
not "infeasible"**. IPOPT solves a nonconvex problem by local search and
reports where *it* got stuck, not a property of your network. A different
starting point can succeed on the same model.

At `deep`, the feasibility relaxation asks the useful question — what would
have to give?

```python
report = opf.diagnose(level="deep")
for issue in report.by_code("RELAXATION_LIMIT_RELAXED"):
    print(issue.element, issue.violation)
```

It clones your model, softens the limits that can physically be relaxed, and
minimises the total violation. Read the answer as *limits in tension*, not as
a cause: it is one relaxation out of many, the weights chose it, and on a
nonconvex problem it is local. Nodal balance is never relaxed, because power
appearing from nowhere explains nothing.

### It fails numerically

A perfectly feasible OPF can fail because its numbers span too many orders of
magnitude. That looks like infeasibility from outside, which is exactly why
it gets its own category:

- `NUMERIC_POOR_SCALING` — coefficients within one constraint family spread
  over more than ten decades.
- `NUMERIC_ILL_CONDITIONED` — the equality Jacobian is close to singular.
- `NUMERIC_NON_FINITE_VALUE` — a NaN or infinity in the model.

These say "the solver could not work with this", not "your network is
impossible".

### It stops without an optimal solution

`SOLVER_ITERATION_LIMIT` and `SOLVER_TIME_LIMIT` are not infeasibility. Raise
the limit, improve the starting point, or look at the scaling findings —
slow convergence is more often conditioning than feasibility.

### The result looks wrong

This is what the `standard` level is for.

- `SOLUTION_CONSTRAINT_BINDING` lists the limits that are shaping the
  answer. A dispatch that looks odd is usually a dispatch pressed against
  something.
- `REPLAY_AGREES` / `REPLAY_MISMATCH` re-solves your optimised dispatch with
  an independent pandapower power flow. On an AC model the two should agree
  to solver tolerance; a gap points at a modelling or result-mapping
  difference, not at a bad network.
- `RESULT_IMPLAUSIBLE_LOSSES`, `RESULT_BALANCE_MISMATCH` — system-level
  sanity in MW and MVAr.

For DC and LPAC the replay is labelled as a plausibility check instead: those
formulations are approximations on purpose, so a difference is expected and
measures the approximation.

---

## Looking around a problem bus

```python
from potpourri.diagnostics.context import DiagnosticContext
from potpourri.diagnostics.solution import explain_bus

ctx = DiagnosticContext.from_model(opf)
explain_bus(ctx, 17)
```

Returns the bus voltage, the branches attached to it and every load,
generator, static generator and storage unit sitting on it, each with its
operating point.

---

## What the diagnostics cannot tell you

Worth being explicit, because a tool like this invites over-reading.

**Nothing here proves feasibility.** The adequacy and island checks are
*necessary* conditions: failing one proves the OPF is infeasible, passing
them all proves nothing, because they ignore the network's impedances
entirely.

**Nothing here establishes causality.** The report says "these limits are in
tension" and "this is the largest violation". It does not say "line 8 caused
the failure", because with a nonconvex problem and a local solver, that is
usually not knowable.

**A local solver gives local answers.** An IPOPT result — optimal or
infeasible — is about the point it reached. It is not a statement about the
global problem.

---

## What each capability needs

| capability | needs | degrades to |
|---|---|---|
| network, bounds, islands, adequacy | pandapower only | — |
| model size, unused variables | Pyomo only | — |
| solution violations, binding limits | a loaded solution | skipped, with the reason |
| power-flow cross-check | a converged pandapower PF | reported as a finding |
| structural singularity | `pyomo.contrib.incidence_analysis` | skipped |
| Jacobian conditioning | SciPy | skipped |
| feasibility relaxation | a solver, one extra solve | skipped |
| Pyomo infeasibility explanation | `pyomo.contrib.iis`, many solves | skipped |

Everything optional degrades gracefully. A missing capability is recorded in
`report.skipped` with the reason, never silently dropped:

```python
>>> report.skipped
{'jacobian': 'SciPy is required for this check: ...'}
```
