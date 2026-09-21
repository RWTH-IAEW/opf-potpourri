# Pyomo conventions

potpourri follows the community
[Pyomo style guide](https://mobook.github.io/MO-book/notebooks/appendix/pyomo-style-guide-update.html).
This page records the parts that apply here, the two places we knowingly
depart from it, and what the linter enforces so none of it decays.

---

## The rules

| | write | not |
|---|---|---|
| import | `import pyomo.environ as pyo` | `from pyomo.environ import *`, `as pe`, explicit-name imports |
| components | `@model.Constraint(model.L)` on the rule | `rule=` plus a separate registration |
| a `Var`'s domain | `domain=pyo.Reals` | `within=pyo.Reals` |
| model kind | `pyo.ConcreteModel` | `AbstractModel` |
| constraints | `pyo.Constraint` | `pyo.ConstraintList` |

`ruff check` fails on a star import (`F403`/`F405`); the rest is
convention, kept by review.

---

## Decorators, and why component names matter

Under `rule=`, the component name is the assignment target and the rule
function's name is free. Under a decorator, **the function name becomes
the component name**:

```python
@self.model.Constraint(self.model.L)
def line_lim_from(model, l):
    return (
        model.pLfrom[l] ** 2 + model.qLfrom[l] ** 2 <= model.SLmax[l] ** 2
    )
```

`self.model.line_lim_from` is now a `Constraint`, and the name
`line_lim_from` in the enclosing scope is bound to that component rather
than to a function.

**Component names are public API.** Tests reach into
`acopf.model.line_lim_from[0]`, scripts call
`ac.model.transf_lim1.deactivate()`, docstrings and the architecture page
name them. Renaming one is a breaking change and belongs in its own
commit with a note in the changelog — never as a side effect of a style
edit.

### The owner expression

The decorator is written on whatever expression already holds the model:

```python
@self.model.Constraint(self.model.L)   # inside a model class
@model.Constraint(model.sGc, model.T)  # inside a device mix-in, which takes `model`
```

Do **not** open a method with `model = self.model` and decorate with
`@model....`. It looks tidier and it is a trap: `create_model()` overrides
call `super().create_model()`, which *replaces* `self.model` partway
through the method, so a binding taken at the top would point at the
previous object.

### Ordering

A decorator constructs the component where the `def` is, not where the
old registration was. When you move a rule, keep it after everything its
index sets and parameters depend on, and keep sibling components in their
original relative order — component order is what the writers hand to the
solver.

Conditional components keep working. A rule defined once per branch of an
`if` carries its own decorator in each branch, which puts the
branch-specific rule next to the condition that selects it:

```python
if thermal_limit == "current":

    @self.model.Constraint(self.model.L)
    def line_lim_from(model, l):
        ...

else:

    @self.model.Constraint(self.model.L)
    def line_lim_from(model, l):
        ...
```

A rule can still return `pyo.Constraint.Skip` to leave an index out.

### When `rule=` has to stay

A decorator can only give a component the name of the function it
decorates. Where the name is computed at run time, `rule=` is correct and
should stay — with a comment saying so. Two places in the package qualify:
`technologies/q_control.py`, which builds `f"{name}_v_link"` from an
argument and registers it with `model.add_component`, and
`Basemodel_multi_period._build_kcl`, which loops over constraint names and
looks each rule up with `getattr`.

---

## Naming

Names in this codebase are transliterated from the power-flow equations —
`qsG`, `pTlv`, `SLmax`, `v`, `delta` — rather than spelled out as
`total_reactive_generation`. That is a deliberate departure from the
guide's naming section: the names match the equations in the
[mathematical modelling](mathematical-modelling.md) page, `pyo_to_net.py`
matches on them, and they appear across `scripts/`, the tutorials and the
published docs.

`E741` stays off in ruff for the same reason: `l` is the standard line
index in power-systems notation.

Sets are ALL_CAPS (`B`, `L`, `TRANSF`, `T`, `WINDc`). A rule's first
parameter is the model the component is being built on, conventionally
named `model` (or `mm` where an outer `m` is already bound).

---

## Model construction

The guide recommends a `build_model(data) -> ConcreteModel` function.
potpourri uses classes plus mix-ins instead — `Basemodel` builds the
model, `AC`/`DC` add power flow, `OPF` adds limits, and device classes
attach their own components to an existing model in `get_all(model)`.
That is a deliberate alternative, not an oversight: it is what supports
attaching a battery or heat pump to a model after it is built. See
[Class Architecture](architecture.md).

---

## Checking your change

```bash
ruff check .                  # star imports, docstrings, line length
ruff format --check .
pytest -m "not integration"   # construction and naming regressions
```

A style change must leave the generated model **numerically identical**.
The test suite catches naming and construction errors; for anything that
touches an equation, also solve a reference case before and after and
compare the objective and `net.res_*` to solver tolerance — for example
with `scripts/acopf_loadcase_analysis.py` and
`scripts/multi_period_acopf.py`.
