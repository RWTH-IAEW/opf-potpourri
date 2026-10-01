# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

r"""Opt-in control of switched shunts (capacitor banks, reactors).

pandapower's `net.shunt` is a constant admittance: `p_mw` and `q_mvar` are
the consumption *per step* at 1 p.u., multiplied by the integer `step` and
by $\rho = (V^{\text{bus}}_n / \texttt{vn\_kv})^2$, so the consumption at
voltage $v$ is $S = (p + jq)\,\text{step}\,\rho\, v^2$. The AC layers carry
that as the fixed parameters `GB` / `BB` in the nodal balance. This module
makes the step of selected shunts a decision variable:

$$
p^{sh}_s = \frac{p_s \rho_s}{S_N}\, m_s\, v_b^2, \qquad
q^{sh}_s = \frac{q_s \rho_s}{S_N}\, m_s\, v_b^2, \qquad
m_s \in [0, m^{\max}_s],
$$

with $m_s$ integer (`mode="discrete"`, a physical bank) or real
(`mode="continuous"`, the relaxation), in the load sign convention of
`net.shunt` (positive `q_mvar` is an inductive reactor, negative a
capacitor). The products $m_s v_b^2$ are bilinear, so a discrete bank makes
the AC OPF a nonconvex MINLP and a continuous one adds a bilinear term to the
NLP; see `docs/research/dso_controllable_equipment.md`.

Movement over a horizon (`shunt_step_up` / `shunt_step_down`), a per-step
change limit, an operation limit and a priced switching term mirror the
on-load tap changer model in [`potpourri.models.oltc`][potpourri.models.oltc],
and `solve_discrete_round_and_fix` rounds both families in one pass.

Nothing here runs unless `enable_shunt_control` is called: without it the
balance keeps its constant shunt terms and this module contributes nothing.
A shunt with `step_dependency_table=True` (table-driven values) is
rejected rather than approximated.

Example:
    ```python
    opf = ACOPF(net)
    opf.add_OPF()
    opf.enable_shunt_control(mode="discrete")   # every eligible shunt
    opf.add_voltage_deviation_objective()
    opf.solve(solver="gurobi_direct_minlp")
    opf.net.res_shunt[["step", "q_mvar", "vm_pu"]]
    ```
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass

import numpy as np
import pandas as pd
import pyomo.environ as pyo
from loguru import logger

from potpourri.models.oltc import (
    OLTC_MODES,
    RoundingGroup,
    _column,
    _finite_integer,
    _is_missing,
    _per_transformer,
    _trafo_of,
)


def shunt_control_terms(model, b, t=None, reactive=False):
    r"""Consumption of the *controllable* shunts at bus `b`, in p.u.

    Returns the sum $\sum_s c_s\, m_s\, v_b^2$ over the controllable
    shunts at the bus, where $c_s$ is `shunt_p_step` (active) or
    `shunt_q_step` (reactive) and $m_s$ the step variable. It is written
    on the consumption side of the balance like the fixed `GB`/`BB` terms it
    replaces for those shunts, so the sign convention is pandapower's: a
    positive `q_mvar` (reactor) consumes reactive power.

    Args:
        model: The Pyomo model being built.
        b: Bus index from `model.B` (a ppc bus number).
        t: Time index for a multi-period model, `None` for single period.
        reactive: `True` for the reactive term, `False` for the active one.

    Returns:
        A Pyomo expression, or the integer `0` when the model has no
        controllable shunt at that bus (so the caller can add it to a sum
        without changing the expression).
    """
    ctrl = getattr(model, "SHUNT_CTRL", None)
    if ctrl is None:
        return 0
    coefficient = model.shunt_q_step if reactive else model.shunt_p_step
    total = 0
    for s in ctrl:
        if (b, s) not in model.SHUNTbs:
            continue
        if t is None:
            step, v = model.shunt_step[s], model.v[b]
        else:
            step, v = model.shunt_step[s, t], model.v[b, t]
        total = total + coefficient[s] * step * v**2
    return total


# ---------------------------------------------------------------------------
# eligibility and data
# ---------------------------------------------------------------------------


@dataclass
class ShuntData:
    r"""Step data of one controllable shunt.

    Attributes:
        index: Shunt index in `net.shunt` and `model.SHUNT`.
        bus_ppc: ppc bus the shunt sits on.
        step_max: Highest admissible step (`max_step`).
        step_init: Step the network was built with.
        p_step_pu: Active consumption per step at 1 p.u., p.u. on the base,
            including pandapower's $(V_n^{bus}/\texttt{vn\_kv})^2$ factor.
        q_step_pu: Reactive consumption per step, same convention.
    """

    index: int
    bus_ppc: int
    step_max: int
    step_init: int
    p_step_pu: float
    q_step_pu: float


@dataclass
class ShuntControlSetup:
    """What `enable_shunt_control` built, kept as `shunt_control_setup`.

    Attributes:
        mode: `"continuous"` or `"discrete"`.
        shunts: Controlled shunt indices, in order.
        data: `ShuntData` per controlled shunt.
        multi_period: Whether the steps carry a time index.
        initial: Reference step per shunt for the first move, or `None`.
        max_change_per_step: Per-step change limit per shunt, if any.
        max_operations: Operation limit per shunt, if any.
        cost_objective: Objective the switching cost was added to, if any.
    """

    mode: str
    shunts: list
    data: dict
    multi_period: bool
    initial: dict | None
    max_change_per_step: dict | None = None
    max_operations: dict | None = None
    cost_objective: str | None = None


def _ineligibility_reason(row, dependency):
    """The first reason a shunt row cannot be step-controlled.

    Args:
        row: A row of `net.shunt`.
        dependency: Its `step_dependency_table` cell.

    Returns:
        A short explanation, or `""` when the shunt is eligible.
    """
    if not bool(row["in_service"]):
        return "out of service"
    if not _is_missing(dependency) and bool(dependency):
        return (
            "step_dependency_table is True (table-driven step values), "
            "which is not supported"
        )
    if not _finite_integer(row.get("max_step")) or int(row["max_step"]) < 1:
        return "max_step is missing or smaller than 1"
    if not _finite_integer(row.get("step")):
        return "step is missing or not an integer"
    if not 0 <= int(row["step"]) <= int(row["max_step"]):
        return "step lies outside [0, max_step]"
    p = row.get("p_mw")
    q = row.get("q_mvar")
    for name, value in (("p_mw", p), ("q_mvar", q)):
        if _is_missing(value) or not math.isfinite(float(value)):
            return f"{name} is missing"
    if float(p) == 0.0 and float(q) == 0.0:
        return "p_mw and q_mvar are both zero, so a step changes nothing"
    vn = row.get("vn_kv")
    if not _is_missing(vn) and float(vn) <= 0:
        return "vn_kv must be positive"
    return ""


def shunt_eligibility(net) -> pd.DataFrame:
    """Which shunts `enable_shunt_control` can control, and why not otherwise.

    A shunt is eligible when it is in service, `max_step >= 1` and
    `0 <= step <= max_step` are integers, `p_mw` and `q_mvar` are present
    and not both zero, and `step_dependency_table` is not set.

    Args:
        net: A pandapower network.

    Returns:
        A DataFrame indexed like `net.shunt` with the columns `eligible`,
        `reason`, `bus`, `step`, `max_step`, `p_mw`, `q_mvar`.
    """
    shunt = net.shunt
    dependency = _column(shunt, "step_dependency_table", False)
    records = []
    for idx, row in shunt.iterrows():
        reason = _ineligibility_reason(row, dependency[idx])
        records.append(
            {
                "eligible": reason == "",
                "reason": reason,
                "bus": row.get("bus"),
                "step": row.get("step"),
                "max_step": row.get("max_step"),
                "p_mw": row.get("p_mw"),
                "q_mvar": row.get("q_mvar"),
            }
        )
    report = pd.DataFrame(
        records,
        index=shunt.index,
        columns=[
            "eligible",
            "reason",
            "bus",
            "step",
            "max_step",
            "p_mw",
            "q_mvar",
        ],
    )
    report.index.name = "shunt"
    return report


def _shunt_data(model_obj, index) -> ShuntData:
    """Collect the step data of one eligible shunt.

    Args:
        model_obj: The model object (its `net`, `bus_lookup`, `baseMVA`).
        index: Shunt index.

    Returns:
        The collected `ShuntData`.
    """
    row = model_obj.net.shunt.loc[index]
    bus = int(row["bus"])
    bus_vn = float(model_obj.net.bus.at[bus, "vn_kv"])
    vn = row.get("vn_kv")
    rho = 1.0 if _is_missing(vn) else (bus_vn / float(vn)) ** 2
    base = float(model_obj.baseMVA)
    return ShuntData(
        index=int(index),
        bus_ppc=int(model_obj.bus_lookup[bus]),
        step_max=int(row["max_step"]),
        step_init=int(row["step"]),
        p_step_pu=float(row["p_mw"]) * rho / base,
        q_step_pu=float(row["q_mvar"]) * rho / base,
    )


def _select_shunts(model_obj, shunts):
    """Resolve the `shunts` argument of `enable_shunt_control`.

    Args:
        model_obj: The model object.
        shunts: `None` for every eligible shunt, or an iterable of indices.

    Returns:
        The chosen indices, in order and without duplicates.

    Raises:
        ValueError: If nothing is eligible or a requested shunt is not.
    """
    report = shunt_eligibility(model_obj.net)
    in_model = set(int(s) for s in model_obj.model.SHUNT)
    if shunts is None:
        chosen = [
            int(s)
            for s in report.index
            if int(s) in in_model and bool(report.at[s, "eligible"])
        ]
        if not chosen:
            reasons = "\n".join(
                f"  shunt {s}: {report.at[s, 'reason']}" for s in report.index
            )
            raise ValueError(
                "No shunt is eligible for step control:\n"
                + (reasons or "  the network has no shunt")
            )
        return chosen
    chosen = list(dict.fromkeys(int(s) for s in shunts))
    if not chosen:
        raise ValueError("shunts must name at least one shunt")
    problems = []
    for s in chosen:
        if s not in report.index:
            problems.append(f"shunt {s}: not a shunt of the network")
        elif s not in in_model:
            problems.append(f"shunt {s}: not in the model (out of service)")
        elif not bool(report.at[s, "eligible"]):
            problems.append(f"shunt {s}: {report.at[s, 'reason']}")
    if problems:
        raise ValueError(
            "Shunt(s) cannot be step-controlled:\n"
            + "\n".join("  " + p for p in problems)
        )
    return chosen


# ---------------------------------------------------------------------------
# attaching the control to a model
# ---------------------------------------------------------------------------


def attach_shunt_control(
    model_obj,
    shunts=None,
    mode="discrete",
    *,
    max_change_per_step=None,
    max_operations=None,
    initial_step="net",
):
    """Make the steps of selected shunts decision variables.

    The implementation behind `ShuntControlMixin.enable_shunt_control`.
    Adds `SHUNT_CTRL`, the per-step consumption parameters, the step and
    movement variables and their constraints, then rebuilds the nodal
    balance so the step variables replace the constant shunt terms.

    Args:
        model_obj: An `OPF`-derived model object whose class sets
            `SHUNT_CONTROL_SUPPORTED`.
        shunts: `None` (every eligible shunt) or indices.
        mode: `"continuous"` or `"discrete"`.
        max_change_per_step: Largest step change between consecutive time
            steps (single period: away from the initial step); scalar or
            per-shunt mapping.
        max_operations: Largest total number of step changes over the
            horizon; scalar or mapping.
        initial_step: `"net"` to measure movement from `net.shunt.step`, a
            mapping of reference steps, or `None` to leave the first move
            untracked.

    Returns:
        The `ShuntControlSetup` describing what was built.

    Raises:
        NotImplementedError: On a formulation without a voltage magnitude.
        ValueError: On a bad `mode` or an ineligible shunt.
        RuntimeError: If already enabled on this model.
    """
    if not getattr(model_obj, "SHUNT_CONTROL_SUPPORTED", False):
        raise NotImplementedError(
            f"{type(model_obj).__name__} has no reactive power or voltage "
            "magnitude to control a shunt with: shunt control needs the "
            "polar AC equations (ACOPF, ACOPF_multi_period)."
        )
    if mode not in OLTC_MODES:
        raise ValueError(f"mode must be one of {OLTC_MODES}, got {mode!r}")
    model = model_obj.model
    if hasattr(model, "SHUNT_CTRL"):
        raise RuntimeError(
            "enable_shunt_control() has already been called on this model"
        )
    time_set = getattr(model, "T", None)
    multi = time_set is not None
    chosen = _select_shunts(model_obj, shunts)
    data = {s: _shunt_data(model_obj, s) for s in chosen}

    if initial_step is None:
        initial = None
    elif isinstance(initial_step, str) and initial_step == "net":
        initial = {s: data[s].step_init for s in chosen}
    else:
        initial = {}
        for s in chosen:
            if s not in initial_step:
                raise ValueError(f"initial_step has no entry for shunt {s}")
            if not _finite_integer(initial_step[s]):
                raise ValueError(
                    f"initial_step for shunt {s} must be an integer, got "
                    f"{initial_step[s]!r}"
                )
            initial[s] = int(initial_step[s])
    change_max = (
        None
        if max_change_per_step is None
        else _per_transformer(
            max_change_per_step, chosen, "max_change_per_step"
        )
    )
    operations_max = (
        None
        if max_operations is None
        else _per_transformer(max_operations, chosen, "max_operations")
    )
    track_first = initial is not None

    model.SHUNT_CTRL = pyo.Set(within=model.SHUNT, initialize=chosen)
    model.shunt_step_max = pyo.Param(
        model.SHUNT_CTRL,
        within=pyo.Integers,
        initialize={s: data[s].step_max for s in chosen},
    )
    model.shunt_step_init = pyo.Param(
        model.SHUNT_CTRL,
        within=pyo.Integers,
        initialize={
            s: (initial[s] if initial is not None else data[s].step_init)
            for s in chosen
        },
    )
    model.shunt_p_step = pyo.Param(
        model.SHUNT_CTRL,
        within=pyo.Reals,
        initialize={s: data[s].p_step_pu for s in chosen},
    )
    model.shunt_q_step = pyo.Param(
        model.SHUNT_CTRL,
        within=pyo.Reals,
        initialize={s: data[s].q_step_pu for s in chosen},
    )
    model.shunt_switching_cost_coeff = pyo.Param(
        model.SHUNT_CTRL,
        within=pyo.NonNegativeReals,
        initialize=0.0,
        mutable=True,
    )
    sets = (model.SHUNT_CTRL,) if not multi else (model.SHUNT_CTRL, time_set)

    def _step_bounds(m, *idx):
        """Step bounds `(0, max_step)` of the shunt in `idx`.

        Args:
            m: The Pyomo model.
            *idx: `s` or `(s, tau)`.

        Returns:
            `(0, max_step)`.
        """
        return (0, m.shunt_step_max[idx[0]])

    def _step_start(m, *idx):
        """Starting value of a step variable: the network's `step`.

        Args:
            m: The Pyomo model.
            *idx: `s` or `(s, tau)`.

        Returns:
            The initial step.
        """
        return data[idx[0]].step_init

    def _move_bounds(m, *idx):
        """Bounds of an up or down move: at most the whole step range.

        Args:
            m: The Pyomo model.
            *idx: `s` or `(s, tau)`.

        Returns:
            `(0, max_step)`.
        """
        return (0.0, float(data[idx[0]].step_max))

    model.shunt_step = pyo.Var(
        *sets,
        domain=pyo.Integers if mode == "discrete" else pyo.Reals,
        bounds=_step_bounds,
        initialize=_step_start,
    )  # number of switched-in steps of the bank
    model.shunt_step_up = pyo.Var(
        *sets, domain=pyo.NonNegativeReals, bounds=_move_bounds, initialize=0.0
    )  # steps switched in since the previous period
    model.shunt_step_down = pyo.Var(
        *sets, domain=pyo.NonNegativeReals, bounds=_move_bounds, initialize=0.0
    )  # steps switched out since the previous period

    @model.Constraint(*sets)
    def shunt_step_movement_def(model, *idx):
        r"""Split the step change since the previous period into moves.

        $m_\tau - m_{\tau^-} = u_\tau - d_\tau$, $u, d \ge 0$, with the
        network's step before the first period. Skipped for the first
        period when the initial state is not tracked.

        Args:
            model: The Pyomo model being built.
            *idx: Shunt index, plus the time index on a multi-period model.

        Returns:
            A Pyomo equality expression, or `Constraint.Skip`.
        """
        s = idx[0]
        if not multi:
            if not track_first:
                return pyo.Constraint.Skip
            previous = model.shunt_step_init[s]
        else:
            tau = idx[1]
            if tau == time_set.first():
                if not track_first:
                    return pyo.Constraint.Skip
                previous = model.shunt_step_init[s]
            else:
                previous = model.shunt_step[s, time_set.prev(tau)]
        return (
            model.shunt_step[idx] - previous
            == model.shunt_step_up[idx] - model.shunt_step_down[idx]
        )

    if not track_first:
        for idx in model.shunt_step_up:
            if not multi or idx[1] == time_set.first():
                model.shunt_step_up[idx].fix(0.0)
                model.shunt_step_down[idx].fix(0.0)

    if change_max is not None:
        model.shunt_step_change_max = pyo.Param(
            model.SHUNT_CTRL,
            within=pyo.NonNegativeReals,
            initialize=change_max,
        )

        @model.Constraint(*sets)
        def shunt_step_change_limit(model, *idx):
            r"""Largest step change between consecutive periods.

            $u_\tau + d_\tau \le \Delta m^{\max}$.

            Args:
                model: The Pyomo model being built.
                *idx: Shunt index, plus the time index on a multi-period
                    model.

            Returns:
                A Pyomo inequality, or `Constraint.Skip`.
            """
            first = (not multi) or idx[1] == time_set.first()
            if first and not track_first:
                return pyo.Constraint.Skip
            return (
                model.shunt_step_up[idx] + model.shunt_step_down[idx]
                <= model.shunt_step_change_max[idx[0]]
            )

    if operations_max is not None:
        model.shunt_step_operations_max = pyo.Param(
            model.SHUNT_CTRL,
            within=pyo.NonNegativeReals,
            initialize=operations_max,
        )

        @model.Constraint(model.SHUNT_CTRL)
        def shunt_step_operations_limit(model, s):
            r"""Largest number of switching operations over the horizon.

            $\sum_\tau (u_\tau + d_\tau) \le N^{\max}$.

            Args:
                model: The Pyomo model being built.
                s: Shunt index from `model.SHUNT_CTRL`.

            Returns:
                A Pyomo inequality expression.
            """
            if not multi:
                moves = model.shunt_step_up[s] + model.shunt_step_down[s]
            else:
                moves = sum(
                    model.shunt_step_up[s, tau] + model.shunt_step_down[s, tau]
                    for tau in time_set
                )
            return moves <= model.shunt_step_operations_max[s]

    @model.Expression()
    def shunt_switching_cost(model):
        r"""Priced shunt switching, $\sum_s c_s \sum_\tau (u + d)$.

        Zero until `penalize_shunt_switching` sets the cost parameters.

        Args:
            model: The Pyomo model being built.

        Returns:
            A Pyomo expression.
        """
        return sum(
            model.shunt_switching_cost_coeff[_trafo_of(idx)]
            * (model.shunt_step_up[idx] + model.shunt_step_down[idx])
            for idx in model.shunt_step_up
        )

    # The balance was built with the constant GB/BB terms of these shunts;
    # rebuild it so the step variables take their place.
    model_obj.rebuild_kcl()

    setup = ShuntControlSetup(
        mode=mode,
        shunts=chosen,
        data=data,
        multi_period=multi,
        initial=initial,
        max_change_per_step=change_max,
        max_operations=operations_max,
    )
    model_obj.shunt_control_setup = setup
    logger.info(
        "enable_shunt_control: {} shunt(s) {} controllable in '{}' mode",
        len(chosen),
        chosen,
        mode,
    )
    return setup


def _require_shunt_control(model_obj) -> ShuntControlSetup:
    """The `ShuntControlSetup` of a model object, or a clear error.

    Args:
        model_obj: The model object.

    Returns:
        Its `shunt_control_setup`.

    Raises:
        RuntimeError: If `enable_shunt_control` has not been called.
    """
    setup = getattr(model_obj, "shunt_control_setup", None)
    if setup is None or not hasattr(model_obj.model, "SHUNT_CTRL"):
        raise RuntimeError("call enable_shunt_control() first")
    return setup


def penalize_shunt_switching(model_obj, cost):
    """Price shunt switching in the active objective.

    Sets `shunt_switching_cost_coeff` and, on the first call, adds
    `shunt_switching_cost` to the single active objective.

    Args:
        model_obj: A model object on which `enable_shunt_control` was called.
        cost: Cost per switched step, in the objective's unit; scalar or
            per-shunt mapping.

    Returns:
        The `shunt_switching_cost` expression.

    Raises:
        RuntimeError: If `enable_shunt_control` has not been called.
        ValueError: If there is no, or more than one, active objective.
    """
    setup = _require_shunt_control(model_obj)
    model = model_obj.model
    for s, c in _per_transformer(cost, setup.shunts, "cost").items():
        model.shunt_switching_cost_coeff[s] = c
    if setup.cost_objective is None:
        objectives = list(
            model.component_data_objects(pyo.Objective, active=True)
        )
        if not objectives:
            raise ValueError(
                "penalize_shunt_switching() needs an active objective; add "
                "one first"
            )
        if len(objectives) > 1:
            raise ValueError(
                "penalize_shunt_switching() found more than one active "
                "objective; deactivate all but one"
            )
        objective = objectives[0]
        sign = 1.0 if objective.sense == pyo.minimize else -1.0
        objective.set_value(objective.expr + sign * model.shunt_switching_cost)
        setup.cost_objective = objective.name
    return model.shunt_switching_cost


# ---------------------------------------------------------------------------
# reading the solution
# ---------------------------------------------------------------------------


def _step_value(model, key, discrete):
    """Solved step at `key`, rounded for a discrete bank.

    Args:
        model: The solved Pyomo model.
        key: `s` or `(s, tau)`.
        discrete: Whether to round to the nearest integer.

    Returns:
        An int (discrete) or float (continuous).
    """
    value = pyo.value(model.shunt_step[key])
    return int(round(value)) if discrete else float(value)


def shunt_schedule(model_obj):
    """Solved steps of the controlled shunts.

    Args:
        model_obj: A solved model object on which `enable_shunt_control`
            was called.

    Returns:
        Multi-period: a DataFrame (time step × shunt); single period: a
        Series indexed by shunt.
    """
    setup = _require_shunt_control(model_obj)
    model = model_obj.model
    discrete = setup.mode == "discrete"
    if not setup.multi_period:
        return pd.Series(
            {s: _step_value(model, s, discrete) for s in setup.shunts},
            name="step",
        )
    steps = list(model.T)
    frame = pd.DataFrame(
        {
            s: [_step_value(model, (s, tau), discrete) for tau in steps]
            for s in setup.shunts
        },
        index=pd.Index(steps, name="t"),
    )
    frame.columns.name = "shunt"
    return frame


def shunt_operations(model_obj):
    """Number of switched steps per controlled shunt, from the solution.

    Args:
        model_obj: A solved model object on which `enable_shunt_control`
            was called.

    Returns:
        A Series indexed by shunt.
    """
    setup = _require_shunt_control(model_obj)
    schedule = shunt_schedule(model_obj)
    counts = {}
    for s in setup.shunts:
        values = [schedule[s]] if not setup.multi_period else list(schedule[s])
        if setup.initial is not None:
            values = [setup.initial[s]] + values
        counts[s] = float(
            sum(abs(b - a) for a, b in zip(values[:-1], values[1:]))
        )
        if setup.mode == "discrete":
            counts[s] = int(round(counts[s]))
    return pd.Series(counts, name="shunt_operations")


def apply_shunt_steps(model_obj, net=None, t=None):
    """Write solved shunt steps into a network's `shunt.step` column.

    Args:
        model_obj: A solved model object on which `enable_shunt_control`
            was called.
        net: Network to write into (default: the model's own copy).
        t: Multi-period time step to apply (default: the last).

    Returns:
        The network written to.

    Warns:
        UserWarning: When a continuous step is rounded.
    """
    setup = _require_shunt_control(model_obj)
    schedule = shunt_schedule(model_obj)
    if setup.multi_period:
        step = model_obj.model.T.last() if t is None else int(t)
        if step not in schedule.index:
            raise ValueError(f"t={t} is not a time step of this model")
        values = schedule.loc[step]
    else:
        values = schedule
    target = model_obj.net if net is None else net
    rounded_any = False
    for s, m in values.items():
        if s not in target.shunt.index:
            raise KeyError(
                f"the target network has no shunt {s}; pass the network the "
                "model was built from (or a copy of it)"
            )
        m_int = int(round(float(m)))
        if abs(float(m) - m_int) > 1e-9:
            rounded_any = True
        target.shunt.at[s, "step"] = m_int
    if rounded_any:
        warnings.warn(
            "continuous shunt steps were rounded to the nearest integer "
            "before writing them to net.shunt.step",
            UserWarning,
            stacklevel=2,
        )
    return target


def shunt_result_columns(net, model, t=None):
    """Solved step per shunt of `net`, for `net.res_shunt["step"]`.

    Controlled shunts report their solved step (rounded when the variable
    is integer), the others the network's `step`.

    Args:
        net: The network the results are written to.
        model: The solved Pyomo model carrying `SHUNT_CTRL`.
        t: Time step for a multi-period model, `None` for single period.

    Returns:
        A Series indexed like `net.shunt`.
    """
    steps = pd.to_numeric(
        _column(net.shunt, "step", np.nan), errors="coerce"
    ).astype(float)
    for s in model.SHUNT_CTRL:
        key = s if t is None else (s, t)
        var = model.shunt_step[key]
        value = pyo.value(var)
        steps[s] = float(round(value)) if var.is_integer() else float(value)
    return steps


def shunt_rounding_group(model_obj) -> RoundingGroup:
    """The `RoundingGroup` of the shunt steps of `model_obj`.

    Args:
        model_obj: A model object on which `enable_shunt_control` was called.

    Returns:
        The group used by `solve_discrete_round_and_fix`.
    """
    setup = _require_shunt_control(model_obj)
    return RoundingGroup(
        name="shunt step",
        var=model_obj.model.shunt_step,
        units=list(setup.shunts),
        bounds={s: (0, d.step_max) for s, d in setup.data.items()},
        initial=setup.initial,
        change_limit=setup.max_change_per_step,
        operations_limit=setup.max_operations,
    )


# ---------------------------------------------------------------------------
# the mix-in that OPF and OPF_multi_period expose
# ---------------------------------------------------------------------------


class ShuntControlMixin:
    """Methods that make selected shunt steps OPF decision variables.

    Mixed into [`OPF`][potpourri.models.OPF.OPF] and `OPF_multi_period`;
    works on the AC formulations and raises `NotImplementedError` on DC
    and LPAC. Nothing changes until `enable_shunt_control` is called.
    """

    def enable_shunt_control(
        self,
        shunts=None,
        mode="discrete",
        *,
        max_change_per_step=None,
        max_operations=None,
        initial_step="net",
    ):
        """Make the steps of selected shunts decision variables.

        Adds `shunt_step` (integer in `"discrete"` mode, real in
        `"continuous"` mode, in `[0, max_step]`) and the up/down movement
        variables, and rebuilds the nodal balance so the controlled shunts
        consume `(p_mw + j q_mvar) * step * (v/vn)²` with the variable step.

        Args:
            shunts: `None` (default) controls every eligible shunt — see
                `shunt_eligibility(net)` — or an iterable of `net.shunt`
                indices.
            mode: `"discrete"` (physical bank; MINLP) or `"continuous"`
                (relaxation; NLP).
            max_change_per_step: Largest step change between consecutive
                time steps. Scalar or per-shunt mapping.
            max_operations: Largest total number of switched steps over
                the horizon. Scalar or mapping.
            initial_step: `"net"` (default) measures movement from
                `net.shunt.step`; a mapping gives other references; `None`
                leaves the first move untracked.

        Returns:
            The `ShuntControlSetup` describing the controlled shunts.

        Raises:
            NotImplementedError: On a DC or LPAC model.
            ValueError: If a requested shunt is not eligible.
            RuntimeError: If already enabled on this model.
        """
        return attach_shunt_control(
            self,
            shunts,
            mode,
            max_change_per_step=max_change_per_step,
            max_operations=max_operations,
            initial_step=initial_step,
        )

    def penalize_shunt_switching(self, cost):
        """Add a cost per switched shunt step to the active objective.

        Args:
            cost: Cost per step change, in the objective's unit; scalar or
                per-shunt mapping.

        Returns:
            The `shunt_switching_cost` expression.
        """
        return penalize_shunt_switching(self, cost)

    def shunt_schedule(self):
        """Solved shunt steps.

        Returns:
            A DataFrame (time step × shunt) on a multi-period model, a
            Series on a single-period one.
        """
        return shunt_schedule(self)

    def shunt_operations(self):
        """Number of switched steps per controlled shunt.

        Returns:
            A Series indexed by shunt.
        """
        return shunt_operations(self)

    def apply_shunt_steps(self, net=None, t=None):
        """Write the solved steps into `net.shunt.step`.

        Args:
            net: Network to write into (default: the model's own copy).
            t: Multi-period time step to apply (default: the last).

        Returns:
            The network written to.
        """
        return apply_shunt_steps(self, net=net, t=t)

    def solve_shunt_round_and_fix(self, solver="ipopt", **solve_kwargs):
        """Relax, round, fix and re-solve the discrete controls.

        The same heuristic as `solve_oltc_round_and_fix`; tap positions
        enabled on the same model are rounded in the same pass.

        Args:
            solver: NLP solver for both stages.
            **solve_kwargs: Forwarded to `solve()`.

        Returns:
            The solver results of the final solve.
        """
        from potpourri.models.oltc import solve_discrete_round_and_fix

        _require_shunt_control(self)
        return solve_discrete_round_and_fix(
            self, solver=solver, **solve_kwargs
        )
