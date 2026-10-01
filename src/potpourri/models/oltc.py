# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

r"""Opt-in on-load tap changer (OLTC) control for the polar AC formulations.

By default every transformer keeps the tap position pandapower built the
network with: `Tap` (HV-side ratio) and `Tap_lv` (LV-side ratio) are fixed
variables and the model is the same NLP as before. `enable_oltc` turns the
tap position of *selected* transformers into a decision variable that the
OPF chooses together with the dispatch, storage and reactive power.

## Physics

The transformer equations in [`AC`][potpourri.models.AC.AC] are the
MATPOWER/PowerModels branch model with a second ideal transformer on the LV
side: $Y_{ff} = (y_s + y_c/2)/a_{hv}^2$, $Y_{tt} = (y_s + y_c/2)/a_{lv}^2$,
$Y_{ft} = -y_s e^{j\varphi}/(a_{hv} a_{lv})$. pandapower's longitudinal
("Ratio") tap changer scales a rated winding voltage by the tap factor

$$
n(k) = 1 + (k - k_{\text{neutral}})\,\frac{\texttt{tap\_step\_percent}}{100},
$$

and refers the short-circuit impedance to the *tapped* LV voltage, so for an
LV-side tap every per-unit admittance of the branch scales with $1/n^2$. In
the two-sided model that is simply the LV-side ratio: with $\tau_0$ the ratio
and $n_0$ the tap factor pandapower built the network with,

| `tap_side` | $a_{hv}$ | $a_{lv}$ |
|---|---|---|
| `"hv"` | $\tau_0\, n(k)/n_0 = r_0\, n(k)$ | $1$ |
| `"lv"` | $\tau_0$ | $n(k)/n_0$ |

where $r_0$ is the nominal mismatch, the rated-voltage ratio
$v_{n,hv}/v_{n,lv}$ over the bus-voltage ratio $V^{bus}_{hv}/V^{bus}_{lv}$.
Both rows reduce to the stored network at $k = k_0$, and both reproduce
`pp.runpp` at every position (`tests/unit_tests/test_oltc.py`). The
derivation and the literature are in
`docs/research/dso_controllable_equipment.md`.

## Modes

* `mode="continuous"`: `trafo_tap_position` is real in
  $[k_{\min}, k_{\max}]$. The model stays an NLP; the result is a
  **relaxation**, not an implementable position.
* `mode="discrete"`: the position is an integer in the same bounds. The
  model is a nonconvex MINLP (`gurobi_direct_minlp`, MindtPy, NEOS
  Bonmin/Couenne).

`solve_oltc_round_and_fix` solves the relaxation, rounds, fixes and re-solves
with any NLP solver, which gives an integer schedule without a MINLP solver.
`solve()` refuses to hand free integer positions to IPOPT.

## Multi-period scheduling

The position of the previous step (the network's `tap_pos` before the first
one) is linked through non-negative up/down moves, which carry an optional
per-step change limit, an optional limit on the number of operations over
the horizon, and an optional cost that `penalize_tap_movement` adds to the
objective.

## Which transformers

`oltc_eligibility(net)` reports, per transformer, whether `enable_oltc` can
control it and why not otherwise. Only in-service two-winding transformers
with `tap_changer_type == "Ratio"`, no `tap_step_degree`, no
`tap_dependency_table`, no second tap changer and complete integer tap data
qualify. A transformer whose `tap_changer_type` is `None` — every SimBench
network as delivered — is *not* eligible, because pandapower itself ignores
its `tap_pos`; set the type to `"Ratio"` first. The pandapower `oltc`
column is a short-circuit flag and plays no part.

Example:
    ```python
    opf = ACOPF(net)
    opf.add_OPF()
    opf.enable_oltc(mode="discrete")            # every eligible transformer
    opf.add_voltage_deviation_objective()
    opf.penalize_tap_movement(cost=1e-4)        # optional
    opf.solve(solver="gurobi_direct_minlp")     # or solve_oltc_round_and_fix()
    opf.net.res_trafo[["tap_pos", "tap_factor", "loading_percent"]]
    ```

See Also:
    - [`OLTCControlMixin`][potpourri.models.oltc.OLTCControlMixin]: the
      methods `OPF` and `OPF_multi_period` expose.
    - [`shunt_control`][potpourri.models.shunt_control]: the same pattern
      for switched capacitor banks and reactors.
"""

from __future__ import annotations

import inspect
import math
import warnings
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import pyomo.environ as pyo
from loguru import logger

#: Accepted values of the `mode` argument of `enable_oltc`.
OLTC_MODES = ("continuous", "discrete")

#: Relative tolerance of the self-consistency check between the ratio
#: pandapower built into the ppc and the one the "Ratio" formula predicts.
_RATIO_TOL = 1e-9


@dataclass
class TapChangerData:
    r"""Tap-changer parameters of one controllable transformer.

    Positions are integers in pandapower's own numbering; `step` is the
    tap step as a fraction (`tap_step_percent / 100`), `ratio_nominal` the
    nominal mismatch $r_0$, `tap_ppc` the ratio pandapower built into the
    ppc and `factor_base` the tap factor $n_0$ that ratio corresponds to.

    Attributes:
        index: Transformer index in `net.trafo` and `model.TRANSF`.
        side: `"hv"` or `"lv"`, the tapped winding.
        pos_min: Lowest admissible position.
        pos_max: Highest admissible position.
        pos_neutral: Position at which the ratio equals the rated one.
        pos_init: Position the network was built with (`tap_pos`).
        step: Tap step as a fraction of the rated voltage per position.
        ratio_nominal: $r_0$, the rated-voltage ratio over the bus-voltage
            ratio; 1 when rated and bus voltages agree.
        tap_ppc: $\tau_0$, the ppc `TAP` value at `pos_init`.
        factor_base: $n_0 = n(\text{pos\_init})$ recovered from the ppc.
    """

    index: int
    side: str
    pos_min: int
    pos_max: int
    pos_neutral: int
    pos_init: int
    step: float
    ratio_nominal: float
    tap_ppc: float
    factor_base: float

    def factor(self, position):
        r"""Tap factor $n(k)$ at a (possibly fractional) position.

        Args:
            position: Tap position $k$.

        Returns:
            $1 + (k - k_{\text{neutral}})\,s$.
        """
        return 1.0 + (position - self.pos_neutral) * self.step

    @property
    def factor_bounds(self):
        """Lower and upper bound of the tap factor over the position range.

        Sorted, so a negative `tap_step_percent` is handled.

        Returns:
            `(n_lo, n_hi)`.
        """
        values = (self.factor(self.pos_min), self.factor(self.pos_max))
        return min(values), max(values)


@dataclass
class OLTCSetup:
    """What `enable_oltc` built, kept on the model object as `oltc_setup`.

    Attributes:
        mode: `"continuous"` or `"discrete"`.
        transformers: Controlled transformer indices, in order.
        data: `TapChangerData` per controlled transformer.
        multi_period: Whether the positions carry a time index.
        initial: Reference position per transformer for the first move, or
            `None` when the move away from the initial state is not tracked.
        max_change_per_step: Per-step change limit per transformer, if any.
        max_operations: Operation limit per transformer, if any.
        cost_objective: Name of the objective the switching cost was added
            to, or `None` until `penalize_tap_movement` has been called.
    """

    mode: str
    transformers: list
    data: dict
    multi_period: bool
    initial: dict | None
    max_change_per_step: dict | None = None
    max_operations: dict | None = None
    cost_objective: str | None = None
    extra: dict = field(default_factory=dict)


# ---------------------------------------------------------------------------
# eligibility
# ---------------------------------------------------------------------------


def _column(table, name, default):
    """Column `name` of `table`, or a Series of `default` when absent.

    Args:
        table: A pandapower element table.
        name: Column name.
        default: Fill value for a missing column.

    Returns:
        A pandas Series indexed like `table`.
    """
    if name in table.columns:
        return table[name]
    return pd.Series([default] * len(table), index=table.index, dtype=object)


def _is_missing(value):
    """Whether a table cell is empty (`None`, NaN or pandas NA).

    Args:
        value: Any cell value.

    Returns:
        True for a missing value.
    """
    if value is None:
        return True
    try:
        return bool(pd.isna(value))
    except (TypeError, ValueError):
        return False


def _finite_integer(value):
    """Whether a cell holds a finite integral number.

    Args:
        value: Any cell value.

    Returns:
        True for e.g. `3`, `3.0`, `-9`; False for NaN, None, `2.5` or text.
    """
    if _is_missing(value):
        return False
    try:
        number = float(value)
    except (TypeError, ValueError):
        return False
    return math.isfinite(number) and float(number).is_integer()


def _ineligibility_reason(row, tap_type, dep_table, tap2_pos, step_degree):
    """The first reason a transformer row cannot be OLTC-controlled.

    Implements rules 1–7 of the eligibility list in
    `docs/research/dso_controllable_equipment.md` (rule 8, the ppc
    consistency check, needs the built model and lives in `_tap_data`).

    Args:
        row: A row of `net.trafo`.
        tap_type: Its `tap_changer_type` cell.
        dep_table: Its `tap_dependency_table` cell.
        tap2_pos: Its `tap2_pos` cell.
        step_degree: Its `tap_step_degree` cell.

    Returns:
        A short explanation, or `""` when the transformer is eligible.
    """
    if not bool(row["in_service"]):
        return "out of service"
    side = row.get("tap_side")
    if side not in ("hv", "lv"):
        return f"tap_side is {side!r}; expected 'hv' or 'lv'"
    if _is_missing(tap_type):
        return (
            "tap_changer_type is None: pandapower treats the transformer as "
            "having no tap changer and ignores tap_pos; set "
            "net.trafo.tap_changer_type = 'Ratio' to make the tap effective"
        )
    if tap_type != "Ratio":
        return (
            f"tap_changer_type {tap_type!r} is not supported (only the "
            "longitudinal 'Ratio' changer is)"
        )
    if not _is_missing(step_degree) and float(step_degree) != 0.0:
        return (
            "tap_step_degree != 0 makes this a cross regulator with a "
            "tap-dependent phase shift, which is not supported"
        )
    if not _is_missing(dep_table) and bool(dep_table):
        return (
            "tap_dependency_table is True (tabular characteristic with "
            "tap-dependent ratio, angle or impedance), which is not supported"
        )
    if not _is_missing(tap2_pos):
        return (
            "a second tap changer (tap2_*) is present, which is not supported"
        )
    step = row.get("tap_step_percent")
    if _is_missing(step) or not math.isfinite(float(step)):
        return "tap_step_percent is missing"
    if float(step) == 0.0:
        return "tap_step_percent is zero"
    for name in ("tap_neutral", "tap_min", "tap_max"):
        if not _finite_integer(row.get(name)):
            return f"{name} is missing or not an integer"
    if int(row["tap_min"]) >= int(row["tap_max"]):
        return "tap_min must be smaller than tap_max"
    if not _finite_integer(row.get("tap_pos")):
        return "tap_pos is missing or not an integer position"
    if not int(row["tap_min"]) <= int(row["tap_pos"]) <= int(row["tap_max"]):
        return "tap_pos lies outside [tap_min, tap_max]"
    return ""


def oltc_eligibility(net) -> pd.DataFrame:
    """Which transformers `enable_oltc` can control, and why not otherwise.

    A transformer is eligible when it is in service, its tap changer is a
    longitudinal `"Ratio"` changer on the `"hv"` or `"lv"` side without
    `tap_step_degree`, without `tap_dependency_table` and without a second
    changer, and `tap_neutral`, `tap_min < tap_max`, `tap_step_percent != 0`
    and `tap_min <= tap_pos <= tap_max` are all present as integers. The
    pandapower `oltc` column (a short-circuit flag) and any controllers in
    `net.controller` are deliberately ignored.

    Args:
        net: A pandapower network (the caller's or a model's copy).

    Returns:
        A DataFrame indexed like `net.trafo` with the columns `eligible`
        (bool), `reason` (empty when eligible), `tap_side`, `tap_pos`,
        `tap_neutral`, `tap_min`, `tap_max`, `tap_step_percent` and
        `n_positions`.

    Examples:
        The four-bus example network carries a 20/0.4 kV unit with a
        complete ±2 × 2.5 % "Ratio" tap changer, so it qualifies; dropping
        the type is enough to disqualify it, because pandapower would then
        ignore its tap position.

        >>> import pandapower as pp
        >>> from potpourri.models.oltc import oltc_eligibility
        >>> net = pp.networks.simple_four_bus_system()
        >>> bool(oltc_eligibility(net)["eligible"].iloc[0])
        True
        >>> net.trafo["tap_changer_type"] = None
        >>> report = oltc_eligibility(net)
        >>> bool(report["eligible"].iloc[0]), report["reason"].iloc[0][:30]
        (False, 'tap_changer_type is None: pand')
    """
    trafo = net.trafo
    if "tap_changer_type" in trafo.columns:
        tap_type = trafo["tap_changer_type"]
    elif "tap_phase_shifter" in trafo.columns:
        # pandapower < 3.0 data: a non-phase-shifting changer with tap data
        # is what 3.x calls "Ratio" (its legacy code path applies it the
        # same way); the ppc consistency check in `_tap_data` guards the
        # inference.
        tap_type = pd.Series(
            [
                None
                if _is_missing(step)
                else ("Ideal" if bool(shifter) else "Ratio")
                for shifter, step in zip(
                    trafo["tap_phase_shifter"].fillna(False),
                    _column(trafo, "tap_step_percent", np.nan),
                )
            ],
            index=trafo.index,
            dtype=object,
        )
    else:
        tap_type = _column(trafo, "tap_changer_type", None)
    dep_table = _column(trafo, "tap_dependency_table", False)
    tap2_pos = _column(trafo, "tap2_pos", np.nan)
    step_degree = _column(trafo, "tap_step_degree", np.nan)
    records = []
    for idx, row in trafo.iterrows():
        reason = _ineligibility_reason(
            row, tap_type[idx], dep_table[idx], tap2_pos[idx], step_degree[idx]
        )
        n_positions = (
            int(row["tap_max"]) - int(row["tap_min"]) + 1
            if _finite_integer(row.get("tap_min"))
            and _finite_integer(row.get("tap_max"))
            else np.nan
        )
        records.append(
            {
                "eligible": reason == "",
                "reason": reason,
                "tap_side": row.get("tap_side"),
                "tap_pos": row.get("tap_pos"),
                "tap_neutral": row.get("tap_neutral"),
                "tap_min": row.get("tap_min"),
                "tap_max": row.get("tap_max"),
                "tap_step_percent": row.get("tap_step_percent"),
                "n_positions": n_positions,
            }
        )
    report = pd.DataFrame(records, index=trafo.index)
    report.index.name = "trafo"
    return report


def _tap_data(model_obj, index) -> TapChangerData:
    """Collect the tap-changer data of one eligible transformer.

    Also performs eligibility rule 8: the ratio pandapower built into the
    ppc must equal the one the "Ratio" formula predicts for `tap_pos`.

    Args:
        model_obj: The model object (its `net` and `trafo_data`).
        index: Transformer index.

    Returns:
        The collected `TapChangerData`.

    Raises:
        ValueError: If the ppc ratio does not match the formula.
    """
    row = model_obj.net.trafo.loc[index]
    bus = model_obj.net.bus
    side = str(row["tap_side"])
    step = float(row["tap_step_percent"]) / 100.0
    pos_neutral = int(row["tap_neutral"])
    pos_init = int(row["tap_pos"])
    ratio_nominal = (float(row["vn_hv_kv"]) / float(row["vn_lv_kv"])) / (
        float(bus.at[int(row["hv_bus"]), "vn_kv"])
        / float(bus.at[int(row["lv_bus"]), "vn_kv"])
    )
    tap_ppc = float(model_obj.trafo_data.at[index, "tap"])
    factor_formula = 1.0 + (pos_init - pos_neutral) * step
    if side == "hv":
        factor_base = tap_ppc / ratio_nominal
    else:
        factor_base = ratio_nominal / tap_ppc
    if abs(factor_base - factor_formula) > _RATIO_TOL * max(
        1.0, abs(factor_formula)
    ):
        raise ValueError(
            f"transformer {index}: pandapower built the ratio {tap_ppc:.9g} "
            f"(tap factor {factor_base:.9g}) but the 'Ratio' tap-changer "
            f"formula gives {factor_formula:.9g} at tap_pos={pos_init}; the "
            "configuration is not one this model describes"
        )
    return TapChangerData(
        index=int(index),
        side=side,
        pos_min=int(row["tap_min"]),
        pos_max=int(row["tap_max"]),
        pos_neutral=pos_neutral,
        pos_init=pos_init,
        step=step,
        ratio_nominal=ratio_nominal,
        tap_ppc=tap_ppc,
        factor_base=factor_formula,
    )


# ---------------------------------------------------------------------------
# attaching the control to a model
# ---------------------------------------------------------------------------


def _per_transformer(value, transformers, what):
    """Broadcast a scalar or mapping to a `{transformer: value}` dict.

    Args:
        value: A number, or a mapping / Series keyed by transformer index.
        transformers: The controlled transformer indices.
        what: Name of the quantity, for error messages.

    Returns:
        A dict with one entry per controlled transformer.

    Raises:
        ValueError: If a mapping misses a controlled transformer or holds a
            negative value.
    """
    if isinstance(value, (pd.Series, dict)):
        out = {}
        for t in transformers:
            if t not in value:
                raise ValueError(f"{what} has no entry for transformer {t}")
            out[t] = float(value[t])
    else:
        out = {t: float(value) for t in transformers}
    for t, v in out.items():
        if not math.isfinite(v) or v < 0:
            raise ValueError(
                f"{what} for transformer {t} must be a finite non-negative "
                f"number, got {v!r}"
            )
    return out


def _trafo_of(idx):
    """The transformer part of a (possibly time-indexed) variable index.

    Args:
        idx: `t` or `(t, tau)`.

    Returns:
        `t`.
    """
    return idx[0] if isinstance(idx, tuple) else idx


def _select_transformers(model_obj, transformers):
    """Resolve the `transformers` argument of `enable_oltc`.

    Args:
        model_obj: The model object.
        transformers: `None` for every eligible transformer, or an iterable
            of `net.trafo` indices.

    Returns:
        The chosen indices, in order and without duplicates.

    Raises:
        ValueError: If nothing is eligible, or an explicitly requested
            transformer is not.
    """
    report = oltc_eligibility(model_obj.net)
    in_model = set(int(t) for t in model_obj.model.TRANSF)
    if transformers is None:
        chosen = [
            int(t)
            for t in report.index
            if int(t) in in_model and bool(report.at[t, "eligible"])
        ]
        if not chosen:
            reasons = "\n".join(
                f"  trafo {t}: {report.at[t, 'reason']}" for t in report.index
            )
            raise ValueError(
                "No transformer is eligible for OLTC control:\n" + reasons
            )
        skipped = [t for t in report.index if int(t) not in chosen]
        if skipped:
            logger.info(
                "enable_oltc: skipping transformer(s) {} -- see "
                "oltc_eligibility(net) for the reasons",
                skipped,
            )
        return chosen
    chosen = list(dict.fromkeys(int(t) for t in transformers))
    if not chosen:
        raise ValueError("transformers must name at least one transformer")
    problems = []
    for t in chosen:
        if t not in report.index:
            problems.append(f"trafo {t}: not a transformer of the network")
        elif t not in in_model:
            problems.append(
                f"trafo {t}: not in the model (out of service, or removed "
                "while fusing bus-bus switches)"
            )
        elif not bool(report.at[t, "eligible"]):
            problems.append(f"trafo {t}: {report.at[t, 'reason']}")
    if problems:
        raise ValueError(
            "Transformer(s) cannot be OLTC-controlled:\n"
            + "\n".join("  " + p for p in problems)
        )
    return chosen


def attach_oltc(
    model_obj,
    transformers=None,
    mode="discrete",
    *,
    max_change_per_step=None,
    max_operations=None,
    initial_tap="net",
):
    """Make the tap positions of selected transformers decision variables.

    The implementation behind `OLTCControlMixin.enable_oltc`; see there
    for the user-facing description. Adds the components listed in
    `docs/research/dso_controllable_equipment.md` (Section 6.2) to
    `model_obj.model`, unfixes `Tap` for HV-side and `Tap_lv` for LV-side
    changers, and records an `OLTCSetup` as `model_obj.oltc_setup`.

    Args:
        model_obj: An `OPF`-derived single- or multi-period model object
            whose class sets `OLTC_SUPPORTED`.
        transformers: `None` (every eligible transformer) or indices.
        mode: `"continuous"` or `"discrete"`.
        max_change_per_step: Largest position change between consecutive
            steps (single period: away from the initial position); scalar
            or per-transformer mapping; `None` for no limit.
        max_operations: Largest total number of position changes over the
            horizon (single period: away from the initial position);
            scalar or mapping; `None` for no limit.
        initial_tap: `"net"` to measure movement from `net.trafo.tap_pos`,
            a mapping of reference positions, or `None` to leave the first
            move untracked.

    Returns:
        The `OLTCSetup` describing what was built.

    Raises:
        NotImplementedError: On a formulation without a variable tap (DC,
            LPAC).
        ValueError: On a bad `mode`, an ineligible transformer or
            inconsistent limits.
        RuntimeError: If `enable_oltc` was already called on the model.
    """
    if not getattr(model_obj, "OLTC_SUPPORTED", False):
        raise NotImplementedError(
            f"{type(model_obj).__name__} has no variable transformer tap: "
            "OLTC control needs the polar AC equations (ACOPF, "
            "ACOPF_multi_period). The DC model has no voltage magnitude and "
            "the LPAC model is linear only for a fixed tap."
        )
    if mode not in OLTC_MODES:
        raise ValueError(f"mode must be one of {OLTC_MODES}, got {mode!r}")
    model = model_obj.model
    if hasattr(model, "TRANSF_OLTC"):
        raise RuntimeError(
            "enable_oltc() has already been called on this model"
        )

    time_set = getattr(model, "T", None)
    multi = time_set is not None
    chosen = _select_transformers(model_obj, transformers)
    data = {t: _tap_data(model_obj, t) for t in chosen}

    if initial_tap is None:
        initial = None
    elif isinstance(initial_tap, str) and initial_tap == "net":
        initial = {t: data[t].pos_init for t in chosen}
    else:
        initial = {}
        for t in chosen:
            if t not in initial_tap:
                raise ValueError(
                    f"initial_tap has no entry for transformer {t}"
                )
            k = initial_tap[t]
            if not _finite_integer(k):
                raise ValueError(
                    f"initial_tap for transformer {t} must be an integer "
                    f"position, got {k!r}"
                )
            initial[t] = int(k)
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

    # --- sets and parameters -------------------------------------------
    model.TRANSF_OLTC = pyo.Set(within=model.TRANSF, initialize=chosen)
    model.TRANSF_OLTC_HV = pyo.Set(
        within=model.TRANSF_OLTC,
        initialize=[t for t in chosen if data[t].side == "hv"],
    )
    model.TRANSF_OLTC_LV = pyo.Set(
        within=model.TRANSF_OLTC,
        initialize=[t for t in chosen if data[t].side == "lv"],
    )
    model.trafo_tap_pos_min = pyo.Param(
        model.TRANSF_OLTC,
        within=pyo.Integers,
        initialize={t: data[t].pos_min for t in chosen},
    )
    model.trafo_tap_pos_max = pyo.Param(
        model.TRANSF_OLTC,
        within=pyo.Integers,
        initialize={t: data[t].pos_max for t in chosen},
    )
    model.trafo_tap_pos_neutral = pyo.Param(
        model.TRANSF_OLTC,
        within=pyo.Integers,
        initialize={t: data[t].pos_neutral for t in chosen},
    )
    model.trafo_tap_pos_init = pyo.Param(
        model.TRANSF_OLTC,
        within=pyo.Integers,
        initialize={
            t: (initial[t] if initial is not None else data[t].pos_init)
            for t in chosen
        },
    )
    model.trafo_tap_step = pyo.Param(
        model.TRANSF_OLTC,
        within=pyo.Reals,
        initialize={t: data[t].step for t in chosen},
    )
    model.trafo_tap_factor_base = pyo.Param(
        model.TRANSF_OLTC,
        within=pyo.PositiveReals,
        initialize={t: data[t].factor_base for t in chosen},
    )
    model.trafo_tap_ratio_nominal = pyo.Param(
        model.TRANSF_OLTC,
        within=pyo.PositiveReals,
        initialize={t: data[t].ratio_nominal for t in chosen},
    )
    model.trafo_tap_switching_cost = pyo.Param(
        model.TRANSF_OLTC,
        within=pyo.NonNegativeReals,
        initialize=0.0,
        mutable=True,
    )

    sets = (model.TRANSF_OLTC,) if not multi else (model.TRANSF_OLTC, time_set)
    hv_sets = (
        (model.TRANSF_OLTC_HV,)
        if not multi
        else (model.TRANSF_OLTC_HV, time_set)
    )
    lv_sets = (
        (model.TRANSF_OLTC_LV,)
        if not multi
        else (model.TRANSF_OLTC_LV, time_set)
    )

    # --- variables ------------------------------------------------------
    def _position_bounds(m, *idx):
        """Position bounds of the transformer in `idx`.

        Args:
            m: The Pyomo model.
            *idx: `t` or `(t, tau)`.

        Returns:
            `(tap_min, tap_max)`.
        """
        t = idx[0]
        return (m.trafo_tap_pos_min[t], m.trafo_tap_pos_max[t])

    def _position_start(m, *idx):
        """Starting value of a position variable: the network's `tap_pos`.

        Args:
            m: The Pyomo model.
            *idx: `t` or `(t, tau)`.

        Returns:
            The initial position.
        """
        return data[idx[0]].pos_init

    def _factor_bounds(m, *idx):
        """Tap-factor bounds of the transformer in `idx`.

        Args:
            m: The Pyomo model.
            *idx: `t` or `(t, tau)`.

        Returns:
            `(n_lo, n_hi)`.
        """
        return data[idx[0]].factor_bounds

    def _factor_start(m, *idx):
        """Starting value of a tap-factor variable.

        Args:
            m: The Pyomo model.
            *idx: `t` or `(t, tau)`.

        Returns:
            $n_0$.
        """
        return data[idx[0]].factor_base

    def _move_bounds(m, *idx):
        """Bounds of an up or down move: at most the whole position range.

        Args:
            m: The Pyomo model.
            *idx: `t` or `(t, tau)`.

        Returns:
            `(0, tap_max - tap_min)`.
        """
        d = data[idx[0]]
        return (0.0, float(d.pos_max - d.pos_min))

    model.trafo_tap_position = pyo.Var(
        *sets,
        domain=pyo.Integers if mode == "discrete" else pyo.Reals,
        bounds=_position_bounds,
        initialize=_position_start,
    )  # tap position k, in pandapower's numbering
    model.trafo_tap_factor = pyo.Var(
        *sets,
        domain=pyo.PositiveReals,
        bounds=_factor_bounds,
        initialize=_factor_start,
    )  # tap factor n(k) = 1 + (k - k_neutral) * step
    model.trafo_tap_up = pyo.Var(
        *sets,
        domain=pyo.NonNegativeReals,
        bounds=_move_bounds,
        initialize=0.0,
    )  # upward position change since the previous step
    model.trafo_tap_down = pyo.Var(
        *sets,
        domain=pyo.NonNegativeReals,
        bounds=_move_bounds,
        initialize=0.0,
    )  # downward position change since the previous step

    # Free the ratio that the tapped winding sits on, with bounds that
    # follow from the position range; the other ratio stays fixed.
    for t in chosen:
        d = data[t]
        n_lo, n_hi = d.factor_bounds
        keys = [t] if not multi else [(t, tau) for tau in time_set]
        for key in keys:
            if d.side == "hv":
                var = model.Tap[key]
                lo, hi = sorted(
                    (d.ratio_nominal * n_lo, d.ratio_nominal * n_hi)
                )
            else:
                var = model.Tap_lv[key]
                lo, hi = sorted((n_lo / d.factor_base, n_hi / d.factor_base))
            var.unfix()
            var.setlb(lo)
            var.setub(hi)

    # --- constraints ----------------------------------------------------
    @model.Constraint(*sets)
    def trafo_tap_factor_def(model, *idx):
        r"""Tap factor of a controlled transformer at its position.

        $n = 1 + (k - k_{\text{neutral}})\,s$ with $s$ =
        `tap_step_percent / 100`: pandapower's longitudinal "Ratio"
        tap changer. Affine in $k$, so the relaxation of an integer $k$
        is exactly the continuous mode.

        Args:
            model: The Pyomo model being built.
            *idx: Transformer index, plus the time index on a multi-period
                model.

        Returns:
            A Pyomo equality expression defining `trafo_tap_factor`.
        """
        t = idx[0]
        return (
            model.trafo_tap_factor[idx]
            == 1.0
            + (model.trafo_tap_position[idx] - model.trafo_tap_pos_neutral[t])
            * model.trafo_tap_step[t]
        )

    @model.Constraint(*hv_sets)
    def trafo_tap_ratio_hv_def(model, *idx):
        r"""HV-side ratio of a transformer tapped on its HV winding.

        $a_{hv} = r_0\, n(k)$: pandapower scales the rated HV voltage by
        the tap factor and the admittances do not depend on the tap, so
        the whole effect is in the from-side ratio `Tap`.

        Args:
            model: The Pyomo model being built.
            *idx: Transformer index, plus the time index on a multi-period
                model.

        Returns:
            A Pyomo equality expression tying `Tap` to the tap factor.
        """
        t = idx[0]
        return (
            model.Tap[idx]
            == model.trafo_tap_ratio_nominal[t] * model.trafo_tap_factor[idx]
        )

    @model.Constraint(*lv_sets)
    def trafo_tap_ratio_lv_def(model, *idx):
        r"""LV-side ratio of a transformer tapped on its LV winding.

        $a_{lv} = n(k)/n_0$: pandapower refers the impedance to the tapped
        LV voltage, so every admittance scales with $1/n^2$. Written with
        the admittances the network was built with (at $n_0$), that is a
        to-side ideal transformer of ratio $n/n_0$ while `Tap` keeps the
        built ratio $\tau_0$.

        Args:
            model: The Pyomo model being built.
            *idx: Transformer index, plus the time index on a multi-period
                model.

        Returns:
            A Pyomo equality expression tying `Tap_lv` to the tap factor.
        """
        t = idx[0]
        return (
            model.Tap_lv[idx]
            == model.trafo_tap_factor[idx] / model.trafo_tap_factor_base[t]
        )

    @model.Constraint(*sets)
    def trafo_tap_movement_def(model, *idx):
        r"""Split the position change since the previous step into moves.

        $k_\tau - k_{\tau^-} = u_\tau - d_\tau$ with $u, d \ge 0$; the
        previous position of the first step (or of the single period) is
        `trafo_tap_pos_init`. Skipped for the first step when the initial
        state is not tracked (`initial_tap=None`).

        Args:
            model: The Pyomo model being built.
            *idx: Transformer index, plus the time index on a multi-period
                model.

        Returns:
            A Pyomo equality expression, or `Constraint.Skip`.
        """
        t = idx[0]
        if not multi:
            if not track_first:
                return pyo.Constraint.Skip
            previous = model.trafo_tap_pos_init[t]
        else:
            tau = idx[1]
            if tau == time_set.first():
                if not track_first:
                    return pyo.Constraint.Skip
                previous = model.trafo_tap_pos_init[t]
            else:
                previous = model.trafo_tap_position[t, time_set.prev(tau)]
        return (
            model.trafo_tap_position[idx] - previous
            == model.trafo_tap_up[idx] - model.trafo_tap_down[idx]
        )

    if not track_first:
        # Nothing to move away from: the first step's moves are zero.
        for idx in model.trafo_tap_up:
            if not multi or idx[1] == time_set.first():
                model.trafo_tap_up[idx].fix(0.0)
                model.trafo_tap_down[idx].fix(0.0)

    if change_max is not None:
        model.trafo_tap_change_max = pyo.Param(
            model.TRANSF_OLTC,
            within=pyo.NonNegativeReals,
            initialize=change_max,
        )

        @model.Constraint(*sets)
        def trafo_tap_change_limit(model, *idx):
            r"""Largest position change between consecutive steps.

            $u_\tau + d_\tau \le \Delta k^{\max}$. Skipped where the move
            itself is not tracked.

            Args:
                model: The Pyomo model being built.
                *idx: Transformer index, plus the time index on a
                    multi-period model.

            Returns:
                A Pyomo inequality, or `Constraint.Skip`.
            """
            t = idx[0]
            first = (not multi) or idx[1] == time_set.first()
            if first and not track_first:
                return pyo.Constraint.Skip
            return (
                model.trafo_tap_up[idx] + model.trafo_tap_down[idx]
                <= model.trafo_tap_change_max[t]
            )

    if operations_max is not None:
        model.trafo_tap_operations_max = pyo.Param(
            model.TRANSF_OLTC,
            within=pyo.NonNegativeReals,
            initialize=operations_max,
        )

        @model.Constraint(model.TRANSF_OLTC)
        def trafo_tap_operations_limit(model, t):
            r"""Largest number of tap operations over the horizon.

            $\sum_\tau (u_\tau + d_\tau) \le N^{\max}$; on a single-period
            model the one move away from the initial position.

            Args:
                model: The Pyomo model being built.
                t: Transformer index from `model.TRANSF_OLTC`.

            Returns:
                A Pyomo inequality expression.
            """
            if not multi:
                moves = model.trafo_tap_up[t] + model.trafo_tap_down[t]
            else:
                moves = sum(
                    model.trafo_tap_up[t, tau] + model.trafo_tap_down[t, tau]
                    for tau in time_set
                )
            return moves <= model.trafo_tap_operations_max[t]

    @model.Expression()
    def trafo_tap_movement_cost(model):
        r"""Priced tap movement, $\sum_t c_t \sum_\tau (u + d)$.

        Zero until `penalize_tap_movement` sets the mutable cost
        parameters and adds this expression to the objective.

        Args:
            model: The Pyomo model being built.

        Returns:
            A Pyomo expression.
        """
        return sum(
            model.trafo_tap_switching_cost[_trafo_of(idx)]
            * (model.trafo_tap_up[idx] + model.trafo_tap_down[idx])
            for idx in model.trafo_tap_up
        )

    setup = OLTCSetup(
        mode=mode,
        transformers=chosen,
        data=data,
        multi_period=multi,
        initial=initial,
        max_change_per_step=change_max,
        max_operations=operations_max,
    )
    model_obj.oltc_setup = setup
    logger.info(
        "enable_oltc: {} tap changer(s) {} controllable in '{}' mode ({})",
        len(chosen),
        chosen,
        mode,
        "multi-period" if multi else "single period",
    )
    return setup


def _require_oltc(model_obj) -> OLTCSetup:
    """The `OLTCSetup` of a model object, or a clear error.

    Args:
        model_obj: The model object.

    Returns:
        Its `oltc_setup`.

    Raises:
        RuntimeError: If `enable_oltc` has not been called.
    """
    setup = getattr(model_obj, "oltc_setup", None)
    if setup is None or not hasattr(model_obj.model, "TRANSF_OLTC"):
        raise RuntimeError("call enable_oltc() first")
    return setup


def penalize_tap_movement(model_obj, cost):
    """Price tap movement in the active objective.

    Sets `trafo_tap_switching_cost` and, on the first call, adds
    `trafo_tap_movement_cost` to the model's single active objective (added
    for a minimisation, subtracted for a maximisation). Later calls only
    update the cost, so a model can be re-solved at several prices without
    rebuilding.

    Args:
        model_obj: A model object on which `enable_oltc` was called.
        cost: Cost per tap operation, in the objective's unit; a scalar or
            a per-transformer mapping.

    Returns:
        The `trafo_tap_movement_cost` expression.

    Raises:
        RuntimeError: If `enable_oltc` has not been called.
        ValueError: If there is no, or more than one, active objective.
    """
    setup = _require_oltc(model_obj)
    model = model_obj.model
    costs = _per_transformer(cost, setup.transformers, "cost")
    for t, c in costs.items():
        model.trafo_tap_switching_cost[t] = c
    if setup.cost_objective is None:
        objectives = list(
            model.component_data_objects(pyo.Objective, active=True)
        )
        if not objectives:
            raise ValueError(
                "penalize_tap_movement() needs an active objective; add one "
                "(e.g. add_voltage_deviation_objective()) first"
            )
        if len(objectives) > 1:
            raise ValueError(
                "penalize_tap_movement() found more than one active "
                "objective; deactivate all but one"
            )
        objective = objectives[0]
        sign = 1.0 if objective.sense == pyo.minimize else -1.0
        objective.set_value(
            objective.expr + sign * model.trafo_tap_movement_cost
        )
        setup.cost_objective = objective.name
        logger.info(
            "penalize_tap_movement: switching cost added to objective '{}'",
            objective.name,
        )
    return model.trafo_tap_movement_cost


# ---------------------------------------------------------------------------
# reading the solution
# ---------------------------------------------------------------------------


def _position_value(model, key, discrete):
    """Solved tap position at `key`, rounded for a discrete changer.

    Args:
        model: The solved Pyomo model.
        key: `t` or `(t, tau)`.
        discrete: Whether to round to the nearest integer.

    Returns:
        A float (continuous mode) or int (discrete mode).
    """
    value = pyo.value(model.trafo_tap_position[key])
    return int(round(value)) if discrete else float(value)


def tap_schedule(model_obj):
    """Solved tap positions of the controlled transformers.

    Args:
        model_obj: A solved model object on which `enable_oltc` was called.

    Returns:
        Multi-period: a DataFrame indexed by time step with one column per
        controlled transformer. Single period: a Series indexed by
        transformer. Values are integers in `"discrete"` mode and floats in
        `"continuous"` mode.
    """
    setup = _require_oltc(model_obj)
    model = model_obj.model
    discrete = setup.mode == "discrete"
    if not setup.multi_period:
        return pd.Series(
            {
                t: _position_value(model, t, discrete)
                for t in setup.transformers
            },
            name="tap_pos",
        )
    steps = list(model.T)
    frame = pd.DataFrame(
        {
            t: [_position_value(model, (t, tau), discrete) for tau in steps]
            for t in setup.transformers
        },
        index=pd.Index(steps, name="t"),
    )
    frame.columns.name = "trafo"
    return frame


def tap_operations(model_obj):
    r"""Number of tap operations per controlled transformer.

    Counted from the solved positions as $\sum |k_\tau - k_{\tau^-}|$,
    including the move away from the initial position when it is tracked —
    never from the auxiliary up/down variables, which can exceed the real
    move when nothing prices or limits them.

    Args:
        model_obj: A solved model object on which `enable_oltc` was called.

    Returns:
        A Series indexed by transformer (floats in continuous mode).
    """
    setup = _require_oltc(model_obj)
    schedule = tap_schedule(model_obj)
    counts = {}
    for t in setup.transformers:
        positions = (
            [schedule[t]] if not setup.multi_period else list(schedule[t])
        )
        if setup.initial is not None:
            positions = [setup.initial[t]] + positions
        counts[t] = float(
            sum(abs(b - a) for a, b in zip(positions[:-1], positions[1:]))
        )
        if setup.mode == "discrete":
            counts[t] = int(round(counts[t]))
    return pd.Series(counts, name="tap_operations")


def apply_tap_positions(model_obj, net=None, t=None):
    """Write solved tap positions into a network's `trafo` table.

    Nothing writes `net.trafo.tap_pos` automatically: `solve()` reports the
    positions in `net.res_trafo["tap_pos"]` and leaves the input data
    alone. Call this to carry a solution over, e.g. before `pp.runpp` or a
    time-series run.

    Args:
        model_obj: A solved model object on which `enable_oltc` was called.
        net: The network to write into. Defaults to the model's own copy
            (`model_obj.net`); pass your original network to update it.
        t: Time step of a multi-period schedule to apply; defaults to the
            last one. Ignored on a single-period model.

    Returns:
        The network written to.

    Warns:
        UserWarning: When a continuous position is rounded, or when a
            transformer's `tap_changer_type` is `None` so pandapower would
            ignore the written position.
    """
    setup = _require_oltc(model_obj)
    schedule = tap_schedule(model_obj)
    if setup.multi_period:
        step = model_obj.model.T.last() if t is None else int(t)
        if step not in schedule.index:
            raise ValueError(f"t={t} is not a time step of this model")
        positions = schedule.loc[step]
    else:
        positions = schedule
    target = model_obj.net if net is None else net
    rounded_any = False
    for trafo, k in positions.items():
        if trafo not in target.trafo.index:
            raise KeyError(
                f"the target network has no transformer {trafo}; pass the "
                "network the model was built from (or a copy of it)"
            )
        k_int = int(round(float(k)))
        if abs(float(k) - k_int) > 1e-9:
            rounded_any = True
        target.trafo.at[trafo, "tap_pos"] = float(k_int)
        if "tap_changer_type" not in target.trafo.columns or _is_missing(
            target.trafo.at[trafo, "tap_changer_type"]
        ):
            warnings.warn(
                f"transformer {trafo}: tap_changer_type is None, so "
                "pandapower ignores the tap_pos just written; set it to "
                "'Ratio'",
                UserWarning,
                stacklevel=2,
            )
    if rounded_any:
        warnings.warn(
            "continuous tap positions were rounded to the nearest integer "
            "before writing them to net.trafo.tap_pos",
            UserWarning,
            stacklevel=2,
        )
    return target


def tap_result_columns(net, model, t=None):
    """`tap_pos` and `tap_factor` for every transformer of `net`.

    Controlled transformers report their solved position and tap factor;
    the others report the network's `tap_pos` and the factor pandapower
    applied to it (1 + (k − k_neutral) s for a "Ratio" changer, NaN for a
    transformer whose factor is not defined that way).

    Args:
        net: The network the results are written to.
        model: The solved Pyomo model carrying `TRANSF_OLTC`.
        t: Time step for a multi-period model, `None` for single period.

    Returns:
        Two Series indexed like `net.trafo`: positions and tap factors.
    """
    trafo = net.trafo
    tap_pos = pd.to_numeric(
        _column(trafo, "tap_pos", np.nan), errors="coerce"
    ).astype(float)
    neutral = pd.to_numeric(
        _column(trafo, "tap_neutral", np.nan), errors="coerce"
    ).astype(float)
    step = pd.to_numeric(
        _column(trafo, "tap_step_percent", np.nan), errors="coerce"
    ).astype(float)
    tap_type = _column(trafo, "tap_changer_type", None)
    is_ratio = pd.Series(
        [(not _is_missing(x)) and x == "Ratio" for x in tap_type],
        index=trafo.index,
    )
    factor = pd.Series(np.nan, index=trafo.index, dtype=float)
    factor[is_ratio] = (
        1.0 + (tap_pos - neutral)[is_ratio] * step[is_ratio] / 100.0
    )
    positions = tap_pos.copy()
    for tr in model.TRANSF_OLTC:
        key = tr if t is None else (tr, t)
        var = model.trafo_tap_position[key]
        value = pyo.value(var)
        positions[tr] = float(round(value)) if var.is_integer() else value
        factor[tr] = pyo.value(model.trafo_tap_factor[key])
    return positions, factor


# ---------------------------------------------------------------------------
# relax - round - fix - resolve
# ---------------------------------------------------------------------------


@dataclass
class RoundingGroup:
    """One family of discrete control variables for the rounding heuristic.

    Built by `_oltc_rounding_group` (tap positions) and by
    `potpourri.models.shunt_control` (shunt steps), so the heuristic can
    round every discrete control of a model in one pass.

    Attributes:
        name: Label used in messages, e.g. `"tap position"`.
        var: The indexed Pyomo variable holding the control.
        units: Controlled element indices, in order.
        bounds: `{unit: (lo, hi)}` integer bounds.
        initial: `{unit: reference value}` for the first move, or `None`.
        change_limit: `{unit: Δmax}` per step, or `None`.
        operations_limit: `{unit: Nmax}` over the horizon, or `None`.
        after_round: Callback `(key, value)` invoked once a control is
            fixed, e.g. to seed a dependent variable's start value.
    """

    name: str
    var: object
    units: list
    bounds: dict
    initial: dict | None = None
    change_limit: dict | None = None
    operations_limit: dict | None = None
    after_round: object = None


def _round_group(group, relaxed, time_steps):
    """Round one family of relaxed controls to integers.

    Sequential in time: each rounded value is clipped to the per-step
    change limit relative to the previous rounded one (the initial value
    first, when tracked) and to the bounds. The operation limit is checked
    afterwards and reported, not enforced.

    Args:
        group: The `RoundingGroup`.
        relaxed: `{key: relaxed value}` for every key of `group.var`.
        time_steps: Ordered time steps, or `[None]` for single period.

    Returns:
        `{key: integer value}` with the same keys as `relaxed`.
    """
    rounded = {}
    for unit in group.units:
        lo, hi = group.bounds[unit]
        previous = group.initial[unit] if group.initial is not None else None
        limit = (
            group.change_limit[unit]
            if group.change_limit is not None
            else None
        )
        moves = 0
        for tau in time_steps:
            key = unit if tau is None else (unit, tau)
            value = int(round(relaxed[key]))
            if previous is not None and limit is not None:
                reach = math.floor(limit)
                value = int(
                    min(max(value, previous - reach), previous + reach)
                )
            value = int(min(max(value, lo), hi))
            if previous is not None:
                moves += abs(value - previous)
            rounded[key] = value
            previous = value
        if (
            group.operations_limit is not None
            and moves > group.operations_limit[unit]
        ):
            warnings.warn(
                f"{group.name} of element {unit}: the rounded schedule needs "
                f"{moves} operations but the limit is "
                f"{group.operations_limit[unit]}; the fixed re-solve may be "
                "infeasible",
                UserWarning,
                stacklevel=3,
            )
    return rounded


def _oltc_rounding_group(model_obj) -> RoundingGroup:
    """The `RoundingGroup` of the tap positions of `model_obj`.

    Args:
        model_obj: A model object on which `enable_oltc` was called.

    Returns:
        The group, with a callback that seeds `trafo_tap_factor`.
    """
    setup = _require_oltc(model_obj)
    model = model_obj.model

    def _seed_factor(key, value):
        """Start the tap factor at the value the fixed position implies.

        Args:
            key: Variable index, `t` or `(t, tau)`.
            value: The fixed integer position.
        """
        model.trafo_tap_factor[key].set_value(
            setup.data[_trafo_of(key)].factor(value)
        )

    return RoundingGroup(
        name="tap position",
        var=model.trafo_tap_position,
        units=list(setup.transformers),
        bounds={t: (d.pos_min, d.pos_max) for t, d in setup.data.items()},
        initial=setup.initial,
        change_limit=setup.max_change_per_step,
        operations_limit=setup.max_operations,
        after_round=_seed_factor,
    )


def solve_discrete_round_and_fix(model_obj, solver="ipopt", **solve_kwargs):
    """Relax, round, fix and re-solve every discrete control of the model.

    1. Solve with the tap positions (and shunt steps) relaxed to
       continuous variables.
    2. Round each control to the nearest admissible integer, in time order
       and within its per-step change limit.
    3. Fix the controls and solve the resulting NLP again.

    A heuristic: the result is a feasible point of the discrete model
    (when the re-solve succeeds), not its global optimum, and the relaxed
    objective is a bound on it. The controls are left **fixed** so a
    further `solve()` re-uses them; `unfix_vars("trafo_tap_position")` /
    `unfix_vars("shunt_step")` release them. The relaxed values and the
    relaxed objective are kept in `model_obj.rounding_info`.

    Args:
        model_obj: A model object on which `enable_oltc` and/or
            `enable_shunt_control` was called.
        solver: NLP solver for both stages.
        **solve_kwargs: Forwarded to `solve()` in both stages (`to_net` is
            forced off for the first; on multi-period models `warm_start`
            defaults to `False` for the second so it starts from the
            relaxed solution).

    Returns:
        The solver results of the final, fixed-control solve.

    Raises:
        RuntimeError: If no discrete control is enabled or the relaxed
            problem does not solve to optimality.
    """
    model = model_obj.model
    groups = []
    if hasattr(model, "TRANSF_OLTC"):
        groups.append(_oltc_rounding_group(model_obj))
    if hasattr(model, "SHUNT_CTRL"):
        from potpourri.models.shunt_control import shunt_rounding_group

        groups.append(shunt_rounding_group(model_obj))
    if not groups:
        raise RuntimeError(
            "call enable_oltc() or enable_shunt_control() first"
        )
    time_steps = [None] if not hasattr(model, "T") else list(model.T)

    saved = []
    for group in groups:
        for idx in group.var:
            var = group.var[idx]
            if var.fixed:
                continue
            saved.append((var, var.domain))
            var.domain = pyo.Reals
    stage1 = dict(solve_kwargs)
    stage1["to_net"] = False
    try:
        relaxed_results = model_obj.solve(solver=solver, **stage1)
    finally:
        for var, domain in saved:
            var.domain = domain
    if not pyo.check_optimal_termination(relaxed_results):
        raise RuntimeError(
            "the relaxed problem did not solve to optimality "
            f"({relaxed_results.solver.termination_condition}); nothing was "
            "rounded"
        )
    relaxed_objective = None
    for objective in model.component_data_objects(pyo.Objective, active=True):
        relaxed_objective = float(pyo.value(objective))

    info = {
        "objective_relaxed": relaxed_objective,
        "relaxed": {},
        "rounded": {},
    }
    for group in groups:
        relaxed = {idx: float(pyo.value(group.var[idx])) for idx in group.var}
        rounded = _round_group(group, relaxed, time_steps)
        for idx, value in rounded.items():
            group.var[idx].fix(value)
            if group.after_round is not None:
                group.after_round(idx, value)
        info["relaxed"][group.var.name] = relaxed
        info["rounded"][group.var.name] = rounded

    stage2 = dict(solve_kwargs)
    parameters = inspect.signature(model_obj.solve).parameters
    if "warm_start" in parameters and "warm_start" not in stage2:
        stage2["warm_start"] = False
    results = model_obj.solve(solver=solver, **stage2)
    info["relaxed_results"] = relaxed_results
    model_obj.rounding_info = info
    logger.info(
        "round-and-fix: relaxed objective {}, rounded controls {}",
        relaxed_objective,
        info["rounded"],
    )
    return results


def solve_oltc_round_and_fix(model_obj, solver="ipopt", **solve_kwargs):
    """Relax, round, fix and re-solve: integer taps without a MINLP solver.

    The same as
    [`solve_discrete_round_and_fix`][potpourri.models.oltc.solve_discrete_round_and_fix];
    kept under the name the OLTC documentation uses. Shunt steps enabled on
    the same model are rounded in the same pass.

    Args:
        model_obj: A model object on which `enable_oltc` was called.
        solver: NLP solver for both stages.
        **solve_kwargs: Forwarded to `solve()`.

    Returns:
        The solver results of the final solve.
    """
    _require_oltc(model_obj)
    return solve_discrete_round_and_fix(
        model_obj, solver=solver, **solve_kwargs
    )


# ---------------------------------------------------------------------------
# the mix-in that OPF and OPF_multi_period expose
# ---------------------------------------------------------------------------


class OLTCControlMixin:
    """Methods that make selected transformer taps OPF decision variables.

    Mixed into [`OPF`][potpourri.models.OPF.OPF] and
    `OPF_multi_period`, so every `ACOPF`, `HC_ACOPF` and
    `ACOPF_multi_period` object has them; on DC and LPAC models they
    raise `NotImplementedError`. Nothing changes until `enable_oltc` is
    called. See the module docstring for the physics and the modes.
    """

    def enable_oltc(
        self,
        transformers=None,
        mode="discrete",
        *,
        max_change_per_step=None,
        max_operations=None,
        initial_tap="net",
    ):
        """Make the tap positions of selected transformers decision variables.

        Adds `trafo_tap_position` (integer in `"discrete"` mode, real in
        `"continuous"` mode, bounded by `tap_min`/`tap_max`), the tap
        factor and its link to the HV- or LV-side ratio of the transformer
        equations, and the up/down movement variables. Call it before or
        after `add_OPF()`, but before solving; once per model.

        Args:
            transformers: `None` (default) controls every eligible
                transformer — see `oltc_eligibility(net)` — or an iterable
                of `net.trafo` indices, each of which must be eligible.
            mode: `"discrete"` (physical positions; MINLP) or
                `"continuous"` (relaxation; NLP).
            max_change_per_step: Largest position change between
                consecutive time steps (single period: away from the
                initial position). Scalar or per-transformer mapping.
            max_operations: Largest total number of position changes over
                the horizon. Scalar or per-transformer mapping.
            initial_tap: `"net"` (default) measures movement from
                `net.trafo.tap_pos`; a mapping gives other reference
                positions; `None` leaves the first move untracked.

        Returns:
            The `OLTCSetup` describing the controlled transformers.

        Raises:
            NotImplementedError: On a DC or LPAC model.
            ValueError: If a requested transformer is not eligible, or no
                transformer is when `transformers=None`.
            RuntimeError: If already enabled on this model.
        """
        return attach_oltc(
            self,
            transformers,
            mode,
            max_change_per_step=max_change_per_step,
            max_operations=max_operations,
            initial_tap=initial_tap,
        )

    def penalize_tap_movement(self, cost):
        """Add a cost per tap operation to the active objective.

        Args:
            cost: Cost per position change, in the unit of the objective;
                a scalar or a per-transformer mapping. Later calls update
                the price without rebuilding.

        Returns:
            The `trafo_tap_movement_cost` expression.
        """
        return penalize_tap_movement(self, cost)

    def tap_schedule(self):
        """Solved tap positions.

        Returns:
            A DataFrame (time step × transformer) on a multi-period model,
            a Series (per transformer) on a single-period one.
        """
        return tap_schedule(self)

    def tap_operations(self):
        """Number of tap operations per controlled transformer.

        Returns:
            A Series indexed by transformer, counted from the solved
            positions including the move away from the initial position.
        """
        return tap_operations(self)

    def apply_tap_positions(self, net=None, t=None):
        """Write the solved tap positions into `net.trafo.tap_pos`.

        Args:
            net: Network to write into (default: the model's own copy).
            t: Multi-period time step to apply (default: the last).

        Returns:
            The network written to.
        """
        return apply_tap_positions(self, net=net, t=t)

    def solve_oltc_round_and_fix(self, solver="ipopt", **solve_kwargs):
        """Relax, round, fix and re-solve the tap positions.

        Integer tap schedules with an NLP solver only; a heuristic. See
        [`solve_oltc_round_and_fix`][potpourri.models.oltc.solve_oltc_round_and_fix].

        Args:
            solver: NLP solver for both stages.
            **solve_kwargs: Forwarded to `solve()`.

        Returns:
            The solver results of the final solve.
        """
        return solve_oltc_round_and_fix(self, solver=solver, **solve_kwargs)
