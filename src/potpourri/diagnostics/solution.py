# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Checks on the solution a solver returned.

`termination_condition == optimal` means the solver is satisfied with its
own tolerances. It does not mean the answer is physically sensible, and it
does not mean the numbers in `net.res_*` are what you expect. Everything
here re-evaluates the model independently of the solver's own verdict.

Three layers, in order of how much they assume:

1. **Feasibility** — evaluate every constraint and report what is
   violated, in engineering units rather than as a squared residual.
2. **Margins** — for constraints that are satisfied, how much room is
   left. This is what answers "the result looks wrong": it shows which
   physical limits are shaping the dispatch.
3. **Plausibility** — system-level sanity. Voltage spread, loading,
   losses, and whether the power balance closes.
"""

from __future__ import annotations

import math

from potpourri.diagnostics.context import DiagnosticContext
from potpourri.diagnostics.metadata import constraint_meta
from potpourri.diagnostics.report import (
    DiagnosticCategory,
    DiagnosticIssue,
    DiagnosticReport,
    DiagnosticSeverity,
)

#: Loss fraction above which the network result is worth a second look.
#: Distribution networks rarely exceed a few percent; 25 % usually means a
#: modelling error rather than a very lossy grid.
_IMPLAUSIBLE_LOSS_FRACTION = 0.25


def check_constraint_violations(ctx: DiagnosticContext) -> DiagnosticReport:
    """Re-evaluate every constraint and report what the solution violates.

    What it checks: the body of every active constraint against its
    bounds, using the variable values currently loaded on the model.

    Necessary or sufficient: exact, within `ctx.tol`. A violation here is
    a real violation of the model as written — which may mean the solver
    stopped early, or that the tolerance is tighter than the solver's.

    Cost: cheap; one expression evaluation per constraint.

    Required state: variable values must be loaded. On a model that never
    solved, the values are the initialisation, and the check says so
    rather than reporting nonsense.

    Args:
        ctx: The diagnostic context.

    Returns:
        A report with one `ERROR` per violated constraint, carrying the
        physical meaning where the metadata registry knows the family.
    """
    report = DiagnosticReport()
    if ctx.model is None:
        report.skip(
            "constraint_violations", "the Pyomo model has not been built"
        )
        return report
    if not ctx.has_solution:
        report.skip(
            "constraint_violations",
            "no solver result is loaded, so the values are the initialisation",
        )
        return report

    import pyomo.environ as pyo

    worst = 0.0
    count = 0
    for con in ctx.model.component_data_objects(pyo.Constraint, active=True):
        body = ctx.value(con.body)
        if body is None:
            continue
        low = ctx.value(con.lower)
        high = ctx.value(con.upper)
        violation = 0.0
        if low is not None and body < low - ctx.tol:
            violation = low - body
        elif high is not None and body > high + ctx.tol:
            violation = body - high
        if violation <= 0:
            continue

        count += 1
        worst = max(worst, violation)
        parent = con.parent_component()
        meta = constraint_meta(parent.local_name)
        index = con.index()
        element, time_step = _element_and_time(ctx, parent.local_name, index)
        scale = max(abs(low or 0.0), abs(high or 0.0), 1.0)

        report.add(
            DiagnosticIssue(
                severity=DiagnosticSeverity.ERROR,
                category=meta.category
                if meta
                else DiagnosticCategory.INFEASIBILITY,
                code="SOLUTION_CONSTRAINT_VIOLATED",
                message=(
                    f"{meta.description if meta else parent.local_name} is "
                    f"violated by {violation:.6g} in the model's own units."
                    + (
                        " The constraint is written on squared quantities, so"
                        " this residual is not a power."
                        if meta and meta.squared
                        else ""
                    )
                ),
                recommendation=(
                    "The solver reported success, so either it stopped within "
                    "a looser tolerance than this check uses, or the solution "
                    "was not loaded."
                ),
                element=element,
                pyomo_component=parent.local_name,
                pyomo_index=index,
                value=body,
                lower_bound=low,
                upper_bound=high,
                violation=violation,
                relative_violation=violation / scale,
                time_step=time_step,
            )
        )

    report.summary["solution.violated constraints"] = count
    if count:
        report.summary["solution.worst violation"] = f"{worst:.3e}"
    return report


def check_binding_constraints(ctx: DiagnosticContext) -> DiagnosticReport:
    """Report which limits are shaping the solution.

    This is the answer to "the result looks wrong": a dispatch that seems
    odd is usually a dispatch pressed against a limit. Listing the binding
    and nearly binding constraints shows which ones.

    What it checks: every satisfied inequality, for how close the body sits
    to its bound, relative to the bound's own magnitude.

    Necessary or sufficient: neither; this is information about a feasible
    point.

    Cost: cheap.

    Args:
        ctx: The diagnostic context.

    Returns:
        A report with one `INFO` per binding or nearly binding constraint,
        capped at the most constrained elements so a large network does
        not produce hundreds of lines.
    """
    report = DiagnosticReport()
    if ctx.model is None or not ctx.has_solution:
        report.skip("binding_constraints", "no solution to inspect")
        return report

    import pyomo.environ as pyo

    found = []
    for con in ctx.model.component_data_objects(pyo.Constraint, active=True):
        if con.equality:
            continue
        body = ctx.value(con.body)
        if body is None:
            continue
        low = ctx.value(con.lower)
        high = ctx.value(con.upper)
        margin = None
        bound = None
        if high is not None:
            margin, bound = high - body, high
        if low is not None:
            other = body - low
            if margin is None or other < margin:
                margin, bound = other, low
        if margin is None or margin < -ctx.tol:
            continue

        parent_name = con.parent_component().local_name
        element, _ = _element_and_time(ctx, parent_name, con.index())
        loading = _loading_percent(ctx, element)
        if loading is not None:
            # A squared thermal constraint's raw margin is in p.u. squared,
            # and comparing it against a floor of 1.0 made a line at 18 %
            # loading look "nearly active". Loading is the same fact on a
            # scale that means something, so it decides here.
            utilisation = loading / 100.0
        else:
            scale = max(abs(bound), abs(body), ctx.tol)
            utilisation = 1.0 - margin / scale
        if utilisation < 1.0 - ctx.near_bound_fraction:
            continue
        found.append((1.0 - utilisation, margin, con, bound))

    found.sort(key=lambda entry: entry[0])
    limit = int(ctx.options.get("max_binding", 15))
    for relative, margin, con, bound in found[:limit]:
        parent = con.parent_component()
        meta = constraint_meta(parent.local_name)
        index = con.index()
        element, time_step = _element_and_time(ctx, parent.local_name, index)
        # A thermal limit is written on squared apparent power, so its
        # raw margin is in p.u. squared and means nothing to a reader.
        # Where pandapower has computed a loading for the same element,
        # report that instead: "98.4 % loaded" is the same fact in units
        # an engineer already thinks in.
        loading = _loading_percent(ctx, element)
        if loading is not None:
            detail = f"the element is at {loading:.1f} % of its rating"
        elif meta and meta.squared:
            detail = (
                f"{margin:.4g} of room left, in the squared units the "
                f"constraint is written in"
            )
        else:
            detail = f"{margin:.4g} of room left"

        report.add(
            DiagnosticIssue(
                severity=DiagnosticSeverity.INFO,
                category=meta.category
                if meta
                else DiagnosticCategory.STRUCTURE,
                code="SOLUTION_CONSTRAINT_BINDING",
                message=(
                    f"{meta.description if meta else parent.local_name} is "
                    f"{'active' if margin <= ctx.tol else 'nearly active'}: "
                    f"{detail}."
                ),
                recommendation=(
                    "This limit is shaping the solution. Relaxing it is what "
                    "would change the result here."
                ),
                element=element,
                pyomo_component=parent.local_name,
                pyomo_index=index,
                value=ctx.value(con.body),
                upper_bound=bound,
                relative_violation=-relative,
                time_step=time_step,
            )
        )
    report.summary["solution.binding constraints"] = len(found)
    return report


def check_power_balance(ctx: DiagnosticContext) -> DiagnosticReport:
    """Check that generation, demand and losses add up, in MW and MVAr.

    What it checks: the system totals from `net.res_*` after the solution
    has been written back, and whether the active balance closes.

    Necessary or sufficient: neither. A balance that closes is expected;
    one that does not points at the result mapping rather than at the
    optimisation.

    Cost: negligible.

    Required state: `net.res_bus` and friends must be populated, which
    `solve(to_net=True)` does.

    Args:
        ctx: The diagnostic context.

    Returns:
        A report carrying the totals in `summary`, with a `WARNING` if the
        active balance does not close or the losses look implausible.
    """
    report = DiagnosticReport()
    net = ctx.net
    if net is None or net.get("res_bus") is None or len(net.res_bus) == 0:
        report.skip(
            "power_balance", "no results have been written to the network"
        )
        return report

    def total(table: str, column: str) -> float:
        """Sum one result column, treating a missing table as zero.

        Args:
            table: Result table name, e.g. `"res_line"`.
            column: Column to sum, e.g. `"p_mw"`.

        Returns:
            The column total, or 0.0 when the table or column is absent.
        """
        frame = net.get(table)
        if frame is None or len(frame) == 0 or column not in frame.columns:
            return 0.0
        return float(frame[column].fillna(0.0).sum())

    ext = total("res_ext_grid", "p_mw")
    gen = total("res_gen", "p_mw")
    sgen = total("res_sgen", "p_mw")
    load = total("res_load", "p_mw")
    storage = total("res_storage", "p_mw")
    line_loss = total("res_line", "pl_mw")
    trafo_loss = total("res_trafo", "pl_mw")

    supply = ext + gen + sgen
    losses = line_loss + trafo_loss
    residual = supply - load - storage - losses

    report.summary.update(
        {
            "balance.external grid": f"{ext:.4g} MW",
            "balance.generators": f"{gen:.4g} MW",
            "balance.static generators": f"{sgen:.4g} MW",
            "balance.load": f"{load:.4g} MW",
            "balance.losses": f"{losses:.4g} MW",
            "balance.residual": f"{residual:.3e} MW",
        }
    )

    tolerance = max(1e-6, 1e-4 * max(abs(supply), abs(load), 1.0))
    if abs(residual) > tolerance:
        report.add(
            DiagnosticIssue(
                severity=DiagnosticSeverity.WARNING,
                category=DiagnosticCategory.POWER_BALANCE,
                code="RESULT_BALANCE_MISMATCH",
                message=(
                    f"Active power in and out differ by {residual:.4g} MW "
                    f"once losses are counted."
                ),
                recommendation=(
                    "The optimisation enforces nodal balance, so a system "
                    "residual points at the result mapping rather than at "
                    "the solution."
                ),
                value=residual,
                unit="MW",
            )
        )

    if load > 0:
        fraction = losses / load
        report.summary["balance.loss fraction"] = f"{fraction * 100:.2f} %"
        if fraction > _IMPLAUSIBLE_LOSS_FRACTION or fraction < -1e-6:
            report.add(
                DiagnosticIssue(
                    severity=DiagnosticSeverity.WARNING,
                    category=DiagnosticCategory.POWER_BALANCE,
                    code="RESULT_IMPLAUSIBLE_LOSSES",
                    message=(
                        f"Network losses are {fraction * 100:.1f} % of the "
                        f"served load, which is outside the range a healthy "
                        f"network usually shows."
                    ),
                    recommendation=(
                        "Negative or very large losses usually mean a branch "
                        "model or a sign convention is wrong, not that the "
                        "optimum is."
                    ),
                    value=losses,
                    unit="MW",
                )
            )
    return report


def check_voltage_profile(ctx: DiagnosticContext) -> DiagnosticReport:
    """Report the voltage spread and any bus outside its band.

    What it checks: `net.res_bus.vm_pu` against `min_vm_pu`/`max_vm_pu`.

    Necessary or sufficient: exact against the declared band.

    Cost: negligible.

    Args:
        ctx: The diagnostic context.

    Returns:
        A report with the extremes in `summary` and an `ERROR` per bus
        outside its band.
    """
    report = DiagnosticReport()
    net = ctx.net
    res = net.get("res_bus") if net is not None else None
    if res is None or len(res) == 0 or "vm_pu" not in res.columns:
        report.skip("voltage_profile", "no bus results are present")
        return report

    values = res["vm_pu"].dropna()
    if values.empty:
        report.skip("voltage_profile", "bus voltages are all missing")
        return report

    report.summary["voltage.minimum"] = f"{values.min():.4f} p.u."
    report.summary["voltage.maximum"] = f"{values.max():.4f} p.u."

    if "min_vm_pu" not in net.bus.columns:
        return report

    reverse = {v: k for k, v in ctx.imap.bus.items()}
    for bus in values.index:
        if bus not in net.bus.index:
            continue
        vm = float(values[bus])
        low = net.bus.at[bus, "min_vm_pu"]
        high = net.bus.at[bus, "max_vm_pu"]
        if low != low or high != high or low <= vm <= high:
            continue
        over = vm > high
        report.add(
            DiagnosticIssue(
                severity=DiagnosticSeverity.ERROR,
                category=DiagnosticCategory.VOLTAGE,
                code="BUS_VOLTAGE_HIGH" if over else "BUS_VOLTAGE_LOW",
                message=(
                    f"voltage {vm:.4f} p.u. "
                    + (
                        f"exceeds max_vm_pu={float(high):.4f}"
                        if over
                        else f"is below min_vm_pu={float(low):.4f}"
                    )
                ),
                element=ctx.imap.bus_ref(reverse[int(bus)])
                if int(bus) in reverse
                else None,
                value=vm,
                lower_bound=float(low),
                upper_bound=float(high),
                violation=vm - float(high) if over else float(low) - vm,
                unit="p.u.",
                pyomo_component="v" if ctx.has("v") else None,
                pyomo_index=reverse.get(int(bus)),
            )
        )
    return report


def check_branch_loading(ctx: DiagnosticContext) -> DiagnosticReport:
    """Report branch loading, and anything over its rating.

    What it checks: `loading_percent` on lines and transformers.

    Necessary or sufficient: exact against what pandapower computed from
    the mapped solution.

    Cost: negligible.

    Args:
        ctx: The diagnostic context.

    Returns:
        A report with the worst loading in `summary` and an `ERROR` per
        overloaded branch.
    """
    report = DiagnosticReport()
    net = ctx.net
    if net is None:
        report.skip("branch_loading", "no network on the model")
        return report

    worst = 0.0
    for table, res_table in (("line", "res_line"), ("trafo", "res_trafo")):
        res = net.get(res_table)
        if (
            res is None
            or len(res) == 0
            or "loading_percent" not in res.columns
        ):
            continue
        loading = res["loading_percent"].dropna()
        if loading.empty:
            continue
        worst = max(worst, float(loading.max()))
        for index in loading.index[loading > 100.0]:
            report.add(
                DiagnosticIssue(
                    severity=DiagnosticSeverity.ERROR,
                    category=DiagnosticCategory.THERMAL,
                    code="LINE_THERMAL_OVERLOAD"
                    if table == "line"
                    else "TRAFO_THERMAL_OVERLOAD",
                    message=(
                        f"loading {float(loading[index]):.1f} %, above the "
                        f"element's rating"
                    ),
                    element=ctx.imap.element_ref(table, index),
                    value=float(loading[index]),
                    upper_bound=100.0,
                    violation=float(loading[index]) - 100.0,
                    unit="%",
                    pyomo_component="line_lim_from"
                    if table == "line"
                    else "transf_lim1",
                    pyomo_index=int(index),
                )
            )
    if worst:
        report.summary["thermal.worst loading"] = f"{worst:.1f} %"
    return report


def explain_bus(ctx: DiagnosticContext, bus: int) -> dict:
    """Everything attached to one bus, for reading around a problem.

    Args:
        ctx: The diagnostic context.
        bus: A bus index in the *caller's* numbering.

    Returns:
        A dict with the bus's voltage, its connected branches, and the
        loads, generators, static generators and storage sitting on it,
        each with its operating point where results exist. Empty values
        where the network does not carry that table.
    """
    net = ctx.net
    out: dict = {"bus": bus, "elements": {}, "branches": []}
    if net is None or bus not in net.bus.index:
        return out

    res = net.get("res_bus")
    if res is not None and bus in getattr(res, "index", []):
        out["vm_pu"] = float(res.at[bus, "vm_pu"])
        if "va_degree" in res.columns:
            out["va_degree"] = float(res.at[bus, "va_degree"])

    for table in ("load", "sgen", "gen", "ext_grid", "storage", "shunt"):
        frame = net.get(table)
        if frame is None or len(frame) == 0 or "bus" not in frame.columns:
            continue
        rows = frame.index[frame["bus"] == bus]
        if not len(rows):
            continue
        res_frame = net.get(f"res_{table}")
        entries = []
        for index in rows:
            entry = {
                "index": int(index),
                "ref": str(ctx.imap.element_ref(table, index)),
            }
            if res_frame is not None and index in getattr(
                res_frame, "index", []
            ):
                for column in ("p_mw", "q_mvar"):
                    if column in res_frame.columns:
                        entry[column] = float(res_frame.at[index, column])
            entries.append(entry)
        out["elements"][table] = entries

    for table, from_col, to_col in (
        ("line", "from_bus", "to_bus"),
        ("trafo", "hv_bus", "lv_bus"),
    ):
        frame = net.get(table)
        if frame is None or len(frame) == 0:
            continue
        rows = frame.index[(frame[from_col] == bus) | (frame[to_col] == bus)]
        res_frame = net.get(f"res_{table}")
        for index in rows:
            entry = {
                "index": int(index),
                "ref": str(ctx.imap.element_ref(table, index)),
            }
            if (
                res_frame is not None
                and index in getattr(res_frame, "index", [])
                and "loading_percent" in res_frame.columns
            ):
                entry["loading_percent"] = float(
                    res_frame.at[index, "loading_percent"]
                )
            out["branches"].append(entry)
    return out


# --- helpers ----------------------------------------------------------


def _loading_percent(ctx: DiagnosticContext, element) -> float | None:
    """The loading pandapower computed for this element, if it has one.

    Args:
        ctx: The diagnostic context.
        element: The `ElementRef` a constraint was mapped to.

    Returns:
        `loading_percent` for a line or transformer, or `None` when the
        element is of another kind or no result exists for it.
    """
    if element is None or element.table not in ("line", "trafo"):
        return None
    res = ctx.net.get(f"res_{element.table}") if ctx.net is not None else None
    if res is None or "loading_percent" not in getattr(res, "columns", []):
        return None
    if element.index not in res.index:
        return None
    value = res.at[element.index, "loading_percent"]
    return None if value != value else float(value)


def _element_and_time(ctx: DiagnosticContext, component: str, index):
    """Resolve a constraint index to a network object and a time step.

    Args:
        ctx: The diagnostic context.
        component: Pyomo component local name.
        index: The constraint's index, possibly a tuple whose last entry
            is a time step.

    Returns:
        `(ElementRef | None, time_step | None)`.
    """
    from potpourri.diagnostics.metadata import element_table_for

    table = element_table_for(component)
    time_step = None
    key = index
    if isinstance(index, tuple):
        if ctx.is_multi_period and len(index) >= 2:
            time_step = _as_int(index[-1])
            key = index[0]
        else:
            key = index[0]
    if table is None or key is None:
        return None, time_step
    try:
        if table == "bus":
            return ctx.imap.bus_ref(key), time_step
        if table == "line":
            return ctx.imap.line_ref(key), time_step
        return ctx.imap.element_ref(table, key), time_step
    except (TypeError, ValueError):
        return None, time_step


def _as_int(value):
    """Best-effort int conversion that returns `None` instead of raising."""
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _finite(value) -> bool:
    """Whether a value is a real, finite number."""
    return (
        value is not None and isinstance(value, float) and math.isfinite(value)
    )
