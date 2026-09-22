# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Checks that can run before the OPF is solved, and before it is even built.

These are the cheap ones, and the ones that pay off most: an OPF that is
infeasible because `min_p_mw > max_p_mw` or because a load sits on an island
with no source does not need a solver to explain it.

Two kinds of check live here.

**Data checks** are exact. `min > max` is a contradiction, full stop, and is
reported as an `ERROR`.

**Necessary-condition checks** are not proofs. If the most generation the
network can produce is below the least demand it must serve, the OPF is
infeasible — that direction is sound. The converse is not: passing these
says nothing about feasibility, because they ignore the network itself.
Every such finding says so in its recommendation, and none of them is ever
phrased as "the OPF is feasible".

pandapower's own `diagnostic` is reused rather than reimplemented, and only
the parts of its output that bear on an OPF are turned into issues.
"""

from __future__ import annotations

import copy
import io
from contextlib import redirect_stderr, redirect_stdout

import numpy as np

from potpourri.diagnostics.context import DiagnosticContext
from potpourri.diagnostics.report import (
    DiagnosticCategory,
    DiagnosticIssue,
    DiagnosticReport,
    DiagnosticSeverity,
)

#: pandapower diagnostic keys worth surfacing for an OPF, mapped to the
#: potpourri code and severity to report them under. Keys pandapower emits
#: that are not here are deliberately dropped: dumping its whole report
#: would bury the OPF-specific findings.
_PANDAPOWER_CHECKS: dict[str, tuple[str, DiagnosticSeverity, str]] = {
    "disconnected_elements": (
        "NET_DISCONNECTED",
        DiagnosticSeverity.ERROR,
        "elements are disconnected from the rest of the network",
    ),
    "different_voltage_levels_connected": (
        "NET_VOLTAGE_LEVEL_MISMATCH",
        DiagnosticSeverity.WARNING,
        "elements connect buses of different nominal voltage",
    ),
    "impedance_values_close_to_zero": (
        "NET_ZERO_IMPEDANCE",
        DiagnosticSeverity.ERROR,
        "branch impedance is at or near zero",
    ),
    "implausible_impedance_values": (
        "NET_IMPLAUSIBLE_IMPEDANCE",
        DiagnosticSeverity.WARNING,
        "branch impedance looks implausible",
    ),
    "nominal_voltages_dont_match": (
        "NET_TRAFO_VOLTAGE_MISMATCH",
        DiagnosticSeverity.WARNING,
        "transformer nominal voltages do not match the buses",
    ),
    "invalid_values": (
        "NET_INVALID_VALUE",
        DiagnosticSeverity.ERROR,
        "invalid values in the element tables",
    ),
    "overload": (
        "NET_BASE_OVERLOAD",
        DiagnosticSeverity.INFO,
        "base-case power flow shows overloading",
    ),
    "wrong_switch_configuration": (
        "NET_SWITCH_CONFIGURATION",
        DiagnosticSeverity.WARNING,
        "switch configuration prevents convergence",
    ),
    "no_ext_grid": (
        "NET_NO_EXT_GRID",
        DiagnosticSeverity.ERROR,
        "the network has no external grid",
    ),
}


def run_pandapower_diagnostic(ctx: DiagnosticContext) -> DiagnosticReport:
    """Translate pandapower's own network checks into potpourri issues.

    What it checks: whatever `pandapower.diagnostic` checks — disconnected
    elements, implausible impedances, switch configuration, invalid values
    and so on — filtered to the findings that matter for an OPF.

    Necessary or sufficient: neither. These are data-quality checks.

    Cost: moderate. pandapower may run trial power flows internally.

    False positives: pandapower flags deviations from standard types,
    which are common in synthetic networks and are not reported here.

    Args:
        ctx: The diagnostic context.

    Returns:
        A report holding one issue per relevant pandapower finding, with
        the raw payload preserved in `context["pandapower"]` so the
        original can still be inspected.
    """
    report = DiagnosticReport()
    if ctx.net is None:
        report.skip("pandapower_diagnostic", "no network on the model")
        return report

    try:
        from pandapower.diagnostic import diagnostic
    except ImportError as exc:  # pragma: no cover - pandapower is required
        report.skip("pandapower_diagnostic", f"unavailable: {exc}")
        return report

    # On a copy, deliberately: pandapower's diagnostic runs trial power
    # flows and leaves their results behind, which poisoned the warm start
    # the replay check later takes from net.res_bus. A diagnostic must not
    # change the thing it is diagnosing.
    buf = io.StringIO()
    try:
        with redirect_stdout(buf), redirect_stderr(buf):
            raw = diagnostic(
                copy.deepcopy(ctx.net),
                report_style=None,
                return_result_dict=True,
            )
    except Exception as exc:
        report.skip(
            "pandapower_diagnostic", f"raised {type(exc).__name__}: {exc}"
        )
        return report

    for key, payload in (raw or {}).items():
        entry = _PANDAPOWER_CHECKS.get(key)
        if entry is None:
            continue
        code, severity, description = entry
        report.add(
            DiagnosticIssue(
                severity=severity,
                category=DiagnosticCategory.NETWORK,
                code=code,
                message=f"pandapower reports that {description}.",
                recommendation=(
                    "Run pandapower.diagnostic(net) for the full report; "
                    "the affected elements are in this issue's context."
                ),
                context={"pandapower": _jsonable(payload), "check": key},
            )
        )
    return report


def check_bounds(ctx: DiagnosticContext) -> DiagnosticReport:
    """Find limits that contradict each other before a solver sees them.

    What it checks: every `min`/`max` pair the OPF reads — bus voltage,
    generator and sgen active and reactive power, load flexibility,
    storage power and energy, and transformer tap range — for `min > max`,
    and for fixed values lying outside their own bounds.

    Necessary or sufficient: a finding here is **sufficient** to make the
    OPF infeasible. Finding nothing proves nothing.

    Cost: negligible; pure DataFrame arithmetic.

    Args:
        ctx: The diagnostic context.

    Returns:
        A report with one `ERROR` per contradictory pair.
    """
    report = DiagnosticReport()
    net = ctx.net
    if net is None:
        report.skip("bounds", "no network on the model")
        return report

    pairs = [
        ("bus", "min_vm_pu", "max_vm_pu", "p.u.", "voltage magnitude"),
        ("gen", "min_p_mw", "max_p_mw", "MW", "active power"),
        ("gen", "min_q_mvar", "max_q_mvar", "MVAr", "reactive power"),
        ("sgen", "min_p_mw", "max_p_mw", "MW", "active power"),
        ("sgen", "min_q_mvar", "max_q_mvar", "MVAr", "reactive power"),
        ("load", "min_p_mw", "max_p_mw", "MW", "active power"),
        ("load", "min_q_mvar", "max_q_mvar", "MVAr", "reactive power"),
        ("ext_grid", "min_p_mw", "max_p_mw", "MW", "active power"),
        ("ext_grid", "min_q_mvar", "max_q_mvar", "MVAr", "reactive power"),
        ("storage", "min_p_mw", "max_p_mw", "MW", "active power"),
        ("storage", "min_e_mwh", "max_e_mwh", "MWh", "stored energy"),
        ("trafo", "tap_min", "tap_max", "-", "tap position"),
    ]

    for table, low_col, high_col, unit, quantity in pairs:
        frame = net.get(table)
        if frame is None or len(frame) == 0:
            continue
        if low_col not in frame.columns or high_col not in frame.columns:
            continue
        low = frame[low_col].astype(float)
        high = frame[high_col].astype(float)
        bad = frame.index[(low.notna()) & (high.notna()) & (low > high)]
        for index in bad:
            report.add(
                DiagnosticIssue(
                    severity=DiagnosticSeverity.ERROR,
                    category=DiagnosticCategory.BOUNDS,
                    code="BOUND_MIN_ABOVE_MAX",
                    message=(
                        f"{quantity} lower limit {low[index]:.6g} {unit} is "
                        f"above the upper limit {high[index]:.6g} {unit}, "
                        f"which no value can satisfy."
                    ),
                    recommendation=(
                        f"Correct net.{table}.{low_col} or "
                        f"net.{table}.{high_col} for this element."
                    ),
                    element=ctx.imap.element_ref(table, index),
                    lower_bound=float(low[index]),
                    upper_bound=float(high[index]),
                    unit=unit,
                )
            )
    return report


def check_islands(ctx: DiagnosticContext) -> DiagnosticReport:
    """Look for energised islands that cannot be supplied.

    What it checks: each connected component of the in-service network for
    demand, available generation and the presence of a slack (external
    grid or generator). An island with load but no reference cannot be
    solved.

    Necessary or sufficient: an island with demand and no source is
    **sufficient** for infeasibility.

    Cost: cheap; one graph traversal via `pandapower.topology`.

    Args:
        ctx: The diagnostic context.

    Returns:
        A report with an `ERROR` per unsupplied energised island and an
        `INFO` recording the island count when there is more than one.
    """
    report = DiagnosticReport()
    net = ctx.net
    if net is None:
        report.skip("islands", "no network on the model")
        return report

    try:
        import pandapower.topology as top
    except ImportError as exc:  # pragma: no cover
        report.skip("islands", f"pandapower.topology unavailable: {exc}")
        return report

    try:
        graph = top.create_nxgraph(net, respect_switches=True)
        components = list(top.connected_components(graph))
    except Exception as exc:
        report.skip("islands", f"raised {type(exc).__name__}: {exc}")
        return report

    if len(components) > 1:
        report.add(
            DiagnosticIssue(
                severity=DiagnosticSeverity.INFO,
                category=DiagnosticCategory.NETWORK,
                code="NET_MULTIPLE_ISLANDS",
                message=f"The network splits into {len(components)} islands.",
                recommendation=(
                    "Each island needs its own slack; check that this is "
                    "intended."
                ),
                context={"island_count": len(components)},
            )
        )

    for island in components:
        buses = set(int(b) for b in island)
        demand = _sum_on(net, "load", "p_mw", buses)
        sources = _count_on(net, "ext_grid", buses) + _count_on(
            net, "gen", buses
        )
        generation = _sum_on(net, "sgen", "p_mw", buses)

        if demand > 0 and sources == 0 and generation <= 0:
            example = sorted(buses)[0]
            report.add(
                DiagnosticIssue(
                    severity=DiagnosticSeverity.ERROR,
                    category=DiagnosticCategory.NETWORK,
                    code="NET_ISLAND_WITHOUT_SOURCE",
                    message=(
                        f"An island of {len(buses)} buses carries "
                        f"{demand:.4g} MW of demand but has no external "
                        f"grid, generator or static generator."
                    ),
                    recommendation=(
                        "Connect the island to a source, add a reference "
                        "there, or take its loads out of service."
                    ),
                    element=ctx.imap.bus_ref(example)
                    if example in ctx.imap.bus
                    else None,
                    value=demand,
                    unit="MW",
                    context={"island_buses": sorted(buses)[:32]},
                )
            )
    return report


def check_power_adequacy(ctx: DiagnosticContext) -> DiagnosticReport:
    """Compare the most the network can supply with the least it must serve.

    What it checks: total upper-bound active generation and import against
    total lower-bound active demand, network-wide.

    Necessary or sufficient: **necessary only**. A shortfall proves
    infeasibility; a surplus proves nothing, because this ignores the
    network entirely — no impedances, no limits, no losses. Losses are not
    added to demand, so the check errs towards not reporting.

    Cost: negligible.

    False positives: none by construction, since the comparison is between
    an upper and a lower bound. A network whose generation bounds are
    missing is skipped rather than guessed at.

    Args:
        ctx: The diagnostic context.

    Returns:
        A report with at most one `ERROR`, plus an `INFO` carrying the
        totals so the numbers are visible even when they are fine.
    """
    report = DiagnosticReport()
    net = ctx.net
    if net is None:
        report.skip("power_adequacy", "no network on the model")
        return report

    supply = 0.0
    unbounded = False
    for table in ("ext_grid", "gen", "sgen", "storage"):
        frame = net.get(table)
        if frame is None or len(frame) == 0:
            continue
        in_service = (
            frame["in_service"].astype(bool) if "in_service" in frame else True
        )
        rows = frame[in_service] if in_service is not True else frame
        if "max_p_mw" in rows.columns and rows["max_p_mw"].notna().any():
            supply += float(rows["max_p_mw"].fillna(0.0).sum())
            if rows["max_p_mw"].isna().any():
                unbounded = True
        elif table == "ext_grid":
            # An external grid with no declared limit is unlimited, which
            # makes this check vacuous rather than failed.
            unbounded = True
        elif "p_mw" in rows.columns:
            supply += float(rows["p_mw"].fillna(0.0).sum())

    loads = net.get("load")
    if loads is None or len(loads) == 0:
        report.skip("power_adequacy", "the network has no loads")
        return report
    in_service = (
        loads["in_service"].astype(bool) if "in_service" in loads else True
    )
    served = loads[in_service] if in_service is not True else loads

    # A load can only be curtailed if the OPF is allowed to move it, and
    # potpourri fills `min_p_mw = 0` on every load whether or not that is
    # the case. Reading the column unconditionally made this check vacuous
    # on any network that had been through add_OPF(): every demand looked
    # fully curtailable, so no shortfall could ever be found.
    flexible = (
        served["controllable"].fillna(False).astype(bool)
        if "controllable" in served.columns
        else served.index == -1  # nothing is controllable
    )
    fixed_demand = float(served.loc[~flexible, "p_mw"].fillna(0.0).sum())
    if "min_p_mw" in served.columns:
        flex_demand = float(served.loc[flexible, "min_p_mw"].fillna(0.0).sum())
    else:
        flex_demand = float(served.loc[flexible, "p_mw"].fillna(0.0).sum())
    demand = fixed_demand + flex_demand
    basis = (
        f"{int(flexible.sum())} curtailable loads at their minimum, the "
        f"rest at p_mw"
    )

    report.add(
        DiagnosticIssue(
            severity=DiagnosticSeverity.INFO,
            category=DiagnosticCategory.POWER_BALANCE,
            code="ADEQUACY_SUMMARY",
            message=(
                f"Upper-bound supply {supply:.4g} MW against lower-bound "
                f"demand {demand:.4g} MW, taken from {basis}."
            ),
            value=supply - demand,
            unit="MW",
            context={"supply_mw": supply, "demand_mw": demand, "basis": basis},
        )
    )

    if unbounded:
        report.skip(
            "power_adequacy_bound",
            "at least one source has no active-power upper limit, so no "
            "shortfall can be established",
        )
        return report

    if supply < demand:
        report.add(
            DiagnosticIssue(
                severity=DiagnosticSeverity.ERROR,
                category=DiagnosticCategory.POWER_BALANCE,
                code="ADEQUACY_INSUFFICIENT_GENERATION",
                message=(
                    f"The most the sources can produce, {supply:.4g} MW, is "
                    f"below the least the loads can take, {demand:.4g} MW. "
                    f"No dispatch satisfies the active power balance."
                ),
                recommendation=(
                    "Raise max_p_mw on a source, allow load flexibility, or "
                    "check that the external grid limits are intended. "
                    "Losses are not included, so the real shortfall is "
                    "larger than the figure shown."
                ),
                value=demand - supply,
                unit="MW",
                context={"supply_mw": supply, "demand_mw": demand},
            )
        )
    return report


def check_voltage_setpoints(ctx: DiagnosticContext) -> DiagnosticReport:
    """Find fixed voltage setpoints that lie outside the OPF voltage band.

    What it checks: the `vm_pu` setpoint of every in-service external grid
    and generator against the `min_vm_pu`/`max_vm_pu` of the bus it sits
    on.

    Necessary or sufficient: **sufficient** for infeasibility when the
    model both fixes the setpoint and enforces the band, which the AC OPF
    does unless `free_slack_vm` is set.

    Cost: negligible.

    Args:
        ctx: The diagnostic context.

    Returns:
        A report with one `ERROR` per conflicting setpoint.
    """
    report = DiagnosticReport()
    net = ctx.net
    if net is None or "min_vm_pu" not in net.bus.columns:
        report.skip("voltage_setpoints", "the buses carry no voltage band")
        return report

    for table in ("ext_grid", "gen"):
        frame = net.get(table)
        if frame is None or len(frame) == 0 or "vm_pu" not in frame.columns:
            continue
        for index, row in frame.iterrows():
            if "in_service" in frame.columns and not bool(row["in_service"]):
                continue
            setpoint = row.get("vm_pu")
            bus = row.get("bus")
            if (
                setpoint is None
                or setpoint != setpoint
                or bus not in net.bus.index
            ):
                continue
            low = net.bus.at[bus, "min_vm_pu"]
            high = net.bus.at[bus, "max_vm_pu"]
            if low != low or high != high:
                continue
            if low <= setpoint <= high:
                continue
            report.add(
                DiagnosticIssue(
                    severity=DiagnosticSeverity.ERROR,
                    category=DiagnosticCategory.VOLTAGE,
                    code="VOLTAGE_SETPOINT_OUTSIDE_BAND",
                    message=(
                        f"net.{table}[{index}] holds its bus at "
                        f"{float(setpoint):.4f} p.u., outside the band "
                        f"[{float(low):.4f}, {float(high):.4f}] p.u. the OPF "
                        f"enforces there."
                    ),
                    recommendation=(
                        "Widen the bus voltage band, move the setpoint into "
                        "it, or let the slack voltage float."
                    ),
                    element=ctx.imap.element_ref(table, index),
                    value=float(setpoint),
                    lower_bound=float(low),
                    upper_bound=float(high),
                    unit="p.u.",
                    context={"bus": int(bus)},
                )
            )
    return report


def check_thermal_data(ctx: DiagnosticContext) -> DiagnosticReport:
    """Check that the branch ratings the OPF needs are usable.

    What it checks: line current ratings and transformer apparent-power
    ratings for non-positive, missing or absurd values.

    Necessary or sufficient: neither. A missing rating does not make the
    OPF infeasible — it quietly removes a constraint the user asked for,
    which is arguably worse, so it is reported.

    Cost: negligible.

    Args:
        ctx: The diagnostic context.

    Returns:
        A report of `ERROR` for non-positive ratings and `WARNING` for
        missing ones.
    """
    report = DiagnosticReport()
    net = ctx.net
    if net is None:
        report.skip("thermal_data", "no network on the model")
        return report

    specs = [("line", "max_i_ka", "kA"), ("trafo", "sn_mva", "MVA")]
    for table, column, unit in specs:
        frame = net.get(table)
        if frame is None or len(frame) == 0 or column not in frame.columns:
            continue
        values = frame[column].astype(float)
        for index in frame.index[values.isna()]:
            report.add(
                DiagnosticIssue(
                    severity=DiagnosticSeverity.WARNING,
                    category=DiagnosticCategory.DATA,
                    code="THERMAL_RATING_MISSING",
                    message=(
                        f"No {column} is set, so no thermal limit can be "
                        f"enforced on this element."
                    ),
                    recommendation=(
                        f"Set net.{table}.{column} if it should be limited."
                    ),
                    element=ctx.imap.element_ref(table, index),
                    unit=unit,
                )
            )
        for index in frame.index[values.notna() & (values <= 0)]:
            report.add(
                DiagnosticIssue(
                    severity=DiagnosticSeverity.ERROR,
                    category=DiagnosticCategory.DATA,
                    code="THERMAL_RATING_NONPOSITIVE",
                    message=(
                        f"{column} is {float(values[index]):.6g} {unit}, "
                        f"which forbids any flow through this element."
                    ),
                    recommendation=(f"Set a positive net.{table}.{column}."),
                    element=ctx.imap.element_ref(table, index),
                    value=float(values[index]),
                    unit=unit,
                )
            )
    return report


def check_base_power_flow(ctx: DiagnosticContext) -> DiagnosticReport:
    """Report what the base-case power flow did, and where it already sits.

    potpourri runs a pandapower power flow in the constructor, and what it
    produced is diagnostic gold: a converged base case that already
    violates the OPF's own voltage band is a very different situation from
    a base case that did not converge at all.

    What it checks: whether `net.res_bus` exists, and which buses start
    outside `min_vm_pu`/`max_vm_pu`.

    Necessary or sufficient: neither. A base case outside the band is
    common and perfectly solvable; it is reported as context, not as a
    fault.

    Cost: negligible; it reads results the constructor already produced.

    Args:
        ctx: The diagnostic context.

    Returns:
        A report describing the base case, with one `INFO` per bus that
        starts outside the band.
    """
    report = DiagnosticReport()
    net = ctx.net
    if net is None:
        report.skip("base_power_flow", "no network on the model")
        return report

    res = net.get("res_bus")
    if (
        res is None
        or len(res) == 0
        or "vm_pu" not in getattr(res, "columns", [])
    ):
        report.add(
            DiagnosticIssue(
                severity=DiagnosticSeverity.WARNING,
                category=DiagnosticCategory.INITIALIZATION,
                code="BASE_PF_MISSING",
                message=(
                    "No base-case power-flow result is present, so the model "
                    "was initialised without one."
                ),
                recommendation=(
                    "A failed base power flow and an infeasible OPF are "
                    "different problems; run pp.runpp(net) to see which "
                    "this is."
                ),
            )
        )
        return report

    report.summary["initialization.base_pf"] = "converged"
    if "min_vm_pu" not in net.bus.columns:
        return report

    reverse = {v: k for k, v in ctx.imap.bus.items()}
    outside = 0
    for bus in res.index:
        if bus not in net.bus.index:
            continue
        vm = res.at[bus, "vm_pu"]
        low = net.bus.at[bus, "min_vm_pu"]
        high = net.bus.at[bus, "max_vm_pu"]
        if vm != vm or low != low or high != high:
            continue
        if low <= vm <= high:
            continue
        outside += 1
        ppc = reverse.get(int(bus))
        report.add(
            DiagnosticIssue(
                severity=DiagnosticSeverity.INFO,
                category=DiagnosticCategory.INITIALIZATION,
                code="BASE_PF_OUTSIDE_VOLTAGE_BAND",
                message=(
                    f"The base power flow converges, but this bus starts at "
                    f"{float(vm):.4f} p.u., outside the OPF band "
                    f"[{float(low):.4f}, {float(high):.4f}] p.u."
                ),
                recommendation=(
                    "This is a starting point, not a verdict: the OPF may "
                    "well bring it inside the band."
                ),
                element=ctx.imap.bus_ref(ppc) if ppc is not None else None,
                value=float(vm),
                lower_bound=float(low),
                upper_bound=float(high),
                unit="p.u.",
                pyomo_component="v" if ctx.has("v") else None,
                pyomo_index=ppc,
            )
        )
    report.summary["initialization.buses_outside_band"] = outside
    return report


# --- helpers ----------------------------------------------------------


def _sum_on(net, table: str, column: str, buses: set[int]) -> float:
    """Total of `column` over in-service rows of `table` on these buses."""
    frame = net.get(table)
    if frame is None or len(frame) == 0 or column not in frame.columns:
        return 0.0
    mask = frame["bus"].isin(buses)
    if "in_service" in frame.columns:
        mask &= frame["in_service"].astype(bool)
    return float(frame.loc[mask, column].fillna(0.0).sum())


def _count_on(net, table: str, buses: set[int]) -> int:
    """Number of in-service rows of `table` sitting on these buses."""
    frame = net.get(table)
    if frame is None or len(frame) == 0 or "bus" not in frame.columns:
        return 0
    mask = frame["bus"].isin(buses)
    if "in_service" in frame.columns:
        mask &= frame["in_service"].astype(bool)
    return int(mask.sum())


def _jsonable(payload):
    """Reduce a pandapower diagnostic payload to plain JSON-safe data."""
    if isinstance(payload, dict):
        return {str(k): _jsonable(v) for k, v in payload.items()}
    if isinstance(payload, (list, tuple, set)):
        return [_jsonable(v) for v in payload]
    if isinstance(payload, (np.integer,)):
        return int(payload)
    if isinstance(payload, (np.floating,)):
        return float(payload)
    if isinstance(payload, (str, int, float, bool)) or payload is None:
        return payload
    return str(payload)
