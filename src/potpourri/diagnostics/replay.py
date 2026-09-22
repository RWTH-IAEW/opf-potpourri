# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Check the OPF solution against an independent pandapower power flow.

The OPF says the network is in a certain state. pandapower's power flow is
an independent implementation of the same physics, so running it on the
optimised dispatch and comparing is the strongest end-to-end check
available: it catches sign errors, indexing mistakes, a missing element
model, a transformer tap applied the wrong way round, and bugs in the
result mapping — none of which the optimiser itself can detect, because
they are baked into the equations it is solving.

**What comparison means depends on the formulation.**

For an AC OPF the two should agree to solver tolerance, and a disagreement
is a defect somewhere. For DC and LPAC the model is an approximation on
purpose, so the AC power flow will differ; the comparison is then a
plausibility check on the size of the approximation error, and is labelled
as such rather than reported as a mismatch.
"""

from __future__ import annotations

import copy
import io
from contextlib import redirect_stderr, redirect_stdout

from potpourri.diagnostics.context import DiagnosticContext
from potpourri.diagnostics.report import (
    DiagnosticCategory,
    DiagnosticIssue,
    DiagnosticReport,
    DiagnosticSeverity,
)

#: Above this voltage difference an AC replay counts as disagreeing, in
#: per unit. Comfortably above solver tolerance, well below anything a
#: user would consider the same operating point.
_AC_VOLTAGE_TOL = 1e-4

#: The same for branch loading, in percentage points.
_AC_LOADING_TOL = 0.5


def replay_power_flow(ctx: DiagnosticContext) -> DiagnosticReport:
    """Re-solve the optimised dispatch with pandapower and compare.

    What it checks: bus voltage magnitude and branch loading from a
    pandapower AC power flow run on the network the OPF wrote its results
    into, against those results.

    Necessary or sufficient: neither, but a disagreement on an AC model is
    strong evidence of a defect rather than of a bad network.

    Cost: moderate. One power flow, no optimisation.

    Required state: the solution must have been written to `net.res_*`,
    which `solve(to_net=True)` does.

    Limitations: only the dispatch is replayed. Controls the OPF chose and
    pandapower cannot represent identically — a continuous transformer tap
    in particular — make the comparison approximate, and the report says
    so when such a control is present.

    Args:
        ctx: The diagnostic context.

    Returns:
        A report with the largest deviations in `summary`, and an issue
        whose severity depends on the formulation: `WARNING` for an AC
        mismatch, `INFO` for an approximate formulation's expected
        divergence.
    """
    report = DiagnosticReport()
    net = ctx.net
    if net is None:
        report.skip("replay", "no network on the model")
        return report
    res_bus = net.get("res_bus")
    if res_bus is None or len(res_bus) == 0 or "vm_pu" not in res_bus.columns:
        report.skip(
            "replay", "no OPF results have been written to the network"
        )
        return report
    if ctx.is_multi_period:
        report.skip(
            "replay",
            "the network holds one time step only, so a multi-period "
            "solution cannot be replayed as a whole; use map_to_net(t) and "
            "replay a single step",
        )
        return report

    import pandapower as pp

    approximate = (
        not ctx.has_reactive or "LPAC" in type(ctx.model_obj).__name__
    )
    caveats = []
    if ctx.has("Tap") and _tap_is_free(ctx):
        caveats.append(
            "the OPF optimised a continuous transformer tap, which "
            "pandapower rounds to a discrete position"
        )

    clone = copy.deepcopy(net)
    _apply_dispatch(clone, net)

    # Warm-started from the OPF result first, since that is the point we
    # want to confirm; a flat start is the fallback, because a warm start
    # that fails says more about the starting point than about the
    # dispatch.
    buf = io.StringIO()
    try:
        try:
            with redirect_stdout(buf), redirect_stderr(buf):
                pp.runpp(clone, init="results", calculate_voltage_angles=True)
        except Exception:
            with redirect_stdout(buf), redirect_stderr(buf):
                pp.runpp(clone, init="flat", calculate_voltage_angles=True)
    except Exception as exc:
        report.add(
            DiagnosticIssue(
                severity=DiagnosticSeverity.WARNING,
                category=DiagnosticCategory.CROSS_CHECK,
                code="REPLAY_POWER_FLOW_FAILED",
                message=(
                    f"A pandapower power flow on the optimised dispatch did "
                    f"not converge ({type(exc).__name__})."
                ),
                recommendation=(
                    "An OPF solution that no power flow reproduces is worth "
                    "investigating; it may indicate the dispatch was not "
                    "mapped back correctly."
                ),
            )
        )
        return report

    max_dv, worst_bus = _largest_difference(
        res_bus["vm_pu"], clone.res_bus["vm_pu"]
    )
    report.summary["cross-check.max |dV|"] = f"{max_dv:.3e} p.u."

    max_dl, worst_branch = 0.0, None
    if "loading_percent" in getattr(net.get("res_line"), "columns", []):
        max_dl, worst_branch = _largest_difference(
            net.res_line["loading_percent"], clone.res_line["loading_percent"]
        )
        report.summary["cross-check.max |dloading|"] = f"{max_dl:.3e} %"

    if approximate:
        report.add(
            DiagnosticIssue(
                severity=DiagnosticSeverity.INFO,
                category=DiagnosticCategory.CROSS_CHECK,
                code="REPLAY_APPROXIMATE_FORMULATION",
                message=(
                    f"An AC power flow on this dispatch differs by at most "
                    f"{max_dv:.3e} p.u. in voltage. This formulation is an "
                    f"approximation, so a difference is expected and is a "
                    f"measure of the approximation, not an error."
                ),
                value=max_dv,
                unit="p.u.",
                element=ctx.imap.element_ref("bus", worst_bus)
                if worst_bus is not None
                else None,
            )
        )
        return report

    if max_dv > _AC_VOLTAGE_TOL or max_dl > _AC_LOADING_TOL:
        report.add(
            DiagnosticIssue(
                severity=DiagnosticSeverity.WARNING,
                category=DiagnosticCategory.CROSS_CHECK,
                code="REPLAY_MISMATCH",
                message=(
                    f"The OPF solution and an independent pandapower power "
                    f"flow disagree by up to {max_dv:.3e} p.u. in voltage and "
                    f"{max_dl:.3e} percentage points in loading."
                ),
                recommendation=(
                    "On an AC model the two should agree to solver "
                    "tolerance. A gap this size usually means a modelling or "
                    "result-mapping difference rather than a bad network."
                    + (
                        " Note that " + "; ".join(caveats) + "."
                        if caveats
                        else ""
                    )
                ),
                element=ctx.imap.element_ref("bus", worst_bus)
                if worst_bus is not None
                else None,
                value=max_dv,
                unit="p.u.",
                context={
                    "max_voltage_difference_pu": max_dv,
                    "max_loading_difference_percent": max_dl,
                    "worst_branch": worst_branch,
                    "caveats": caveats,
                },
            )
        )
    else:
        report.add(
            DiagnosticIssue(
                severity=DiagnosticSeverity.INFO,
                category=DiagnosticCategory.CROSS_CHECK,
                code="REPLAY_AGREES",
                message=(
                    f"An independent pandapower power flow reproduces the "
                    f"OPF solution to {max_dv:.3e} p.u. in voltage."
                ),
                value=max_dv,
                unit="p.u.",
            )
        )
    return report


def _apply_dispatch(clone, source) -> None:
    """Copy the optimised injections onto the replay network.

    Args:
        clone: The network the power flow will run on, modified in place.
        source: The network carrying the OPF results.

    Returns:
        None.
    """
    for table in ("sgen", "gen", "storage", "load"):
        res = source.get(f"res_{table}")
        frame = clone.get(table)
        if res is None or frame is None or len(frame) == 0:
            continue
        for column in ("p_mw", "q_mvar"):
            if column in res.columns and column in frame.columns:
                shared = frame.index.intersection(res.index)
                frame.loc[shared, column] = res.loc[shared, column].values
    res_bus = source.get("res_bus")
    if res_bus is not None and "vm_pu" in res_bus.columns:
        for table in ("ext_grid", "gen"):
            frame = clone.get(table)
            if (
                frame is None
                or len(frame) == 0
                or "vm_pu" not in frame.columns
            ):
                continue
            for index, row in frame.iterrows():
                bus = row.get("bus")
                if bus in res_bus.index:
                    frame.at[index, "vm_pu"] = float(res_bus.at[bus, "vm_pu"])


def _largest_difference(left, right):
    """The largest absolute difference between two aligned series.

    Args:
        left: One series, indexed by element.
        right: The other.

    Returns:
        `(magnitude, index)` of the worst difference, or `(0.0, None)`
        when the series do not overlap.
    """
    if left is None or right is None:
        return 0.0, None
    shared = left.index.intersection(right.index)
    if not len(shared):
        return 0.0, None
    delta = (left.loc[shared] - right.loc[shared]).abs().dropna()
    if delta.empty:
        return 0.0, None
    return float(delta.max()), int(delta.idxmax())


def _tap_is_free(ctx: DiagnosticContext) -> bool:
    """Whether the model left any transformer tap as a free variable."""
    tap = getattr(ctx.model, "Tap", None)
    if tap is None:
        return False
    try:
        return any(not var.fixed for var in tap.values())
    except (AttributeError, TypeError):
        return False
