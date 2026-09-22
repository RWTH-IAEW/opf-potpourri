# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Decide which checks to run, run them, and collect the findings.

Three levels, chosen so the default is worth running every time:

| level | what it adds | cost |
|---|---|---|
| `basic` | network data, bounds, islands, adequacy, model size, | no solver, |
|         | solver verdict, solution violations | no power flow |
| `standard` | binding limits, plausibility, scaling, replay | one power flow |
| `deep` | structural singularity, Jacobian conditioning, | extra solver |
|        | feasibility relaxation | calls |

`standard` is the default. `basic` never runs a solver or a power flow, so
it is safe inside a loop. `deep` can take longer than the original solve
and is opt-in.

Every check is wrapped: one raising must not take the report down with it.
A diagnostic tool that fails because the thing it is diagnosing is broken
is not much use, so a crashing check is recorded as a skip and the rest
carry on.
"""

from __future__ import annotations

from typing import Any

from potpourri.diagnostics import infeasibility, model as model_checks
from potpourri.diagnostics import network, numerical, replay, solution, solver
from potpourri.diagnostics.context import DiagnosticContext
from potpourri.diagnostics.report import DiagnosticReport

#: Levels in increasing order of cost.
LEVELS = ("basic", "standard", "deep")


def diagnose(
    model_obj,
    level: str = "standard",
    *,
    include_solver: bool = True,
    numerical_checks: bool | None = None,
    power_system: bool = True,
    cross_check: bool | None = None,
    relaxation: bool | None = None,
    solver_name: str = "ipopt",
    tol: float = 1e-6,
    print_report: bool = False,
    **options: Any,
) -> DiagnosticReport:
    """Run the diagnostic suite against one potpourri model.

    Args:
        model_obj: A constructed potpourri model, e.g. an `ACOPF`. It need
            not have been solved, and `add_OPF()` need not have run.
        level: `"basic"`, `"standard"` or `"deep"`.
        include_solver: Read and interpret the solver's own verdict.
        numerical_checks: Run the scaling and conditioning checks.
            Defaults to on from `standard` upwards.
        power_system: Run the plausibility checks — voltage profile,
            branch loading, power balance.
        cross_check: Replay the solution through a pandapower power flow.
            Defaults to on from `standard` upwards.
        relaxation: Solve the elastic relaxation to quantify what would
            have to give. Defaults to on at `deep`, and only when the
            model looks infeasible.
        solver_name: Solver for the checks that need one.
        tol: Absolute tolerance for calling a constraint violated.
        print_report: Print the report before returning it.
        **options: Passed through to the checks, e.g. `max_binding`.

    Returns:
        A `DiagnosticReport`. Checks that could not run are listed in
        `report.skipped` with the reason.

    Raises:
        ValueError: If `level` is not one of `LEVELS`.
    """
    if level not in LEVELS:
        raise ValueError(f"level must be one of {LEVELS}, got {level!r}")

    depth = LEVELS.index(level)
    if numerical_checks is None:
        numerical_checks = depth >= 1
    if cross_check is None:
        cross_check = depth >= 1

    ctx = DiagnosticContext.from_model(model_obj, tol=tol, options=options)
    report = DiagnosticReport()
    report.summary["model.formulation"] = ctx.formulation
    _add_network_summary(ctx, report)

    # --- always, and never needing a solver ---------------------------
    _run(report, "network data", network.run_pandapower_diagnostic, ctx)
    _run(report, "bounds", network.check_bounds, ctx)
    _run(report, "islands", network.check_islands, ctx)
    _run(report, "power adequacy", network.check_power_adequacy, ctx)
    _run(report, "voltage setpoints", network.check_voltage_setpoints, ctx)
    _run(report, "thermal data", network.check_thermal_data, ctx)
    _run(report, "base power flow", network.check_base_power_flow, ctx)
    _run(report, "model size", model_checks.check_model_size, ctx)

    if include_solver:
        _run(report, "solver result", solver.check_solver_result, ctx)

    if ctx.has_solution:
        _run(
            report,
            "constraint violations",
            solution.check_constraint_violations,
            ctx,
        )
        _run(
            report,
            "infeasible points",
            infeasibility.check_infeasible_points,
            ctx,
        )

    # --- standard -----------------------------------------------------
    if depth >= 1:
        _run(
            report,
            "unused variables",
            model_checks.check_unconstrained_variables,
            ctx,
        )
        if ctx.has_solution:
            _run(
                report,
                "binding limits",
                solution.check_binding_constraints,
                ctx,
            )
        if power_system:
            _run(
                report, "voltage profile", solution.check_voltage_profile, ctx
            )
            _run(report, "branch loading", solution.check_branch_loading, ctx)
            _run(report, "power balance", solution.check_power_balance, ctx)
        if numerical_checks:
            _run(
                report, "variable values", numerical.check_variable_values, ctx
            )
            _run(
                report,
                "constraint scaling",
                numerical.check_constraint_scaling,
                ctx,
            )
        if cross_check:
            _run(report, "power-flow replay", replay.replay_power_flow, ctx)

    # --- deep ---------------------------------------------------------
    if depth >= 2:
        _run(
            report,
            "structural singularity",
            model_checks.check_structural_singularity,
            ctx,
        )
        if numerical_checks:
            _run(
                report, "jacobian", numerical.check_jacobian_conditioning, ctx
            )

    if relaxation is None:
        relaxation = depth >= 2 and _looks_infeasible(report)
    if relaxation:
        _run(
            report,
            "feasibility relaxation",
            infeasibility.relax_for_feasibility,
            ctx,
            solver=solver_name,
        )

    if print_report:
        print(report)
    return report


def _run(report: DiagnosticReport, label: str, check, ctx, **kwargs) -> None:
    """Run one check, folding its findings in and surviving a crash.

    Args:
        report: The report to extend.
        label: Human-readable name of the check, used in `skipped`.
        check: The check function, taking `ctx` and returning a report.
        ctx: The diagnostic context.
        **kwargs: Extra arguments for the check.

    Returns:
        None.
    """
    try:
        report.extend(check(ctx, **kwargs))
    except Exception as exc:  # noqa: BLE001 - a check must never abort the run
        report.skip(
            label, f"the check itself raised {type(exc).__name__}: {exc}"
        )


def _add_network_summary(
    ctx: DiagnosticContext, report: DiagnosticReport
) -> None:
    """Put the network's size at the top of the report."""
    net = ctx.net
    if net is None:
        return
    for table, label in (
        ("bus", "buses"),
        ("line", "lines"),
        ("trafo", "transformers"),
        ("load", "loads"),
        ("sgen", "static generators"),
        ("storage", "storage units"),
    ):
        frame = net.get(table)
        if frame is not None and len(frame):
            report.summary[f"network.{label}"] = len(frame)
    if ctx.is_multi_period:
        report.summary["network.time steps"] = len(ctx.time_steps)


def _looks_infeasible(report: DiagnosticReport) -> bool:
    """Whether anything so far suggests the model did not solve.

    Args:
        report: The report gathered so far.

    Returns:
        True when the solver reported an infeasible or locally infeasible
        termination, or a constraint violation was found.
    """
    codes = {issue.code for issue in report.issues}
    return bool(
        codes
        & {
            "SOLVER_INFEASIBLE",
            "SOLVER_LOCALLY_INFEASIBLE",
            "PYOMO_INFEASIBLE_CONSTRAINT",
            "SOLUTION_CONSTRAINT_VIOLATED",
        }
    )
