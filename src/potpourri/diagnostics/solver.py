# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Turn a solver's verdict into something actionable.

The distinction that matters most: **"locally infeasible" is not
"infeasible"**. IPOPT solves a nonconvex problem by local search, so it
reports where *it* got stuck, not a proof about the problem. Translating
every non-optimal result into "the OPF is infeasible" is the single most
misleading thing a diagnostic tool can do, and this module does not do it.

Each termination condition maps to a plain-English reading, a severity,
and what to try next. Solver-specific detail is read where the interface
exposes it and skipped where it does not.
"""

from __future__ import annotations

from potpourri.diagnostics.context import DiagnosticContext
from potpourri.diagnostics.report import (
    DiagnosticCategory,
    DiagnosticIssue,
    DiagnosticReport,
    DiagnosticSeverity,
)

#: Termination condition -> (code, severity, reading, what to try next).
#: Keys are matched case-insensitively against the string form of Pyomo's
#: `TerminationCondition`.
_TERMINATION: dict[str, tuple[str, DiagnosticSeverity, str, str]] = {
    "optimal": (
        "SOLVER_OPTIMAL",
        DiagnosticSeverity.INFO,
        "The solver converged to a point satisfying its optimality "
        "tolerances.",
        "For a nonconvex AC OPF this is a local optimum; it is not a proof "
        "of global optimality.",
    ),
    "locallyoptimal": (
        "SOLVER_LOCALLY_OPTIMAL",
        DiagnosticSeverity.INFO,
        "The solver converged to a local optimum.",
        "Another starting point may reach a different one.",
    ),
    "feasible": (
        "SOLVER_FEASIBLE",
        DiagnosticSeverity.INFO,
        "The solver found a feasible point but did not prove optimality.",
        "Check the iteration and time limits.",
    ),
    "infeasible": (
        "SOLVER_INFEASIBLE",
        DiagnosticSeverity.ERROR,
        "The solver reports the problem as infeasible.",
        "Start with the bound and adequacy checks above; they find "
        "contradictions that need no solver.",
    ),
    "locallyinfeasible": (
        "SOLVER_LOCALLY_INFEASIBLE",
        DiagnosticSeverity.ERROR,
        "The solver converged to a locally infeasible point. For a "
        "nonconvex problem this is not a proof that the OPF is infeasible.",
        "Try a different starting point, then look at the pre-solve checks "
        "and the feasibility relaxation to see which limits are in conflict.",
    ),
    "maxiterations": (
        "SOLVER_ITERATION_LIMIT",
        DiagnosticSeverity.WARNING,
        "The solver hit its iteration limit before converging.",
        "Raise max_iter, improve the starting point, or look at the "
        "scaling diagnostics: slow convergence is often a conditioning "
        "problem rather than a feasibility one.",
    ),
    "maxtimelimit": (
        "SOLVER_TIME_LIMIT",
        DiagnosticSeverity.WARNING,
        "The solver ran out of time.",
        "Raise the time limit, or reduce the horizon on a multi-period model.",
    ),
    "unbounded": (
        "SOLVER_UNBOUNDED",
        DiagnosticSeverity.ERROR,
        "The objective is unbounded.",
        "Usually a missing limit on a source, or an objective with the "
        "wrong sign.",
    ),
    "error": (
        "SOLVER_ERROR",
        DiagnosticSeverity.ERROR,
        "The solver stopped with an error rather than a result.",
        "Check the solver log; an evaluation error in the model is the "
        "common cause.",
    ),
    "intermediatenoninteger": (
        "SOLVER_NONINTEGER",
        DiagnosticSeverity.WARNING,
        "The solver stopped at a point whose integer variables are not "
        "integral.",
        "A MINLP needs a MINLP solver; check that the right one is selected.",
    ),
}


def check_solver_result(ctx: DiagnosticContext) -> DiagnosticReport:
    """Read the solver's own verdict and say what it means.

    What it checks: the termination condition, solver status and whatever
    numeric detail the interface exposed.

    Necessary or sufficient: neither. This reports what the solver said;
    it does not verify it. `check_constraint_violations` does that
    independently.

    Cost: negligible.

    Required state: `solve()` must have been called. Without a result the
    check says so and stops.

    Args:
        ctx: The diagnostic context.

    Returns:
        A report with one issue for the termination condition plus the
        solver's own numbers in `summary`.
    """
    report = DiagnosticReport()
    results = getattr(ctx.model_obj, "results", None)
    if results is None:
        report.skip("solver_result", "solve() has not been called")
        return report

    solver = _first(results, "solver")
    if solver is None:
        report.skip(
            "solver_result", "the result object carries no solver section"
        )
        return report

    termination = str(getattr(solver, "termination_condition", "unknown"))
    status = str(getattr(solver, "status", "unknown"))
    report.summary["solver.termination"] = termination
    report.summary["solver.status"] = status
    for attribute, label in (
        ("time", "solver.wall time"),
        ("iterations", "solver.iterations"),
        ("message", "solver.message"),
    ):
        value = getattr(solver, attribute, None)
        if value is not None and str(value) not in ("", "None"):
            report.summary[label] = (
                f"{value:.2f} s" if attribute == "time" else str(value)
            )

    key = termination.split(".")[-1].replace("_", "").lower()
    entry = _TERMINATION.get(key)
    if entry is None:
        report.add(
            DiagnosticIssue(
                severity=DiagnosticSeverity.WARNING,
                category=DiagnosticCategory.SOLVER,
                code="SOLVER_UNRECOGNISED_TERMINATION",
                message=(
                    f"The solver terminated with '{termination}', which this "
                    f"tool has no reading for."
                ),
                recommendation="Check the solver's own documentation.",
                context={"termination": termination, "status": status},
            )
        )
        return report

    code, severity, reading, advice = entry
    report.add(
        DiagnosticIssue(
            severity=severity,
            category=DiagnosticCategory.SOLVER,
            code=code,
            message=reading,
            recommendation=advice,
            context={"termination": termination, "status": status},
        )
    )
    return report


def _first(results, section):
    """Fetch a Pyomo results section without raising when it is absent.

    Pyomo exposes `results.solver` as a list-like that raises on an empty
    result, so the access is guarded rather than assumed.

    Args:
        results: A Pyomo results object.
        section: Section name, e.g. `"solver"`.

    Returns:
        The first entry of that section, or `None`.
    """
    try:
        holder = getattr(results, section)
    except (AttributeError, KeyError, IndexError):
        return None
    if holder is None:
        return None
    try:
        return holder[0]
    except (TypeError, KeyError, IndexError):
        return holder
