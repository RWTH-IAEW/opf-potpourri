# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Why an infeasible OPF is infeasible, as far as that can be established.

Two routes, with very different strength of evidence.

**Residuals** (`check_infeasible_points`) evaluate the model at the point
the solver stopped at and report what does not hold there. Cheap, exact
about that point, and says nothing about the problem as a whole.

**Elastic relaxation** (`relax_for_feasibility`) is the useful one. It
clones the model, adds a non-negative slack to each limit that can
physically be relaxed, and minimises a weighted sum of those slacks. The
answer reads as "feasibility could be restored by allowing bus 17 to rise
0.006 p.u., or by 3 % more current on line 4". Three caveats travel with
every such answer and are repeated in the issues themselves:

* it is **one** relaxation, not the only one;
* the weights decide which one you get;
* for a nonconvex AC OPF the solver finds a local answer.

The user's model is never modified. The relaxation runs on a clone.

A classical LP/MILP IIS is deliberately not offered for the AC OPF: an IIS
is defined for linear systems, and presenting one for a nonconvex problem
would imply a guarantee that does not exist. Pyomo's
`compute_infeasibility_explanation` is used where the formulation and
solver make it meaningful.
"""

from __future__ import annotations

import copy

from potpourri.diagnostics.context import DiagnosticContext
from potpourri.diagnostics.metadata import constraint_meta
from potpourri.diagnostics.report import (
    DiagnosticCategory,
    DiagnosticIssue,
    DiagnosticReport,
    DiagnosticSeverity,
)
from potpourri.diagnostics.solution import _element_and_time

#: Constraint families the relaxation is allowed to soften, with the
#: weight applied to their slack. Weights are relative, and chosen so that
#: relaxing a physical rating is "cheaper" than breaking a power balance:
#: a balance residual is not something an operator can grant.
_RELAXABLE: dict[str, float] = {
    "v_pyo": 1.0,
    "v_constraint": 1.0,
    "line_lim_from": 1.0,
    "line_lim_to": 1.0,
    "transf_lim1": 1.0,
    "transf_lim2": 1.0,
    "QsG_pyo": 2.0,
    "QG_pyo": 2.0,
    "PsG_Constraint": 2.0,
    "PG_Constraint": 2.0,
    "QD_pyo": 5.0,
    "PD_Constraint": 5.0,
    "stor_soc_bounds": 5.0,
}


def check_infeasible_points(ctx: DiagnosticContext) -> DiagnosticReport:
    """Report what does not hold at the point the solver stopped at.

    What it checks: Pyomo's own `find_infeasible_constraints` and
    `find_infeasible_bounds` against the loaded variable values.

    Necessary or sufficient: exact about that point only. A solver that
    stopped at a locally infeasible point leaves the violated constraints
    behind, and those are worth reading — but they are where this solve
    got stuck, not a property of the problem.

    Cost: cheap; no solver call.

    Args:
        ctx: The diagnostic context.

    Returns:
        A report with one `ERROR` per violated constraint and per variable
        outside its bounds, mapped back to network objects.
    """
    report = DiagnosticReport()
    if ctx.model is None:
        report.skip("infeasible_points", "the Pyomo model has not been built")
        return report

    try:
        from pyomo.util.infeasible import (
            find_infeasible_bounds,
            find_infeasible_constraints,
        )
    except ImportError as exc:
        report.skip(
            "infeasible_points", f"pyomo.util.infeasible unavailable: {exc}"
        )
        return report

    try:
        constraints = list(find_infeasible_constraints(ctx.model, tol=ctx.tol))
    except Exception as exc:
        report.skip("infeasible_points", f"raised {type(exc).__name__}: {exc}")
        return report

    for entry in constraints:
        con = entry[0] if isinstance(entry, tuple) else entry
        parent = con.parent_component()
        meta = constraint_meta(parent.local_name)
        index = con.index()
        element, time_step = _element_and_time(ctx, parent.local_name, index)
        body = ctx.value(con.body)
        low = ctx.value(con.lower)
        high = ctx.value(con.upper)
        report.add(
            DiagnosticIssue(
                severity=DiagnosticSeverity.ERROR,
                category=meta.category
                if meta
                else DiagnosticCategory.INFEASIBILITY,
                code="PYOMO_INFEASIBLE_CONSTRAINT",
                message=(
                    f"{meta.description if meta else parent.local_name} does "
                    f"not hold at the point the solver stopped at."
                ),
                recommendation=(
                    "This is where this solve ended, not a proof about the "
                    "problem. Run diagnose(level='deep') for a relaxation "
                    "that quantifies what would have to give."
                ),
                element=element,
                pyomo_component=parent.local_name,
                pyomo_index=index,
                value=body,
                lower_bound=low,
                upper_bound=high,
                time_step=time_step,
            )
        )

    try:
        bounds = list(find_infeasible_bounds(ctx.model, tol=ctx.tol))
    except Exception:
        bounds = []
    if bounds:
        report.add(
            DiagnosticIssue(
                severity=DiagnosticSeverity.WARNING,
                category=DiagnosticCategory.BOUNDS,
                code="PYOMO_VARIABLE_OUTSIDE_BOUNDS",
                message=(
                    f"{len(bounds)} variables lie outside their own bounds at "
                    f"the point the solver stopped at."
                ),
                context={
                    "examples": [
                        (entry[0] if isinstance(entry, tuple) else entry).name
                        for entry in bounds[:10]
                    ]
                },
            )
        )
    return report


def relax_for_feasibility(
    ctx: DiagnosticContext, solver: str = "ipopt", max_reported: int = 10
) -> DiagnosticReport:
    """Ask what would have to give for the OPF to become feasible.

    Clones the model, replaces each relaxable limit with a softened
    version carrying a non-negative slack, and minimises the weighted sum
    of those slacks. Where the result puts slack, the corresponding limit
    is one that has to move.

    What it checks: the families in `_RELAXABLE` — voltage bands, branch
    thermal limits, generator and sgen power limits, load flexibility and
    storage energy limits. Nodal balance is **not** relaxed: allowing
    power to appear from nowhere would produce an answer that means
    nothing physically.

    Necessary or sufficient: **neither, and this matters.** The result is
    one feasible relaxation out of many. Different weights give different
    answers, and on a nonconvex AC OPF the solver finds a local one. It
    identifies limits that are in tension; it does not prove that any one
    of them is the cause.

    Cost: expensive. One additional solve. Never run by default.

    Args:
        ctx: The diagnostic context.
        solver: Solver name for the relaxed problem.
        max_reported: How many relaxed limits to report, worst first.

    Returns:
        A report with one `INFO` per limit the relaxation had to soften,
        each carrying the amount in engineering units where the family's
        units are known.
    """
    report = DiagnosticReport()
    if ctx.model is None:
        report.skip("relaxation", "the Pyomo model has not been built")
        return report

    import pyomo.environ as pyo

    try:
        clone = ctx.model.clone()
    except Exception as exc:
        report.skip("relaxation", f"the model could not be cloned: {exc}")
        return report

    slacks = _soften(clone, pyo)
    if not slacks:
        report.skip(
            "relaxation",
            "none of the relaxable constraint families are present in this "
            "model",
        )
        return report

    for objective in clone.component_data_objects(pyo.Objective, active=True):
        objective.deactivate()
    clone._diagnostic_objective = pyo.Objective(
        expr=sum(weight * var for var, weight in slacks.values()),
        sense=pyo.minimize,
    )

    try:
        opt = pyo.SolverFactory(solver)
        if not opt.available(exception_flag=False):
            report.skip("relaxation", f"solver '{solver}' is not available")
            return report
        result = opt.solve(clone, load_solutions=True, tee=False)
    except Exception as exc:
        report.skip("relaxation", f"the relaxed solve failed: {exc}")
        return report

    termination = str(result.solver.termination_condition)
    if "optimal" not in termination.lower():
        report.add(
            DiagnosticIssue(
                severity=DiagnosticSeverity.WARNING,
                category=DiagnosticCategory.INFEASIBILITY,
                code="RELAXATION_NOT_SOLVED",
                message=(
                    f"The relaxed problem itself terminated '{termination}', "
                    f"so no relaxation could be quantified."
                ),
                recommendation=(
                    "With every relaxable limit softened and the problem "
                    "still not solving, look at the network and structural "
                    "checks instead."
                ),
            )
        )
        return report

    active = []
    for key, (var, _weight) in slacks.items():
        amount = ctx.value(var, 0.0) or 0.0
        if amount > ctx.tol:
            active.append((amount, key))
    active.sort(reverse=True)

    if not active:
        report.add(
            DiagnosticIssue(
                severity=DiagnosticSeverity.INFO,
                category=DiagnosticCategory.INFEASIBILITY,
                code="RELAXATION_NO_SLACK_NEEDED",
                message=(
                    "The relaxed problem solved with every slack at zero, so "
                    "the relaxable limits are not what is blocking it."
                ),
                recommendation=(
                    "Look at the nodal balance, the network topology and the "
                    "numerical checks; those are not relaxed here."
                ),
            )
        )
        return report

    report.add(
        DiagnosticIssue(
            severity=DiagnosticSeverity.INFO,
            category=DiagnosticCategory.INFEASIBILITY,
            code="RELAXATION_SUMMARY",
            message=(
                f"Feasibility could be restored by relaxing {len(active)} "
                f"limits. This is one such relaxation, not the only one."
            ),
            recommendation=(
                "The weights chose this answer; different weights give a "
                "different set. Treat these as limits in tension, not as "
                "the cause."
            ),
            context={"relaxed_count": len(active)},
        )
    )

    for amount, (component, index) in active[:max_reported]:
        meta = constraint_meta(component)
        element, time_step = _element_and_time(ctx, component, index)
        report.add(
            DiagnosticIssue(
                severity=DiagnosticSeverity.INFO,
                category=meta.category
                if meta
                else DiagnosticCategory.INFEASIBILITY,
                code="RELAXATION_LIMIT_RELAXED",
                message=(
                    f"{meta.description if meta else component} would have to "
                    f"give by {amount:.6g} in the model's own units."
                ),
                recommendation=(
                    "One of several possible relaxations; the others are "
                    "equally valid."
                ),
                element=element,
                pyomo_component=component,
                pyomo_index=index,
                violation=amount,
                time_step=time_step,
            )
        )
    return report


def _soften(clone, pyo) -> dict:
    """Add a slack variable to every relaxable inequality on the clone.

    Args:
        clone: A cloned Pyomo model, modified in place.
        pyo: The `pyomo.environ` module, passed in to keep the import at
            the caller.

    Returns:
        `{(component_name, index): (slack_var, weight)}` for every
        constraint that was softened. Equalities are left alone.
    """
    slacks: dict = {}
    block = pyo.Block(concrete=True)
    clone._diagnostic_slacks = block
    entries = []

    for con in list(clone.component_data_objects(pyo.Constraint, active=True)):
        name = con.parent_component().local_name
        weight = _RELAXABLE.get(name)
        if weight is None or con.equality:
            continue
        entries.append((con, name, weight))

    block.index = pyo.Set(initialize=range(len(entries)))
    block.slack = pyo.Var(
        block.index, domain=pyo.NonNegativeReals, initialize=0.0
    )

    for position, (con, name, weight) in enumerate(entries):
        var = block.slack[position]
        low, high = con.lower, con.upper
        body = con.body
        con.deactivate()
        relaxed = pyo.Constraint(
            expr=_relaxed_expression(body, low, high, var)
        )
        setattr(block, f"relaxed_{position}", relaxed)
        slacks[(name, con.index())] = (var, weight)
    return slacks


def _relaxed_expression(body, low, high, slack):
    """Build the softened form of one inequality.

    Args:
        body: The constraint body expression.
        low: Lower bound, or `None`.
        high: Upper bound, or `None`.
        slack: The non-negative slack variable to widen the bounds by.

    Returns:
        A Pyomo relational expression with the bounds widened by `slack`
        on whichever side exists.
    """
    if low is not None and high is not None:
        return (low - slack, body, high + slack)
    if high is not None:
        return body <= high + slack
    return body >= low - slack


def explain_with_pyomo(ctx: DiagnosticContext, solver: str = "ipopt"):
    """Run Pyomo's own infeasibility explanation, where it is meaningful.

    Pyomo's `compute_infeasibility_explanation` builds a minimal
    intractable subsystem. It needs a solver and many solves, so it is
    never run by default.

    Necessary or sufficient: the subsystem it returns is jointly
    infeasible, which is a real statement — but for a nonconvex AC OPF
    that statement is local, like everything else the solver reports.

    Cost: very expensive; many solver calls.

    Args:
        ctx: The diagnostic context.
        solver: Solver to use.

    Returns:
        A report with the explanation as a single `INFO`, or a recorded
        skip naming what was missing.
    """
    report = DiagnosticReport()
    if ctx.model is None:
        report.skip("pyomo_explanation", "the Pyomo model has not been built")
        return report
    try:
        from pyomo.contrib.iis import compute_infeasibility_explanation
    except ImportError as exc:
        report.skip(
            "pyomo_explanation", f"pyomo.contrib.iis unavailable: {exc}"
        )
        return report

    import io
    import logging

    stream = io.StringIO()
    handler = logging.StreamHandler(stream)
    logger = logging.getLogger("pyomo.contrib.iis")
    logger.addHandler(handler)
    try:
        compute_infeasibility_explanation(
            copy.deepcopy(ctx.model), solver=solver
        )
    except Exception as exc:
        report.skip("pyomo_explanation", f"raised {type(exc).__name__}: {exc}")
        return report
    finally:
        logger.removeHandler(handler)

    text = stream.getvalue().strip()
    if not text:
        report.skip("pyomo_explanation", "the routine produced no explanation")
        return report
    report.add(
        DiagnosticIssue(
            severity=DiagnosticSeverity.INFO,
            category=DiagnosticCategory.INFEASIBILITY,
            code="PYOMO_INFEASIBILITY_EXPLANATION",
            message="Pyomo identified a jointly infeasible subsystem.",
            recommendation=(
                "These constraints cannot all hold together. For a "
                "nonconvex model the statement is local."
            ),
            context={"explanation": text[:4000]},
        )
    )
    return report
