# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Structural checks on the Pyomo model itself.

These look at the shape of the optimisation problem rather than at its
numbers: how big it is, whether any variable is free but unconstrained,
and whether the equality system is structurally singular.

On degrees of freedom: an optimisation model is *supposed* to have them.
The "degrees of freedom must be zero" rule applies to square equation
systems, not to an OPF, whose whole purpose is to choose among feasible
dispatches. The count is reported as structural information, never as a
fault.
"""

from __future__ import annotations

from potpourri.diagnostics.context import DiagnosticContext
from potpourri.diagnostics.metadata import constraint_meta
from potpourri.diagnostics.report import (
    DiagnosticCategory,
    DiagnosticIssue,
    DiagnosticReport,
    DiagnosticSeverity,
)


def check_model_size(ctx: DiagnosticContext) -> DiagnosticReport:
    """Count what the model contains.

    What it checks: active and free variables, equality and inequality
    constraints, objectives and integer variables, using Pyomo's own
    `build_model_size_report` where it applies and direct component walks
    for the split Pyomo does not give.

    Necessary or sufficient: neither; this is descriptive.

    Cost: cheap, one pass over the model's components.

    Args:
        ctx: The diagnostic context.

    Returns:
        A report carrying the counts in `summary` and no issues, except an
        `ERROR` when the model has no objective and an OPF has been built.
    """
    report = DiagnosticReport()
    if ctx.model is None:
        report.skip("model_size", "the Pyomo model has not been built")
        return report

    import pyomo.environ as pyo

    free = fixed = 0
    integer = 0
    for var in ctx.model.component_data_objects(pyo.Var, active=True):
        if var.fixed:
            fixed += 1
        else:
            free += 1
            if var.is_integer() or var.is_binary():
                integer += 1

    equalities = inequalities = 0
    for con in ctx.model.component_data_objects(pyo.Constraint, active=True):
        if con.equality:
            equalities += 1
        else:
            inequalities += 1

    objectives = [
        obj.name
        for obj in ctx.model.component_data_objects(pyo.Objective, active=True)
    ]

    report.summary.update(
        {
            "model.formulation": ctx.formulation,
            "model.free variables": free,
            "model.fixed variables": fixed,
            "model.integer variables": integer,
            "model.equalities": equalities,
            "model.inequalities": inequalities,
            "model.objectives": len(objectives),
            "model.degrees of freedom": free - equalities,
        }
    )

    if ctx.is_opf and not objectives:
        report.add(
            DiagnosticIssue(
                severity=DiagnosticSeverity.ERROR,
                category=DiagnosticCategory.STRUCTURE,
                code="MODEL_NO_OBJECTIVE",
                message=(
                    "The model has operating limits but no active objective, "
                    "so there is nothing to optimise."
                ),
                recommendation=(
                    "Call one of the add_*_objective() methods before solve()."
                ),
            )
        )
    elif len(objectives) > 1:
        report.add(
            DiagnosticIssue(
                severity=DiagnosticSeverity.ERROR,
                category=DiagnosticCategory.STRUCTURE,
                code="MODEL_MULTIPLE_OBJECTIVES",
                message=(
                    f"{len(objectives)} objectives are active at once: "
                    f"{', '.join(objectives)}. Most solvers reject this."
                ),
                recommendation="Deactivate all but one.",
                context={"objectives": objectives},
            )
        )

    if free - equalities < 0:
        report.add(
            DiagnosticIssue(
                severity=DiagnosticSeverity.WARNING,
                category=DiagnosticCategory.STRUCTURE,
                code="MODEL_OVERDETERMINED",
                message=(
                    f"There are {equalities} equality constraints for {free} "
                    f"free variables, so the equality system is "
                    f"overdetermined on its own."
                ),
                recommendation=(
                    "Some equalities are likely redundant; that is not "
                    "automatically an error, but it is worth checking."
                ),
                context={"free_variables": free, "equalities": equalities},
            )
        )
    return report


def check_unconstrained_variables(ctx: DiagnosticContext) -> DiagnosticReport:
    """Find free variables that no active constraint or objective touches.

    What it checks: every unfixed variable, against the set of variables
    appearing in active constraints and objectives.

    Necessary or sufficient: neither. A variable that appears nowhere is
    harmless to the optimum but usually means a constraint was not added,
    so it is worth a warning.

    Cost: moderate; walks every expression in the model once.

    False positives: auxiliary variables that a formulation declares and
    deliberately leaves loose would show up here.

    Args:
        ctx: The diagnostic context.

    Returns:
        A report with one `WARNING` summarising the unused variables,
        naming a few, rather than one issue per variable.
    """
    report = DiagnosticReport()
    if ctx.model is None:
        report.skip(
            "unconstrained_variables", "the Pyomo model has not been built"
        )
        return report

    import pyomo.environ as pyo
    from pyomo.core.expr.visitor import identify_variables

    used: set[int] = set()
    for con in ctx.model.component_data_objects(pyo.Constraint, active=True):
        if con.body is None:
            continue
        for var in identify_variables(con.body, include_fixed=False):
            used.add(id(var))
    for obj in ctx.model.component_data_objects(pyo.Objective, active=True):
        for var in identify_variables(obj.expr, include_fixed=False):
            used.add(id(var))

    loose = [
        var
        for var in ctx.model.component_data_objects(pyo.Var, active=True)
        if not var.fixed and id(var) not in used
    ]
    if not loose:
        return report

    report.add(
        DiagnosticIssue(
            severity=DiagnosticSeverity.WARNING,
            category=DiagnosticCategory.STRUCTURE,
            code="MODEL_UNUSED_VARIABLE",
            message=(
                f"{len(loose)} free variables appear in no active constraint "
                f"or objective."
            ),
            recommendation=(
                "Usually a constraint that was not added. Harmless to the "
                "optimum, but the model is not the one you think it is."
            ),
            context={"examples": [v.name for v in loose[:10]]},
        )
    )
    return report


def check_structural_singularity(ctx: DiagnosticContext) -> DiagnosticReport:
    """Look for structurally singular blocks in the equality system.

    Uses Pyomo's incidence analysis and the Dulmage-Mendelsohn
    decomposition to split the equality constraints and the variables they
    contain into under-determined, well-determined and over-determined
    blocks. Anything in the first or last block is structurally singular:
    no assignment of numbers can fix it, because the sparsity pattern
    itself is wrong.

    What it checks: equality constraints only. Inequalities do not
    participate in a square-system decomposition.

    Necessary or sufficient: a non-empty under- or over-determined block
    means the equality system is singular regardless of the numbers, which
    is strong evidence but not the same as OPF infeasibility — an OPF has
    more variables than equalities by design, so the *variable* side being
    under-determined is expected and is not reported.

    Cost: moderate to expensive on large models; part of the `deep` level.

    Args:
        ctx: The diagnostic context.

    Returns:
        A report naming the over-determined constraint block, mapped back
        to network objects where the metadata registry knows the family.
    """
    report = DiagnosticReport()
    if ctx.model is None:
        report.skip(
            "structural_singularity", "the Pyomo model has not been built"
        )
        return report

    try:
        from pyomo.contrib.incidence_analysis import IncidenceGraphInterface
    except ImportError as exc:
        report.skip(
            "structural_singularity", f"incidence analysis unavailable: {exc}"
        )
        return report

    try:
        igraph = IncidenceGraphInterface(ctx.model, include_inequality=False)
        var_dm, con_dm = igraph.dulmage_mendelsohn()
    except Exception as exc:
        report.skip(
            "structural_singularity",
            f"raised {type(exc).__name__}: {exc}",
        )
        return report

    over = list(getattr(con_dm, "unmatched", []) or [])
    if not over:
        report.summary["structure.dulmage_mendelsohn"] = "no singular block"
        return report

    described = []
    for con in over[:10]:
        parent = con.parent_component()
        meta = constraint_meta(parent.local_name)
        described.append(
            {
                "pyomo": con.name,
                "meaning": meta.description if meta else None,
            }
        )

    report.add(
        DiagnosticIssue(
            severity=DiagnosticSeverity.WARNING,
            category=DiagnosticCategory.STRUCTURE,
            code="MODEL_STRUCTURALLY_SINGULAR",
            message=(
                f"{len(over)} equality constraints could not be matched to a "
                f"variable, so the equality system is structurally singular."
            ),
            recommendation=(
                "These constraints form the unmatched block of the "
                "Dulmage-Mendelsohn decomposition. That identifies where the "
                "sparsity pattern is wrong, not which constraint is at fault."
            ),
            context={"unmatched": described, "unmatched_count": len(over)},
        )
    )
    unmatched_vars = list(getattr(var_dm, "unmatched", []) or [])
    report.summary["structure.unmatched constraints"] = len(over)
    report.summary["structure.unmatched variables"] = len(unmatched_vars)
    return report
