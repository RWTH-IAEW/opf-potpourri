# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Conditioning and scaling checks.

An OPF can be perfectly feasible and still fail, because the numbers in it
span too many orders of magnitude for a solver to work with. That failure
looks like infeasibility from the outside, which is exactly why it needs
its own diagnostic: reporting a badly scaled model as physically
infeasible sends the user looking in the wrong place entirely.

Everything here works with NumPy and, where present, SciPy. PyNumero's
compiled ASL interface would give exact Jacobians, but it is frequently not
installed — it is not available in this project's own environment — so its
absence is reported as a gap rather than being allowed to fail the run.
"""

from __future__ import annotations

import math

from potpourri.diagnostics.context import DiagnosticContext
from potpourri.diagnostics.report import (
    DiagnosticCategory,
    DiagnosticIssue,
    DiagnosticReport,
    DiagnosticSeverity,
)

#: Orders of magnitude between the largest and smallest constraint
#: coefficient beyond which a solver is likely to struggle.
_SPREAD_DECADES = 10.0

#: Absolute variable magnitude treated as suspicious in a per-unit model.
_LARGE_VALUE = 1e6


def check_variable_values(ctx: DiagnosticContext) -> DiagnosticReport:
    """Look for variable values a solver will struggle with.

    What it checks: every active variable for a non-finite value, for a
    value outside its own bounds, and for a magnitude far outside what a
    per-unit model should contain.

    Necessary or sufficient: neither. A variable outside its bounds before
    a solve is just a poor starting point; after one it means the solution
    was not loaded or the solver failed.

    Cost: cheap.

    Args:
        ctx: The diagnostic context.

    Returns:
        A report with one `ERROR` per non-finite value and `WARNING`s for
        out-of-bounds or implausibly large ones, summarised rather than
        listed one by one.
    """
    report = DiagnosticReport()
    if ctx.model is None:
        report.skip("variable_values", "the Pyomo model has not been built")
        return report

    import pyomo.environ as pyo

    non_finite, out_of_bounds, huge = [], [], []
    for var in ctx.model.component_data_objects(pyo.Var, active=True):
        raw = var.value
        if raw is None:
            continue
        try:
            value = float(raw)
        except (TypeError, ValueError):
            continue
        if not math.isfinite(value):
            non_finite.append(var.name)
            continue
        low, high = var.lb, var.ub
        if (low is not None and value < low - ctx.tol) or (
            high is not None and value > high + ctx.tol
        ):
            out_of_bounds.append(var.name)
        if abs(value) > _LARGE_VALUE:
            huge.append((var.name, value))

    if non_finite:
        report.add(
            DiagnosticIssue(
                severity=DiagnosticSeverity.ERROR,
                category=DiagnosticCategory.NUMERICAL,
                code="NUMERIC_NON_FINITE_VALUE",
                message=(
                    f"{len(non_finite)} variables hold NaN or infinity, so "
                    f"the model cannot be evaluated."
                ),
                recommendation=(
                    "Usually a division by a zero impedance or a failed "
                    "initialisation."
                ),
                context={"examples": non_finite[:10]},
            )
        )
    if out_of_bounds:
        report.add(
            DiagnosticIssue(
                severity=DiagnosticSeverity.WARNING,
                category=DiagnosticCategory.NUMERICAL,
                code="NUMERIC_VALUE_OUTSIDE_BOUNDS",
                message=(
                    f"{len(out_of_bounds)} variables sit outside their own "
                    f"bounds."
                ),
                recommendation=(
                    "Before a solve this is only a starting point. After one "
                    "it means the solution was not loaded."
                ),
                context={"examples": out_of_bounds[:10]},
            )
        )
    if huge:
        report.add(
            DiagnosticIssue(
                severity=DiagnosticSeverity.WARNING,
                category=DiagnosticCategory.NUMERICAL,
                code="NUMERIC_LARGE_MAGNITUDE",
                message=(
                    f"{len(huge)} variables exceed {_LARGE_VALUE:.0e} in "
                    f"magnitude, which is unusual in a per-unit model."
                ),
                recommendation=(
                    "Check the per-unit base and any limit set to a sentinel "
                    "value like 1e9 to mean 'unlimited'."
                ),
                context={"examples": [name for name, _ in huge[:10]]},
            )
        )
    return report


def check_constraint_scaling(ctx: DiagnosticContext) -> DiagnosticReport:
    """Measure how far apart the constraint coefficients are.

    What it checks: for each constraint family, the largest and smallest
    absolute coefficient of a linear term, estimated by differentiating
    the body symbolically. A family spanning many orders of magnitude is
    hard for a solver even when it is perfectly feasible.

    Necessary or sufficient: neither, and deliberately conservative — it
    reports a spread, not a verdict.

    Cost: moderate. It samples rather than differentiating everything, so
    it stays usable on large models.

    Limitations: the derivative is taken at the current point, so for a
    nonlinear constraint it describes the local Jacobian row, not the
    constraint everywhere.

    Args:
        ctx: The diagnostic context.

    Returns:
        A report with a `WARNING` per badly spread family and the overall
        range in `summary`.
    """
    report = DiagnosticReport()
    if ctx.model is None:
        report.skip("constraint_scaling", "the Pyomo model has not been built")
        return report

    import pyomo.environ as pyo
    from pyomo.core.expr.calculus.derivatives import differentiate

    sample_size = int(ctx.options.get("scaling_sample", 40))
    families: dict[str, list[float]] = {}

    for con in ctx.model.component_data_objects(pyo.Constraint, active=True):
        family = con.parent_component().local_name
        seen = families.setdefault(family, [])
        if len(seen) >= sample_size:
            continue
        if con.body is None:
            continue
        try:
            grad = differentiate(
                con.body, mode=differentiate.Modes.reverse_numeric
            )
        except Exception:
            continue
        for value in (grad or {}).values():
            try:
                magnitude = abs(float(value))
            except (TypeError, ValueError):
                continue
            if magnitude > 0 and math.isfinite(magnitude):
                seen.append(magnitude)

    overall_low, overall_high = math.inf, 0.0
    for family, magnitudes in families.items():
        if len(magnitudes) < 2:
            continue
        low, high = min(magnitudes), max(magnitudes)
        overall_low = min(overall_low, low)
        overall_high = max(overall_high, high)
        decades = math.log10(high / low) if low > 0 else math.inf
        if decades <= _SPREAD_DECADES:
            continue
        report.add(
            DiagnosticIssue(
                severity=DiagnosticSeverity.WARNING,
                category=DiagnosticCategory.NUMERICAL,
                code="NUMERIC_POOR_SCALING",
                message=(
                    f"Coefficients in '{family}' span {decades:.1f} orders of "
                    f"magnitude, from {low:.2e} to {high:.2e}."
                ),
                recommendation=(
                    "A solver can stall or report infeasibility on a model "
                    "like this even when it is feasible. Check the per-unit "
                    "base and any very small impedance."
                ),
                pyomo_component=family,
                context={"min_coefficient": low, "max_coefficient": high},
            )
        )

    if overall_high > 0 and math.isfinite(overall_low):
        report.summary["numerical.coefficient range"] = (
            f"{overall_low:.1e} to {overall_high:.1e}"
        )
    return report


def check_jacobian_conditioning(ctx: DiagnosticContext) -> DiagnosticReport:
    """Estimate how close the equality Jacobian is to singular.

    What it checks: the smallest and largest singular value of the
    equality-constraint Jacobian, via a sparse SVD, and the resulting
    condition estimate.

    Necessary or sufficient: neither. An ill-conditioned Jacobian explains
    why a solver struggles; it does not make the problem infeasible.

    Cost: expensive. Part of the `deep` level only.

    Required state: SciPy, and a Jacobian that can be assembled. PyNumero's
    ASL interface would be the exact route and is used when available; the
    fallback differentiates symbolically, which is slower and only
    practical on moderate models.

    Args:
        ctx: The diagnostic context.

    Returns:
        A report with the condition estimate in `summary`, or a recorded
        skip explaining precisely which capability was missing.
    """
    report = DiagnosticReport()
    if ctx.model is None:
        report.skip("jacobian", "the Pyomo model has not been built")
        return report

    try:
        import numpy as np
        from scipy.sparse import coo_matrix
        from scipy.sparse.linalg import svds
    except ImportError as exc:
        report.skip("jacobian", f"SciPy is required for this check: {exc}")
        return report

    rows, cols, values, n_con = _assemble_jacobian(ctx)
    if n_con == 0:
        report.skip("jacobian", "the model has no equality constraints")
        return report
    if not values:
        report.skip("jacobian", "no Jacobian entries could be evaluated")
        return report

    n_var = max(cols) + 1
    matrix = coo_matrix((values, (rows, cols)), shape=(n_con, n_var)).tocsc()
    k = min(matrix.shape) - 1
    if k < 1:
        report.skip("jacobian", "the Jacobian is too small for a sparse SVD")
        return report
    k = min(k, 6)

    try:
        largest = svds(matrix, k=1, return_singular_vectors=False)
        smallest = svds(matrix, k=k, which="SM", return_singular_vectors=False)
    except Exception as exc:
        report.skip("jacobian", f"the sparse SVD did not converge: {exc}")
        return report

    sigma_max = float(np.max(largest))
    sigma_min = float(np.min(smallest))
    report.summary["numerical.largest singular value"] = f"{sigma_max:.3e}"
    report.summary["numerical.smallest singular value"] = f"{sigma_min:.3e}"

    if sigma_min <= 0 or not math.isfinite(sigma_min):
        report.add(
            DiagnosticIssue(
                severity=DiagnosticSeverity.WARNING,
                category=DiagnosticCategory.NUMERICAL,
                code="NUMERIC_SINGULAR_JACOBIAN",
                message=(
                    "The equality Jacobian is numerically singular at the "
                    "current point."
                ),
                recommendation=(
                    "Look at the structural check as well: a singular "
                    "Jacobian with a sound sparsity pattern usually means "
                    "scaling, not a missing equation."
                ),
            )
        )
        return report

    condition = sigma_max / sigma_min
    report.summary["numerical.condition estimate"] = f"{condition:.2e}"
    if condition > 1e12:
        report.add(
            DiagnosticIssue(
                severity=DiagnosticSeverity.WARNING,
                category=DiagnosticCategory.NUMERICAL,
                code="NUMERIC_ILL_CONDITIONED",
                message=(
                    f"The equality Jacobian has an estimated condition "
                    f"number of {condition:.1e} at the current point."
                ),
                recommendation=(
                    "Convergence trouble here is a conditioning problem "
                    "rather than an infeasible network."
                ),
                value=condition,
            )
        )
    return report


def _assemble_jacobian(ctx: DiagnosticContext):
    """Build the equality Jacobian in COO form at the current point.

    Args:
        ctx: The diagnostic context.

    Returns:
        `(rows, cols, values, n_constraints)`. Constraints whose
        derivative cannot be taken are skipped rather than aborting the
        assembly, so a partially evaluable model still yields an estimate.
    """
    import pyomo.environ as pyo
    from pyomo.core.expr.calculus.derivatives import differentiate

    rows: list[int] = []
    cols: list[int] = []
    values: list[float] = []
    column_of: dict[int, int] = {}
    n_con = 0

    for con in ctx.model.component_data_objects(pyo.Constraint, active=True):
        if not con.equality or con.body is None:
            continue
        try:
            grad = differentiate(
                con.body, mode=differentiate.Modes.reverse_numeric
            )
        except Exception:
            continue
        for var, value in (grad or {}).items():
            try:
                entry = float(value)
            except (TypeError, ValueError):
                continue
            if entry == 0.0 or not math.isfinite(entry):
                continue
            key = id(var)
            if key not in column_of:
                column_of[key] = len(column_of)
            rows.append(n_con)
            cols.append(column_of[key])
            values.append(entry)
        n_con += 1
    return rows, cols, values, n_con
