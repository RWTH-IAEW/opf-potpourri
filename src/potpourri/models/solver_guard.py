# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Refuse to hand integer variables to a solver that would ignore them.

Pyomo's NL writer flags integer and binary variables as such, but IPOPT and
the other continuous NLP codes simply treat them as continuous, so a
discrete tap position, shunt step or placement binary comes back fractional
and is reported as if it were a solution of the discrete model. Both
`Basemodel.solve` and `Basemodel_multi_period.solve` call
`check_integrality_support` before invoking the solver, so that combination
raises a `ValueError` naming the offending variables unless the caller asks
for the continuous relaxation explicitly (`relax_integrality=True`).

The module is deliberately separate from the model classes: it is about
solver capability, not about reading a network, and it has no dependency on
anything else in the package.
"""

import pyomo.environ as pyo
from loguru import logger

#: Solver names (prefixes, lower case) that solve continuous problems only.
#: Pyomo's NL writer flags integer variables as integer, and these solvers
#: ignore the flag, so a model with free integer variables would be solved
#: as its continuous relaxation without anyone noticing. `solve` refuses
#: that combination unless asked for the relaxation explicitly.
CONTINUOUS_ONLY_SOLVERS = ("ipopt", "conopt", "snopt", "minos")


def free_integer_variables(model):
    """Names of the integer and binary variables of `model` that are free.

    A fixed variable is a constant to every writer and solver, so only the
    unfixed ones count. Walks every active block of the model.

    Args:
        model: A Pyomo `ConcreteModel`.

    Returns:
        A sorted list of fully qualified variable names, e.g.
        `["trafo_tap_position[0]", "y[3]"]`. Empty when the model is
        continuous.
    """
    names = []
    for var in model.component_data_objects(
        pyo.Var, active=True, descend_into=True
    ):
        if var.fixed:
            continue
        if var.is_integer() or var.is_binary():
            names.append(var.name)
    return sorted(names)


def check_integrality_support(
    model, solver, relax_integrality=False, neos_opt=None
):
    """Refuse to hand free integer variables to a continuous-only solver.

    IPOPT and the other NLP codes in `CONTINUOUS_ONLY_SOLVERS` ignore the
    integrality of a variable, so a discrete tap position or a binary
    placement variable would come back fractional and be reported as if it
    were a solution of the discrete model. That used to happen silently.

    Args:
        model: The Pyomo model about to be solved.
        solver: Solver name as passed to `solve`.
        relax_integrality: Pass `True` to solve the continuous relaxation
            knowingly; the names of the relaxed variables are then logged
            and returned instead of raising.
        neos_opt: The remote solver name when `solver == "neos"`.

    Returns:
        The names of the integer variables that are being relaxed (empty
        when the solver can handle them or the model has none).

    Raises:
        ValueError: If the solver is continuous-only, the model has free
            integer variables and `relax_integrality` is `False`.
    """
    name = (neos_opt if solver == "neos" else solver) or ""
    if not str(name).lower().startswith(CONTINUOUS_ONLY_SOLVERS):
        return []
    free = free_integer_variables(model)
    if not free:
        return []
    families = sorted({n.split("[", 1)[0] for n in free})
    if relax_integrality:
        logger.warning(
            "Solving the continuous relaxation: {} integer variable(s) "
            "({}) are treated as continuous by solver '{}'.",
            len(free),
            ", ".join(families),
            name,
        )
        return free
    raise ValueError(
        f"The model has {len(free)} free integer variable(s) "
        f"({', '.join(families)}) but solver {str(name)!r} solves "
        "continuous problems only and would silently ignore their "
        "integrality. Use a MINLP solver ('gurobi_direct_minlp', "
        "'mindtpy', or NEOS with 'bonmin'/'couenne'), pass "
        "relax_integrality=True to solve the continuous relaxation "
        "knowingly, or, for tap changers, call "
        "solve_oltc_round_and_fix()."
    )
