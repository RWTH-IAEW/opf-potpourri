"""Computation and plotting of the Feasible Operation Region (FOR).

The Feasible Operation Region is the set of all active / reactive power
(P, Q) combinations that the external grid connection points can realise
while satisfying all network constraints (voltage limits, line loading
limits, generator capability curves).

This module provides several sampling strategies for approximating the FOR
boundary:

* ``run_feasible_operation_region`` – angle-based sampling (36 directions).
* ``for_angle_based_sampling`` – refined angle-based sampling split into four
  quadrants.
* ``for_setpoint_based`` – eight cardinal direction setpoint optimisations.
* ``for_setpoint_based_with_directions`` – adaptive refinement that starts
  from eight directions and inserts additional samples between points that
  are too far apart on the boundary.
* ``node_for_setpoint_based_with_directions`` – same as above but for a
  single static generator node instead of the external grid.

Plotting
--------
``get_pq_to_plot`` flattens the nested boundary lists into plottable arrays.
``plot_hull`` visualises individual grid connection points plus their convex
or concave hull (requires *shapely* and *geopandas*).

Usage example
-------------
::

    net = sb.get_simbench_net("1-HV-mixed--0-no_sw")
    acopf = ACOPF(net)
    acopf.add_OPF()
    p, q, u, nets = for_setpoint_based_with_directions(acopf, stepsize=40)

Institut für Elektrische Anlagen und Netze, Digitalisierung und
Energiewirtschaft (IAEW)
(c) 2023, Steffen Kortmann
"""

from __future__ import division

import copy
import math

import matplotlib.pyplot as plt
import numpy as np
import pandapower as pp
import simbench as sb
from pyomo.environ import (
    Constraint,
    Objective,
    Param,
    check_optimal_termination,
    maximize,
    minimize,
    value,
)

# from shapely import concave_hull, MultiPoint, convex_hull
# import geopandas as gpd

from potpourri.models.ACOPF_base import ACOPF
from potpourri.models.HC_ACOPF import HC_ACOPF


def run_feasible_operation_region(opf):
    """Sample the FOR boundary by sweeping 36 power-factor angles.

    For each angle θ ∈ [0, 2π) a power-factor constraint
    Q = tan(θ)·P is imposed on every external-grid bus and the
    squared apparent power is maximised.  The resulting (P, Q) pairs
    trace the FOR boundary.

    Args:
        opf: A solved ``ACOPF`` (or compatible) instance whose Pyomo model
            contains ``eG``, ``pG``, and ``qG`` components.

    Returns:
        A tuple ``(boundary_P_values, boundary_Q_values)`` where each
        element is a list of lists indexed by external-grid index.
    """
    theta_values = np.linspace(0, 2 * np.pi, 36)
    n_ext_grids = len(list(opf.model.eG))
    boundary_P_values = [[] for _ in range(n_ext_grids)]
    boundary_Q_values = [[] for _ in range(n_ext_grids)]

    opf.model.obj = Objective(
        expr=sum(
            opf.model.pG[b0] ** 2 + opf.model.qG[b0] ** 2
            for b0 in opf.model.eG
        ),
        sense=maximize,
    )

    opf.model.tan_theta = Param(mutable=True, initialize=0.0)

    def constrain_pf(model, b0):
        return model.qG[b0] == model.tan_theta * model.pG[b0]

    opf.model.power_factor_constraint = Constraint(
        opf.model.eG, rule=constrain_pf
    )

    for theta in theta_values:
        opf.model.tan_theta = np.tan(theta)
        opf.solve()

        for g in opf.model.eG:
            print(value(opf.model.pG[g]))

        for b0 in opf.model.eG:
            boundary_P_values[b0].append(value(opf.model.pG[b0]))
            boundary_Q_values[b0].append(value(opf.model.qG[b0]))

    return boundary_P_values, boundary_Q_values


def for_angle_based_sampling(opf, n=36):
    """Sample the FOR boundary using angle-based optimisation in four quadrants.

    Splits *n* evenly spaced angles into four quadrant groups.  Each quadrant
    uses a sign combination (a, b) so that the linear objective
    ``a·P + b·|tan(θ)|·Q`` points into the correct half-plane, avoiding
    degenerate solutions near the axes.

    Args:
        opf: A solved ``ACOPF`` instance.
        n: Number of angle samples.  Must be a multiple of 4.

    Returns:
        A tuple ``(boundary_P_values, boundary_Q_values, boundary_U_values)``
        where each element is a list of lists indexed by external-grid index.
    """
    theta_values = np.linspace(0, 2 * np.pi, n)
    if n % 4 == 0:
        theta_splits = np.split(theta_values, 4)
    else:
        print("n has to be a multiple of 4")
        return [], [], []

    a = [-1, 1, 1, -1]
    b = [-1, -1, 1, 1]

    n_ext_grids = len(list(opf.model.eG))
    boundary_P_values = [[] for _ in range(n_ext_grids)]
    boundary_Q_values = [[] for _ in range(n_ext_grids)]
    boundary_U_values = [[] for _ in range(n_ext_grids)]

    opf.model.tan_theta = Param(mutable=True, initialize=0.0)

    def constrain_pf(model, b0):
        return model.qG[b0] == model.tan_theta * model.pG[b0]

    opf.model.power_factor_constraint = Constraint(
        opf.model.eG, rule=constrain_pf
    )

    opf.model.a = Param(initialize=-1, mutable=True)
    opf.model.b = Param(initialize=-1, mutable=True)
    opf.model.obj = Objective(
        expr=sum(
            opf.model.a * opf.model.pG[b0]
            + opf.model.b * abs(opf.model.tan_theta) * opf.model.qG[b0]
            for b0 in opf.model.eG
        ),
        sense=minimize,
    )

    for quad_idx, theta_i in enumerate(theta_splits):
        opf.model.a = a[quad_idx]
        opf.model.b = b[quad_idx]

        for theta in theta_i:
            opf.model.tan_theta = np.tan(theta)
            opf.solve()

            for g in opf.model.eG:
                print(value(opf.model.pG[g]))

            for b0 in opf.model.eG:
                boundary_P_values[b0].append(value(opf.model.pG[b0]))
                boundary_Q_values[b0].append(value(opf.model.qG[b0]))
                boundary_U_values[b0].append(value(opf.model.v[b0]))

    return boundary_P_values, boundary_Q_values, boundary_U_values


def for_setpoint_based(opf, n=36):
    """Approximate the FOR using eight cardinal setpoint directions.

    Minimises ``-alpha·P - beta·Q`` for each of the eight (alpha, beta)
    combinations that point into the eight cardinal directions of the (P, Q)
    plane.

    Args:
        opf: A solved ``ACOPF`` instance.
        n: Unused parameter kept for API compatibility.

    Returns:
        A tuple ``(boundary_P_values, boundary_Q_values)`` where each
        element is a list of lists indexed by external-grid index.
    """
    alpha_beta = [
        (1, 0),
        (1, 1),
        (0, 1),
        (-1, 1),
        (-1, 0),
        (-1, -1),
        (0, -1),
        (1, -1),
    ]
    opf.model.alpha = Param(initialize=0, mutable=True)
    opf.model.beta = Param(initialize=0, mutable=True)

    def setpoint_based(model):
        return sum(
            -model.pG[g] * model.alpha + -model.qG[g] * model.beta
            for g in model.eG
        )

    opf.model.obj = Objective(expr=setpoint_based, sense=minimize)

    n_ext_grids = len(list(opf.model.eG))
    boundary_P_values = [[] for _ in range(n_ext_grids)]
    boundary_Q_values = [[] for _ in range(n_ext_grids)]

    for alpha, beta in alpha_beta:
        opf.model.alpha = alpha
        opf.model.beta = beta
        opf.solve()
        print(value(opf.model.obj))

        for g in opf.model.eG:
            print(value(opf.model.pG[g]))

        for b0 in opf.model.eG:
            boundary_P_values[b0].append(value(opf.model.pG[b0]))
            boundary_Q_values[b0].append(value(opf.model.qG[b0]))

    return boundary_P_values, boundary_Q_values


def for_setpoint_based_with_directions(opf, stepsize=100, solver="ipopt"):
    """Compute an adaptive FOR boundary with automatic gap refinement.

    Algorithm
    ---------
    1. Solve the OPF for eight cardinal (alpha, beta) directions to obtain
       coarse boundary samples.
    2. Close the polygon by appending the first point at the end.
    3. Compute the normalised arc-length step ``d`` between consecutive
       points.  Insert additional samples wherever ``d >= stepsize / delta_max``.
    4. For gaps where ΔP dominates, sweep P linearly between the bounding
       points while optimising Q freely.
    5. For gaps where ΔQ dominates, sweep Q linearly while optimising P freely.

    Args:
        opf: An ``ACOPF`` instance with OPF constraints already added.
        stepsize: Target normalised step size.  Smaller values → denser
            sampling.
        solver: Pyomo solver name passed to ``opf.solve``.

    Returns:
        A tuple ``(boundary_P_values, boundary_Q_values, boundary_V_values,
        nets)`` where boundary lists are indexed as
        ``[ext_grid_index][refinement_pass]`` and *nets* contains deep-copied
        network states at each coarse corner point.
    """
    alpha_beta = [
        (1, 0),
        (1, 1),
        (0, 1),
        (-1, 1),
        (-1, 0),
        (-1, -1),
        (0, -1),
        (1, -1),
    ]

    n_ext_grids = len(opf.model.eG)
    range_ext_grids = range(n_ext_grids)

    def _set_objective_function(model):
        model.alpha = Param(initialize=0, mutable=True)
        model.beta = Param(initialize=0, mutable=True)

        def setpoint_based(model):
            return sum(
                model.pG[g] * model.alpha + model.qG[g] * model.beta
                for g in model.eG
            )

        model.obj = Objective(expr=setpoint_based, sense=maximize)

    _set_objective_function(opf.model)

    boundary_P_values = [[] for _ in range_ext_grids]
    boundary_Q_values = [[] for _ in range_ext_grids]
    boundary_V_values = [[] for _ in range_ext_grids]
    nets = []
    p = [[] for _ in range_ext_grids]
    q = [[] for _ in range_ext_grids]
    v = [[] for _ in range_ext_grids]

    # ---- coarse pass: eight cardinal directions ---- #
    for alpha, beta in alpha_beta:
        opf.model.alpha = alpha
        opf.model.beta = beta
        opf.solve(solver=solver)

        for b0 in opf.model.eG:
            p[b0].append(value(opf.model.pG[b0]))
            q[b0].append(value(opf.model.qG[b0]))
            v[b0].append(value(opf.model.v[b0]))

    for g in range_ext_grids:
        boundary_P_values[g].append(p[g])
        boundary_Q_values[g].append(q[g])
        boundary_V_values[g].append(v[g])

        # Close the polygon for gap detection.
        boundary_P_values[g][-1].append(boundary_P_values[g][-1][0])
        boundary_Q_values[g][-1].append(boundary_Q_values[g][-1][0])
        boundary_V_values[g][-1].append(boundary_V_values[g][-1][0])

    p_max = np.array([max(boundary_P_values[i][0]) for i in range_ext_grids])
    p_min = np.array([min(boundary_P_values[i][0]) for i in range_ext_grids])
    q_max = np.array([max(boundary_Q_values[i][0]) for i in range_ext_grids])
    q_min = np.array([min(boundary_Q_values[i][0]) for i in range_ext_grids])

    delta_p = abs(p_max - p_min)
    delta_q = abs(q_max - q_min)
    delta_max = np.max((delta_p, delta_q))
    d_max = stepsize / delta_max

    p_diff = [np.diff(boundary_P_values[i][0]) for i in range_ext_grids]
    q_diff = [np.diff(boundary_Q_values[i][0]) for i in range_ext_grids]
    ds = [
        np.sqrt((p_diff[g] / delta_p[g]) ** 2 + (q_diff[g] / delta_q[g]) ** 2)
        for g in range_ext_grids
    ]
    ind_next = [np.argwhere(ds[g] >= d_max) for g in range_ext_grids]
    ind_next = np.unique(np.concatenate(ind_next))

    tol = 0.1

    def _add_setpoint_constraints(model):
        """Add mutable P and Q setpoint band constraints to the model."""
        model.p_sp = Param(model.eG, mutable=True)

        def p_eg_upper(model, g):
            return model.pG[g] <= model.p_sp[g] + tol

        model.p_eg_max = Constraint(model.eG, rule=p_eg_upper)

        def p_eg_lower(model, g):
            return model.pG[g] >= model.p_sp[g] - tol

        model.p_eg_min = Constraint(model.eG, rule=p_eg_lower)

        model.q_sp = Param(model.eG, mutable=True)

        def q_eg_upper(model, g):
            return model.qG[g] <= model.q_sp[g] + tol

        model.q_eg_max = Constraint(model.eG, rule=q_eg_upper)

        def q_eg_lower(model, g):
            return model.qG[g] >= model.q_sp[g] - tol

        model.q_eg_min = Constraint(model.eG, rule=q_eg_lower)

    _add_setpoint_constraints(opf.model)

    opf.model.q_eg_min.deactivate()
    opf.model.q_eg_max.deactivate()

    # ---- refinement pass: P-dominated gaps ---- #
    p = [[] for _ in range_ext_grids]
    q = [[] for _ in range_ext_grids]
    v = [[] for _ in range_ext_grids]

    ind_alpha_0 = np.argwhere(
        abs(p_diff[0][ind_next]) >= abs(q_diff[0][ind_next])
    )

    for row in ind_alpha_0:
        ind = ind_next[row[0]]
        alpha = alpha_beta[ind][0]
        beta = alpha_beta[ind][1]

        opf.model.alpha = alpha
        opf.model.beta = beta

        opf.model.p_eg_min.deactivate()
        opf.model.p_eg_max.deactivate()
        opf.solve(solver=solver)

        nets.append(copy.deepcopy(opf.net))

        opf.model.p_eg_min.activate()
        opf.model.p_eg_max.activate()

        n_steps = math.ceil(
            max(abs(p_diff[g][ind]) for g in range(len(p_diff))) / stepsize
        )
        p_sp = np.array(
            [
                np.linspace(
                    boundary_P_values[g][0][ind],
                    boundary_P_values[g][0][ind + 1],
                    n_steps,
                )
                for g in range(len(boundary_P_values))
            ]
        )

        opf.model.alpha = 0

        for step in range(len(p_sp[0])):
            for g in opf.model.eG:
                opf.model.p_sp[g] = p_sp[g, step]
            opf.solve(solver=solver)

            if check_optimal_termination(opf.results):
                for g in opf.model.eG:
                    p[g].append(value(opf.model.pG[g]))
                    q[g].append(value(opf.model.qG[g]))
                for b in opf.model.b0:
                    v[b].append(value(opf.model.v[b]))

    for i in range_ext_grids:
        boundary_P_values[i].append(p[i])
        boundary_Q_values[i].append(q[i])
        boundary_V_values[i].append(v[i])

    opf.model.p_eg_max.deactivate()
    opf.model.p_eg_min.deactivate()

    # ---- refinement pass: Q-dominated gaps ---- #
    p = [[] for _ in range_ext_grids]
    q = [[] for _ in range_ext_grids]
    v = [[] for _ in range_ext_grids]

    ind_beta_0 = np.argwhere(
        abs(p_diff[0][ind_next]) <= abs(q_diff[0][ind_next])
    )
    for row in ind_beta_0:
        ind = ind_next[row[0]]
        beta = alpha_beta[ind][1]
        alpha = alpha_beta[ind][0]

        n_steps = math.ceil(
            max(abs(q_diff[g][ind]) for g in range(len(q_diff))) / stepsize
        )

        q_sp = np.array(
            [
                np.linspace(
                    boundary_Q_values[g][0][ind],
                    boundary_Q_values[g][0][ind + 1],
                    n_steps,
                )
                for g in range(len(boundary_Q_values))
            ]
        )

        opf.model.alpha = alpha
        opf.model.beta = beta

        opf.model.q_eg_min.deactivate()
        opf.model.q_eg_max.deactivate()
        opf.solve(solver=solver)

        nets.append(copy.deepcopy(opf.net))

        opf.model.q_eg_min.activate()
        opf.model.q_eg_max.activate()

        opf.model.beta = 0

        for step in range(len(q_sp[0])):
            for g in opf.model.eG:
                opf.model.q_sp[g] = q_sp[g, step]
            opf.solve(solver=solver)
            for g in opf.model.eG:
                print(value(opf.model.pG[g]))
            if check_optimal_termination(opf.results):
                for g in opf.model.eG:
                    p[g].append(value(opf.model.pG[g]))
                    q[g].append(value(opf.model.qG[g]))
                for b in opf.model.b0:
                    v[b].append(value(opf.model.v[b]))

    for i in range_ext_grids:
        boundary_P_values[i].append(p[i])
        boundary_Q_values[i].append(q[i])
        boundary_V_values[i].append(v[i])

    return boundary_P_values, boundary_Q_values, boundary_V_values, nets


def node_for_setpoint_based_with_directions(opf, w, stepsize=100):
    """Compute the FOR for a single static generator node *w*.

    Mirrors ``for_setpoint_based_with_directions`` but sweeps the capability
    region of one static generator (indexed *w* in ``opf.model.sG``) rather
    than the aggregated external-grid injection.

    Args:
        opf: An ``ACOPF`` instance with OPF constraints already added.
            Must already have ``p_eg_min``/``p_eg_max`` and
            ``q_eg_min``/``q_eg_max`` constraints (added externally before
            calling this function).
        w: Integer index of the static generator in ``opf.model.sG``.
        stepsize: Target normalised step size for gap refinement.

    Returns:
        A tuple ``(boundary_P_values, boundary_Q_values, boundary_V_values,
        nets)`` for the single generator node.
    """
    alpha_beta = [
        (1, 0),
        (1, 1),
        (0, 1),
        (-1, 1),
        (-1, 0),
        (-1, -1),
        (0, -1),
        (1, -1),
    ]

    n_ext_grids = 1
    range_ext_grids = range(n_ext_grids)

    def _set_objective_function(model):
        model.alpha = Param(initialize=0, mutable=True)
        model.beta = Param(initialize=0, mutable=True)

        def setpoint_based(model):
            return model.psG[w] * model.alpha + model.qsG[w] * model.beta

        model.obj = Objective(expr=setpoint_based, sense=maximize)

    _set_objective_function(opf.model)

    boundary_P_values = [[] for _ in range_ext_grids]
    boundary_Q_values = [[] for _ in range_ext_grids]
    boundary_V_values = [[] for _ in range_ext_grids]
    nets = []
    p = [[] for _ in range_ext_grids]
    q = [[] for _ in range_ext_grids]
    v = [[] for _ in range_ext_grids]

    # ---- coarse pass ---- #
    for alpha, beta in alpha_beta:
        opf.model.alpha = alpha
        opf.model.beta = beta
        opf.solve()

        for b0 in range_ext_grids:
            p[b0].append(value(opf.model.psG[w]))
            q[b0].append(value(opf.model.qsG[w]))

    for g in range_ext_grids:
        boundary_P_values[g].append(p[g])
        boundary_Q_values[g].append(q[g])
        boundary_V_values[g].append(v[g])

        boundary_P_values[g][-1].append(boundary_P_values[g][-1][0])
        boundary_Q_values[g][-1].append(boundary_Q_values[g][-1][0])

    p_max = np.array([max(boundary_P_values[i][0]) for i in range_ext_grids])
    p_min = np.array([min(boundary_P_values[i][0]) for i in range_ext_grids])
    q_max = np.array([max(boundary_Q_values[i][0]) for i in range_ext_grids])
    q_min = np.array([min(boundary_Q_values[i][0]) for i in range_ext_grids])

    delta_p = abs(p_max - p_min)
    delta_q = abs(q_max - q_min)
    delta_max = np.max((delta_p, delta_q))
    d_max = stepsize / delta_max

    p_diff = [np.diff(boundary_P_values[i][0]) for i in range_ext_grids]
    q_diff = [np.diff(boundary_Q_values[i][0]) for i in range_ext_grids]
    ds = [
        np.sqrt((p_diff[g] / delta_p[g]) ** 2 + (q_diff[g] / delta_q[g]) ** 2)
        for g in range_ext_grids
    ]
    ind_next = [np.argwhere(ds[g] >= d_max) for g in range_ext_grids]
    ind_next = np.unique(np.concatenate(ind_next))

    opf.model.q_eg_min.deactivate()
    opf.model.q_eg_max.deactivate()

    # ---- refinement: P-dominated gaps ---- #
    p = [[], [], []]
    q = [[], [], []]
    v = [[], [], []]

    ind_alpha_0 = np.argwhere(
        abs(p_diff[0][ind_next]) >= abs(q_diff[0][ind_next])
    )

    for row in ind_alpha_0:
        ind = ind_next[row[0]]
        alpha = alpha_beta[ind][0]
        beta = alpha_beta[ind][1]

        opf.model.alpha = alpha
        opf.model.beta = beta

        opf.model.p_eg_min.deactivate()
        opf.model.p_eg_max.deactivate()
        opf.solve()

        nets.append(copy.deepcopy(opf.net))

        opf.model.p_eg_min[w].activate()
        opf.model.p_eg_max[w].activate()

        n_steps = math.ceil(
            max(abs(p_diff[g][ind]) for g in range(len(p_diff))) / stepsize
        )
        p_sp = np.array(
            [
                np.linspace(
                    boundary_P_values[g][0][ind],
                    boundary_P_values[g][0][ind + 1],
                    n_steps,
                )
                for g in range(len(boundary_P_values))
            ]
        )

        opf.model.alpha = 0

        for step in range(len(p_sp[0])):
            for g in range_ext_grids:
                opf.model.p_sp[w] = p_sp[g, step]
            opf.solve()

            if check_optimal_termination(opf.results):
                for g in range_ext_grids:
                    p[g].append(value(opf.model.psG[w]))
                    q[g].append(value(opf.model.qsG[w]))

    for i in range_ext_grids:
        boundary_P_values[i].append(p[i])
        boundary_Q_values[i].append(q[i])
        boundary_V_values[i].append(v[i])

    opf.model.p_eg_max.deactivate()
    opf.model.p_eg_min.deactivate()

    # ---- refinement: Q-dominated gaps ---- #
    ind_beta_0 = np.argwhere(
        abs(p_diff[0][ind_next]) <= abs(q_diff[0][ind_next])
    )
    for row in ind_beta_0:
        ind = ind_next[row[0]]
        beta = alpha_beta[ind][1]
        alpha = alpha_beta[ind][0]

        n_steps = math.ceil(
            max(abs(q_diff[g][ind]) for g in range(len(q_diff))) / stepsize
        )

        q_sp = np.array(
            [
                np.linspace(
                    boundary_Q_values[g][0][ind],
                    boundary_Q_values[g][0][ind + 1],
                    n_steps,
                )
                for g in range(len(boundary_Q_values))
            ]
        )

        opf.model.alpha = alpha
        opf.model.beta = beta

        opf.model.q_eg_min.deactivate()
        opf.model.q_eg_max.deactivate()
        opf.solve()

        nets.append(copy.deepcopy(opf.net))

        opf.model.q_eg_min[w].activate()
        opf.model.q_eg_max[w].activate()

        opf.model.beta = 0

        for step in range(len(q_sp[0])):
            for g in range_ext_grids:
                opf.model.q_sp[w] = q_sp[g, step]
            opf.solve()

            if check_optimal_termination(opf.results):
                for g in range_ext_grids:
                    p[g].append(value(opf.model.psG[w]))
                    q[g].append(value(opf.model.qsG[w]))

    for i in range_ext_grids:
        boundary_P_values[i].append(p[i])
        boundary_Q_values[i].append(q[i])
        boundary_V_values[i].append(v[i])

    return boundary_P_values, boundary_Q_values, boundary_V_values, nets


def get_pq_to_plot(p, q):
    """Flatten FOR boundary lists into plottable (P, Q) arrays.

    Accepts both the flat ``[v1, v2, ...]`` structure returned by
    ``for_setpoint_based`` and the nested ``[[v1, ...], [v2, ...]]``
    structure returned by ``for_setpoint_based_with_directions``.

    Args:
        p: List indexed by ext_grid, each element either a flat list of
            P-values or a list of lists (one per refinement pass).
        q: Corresponding structure for Q values.

    Returns:
        A tuple ``(pq, pq_tot)`` where *pq* is a list of Nx2 arrays (one per
        external grid) and *pq_tot* is the element-wise sum of all grids.
    """
    pq = []
    for i in range(len(p)):
        # Flat list when elements are scalars; nested list otherwise.
        if p[i] and isinstance(p[i][0], (int, float, np.floating)):
            p_i = list(p[i])
            q_i = list(q[i])
        else:
            p_i = [v for p_list in p[i] for v in p_list]
            q_i = [v for q_list in q[i] for v in q_list]
        pq.append(np.array((p_i, q_i)).T)

    pq_tot = sum(pq_i for pq_i in pq)
    return pq, pq_tot


def plot_hull(p, q, ratio=0.1):
    """Plot the FOR boundary points and their convex / concave hulls.

    Requires *shapely* and *geopandas* (commented-out imports at the top of
    this module).  Falls back to convex hull if a concave hull cannot be
    computed.

    Args:
        p: Nested boundary P list as returned by the sampling functions.
        q: Corresponding nested boundary Q list.
        ratio: Concave-hull ratio passed to ``shapely.concave_hull``.  A
            value of 0 gives a convex hull; 1 gives the tightest fit.

    Returns:
        A tuple ``(figs, axs)`` – lists of Matplotlib Figure and Axes objects,
        one entry per external grid plus one combined figure.
    """
    pq = []
    for i in range(len(p)):
        p_i = [p_mw for p_list in p[i] for p_mw in p_list]
        q_i = [q_mw for q_list in q[i] for q_mw in q_list]
        pq.append(np.array((p_i, q_i)).T)

    clrs = [
        "#00549F",
        "#000000",
        "#E30066",
        "#FFED00",
        "#006165",
        "#0098A1",
        "#57AB27",
        "#BDCD00",
        "#F6A800",
        "#CC071E",
        "#A11035",
        "#612158",
        "#7A6FAC",
    ]

    fig, ax = plt.subplots()
    figs = [fig]
    axs = [ax]

    for i, points in enumerate(pq):
        fig_i, ax_i = plt.subplots()
        clr = clrs[i + 6]

        ax.plot(
            points[:, 0],
            points[:, 1],
            ".",
            label="Slack #" + str(i),
            color=clr,
        )
        ax.legend()

        ax_i.plot(
            points[:, 0],
            points[:, 1],
            ".",
            label="Slack #" + str(i),
            color=clr,
        )

        try:
            hull = concave_hull(MultiPoint(points), ratio)  # noqa: F821
        except Exception as e:
            print("Creating convex hull instead of concave hull: " + str(e))
            hull = convex_hull(MultiPoint(points))  # noqa: F821

        polygon = gpd.GeoSeries(hull)  # noqa: F821
        polygon.plot(ax=ax, alpha=0.2, color=clr)
        polygon.boundary.plot(ax=ax, edgecolor=clr, linewidth=2)

        polygon.plot(ax=ax_i, alpha=0.2, color=clr)
        polygon.boundary.plot(ax=ax_i, edgecolor=clr, linewidth=2)
        ax_i.set_xlabel("P [MW]")
        ax_i.set_ylabel("Q [MVar]")
        figs.append(fig_i)
        axs.append(ax_i)

        ax.set_xlabel("P [MW]")
        ax.set_ylabel("Q [MVar]")

    for axes in axs:
        axes.grid()
        axes.set_aspect("equal")

    return figs, axs


if __name__ == "__main__":
    # ------------------------------------------------------------------ #
    # Network setup                                                        #
    # ------------------------------------------------------------------ #
    # Note: the __main__ block uses the HV-mixed simbench network which has
    # three external-grid connection points and includes Wind generators.
    # The HC_ACOPF hosting-capacity step requires custom columns
    # ``net.sgen["wind_hc"]`` and ``net.sgen["var_q"]`` that must be added
    # via ``net_augmentation/prepare_net.py`` before running.

    net = sb.get_simbench_net("1-HV-mixed--0-no_sw")

    case = "lW"
    factors = net.loadcases.loc[case]
    net.load.p_mw *= factors["pload"]
    net.load.q_mvar *= factors["qload"]
    net.sgen.loc[net.sgen.type == "Wind", "scaling"] = factors["Wind_p"]
    net.sgen.loc[net.sgen.type == "PV", "scaling"] = factors["PV_p"]
    net.sgen.loc[
        (net.sgen.type != "Wind") & (net.sgen.type != "Solar"), "scaling"
    ] = factors["RES_p"]
    net.ext_grid.vm_pu = factors["Slack_vm"]

    # ---- Step 1: Hosting-Capacity OPF to determine wind placement ---- #
    hc = HC_ACOPF(net)
    hc.solve()
    hc.add_OPF(SWmin=10)
    hc.add_tap_changer_linear()
    hc.solve(solver="neos", neos_opt="bonmin")

    net_hc = hc.net.copy()

    # ---- Step 2: Add HC wind generators to a fresh copy of the net ---- #
    # ``res_sgen.y_wind`` contains the binary HC decision variable written
    # back by pyo_to_net.  Requires prepare_net.py to have added the
    # ``wind_hc`` column first.
    input_hc_net_dir = sb.get_simbench_net("1-HV-mixed--0-no_sw")
    input_hc_net_dir.ext_grid.vm_pu = factors["Slack_vm"]
    wind_hc_index = net_hc.sgen.index[net_hc.res_sgen.y_wind == 1]

    pp.create_sgens(
        net_hc,
        net_hc.sgen.bus[wind_hc_index],
        p_mw=net_hc.sgen.p_mw[wind_hc_index],
        var_q=0,
        type="Wind",
        wind_hc=True,
    )

    # ---- Step 3: FOR computation on the wind-augmented network ---- #
    net_wind = copy.deepcopy(net)
    net_wind.sgen["wind_hc"].fillna(False, inplace=True)

    net_wind.sgen["controllable"] = True
    net_wind.sgen.loc[net_wind.sgen.type == "Wind", "controllable"] = True
    net_wind.sgen["p_inst_mw"] = net_wind.sgen["p_mw"]
    net_wind.sgen.loc[net_wind.sgen.type == "Wind", "var_q"] = 1
    net_wind.sgen.loc[net_wind.sgen.wind_hc, "var_q"] = 0
    net_wind.load["controllable"] = True

    acopf = ACOPF(net_wind)
    acopf.add_OPF()
    acopf.add_tap_changer_linear()

    p, q, u, nets = for_setpoint_based_with_directions(acopf, stepsize=40)
