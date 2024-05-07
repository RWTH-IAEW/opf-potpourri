#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Calculation of the feasible operation region

example of usage:


Institut für Elektrische Anlagen und Netze, Digitalisierung und Energiewirtschaft (IAEW)
(c) 2023, Steffen Kortmann
"""

# ==========Import==========
from __future__ import division

import copy
import math
from pyomo.environ import *
import pickle
from shapely import concave_hull, MultiPoint, convex_hull
import geopandas as gpd
import simbench as sb

from potpourri.models.HC_ACOPF import HC_ACOPF
from potpourri.models.ACOPF_base import ACOPF
# from scripts.plot_functions import *

# from scripts.plot_functions import set_plt_config

# ==========================
def run_feasible_operation_region(self):
    print("Run feasible operation region")
    import numpy as np

    theta_values = np.linspace(0, 2 * np.pi, 36)  # One degree resolution
    boundary_P_values = [[], [], []]
    boundary_Q_values = [[], [], []]
    boundary_U_values = []

    # Set objective to maximize absolute power for the given angle
    self.model.obj = Objective(
        expr=sum(self.model.pG[b0] ** 2 + self.model.qG[b0] ** 2 for b0 in self.model.eG),
        sense=maximize)

    # Constraint for the given power factor angle
    self.model.tan_theta = Param(mutable=True, initialize=0.0)

    def constrain_pf(model, b0):
        return model.qG[b0] == model.tan_theta * model.pG[b0]

    self.model.power_factor_constraint = Constraint(self.model.eG, rule=constrain_pf)

    for theta in theta_values:
        # Constraint for the given power factor angle
        tan_theta = np.tan(theta)
        self.model.tan_theta = tan_theta

        # Solve
        self.solve()

        for i in self.model.eG:
            print(value(self.model.pG[i]))

        for b0 in self.model.eG:
            boundary_P_values[b0].append(value(self.model.pG[b0]))
            boundary_Q_values[b0].append(value(self.model.qG[b0]))

    return boundary_P_values, boundary_Q_values


def for_angle_based_sampling(opf, n=36):
    theta_values = np.linspace(0, 2 * np.pi, n)  # One degree resolution
    if n % 4 == 0:
        theta_splits = np.split(theta_values, 4)
    else:
        print('n has to be a multiple of 4')
    a = [-1, 1, 1, -1]
    b = [-1, -1, 1, 1]

    boundary_P_values = [[], [], []]
    boundary_Q_values = [[], [], []]
    boundary_U_values = [[], [], []]

    # Constraint for the given power factor angle
    opf.model.tan_theta = Param(mutable=True, initialize=0.0)

    def constrain_pf(model, b0):
        return model.qG[b0] == model.tan_theta * model.pG[b0]

    opf.model.power_factor_constraint = Constraint(opf.model.eG, rule=constrain_pf)

    # Set objective to maximize absolute power for the given angle
    opf.model.a = Param(initialize=-1, mutable=True)
    opf.model.b = Param(initialize=-1, mutable=True)
    opf.model.obj = Objective(
        expr=sum(opf.model.a * opf.model.pG[b0] + opf.model.b * abs(opf.model.tan_theta) * opf.model.qG[b0] for b0 in
                 opf.model.eG),
        sense=minimize)

    for i, theta_i in enumerate(theta_splits):
        opf.model.a = a[i]
        opf.model.b = b[i]

        for theta in theta_i:
            # Constraint for the given power factor angle
            tan_theta = np.tan(theta)
            opf.model.tan_theta = tan_theta

            # Solve
            opf.solve()

            for i in opf.model.eG:
                print(value(opf.model.pG[i]))

            for b0 in opf.model.eG:
                boundary_P_values[b0].append(value(opf.model.pG[b0]))
                boundary_Q_values[b0].append(value(opf.model.qG[b0]))
                boundary_U_values[b0].append(value(opf.model.v[b0]))

    return boundary_P_values, boundary_Q_values, boundary_U_values


def for_setpoint_based(opf, n=36):
    alpha_beta = [(1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1)]
    opf.model.alpha = Param(initialize=0, mutable=True)
    opf.model.beta = Param(initialize=0, mutable=True)

    def setpoint_based(model):
        return sum(-model.pG[g] * model.alpha + -model.qG[g] * model.beta for g in model.eG)

    opf.model.obj = Objective(expr=setpoint_based, sense=minimize)

    boundary_P_values = [[], [], []]
    boundary_Q_values = [[], [], []]

    for (alpha, beta) in alpha_beta:
        opf.model.alpha = alpha
        opf.model.beta = beta
        opf.solve()
        print(value(opf.model.obj))

        for i in opf.model.eG:
            print(value(opf.model.pG[i]))

        for b0 in opf.model.eG:
            boundary_P_values[b0].append(value(opf.model.peG[b0]))
            boundary_Q_values[b0].append(value(opf.model.qeG[b0]))

    return boundary_P_values, boundary_Q_values


def for_setpoint_based_with_directions(opf, stepsize=100, solver='ipopt'):
    alpha_beta = [(1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1)]

    n_ext_grids = len(opf.model.eG)
    range_ext_grids = range(n_ext_grids)

    def _set_objective_function(model):
        model.alpha = Param(initialize=0, mutable=True)
        model.beta = Param(initialize=0, mutable=True)

        def setpoint_based(model):
            return sum(model.pG[g] * model.alpha + model.qG[g] * model.beta for g in model.eG)

        model.obj = Objective(expr=setpoint_based, sense=maximize)

    _set_objective_function(opf.model)

    boundary_P_values = [[] for _ in range_ext_grids]
    boundary_Q_values = [[] for _ in range_ext_grids]
    boundary_V_values = [[] for _ in range_ext_grids]
    nets = []
    p = [[] for _ in range_ext_grids]
    q = [[] for _ in range_ext_grids]
    v = [[] for _ in range_ext_grids]

    for (alpha, beta) in alpha_beta:
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

        boundary_P_values[g][-1].append(boundary_P_values[g][-1][0])
        boundary_Q_values[g][-1].append(boundary_Q_values[g][-1][0])
        boundary_V_values[g][-1].append(boundary_V_values[g][-1][0])

    p_max = np.array([max(boundary_P_values[i][0]) for i in range_ext_grids])
    p_min = np.array([min(boundary_P_values[i][0]) for i in range_ext_grids])

    q_max = np.array([max(boundary_Q_values[i][0]) for i in range_ext_grids])
    q_min = np.array([min(boundary_Q_values[i][0]) for i in range_ext_grids])

    delta_p = abs(p_max - p_min)
    delta_q = abs(q_max - q_min)

    # ----
    delta_max = np.max((delta_p, delta_q))
    d_max = stepsize / delta_max

    p_diff = [np.diff(boundary_P_values[i][0]) for i in range_ext_grids]
    q_diff = [np.diff(boundary_Q_values[i][0]) for i in range_ext_grids]
    ds = [np.sqrt((p_diff[g] / delta_p[g]) ** 2 + (q_diff[g] / delta_q[g]) ** 2) for g in range_ext_grids]
    ind_next = [np.argwhere(ds[g] >= d_max) for g in range_ext_grids]
    ind_next = np.unique(np.concatenate(ind_next))
    # p_next = []
    # q_next = []
    # for g in range(len(boundary_P_values)):
    #     p_next.append([(boundary_P_values[g][-1][i] + boundary_P_values[g][-1][i + 1]) / 2 for i in ind_next])
    #     q_next.append([(boundary_Q_values[g][-1][i] + boundary_Q_values[g][-1][i + 1]) / 2 for i in ind_next])

    tol = 0.1

    def _add_setpoint_constraints(model):
        model.p_sp = Param(model.eG, mutable=True)

        def p_eg_upper(model, g):
            return (model.pG[g]) <= (model.p_sp[g]) + tol

        model.p_eg_max = Constraint(model.eG, rule=p_eg_upper)

        def p_eg_lower(model, g):
            return (model.pG[g]) >= (model.p_sp[g]) - tol

        model.p_eg_min = Constraint(model.eG, rule=p_eg_lower)

        model.q_sp = Param(model.eG, mutable=True)

        def q_eg_upper(model, g):
            return (model.qG[g]) <= (model.q_sp[g]) + tol

        model.q_eg_max = Constraint(model.eG, rule=q_eg_upper)

        def q_eg_lower(model, g):
            return (model.qG[g]) >= (model.q_sp[g]) - tol

        model.q_eg_min = Constraint(model.eG, rule=q_eg_lower)

    _add_setpoint_constraints(opf.model)

    opf.model.q_eg_min.deactivate()
    opf.model.q_eg_max.deactivate()

    ind_alpha_0 = np.argwhere(abs(p_diff[0][ind_next]) >= abs(q_diff[0][ind_next]))

    p = [[], [], []]
    q = [[], [], []]
    v = [[], [], []]

    for i in ind_alpha_0:
        ind = ind_next[i[0]]
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

        n_steps = math.ceil(max(abs(p_diff[g][ind]) for g in range(len(p_diff))) / stepsize)
        p_sp = np.array([np.linspace(boundary_P_values[g][0][ind], boundary_P_values[g][0][ind + 1], n_steps) for g in
                         range(len(boundary_P_values))])

        opf.model.alpha = 0

        for i in range(len(p_sp[0])):
            for g in opf.model.eG:
                opf.model.p_sp[g] = p_sp[g, i]
            opf.solve(solver=solver)

            if check_optimal_termination(opf.results):
                for g in opf.model.eG:
                    p[g].append(value(opf.model.pG[g]))
                    q[g].append(value(opf.model.qG[g]))
                for b in opf.model.b0:
                    v[b].append(value(opf.model.v[b]))
            else:
                pass

    for i in range_ext_grids:
        boundary_P_values[i].append(p[i])
        boundary_Q_values[i].append(q[i])
        boundary_V_values[i].append(v[i])

    opf.model.p_eg_max.deactivate()
    opf.model.p_eg_min.deactivate()

    ind_beta_0 = np.argwhere(abs(p_diff[0][ind_next]) <= abs(q_diff[0][ind_next]))
    for i in ind_beta_0:
        ind = ind_next[i[0]]
        beta = alpha_beta[ind][1]
        alpha = alpha_beta[ind][0]

        n_steps = math.ceil(max(abs(q_diff[g][ind]) for g in range(len(q_diff))) / stepsize)

        q_sp = np.array([np.linspace(boundary_Q_values[g][0][ind], boundary_Q_values[g][0][ind + 1], n_steps) for g in
                         range(len(boundary_Q_values))])

        opf.model.alpha = alpha
        opf.model.beta = beta

        opf.model.q_eg_min.deactivate()
        opf.model.q_eg_max.deactivate()
        opf.solve(solver=solver)

        nets.append(copy.deepcopy(opf.net))

        opf.model.q_eg_min.activate()
        opf.model.q_eg_max.activate()

        opf.model.beta = 0

        for i in range(len(q_sp[0])):
            for g in opf.model.eG:
                opf.model.q_sp[g] = q_sp[g, i]
            opf.solve(solver=solver)
            for i in opf.model.eG:
                print(value(opf.model.pG[i]))
            if check_optimal_termination(opf.results):
                for g in opf.model.eG:
                    p[g].append(value(opf.model.pG[g]))
                    q[g].append(value(opf.model.qG[g]))
                for b in opf.model.b0:
                    v[b].append(value(opf.model.v[b]))
            else:
                pass

    for i in range(len(p)):
        boundary_P_values[i].append(p[i])
        boundary_Q_values[i].append(q[i])
        boundary_V_values[i].append(v[i])

    return boundary_P_values, boundary_Q_values, boundary_V_values, nets


def node_for_setpoint_based_with_directions(opf, w, stepsize=100):
    alpha_beta = [(1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1)]

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

    for (alpha, beta) in alpha_beta:
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
        # boundary_V_values[g][-1].append(boundary_V_values[g][-1][0])

    p_max = np.array([max(boundary_P_values[i][0]) for i in range_ext_grids])
    p_min = np.array([min(boundary_P_values[i][0]) for i in range_ext_grids])

    q_max = np.array([max(boundary_Q_values[i][0]) for i in range_ext_grids])
    q_min = np.array([min(boundary_Q_values[i][0]) for i in range_ext_grids])

    delta_p = abs(p_max - p_min)
    delta_q = abs(q_max - q_min)

    # ----
    delta_max = np.max((delta_p, delta_q))
    d_max = stepsize / delta_max

    p_diff = [np.diff(boundary_P_values[i][0]) for i in range_ext_grids]
    q_diff = [np.diff(boundary_Q_values[i][0]) for i in range_ext_grids]
    ds = [np.sqrt((p_diff[g] / delta_p[g]) ** 2 + (q_diff[g] / delta_q[g]) ** 2) for g in range_ext_grids]
    ind_next = [np.argwhere(ds[g] >= d_max) for g in range_ext_grids]
    ind_next = np.unique(np.concatenate(ind_next))
    # p_next = []
    # q_next = []
    # for g in range(len(boundary_P_values)):
    #     p_next.append([(boundary_P_values[g][-1][i] + boundary_P_values[g][-1][i + 1]) / 2 for i in ind_next])
    #     q_next.append([(boundary_Q_values[g][-1][i] + boundary_Q_values[g][-1][i + 1]) / 2 for i in ind_next])

    tol = 0.1

    # def _add_setpoint_constraints(model):
    #     model.p_sp = Param([w], mutable=True)
    #
    #     def p_eg_upper(model, g):
    #         return (model.psG[g]) <= (model.p_sp[g]) + tol
    #
    #     model.p_eg_max = Constraint([w], rule=p_eg_upper)
    #
    #     def p_eg_lower(model, g):
    #         return (model.psG[g]) >= (model.p_sp[g]) - tol
    #
    #     model.p_eg_min = Constraint([w], rule=p_eg_lower)
    #
    #     model.q_sp = Param([w], mutable=True)
    #
    #     def q_eg_upper(model, g):
    #         return (model.qsG[g]) <= (model.q_sp[g]) + tol
    #
    #     model.q_eg_max = Constraint([w], rule=q_eg_upper)
    #
    #     def q_eg_lower(model, g):
    #         return (model.qsG[g]) >= (model.q_sp[g]) - tol
    #
    #     model.q_eg_min = Constraint([w], rule=q_eg_lower)
    #
    # _add_setpoint_constraints(opf.model)

    opf.model.q_eg_min.deactivate()
    opf.model.q_eg_max.deactivate()

    ind_alpha_0 = np.argwhere(abs(p_diff[0][ind_next]) >= abs(q_diff[0][ind_next]))

    p = [[], [], []]
    q = [[], [], []]
    v = [[], [], []]

    for i in ind_alpha_0:
        ind = ind_next[i[0]]
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

        n_steps = math.ceil(max(abs(p_diff[g][ind]) for g in range(len(p_diff))) / stepsize)
        p_sp = np.array([np.linspace(boundary_P_values[g][0][ind], boundary_P_values[g][0][ind + 1], n_steps) for g in
                         range(len(boundary_P_values))])

        opf.model.alpha = 0

        for i in range(len(p_sp[0])):
            for g in range_ext_grids:
                opf.model.p_sp[w] = p_sp[g, i]
            opf.solve()

            if check_optimal_termination(opf.results):
                for g in range_ext_grids:
                    p[g].append(value(opf.model.psG[w]))
                    q[g].append(value(opf.model.qsG[w]))
                # for b in opf.model.b0:
                #     v[b].append(value(opf.model.v[b]))
            else:
                pass

    for i in range_ext_grids:
        boundary_P_values[i].append(p[i])
        boundary_Q_values[i].append(q[i])
        boundary_V_values[i].append(v[i])

    opf.model.p_eg_max.deactivate()
    opf.model.p_eg_min.deactivate()

    ind_beta_0 = np.argwhere(abs(p_diff[0][ind_next]) <= abs(q_diff[0][ind_next]))
    for i in ind_beta_0:
        ind = ind_next[i[0]]
        beta = alpha_beta[ind][1]
        alpha = alpha_beta[ind][0]

        n_steps = math.ceil(max(abs(q_diff[g][ind]) for g in range(len(q_diff))) / stepsize)

        q_sp = np.array([np.linspace(boundary_Q_values[g][0][ind], boundary_Q_values[g][0][ind + 1], n_steps) for g in
                         range(len(boundary_Q_values))])

        opf.model.alpha = alpha
        opf.model.beta = beta

        opf.model.q_eg_min.deactivate()
        opf.model.q_eg_max.deactivate()
        opf.solve()

        nets.append(copy.deepcopy(opf.net))

        opf.model.q_eg_min[w].activate()
        opf.model.q_eg_max[w].activate()

        opf.model.beta = 0

        for i in range(len(q_sp[0])):
            for g in range_ext_grids:
                opf.model.q_sp[w] = q_sp[g, i]
            opf.solve()

            if check_optimal_termination(opf.results):
                for g in range_ext_grids:
                    p[g].append(value(opf.model.psG[w]))
                    q[g].append(value(opf.model.qsG[w]))
                # for b in opf.model.b0:
                #     v[b].append(value(opf.model.v[b]))
            else:
                pass

    for i in range_ext_grids:
        boundary_P_values[i].append(p[i])
        boundary_Q_values[i].append(q[i])
        boundary_V_values[i].append(v[i])

    return boundary_P_values, boundary_Q_values, boundary_V_values, nets


def get_pq_to_plot(p, q):
    pq = []
    for i in range(len(p)):
        p_i = [p_mw for p_list in p[i] for p_mw in p_list]
        q_i = [q_mw for q_list in q[i] for q_mw in q_list]
        pq.append(np.array((p_i, q_i)).T)

    pq_tot = sum(pq_i for pq_i in pq)

    return pq, pq_tot


def plot_hull(p, q, ratio=0.1):
    pq = []
    for i in range(len(p)):
        p_i = [p_mw for p_list in p[i] for p_mw in p_list]
        q_i = [q_mw for q_list in q[i] for q_mw in q_list]
        pq.append(np.array((p_i, q_i)).T)

    pq_tot = sum(pq_i for pq_i in pq)

    clrs = ['#00549F', '#000000', '#E30066', '#FFED00', '#006165',
            '#0098A1', '#57AB27', '#BDCD00', '#F6A800', '#CC071E',
            '#A11035', '#612158', '#7A6FAC']

    fig, ax = plt.subplots()

    clr = clrs[5]
    # ax.plot(pq_tot[:, 0], pq_tot[:, 1], '.', label='Total', color=clr)
    # try:
    #     hull = concave_hull(MultiPoint(pq_tot), ratio)
    # except Exception as e:
    #     print('Creating convex hull instead of concave hull: ' + str(e))
    #     hull = convex_hull(MultiPoint(pq_tot))
    # polygon = gpd.GeoSeries(hull)
    # polygon.plot(ax=ax, alpha=0.2, color=clr)
    # polygon.boundary.plot(ax=ax, edgecolor=clr, linewidth=2)

    # alpha_shape = alphashape.alphashape(pq_tot, 0.001)
    # ax.add_patch(PolygonPatch(alpha_shape, alpha=1, fill=False, edgecolor=clr, linewidth=2, color=clr))
    # ax.add_patch(PolygonPatch(alpha_shape, alpha=0.2, edgecolor=clr, linewidth=2, color=clr))
    figs = [fig]
    axs = [ax]
    for i, points in enumerate(pq):
        fig_i, ax_i = plt.subplots()

        clr = clrs[i + 6]
        ax.plot(points[:, 0], points[:, 1], '.', label='Slack #' + str(i), color=clr)
        ax.legend()

        # plot convex hull for slack alone
        ax_i.plot(points[:, 0], points[:, 1], '.', label='Slack #' + str(i), color=clr)

        try:
            hull = concave_hull(MultiPoint(points), ratio)
        except Exception as e:
            print('Creating convex hull instead of concave hull: ' + str(e))
            hull = convex_hull(MultiPoint(points))

        polygon = gpd.GeoSeries(hull)
        polygon.plot(ax=ax, alpha=0.2, color=clr)
        polygon.boundary.plot(ax=ax, edgecolor=clr, linewidth=2)

        polygon.plot(ax=ax_i, alpha=0.2, color=clr)
        polygon.boundary.plot(ax=ax_i, edgecolor=clr, linewidth=2)
        ax_i.set_xlabel('P [MW]')
        ax_i.set_ylabel('Q [MVar]')
        figs.append(fig_i)
        axs.append(ax_i)
        # hull = ConvexHull(pq[i])
        # for simplex in hull.simplices:
        #     plt.plot(points[simplex, 0], points[simplex, 1], color=clr)
        # plt.legend()

        # alpha_shape = alphashape.alphashape(pq[i], 0.001)
        # ax.add_patch(PolygonPatch(alpha_shape, alpha=1, fill=False, edgecolor=clr, linewidth=2, color=clr))
        # ax.add_patch(PolygonPatch(alpha_shape, alpha=0.2, edgecolor=clr, linewidth=2, color=clr))

        ax.set_xlabel('P [MW]')
        ax.set_ylabel('Q [MVar]')

    config = set_plt_config()
    fig.set_size_inches((config['textbreite'] * 0.8, 0.6 * config['textbreite']))
    ax.set_aspect('equal')
    # plt.legend(bbox_to_anchor=(1.03, 1.0), loc='upper left')
    # plt.tight_layout()

    for axes in axs:
        axes.grid()
        axes.set_aspect('equal')
        # axes.set_ylim(-240, 415)
        # axes.set_xlim(-750, 30)

    # axs[0].set_ylim(-700, 1200)
    # axs[0].set_xlim(-2140, 80)

    # dir = 'C:\\Users\\f.lohse\PycharmProjects\potpourri\potpourri\\results\\for\\sb_hv_grid_with_potential_3MW_230m\\'
    # for i, figure in enumerate(figs):
    #     figure.set_size_inches((config['textbreite'] * 0.6, 0.4 * config['textbreite']))
    #     figure.savefig(dir + 'for_'+ str(i) +'.pdf', format='pdf', bbox_inches='tight')

    return figs, axs


if __name__ == "__main__":
     net = pp.networks.create_cigre_network_mv()

     with open('potpourri/data/windpot/sb_hv_grid_with_potential_3MW_230m.pkl',
               'rb') as f:
         net = pickle.load(f)

    net = sb.get_simbench_net("1-HV-mixed--0-no_sw")

    case = 'lW'
    factors = net.loadcases.loc[case]
    net.load.p_mw *= factors['pload']
    net.load.q_mvar *= factors['qload']
    net.sgen.scaling[net.sgen.type == 'Wind'] = factors['Wind_p']
    net.sgen.scaling[net.sgen.type == 'PV'] = factors['PV_p']
    net.sgen.scaling[(net.sgen.type != 'Wind') & (net.sgen.type != 'Solar')] = factors['RES_p']
    net.ext_grid.vm_pu = factors['Slack_vm']

    hc = HC_ACOPF(net)
    hc.solve()
    hc.add_OPF(SWmin=10)
    hc.add_tap_changer_linear()
    hc.solve(solver='mindtpy')

    with open('misc/LM_Sans_10/lmsans10-regular.otf',
              'rb') as f:
        net_wind = pickle.load(f)
    net_wind.ext_grid.vm_pu = factors['Slack_vm']

    input_hc_net_dir ='potpourri/data/windpot/sb_hv_grid_with_potential_3MW_230m.pkl'
    #'input_hc_net_dir = sb.get_simbench_net("1-HV-mixed--0-no_sw") #test input
    with open (input_hc_net_dir, 'rb') as f:
        net_hc = pickle.load(f)
    wind_hc_index = net_hc.sgen.index[net_hc.res_sgen.y_wind == 1]
    pp.create_sgens(net_wind, net_hc.sgen.bus[wind_hc_index], p_mw=net_hc.sgen.p_mw[wind_hc_index], var_q=0, type='Wind', wind_hc=True)

    net_wind = copy.deepcopy(net)
    # create wind generators in original net
    # wind_hc_index = hc.net.sgen.index[hc.net.res_sgen.y_wind == 1]
    # pp.create_sgens(net_wind, hc.net.sgen.bus[wind_hc_index], p_mw=hc.net.sgen.p_mw[wind_hc_index], var_q=0,
    #                 type='Wind', wind_hc=True)

    net_wind.sgen['wind_hc'].fillna(False, inplace=True)

    # -- ext grids for
    net_wind.sgen['controllable'] = True
    net_wind.sgen['controllable'][net_wind.sgen.type == 'Wind'] = True
    net_wind.sgen['p_inst_mw'] = net_wind.sgen['p_mw']
    net_wind.sgen['var_q'][net_wind.sgen.type == 'Wind'] = 1
    net_wind.sgen['var_q'][net_wind.sgen.wind_hc] = 0
    net_wind.load['controllable'] = True

    acopf = ACOPF(net_wind)
    acopf.add_OPF()
    acopf.add_tap_changer_linear()

    p, q, u, nets = for_setpoint_based_with_directions(acopf, stepsize=40)

    #
    # # -- node for
    # net_wind.sgen['controllable'] = False
    # net_wind.sgen['controllable'][(net_wind.sgen.type == 'Wind') & ~net_wind.sgen.wind_hc] = True
    # net_wind.sgen['p_inst_mw'] = net_wind.sgen['p_mw']
    # net_wind.sgen['var_q'][(net_wind.sgen.type == 'Wind') & ~net_wind.sgen.wind_hc] = 1
    #
    # hc_node = HC_ACOPF(net_wind)
    # hc_node.add_OPF()
    #
    # hc_node.fix_vars('y', 0)
    # hc_node.model.obj.deactivate()
    #
    # p_nodes = []
    # q_nodes = []
    # u_nodes = []
    # nets_nodes = []
    #
    # def _add_setpoint_constraints(model):
    #     model.p_sp = Param(model.WIND_HC, mutable=True)
    #
    #     tol = 0.1
    #     def p_eg_upper(model, g):
    #         return (model.psG[g]) <= (model.p_sp[g]) + tol
    #
    #     model.p_eg_max = Constraint(model.WIND_HC, rule=p_eg_upper)
    #
    #     def p_eg_lower(model, g):
    #         return (model.psG[g]) >= (model.p_sp[g]) - tol
    #
    #     model.p_eg_min = Constraint(model.WIND_HC, rule=p_eg_lower)
    #
    #     model.q_sp = Param(model.WIND_HC, mutable=True)
    #
    #     def q_eg_upper(model, g):
    #         return (model.qsG[g]) <= (model.q_sp[g]) + tol
    #
    #     model.q_eg_max = Constraint(model.WIND_HC, rule=q_eg_upper)
    #
    #     def q_eg_lower(model, g):
    #         return (model.qsG[g]) >= (model.q_sp[g]) - tol
    #
    #     model.q_eg_min = Constraint(model.WIND_HC, rule=q_eg_lower)
    #
    #
    #
    # _add_setpoint_constraints(hc_node.model)
    #
    # for w in hc_node.model.WIND_HC:
    #     hc_node.model.p_eg_max.deactivate()
    #     hc_node.model.p_eg_min.deactivate()
    #     hc_node.model.q_eg_max.deactivate()
    #     hc_node.model.q_eg_min.deactivate()
    #     hc_node.model.y[w].fix(1)
    #
    #     p, q, u, nets = node_for_setpoint_based_with_directions(hc_node, w, stepsize=60)
    #     plot_hull(p, q)
    #
    #     p_nodes.append(p)
    #     q_nodes.append(q)
    #     u_nodes.append(u)
    #     nets_nodes.append(nets)
    #
    #     hc_node.model.y[w].fix(0)
    #     hc_node.model.psG[w] = 0
    #     hc_node.model.qsG[w] = 0
    #



    # p, q = for_setpoint_based(hc_for, n=9)
#    p, q = run_feasible_operation_region(hc_for)

    # p_tot_slack = value(sum(hc.model.pG[g] for g in hc.model.eG))
    # q_tot_slack = value(sum(hc.model.qG[g] for g in hc.model.eG))
    #
    # p_cases = []
    # q_cases = []
    # u_cases = []
    # nets_cases = []
    # for i in ['hL', 'hW', 'hPV', 'lW', 'lPV']:
    #
    #     factors = net.loadcases.loc[i]
    #
    #     net_case = copy.deepcopy(net_wind)
    #
    #     net_case.load.p_mw *= factors['pload']
    #     net_case.load.q_mvar *= factors['qload']
    #     net_case.sgen.p_mw[net_case.sgen.type == 'Wind'] *= factors['Wind_p']
    #     net_case.sgen.p_mw[net_case.sgen.wind_hc] *= factors['Wind_p']
    #     net_case.sgen.p_mw[net_case.sgen.type == 'PV'] *= factors['PV_p']
    #     net_case.sgen.p_mw[(net_case.sgen.type != 'Wind') & (net_case.sgen.type != 'Solar') & (net_case.sgen.wind_hc == False)] *= factors['RES_p']
    #
    #     acopf = ACOPF(net_case)
    #     acopf.add_OPF()
    #
    #     p, q, u, nets = for_setpoint_based_iterative(acopf, stepsize=120)
    #     p_cases.append(p)
    #     q_cases.append(q)
    #     u_cases.append(u)
    #     nets_cases.append(nets)
    #
    # for i in range(len(p_cases)):
    #     plot_hull(p_cases[i], q_cases[i])
    #     plot_qu_res(nets_cases[i])
    #     plot_pq_res(nets_cases[i])

# p, q, u = for_angle_based_sampling(hc_for, n=360)
#  model = feasible_operation_region(net)
# fig, ax = plt.subplots()
# for g in hc_for.model.eG:
#     ax.plot(p[g], q[g], '.-', label = 'Slack #' + str(g))
# ax.legend()
# ax.set_xlabel('P [MW]')
# ax.set_ylabel('Q [MVAr]')
