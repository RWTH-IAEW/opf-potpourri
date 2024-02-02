import copy

import pandas as pd
import simbench as sb
import numpy as np
import pyomo.environ as pe
import matplotlib.pyplot as plt
from potpourri.models.class_based.HC_ACOPF import HC_ACOPF
from potpourri.scripts.classbased.plot_functions import plot_wind_hc_results


def max_wind_min_loss(net, **kwargs):
    hc = HC_ACOPF(net, SWmax=10000, SWmin=0, peGmax=1000000)
    hc.solve()
    hc.add_OPF()
    if kwargs.get('tap_linear', False):
        hc.add_tap_changer_linear()

    if kwargs.get('tap_discrete', False):
        hc.add_tap_changer_discrete()

    hc.add_loss_obj()

    obj = []
    p_wind = []
    p_loss = []

    for i in range(11):
        hc.model.eps = 1 - i / 10
        hc.solve(solver='mindtpy', mip_solver='gurobi')

        if pe.check_optimal_termination(hc.results):
            obj.append(pe.value(hc.model.OBJ_with_loss))
            p_wind.append(pe.value(sum(hc.model.pG[w] for w in hc.model.WIND)))
            p_loss.append(pe.value(sum(hc.model.pLfrom[l] + hc.model.pLto[l] for l in hc.model.L)))
        else:
            obj.append(None)
            p_wind.append(None)
            p_loss.append(None)

    fig, ax = plt.subplots()
    ax.scatter(p_loss, p_wind)

    return obj, p_wind, p_loss


def weighted_sum(net, **kwargs):
    hc = HC_ACOPF(net, SWmax=10000, SWmin=0, peGmax=1000000)
    hc.solve()
    hc.add_OPF()

    if kwargs.get('tap_linear', False):
        hc.add_tap_changer_linear()

    if kwargs.get('tap_discrete', False):
        hc.add_tap_changer_discrete()

    if kwargs.get('SWmin', False):
        hc.change_vals('SWmin', kwargs.get('SWmin'))

    hc.solve(solver='mindtpy', mip_solver='gurobi')

    p_wind_max = pe.value(hc.model.OBJ)

    hc.model.w_wind = pe.Param(within=pe.Reals, initialize=0., mutable=True)
    hc.model.w_loss = pe.Param(within=pe.Reals, initialize=1., mutable=True)

    def objective_pwind_loss(model):
        return model.w_wind * sum(model.pG[w] for w in model.WIND) + model.w_loss * (
            - sum(model.pLfrom[l] + model.pLto[l] for l in model.L))

    hc.model.OBJ.deactivate()
    hc.model.OBJ_with_loss = pe.Objective(rule=objective_pwind_loss, sense=pe.maximize)

    hc.solve(solver='mindtpy', mip_solver='gurobi')
    p_loss_min = -pe.value(hc.model.OBJ_with_loss)

    obj = []
    p_wind = []
    p_loss = []
    hcs = []

    for i in range(21):
        hc.model.w_wind = i / 20 / p_wind_max
        hc.model.w_loss = (1 - i / 20) / p_loss_min

        hc.solve(solver='mindtpy', mip_solver='gurobi')

        if pe.check_optimal_termination(hc.results):
            obj.append(pe.value(hc.model.OBJ_with_loss))
            p_wind.append(pe.value(sum(hc.model.pG[w] for w in hc.model.WIND)))
            p_loss.append(pe.value(sum(hc.model.pLfrom[l] + hc.model.pLto[l] for l in hc.model.L)))
            hcs.append(copy.deepcopy(hc))
        else:
            obj.append(None)
            p_wind.append(None)
            p_loss.append(None)
            hcs.append(None)

    return obj, p_wind, p_loss


def epsilon_constraint(net, n=10, plot_hc=False, **kwargs):
    hc = HC_ACOPF(net, SWmax=10000, SWmin=0, peGmax=1000000)
    hc.solve()
    hc.add_OPF()

    if kwargs.get('tap_linear', False):
        hc.add_tap_changer_linear()

    if kwargs.get('tap_discrete', False):
        hc.add_tap_changer_discrete()

    if kwargs.get('SWmin', False):
        hc.change_vals('SWmin', kwargs.get('SWmin'))

    hc.solve(solver='mindtpy', mip_solver='gurobi')

    p_wind_max = pe.value(hc.model.OBJ)
    p_loss_max = pe.value(sum(hc.model.pLfrom[l] + hc.model.pLto[l] for l in hc.model.L))
    obj_max = pe.value(hc.model.OBJ)

    hc.model.OBJ.deactivate()

    # get minimum losses
    def obj_min_loss(model):
        return sum(model.pLfrom[l] + model.pLto[l] for l in model.L)
    hc.model.OBJ_loss = pe.Objective(rule=obj_min_loss, sense=pe.minimize)

    hc.solve(solver='mindtpy', mip_solver='gurobi')
    p_loss_min = pe.value(hc.model.OBJ_loss)
    p_wind_min = pe.value(sum(hc.model.pG[w] for w in hc.model.WIND))

    hc.model.OBJ_loss.deactivate()
    hc.model.OBJ.activate()

    # add loss objective function as constraint
    hc.model.eps = pe.Param(within=pe.Reals, initialize=p_loss_min, mutable=True)
    hc.model.constr_loss = pe.Constraint(
        expr=sum(hc.model.pLfrom[l] + hc.model.pLto[l] for l in hc.model.L) <= hc.model.eps)

    # generate pareto front
    obj, p_wind, p_loss, nets = [], [], [], []

    step = (p_loss_max - p_loss_min) / n
    steps = np.arange(p_loss_min, p_loss_max, step)

    while step > 0.1:
        obj_eps, p_wind_eps, p_loss_eps, nets_eps = epsilon_constraint2(hc, step, steps)

        obj += obj_eps
        p_wind += p_wind_eps
        p_loss += p_loss_eps
        nets += nets_eps

        step = (p_loss_eps[-2] - p_loss_eps[-3]) / n
        steps = np.arange(p_loss_eps[-3], p_loss_eps[-2], step)

    p_wind.append(p_wind_max)
    p_loss.append(p_loss_max)
    obj.append(obj_max)

    results = pd.DataFrame({'obj': obj, 'p_wind': p_wind, 'p_loss': p_loss})

    # # plot pareto front
    # fig, ax = plt.subplots()
    # ax.plot(p_loss, p_wind, 'o-.')
    # ax.set_xlabel('$P_{losses}$ [MW]')
    # ax.set_ylabel('$P_{wind}$ [MW]')

    return obj, p_wind, p_loss, nets, results


def epsilon_constraint2(hc, step, steps):
    # step = (p_max - p_min) / n
    # steps = np.arange(p_min, p_max, step)

    obj = []
    p_wind = []
    p_loss = []
    nets = []

    for i in steps:
        hc.model.eps = i
        hc.solve(solver='mindtpy', mip_solver='gurobi')

        if pe.check_optimal_termination(hc.results) | (
                hc.results.solver.termination_condition == pe.TerminationCondition.feasible):
            obj.append(pe.value(hc.model.OBJ))
            p_wind.append(pe.value(sum(hc.model.pG[w] for w in hc.model.WIND)))
            p_loss.append(pe.value(sum(hc.model.pLfrom[l] + hc.model.pLto[l] for l in hc.model.L)))
            nets.append(copy.deepcopy(hc.net))

            # calculate difference between last two values
            if len(p_wind) > 1:
                p_wind_diff = p_wind[-1] - p_wind[-2]
                if abs(p_wind_diff - step) < 1e-3:
                    break

        else:
            hc.model.OBJ.display()
            obj.append(None)
            p_wind.append(None)
            p_loss.append(None)
            nets.append(None)

    return obj, p_wind, p_loss, nets


if __name__ == '__main__':
    net = sb.get_simbench_net("1-HV-mixed--0-no_sw")

    case = 'hL'

    factors = net.loadcases.loc[case]
    net.load.p_mw *= factors['pload']
    net.load.q_mvar *= factors['qload']
    net.sgen.scaling[net.sgen.type == 'Wind'] = factors['Wind_p']
    net.sgen.scaling[net.sgen.type == 'PV'] = factors['PV_p']
    net.sgen.scaling[(net.sgen.type != 'Wind') & (net.sgen.type != 'Solar')] = factors['RES_p']
    net.ext_grid.vm_pu = factors['Slack_vm']

    # obj_eps, p_wind_eps, p_loss_eps, hcs_eps = epsilon_constraint(net)
    # obj_ws, p_wind_ws, p_loss_ws = weighted_sum(net)
    # #
    # fig, ax = plt.subplots()
    # ax.scatter(p_loss_eps, p_wind_eps, label='Epsilon Constraint')
    # ax.scatter(p_loss_ws, p_wind_ws, label='Weighted Sum')
    # ax.legend()
    # ax.set_xlabel('$P_{losses}$ [MW]')
    # ax.set_ylabel('$P_{wind}$ [MW]')

    # p_loss = p_loss_ws + p_loss_eps
    # p_wind = p_wind_ws + p_wind_eps
    # ax.plot(p_loss, p_wind, 'o-.')
