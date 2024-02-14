import copy
import pickle
import pandas as pd
import numpy as np
import pyomo.environ as pe
import matplotlib.pyplot as plt
from tqdm import tqdm

from potpourri.models.class_based.HC_ACOPF import HC_ACOPF
from potpourri.scripts.classbased.plot_functions import plot_wind_hc_results


def max_wind_min_loss(net, n=11, **kwargs):
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

    for i in range(n):
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
    ax.plot(p_loss, p_wind, 'o-.')

    return obj, p_wind, p_loss


def weighted_sum(net, n=11, **kwargs):
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

    for i in range(n):
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


def pareto_front(net, n=10, **kwargs):
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
    p_wind, p_loss, nets = [], [], []

    step = (p_loss_max - p_loss_min) / n
    steps = np.linspace(p_loss_max, p_loss_min, n)

    p_wind, p_loss, nets, p_wind_opt, p_loss_opt = epsilon_constraint(hc, step, steps)

    # p_wind = np.append(p_wind, p_wind_max)
    # p_loss = np.append(p_loss, p_loss_max)

    results = pd.DataFrame({'p_wind': p_wind, 'p_loss': p_loss})

    # # plot pareto front
    # fig, ax = plt.subplots()
    # ax.plot(p_loss, p_wind, 'o-.')
    # ax.set_xlabel('$P_{losses}$ [MW]')
    # ax.set_ylabel('$P_{wind}$ [MW]')

    return p_wind, p_loss, nets, results


def p_wind_loss_opt(net, n=10, **kwargs):
    hc = HC_ACOPF(net, SWmax=10000, SWmin=0, peGmax=1000000)
    hc.solve()
    hc.add_OPF()

    if kwargs.get('tap_linear', False):
        hc.add_tap_changer_linear()

    # if kwargs.get('tap_discrete', False):
    #     hc.add_tap_changer_discrete()

    if kwargs.get('SWmin', False):
        hc.change_vals('SWmin', kwargs.get('SWmin'))

    hc.solve(solver='mindtpy', mip_solver='gurobi')

    p_wind_max = pe.value(hc.model.OBJ)
    p_loss_max = pe.value(sum(hc.model.pLfrom[l] + hc.model.pLto[l] for l in hc.model.L))

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

    if kwargs.get('tap_discrete', False):
        hc.add_tap_changer_discrete()

    # add loss objective function as constraint
    hc.model.eps = pe.Param(within=pe.Reals, initialize=p_loss_min, mutable=True)
    hc.model.constr_loss = pe.Constraint(
        expr=sum(hc.model.pLfrom[l] + hc.model.pLto[l] for l in hc.model.L) <= hc.model.eps)

    # generate pareto front
    p_wind, p_loss, nets = [[p_wind_min]], [[p_loss_min]], [copy.deepcopy(hc.net)]
    p_wind_opt = p_wind_min
    p_loss_opt = p_loss_min

    p_wind_next = p_wind_max
    p_loss_next = p_loss_max

    step = (p_loss_max - p_loss_min) / (n-1)
    steps = np.linspace(p_loss_min+step, p_loss_max-step, (n-2))

    # prev_p_wind_opt = p_wind_min
    # prev_p_loss_opt = p_loss_min
    # tolerance = 1e-5  # Define your tolerance here
    # same_value_counter = 0

    while step > 0.005:
        tqdm.write(str(steps))

        p_w = [p_wind_opt, p_wind_next]
        p_l = [p_loss_opt, p_loss_next]

        p_wind_eps, p_loss_eps, nets_eps, p_wind_opt, p_loss_opt, p_wind_next, p_loss_next = epsilon_constraint(hc, steps, p_w, p_l, mode='opt')

        nets.append(nets_eps)
        p_wind.append(p_wind_eps)
        p_loss.append(p_loss_eps)

        if p_wind_opt is None:
            p_wind_opt = p_wind_eps[-1]
            p_loss_opt = p_loss_eps[-1]

        # # If the new values are the same as the previous ones within the tolerance, increment the counter
        # if np.isclose(p_wind_opt, prev_p_wind_opt, atol=tolerance) and np.isclose(p_loss_opt, prev_p_loss_opt, atol=tolerance):
        #     same_value_counter += 1
        #     # If the counter reaches 3, break the loop
        #     if same_value_counter == 5:
        #         break
        # else:
        #     same_value_counter = 0  # Reset counter if values are different
        #
        # # Update the previous values
        # prev_p_wind_opt = p_wind_opt
        # prev_p_loss_opt = p_loss_opt

        # # Get the next smaller value from p_loss_eps
        # smaller_elements = p_loss_eps[p_loss_eps < p_loss_opt]
        # next_smaller_value = np.max(smaller_elements) if smaller_elements.size > 0 else p_loss_opt

        # Calculate the new steps
        step = abs(p_loss_opt - p_loss_next) / (n-1)
        steps = np.linspace(p_loss_opt - step, p_loss_next + step, n - 2)

    p_wind.append(p_wind_max)
    p_loss.append(p_loss_max)

    results = pd.DataFrame({'p_wind': p_wind, 'p_loss': p_loss})

    # # plot pareto front
    # fig, ax = plt.subplots()
    # ax.plot(p_loss, p_wind, 'o-.')
    # ax.set_xlabel('$P_{losses}$ [MW]')
    # ax.set_ylabel('$P_{wind}$ [MW]')

    return p_wind, p_loss, nets, results, p_wind_opt, p_loss_opt


def epsilon_constraint(hc, steps, p_w, p_l, mode=None):
    p_wind = [p_w[0]]
    p_loss = [p_l[0]]
    nets = []

    p_wind_opt = None
    p_loss_opt = None

    for i in tqdm(steps):
        hc.model.eps = i
        hc.solve(solver='mindtpy', mip_solver='gurobi')

        p_wind.append(pe.value(sum(hc.model.pG[w] for w in hc.model.WIND)))
        p_loss.append(pe.value(sum(hc.model.pLfrom[l] + hc.model.pLto[l] for l in hc.model.L)))
        nets.append(copy.deepcopy(hc.net))

        if not mode:
            continue
        # calculate difference between last two values
        if len(p_wind) > 1:
            p_wind_diff = p_wind[-1] - p_wind[-2]
            p_loss_diff = p_loss[-1] - p_loss[-2]
            if p_wind_diff - p_loss_diff <= -1e-3:
                p_wind_opt = p_wind[-2]
                p_loss_opt = p_loss[-2]

                if p_wind_diff <= 0:
                    p_wind_next = p_wind[-1]
                    p_loss_next = p_loss[-1]
                else:
                    p_wind_next = p_wind[-3]
                    p_loss_next = p_loss[-3]

                p_wind.pop(0)
                p_loss.pop(0)

                return np.array(p_wind), np.array(p_loss), nets, p_wind_opt, p_loss_opt, p_wind_next, p_loss_next

    if mode:
        p_wind_next = None
        p_loss_next = None

        if len(p_wind) > 1:
            p_wind_diff = p_w[-1] - p_wind[-1]
            p_loss_diff = p_l[-1] - p_loss[-1]

            if p_wind_diff - p_loss_diff <= -1e-3:
                p_wind_opt = p_wind[-1]
                p_loss_opt = p_loss[-1]

                if p_wind_diff <= 0:
                    p_wind_next = p_w[-1]
                    p_loss_next = p_l[-1]
                else:
                    p_wind_next = p_wind[-2]
                    p_loss_next = p_loss[-2]

        p_wind.pop(0)
        p_loss.pop(0)

        return np.array(p_wind), np.array(p_loss), nets, p_wind_opt, p_loss_opt, p_wind_next, p_loss_next

    return np.array(p_wind), np.array(p_loss), nets, p_wind_opt, p_loss_opt


if __name__ == '__main__':
    with open('C:/Users/f.lohse/PycharmProjects/potpourri/potpourri/data/simbench_hv_grid_with_potential_pkl.pkl', 'rb') as f:
        net = pickle.load(f)

    # net = sb.get_simbench_net("1-HV-mixed--0-no_sw")

    case = 'lW'

    factors = net.loadcases.loc[case]
    net.load.p_mw *= factors['pload']
    net.load.q_mvar *= factors['qload']
    net.sgen.scaling[net.sgen.type == 'Wind'] = factors['Wind_p']
    net.sgen.scaling[net.sgen.type == 'PV'] = factors['PV_p']
    net.sgen.scaling[(net.sgen.type != 'Wind') & (net.sgen.type != 'Solar')] = factors['RES_p']
    net.ext_grid.vm_pu = factors['Slack_vm']
