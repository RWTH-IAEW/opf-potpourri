import copy
import pickle
import pandas as pd
import numpy as np
import pyomo.environ as pe
import matplotlib.pyplot as plt
from tqdm import tqdm

from potpourri.models.class_based.HC_ACOPF import HC_ACOPF
from potpourri.scripts.classbased.plot_functions import plot_wind_hc_results
from potpourri.scripts.classbased.run_scenarios import create_output_writer


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

    if kwargs.get('tap_discrete', False):
        hc.add_tap_changer_discrete()

    hc.model.OBJ_loss.deactivate()
    hc.model.OBJ.activate()

    # add loss objective function as constraint
    hc.model.eps = pe.Param(within=pe.Reals, initialize=p_loss_min, mutable=True)
    hc.model.constr_loss = pe.Constraint(
        expr=sum(hc.model.pLfrom[l] + hc.model.pLto[l] for l in hc.model.L) <= hc.model.eps)

    # generate pareto front
    p_w = [p_wind_min, p_wind_max]
    p_l = [p_loss_min, p_loss_max]

    step = (p_loss_max - p_loss_min) / (n - 1)
    steps = np.linspace(p_loss_min + step, p_loss_max - step, (n - 2))

    p_wind, p_loss, nets, p_wind_opt, p_loss_opt = epsilon_constraint(hc, steps, p_w, p_l)

    results = pd.DataFrame({'p_wind': p_wind, 'p_loss': p_loss})

    return p_wind, p_loss, nets, results


def p_wind_loss_opt(net, n=10, **kwargs):
    if kwargs.get('hc', False):
        hc = kwargs.get('hc')
    else:
        hc = HC_ACOPF(net, SWmax=10000, SWmin=0, peGmax=1000000)
        hc.solve()
        hc.add_OPF()

    if kwargs.get('tap_linear', False):
        hc.add_tap_changer_linear()

    if kwargs.get('SWmin', False):
        hc.change_vals('SWmin', kwargs.get('SWmin'))

    try:
        hc.model.constr_loss.deactivate()
    except AttributeError:
        pass

    hc.solve(solver='mindtpy', mip_solver='gurobi')

    p_wind_max = pe.value(hc.model.OBJ)
    p_loss_max = pe.value(sum(hc.model.pLfrom[l] + hc.model.pLto[l] for l in hc.model.L))

    hc.model.OBJ.deactivate()

    try:
        hc.model.OBJ_loss.activate()
    except AttributeError:
        def obj_min_loss(model):
            return sum(model.pLfrom[l] + model.pLto[l] for l in model.L)

        hc.model.OBJ_loss = pe.Objective(rule=obj_min_loss, sense=pe.minimize)

    # get minimum losses
    hc.solve(solver='mindtpy', mip_solver='gurobi')
    p_loss_min = pe.value(hc.model.OBJ_loss)
    p_wind_min = pe.value(sum(hc.model.pG[w] for w in hc.model.WIND))

    hc.model.OBJ_loss.deactivate()
    hc.model.OBJ.activate()

    if kwargs.get('tap_discrete', False):
        hc.add_tap_changer_discrete()

    # add loss objective function as constraint
    try:
        hc.model.constr_loss.activate()
    except AttributeError:
        hc.model.eps = pe.Param(within=pe.Reals, initialize=p_loss_min, mutable=True)
        hc.model.constr_loss = pe.Constraint(
            expr=sum(hc.model.pLfrom[l] + hc.model.pLto[l] for l in hc.model.L) <= hc.model.eps)

    # generate pareto front
    p_wind, p_loss, nets = [[p_wind_min]], [[p_loss_min]], [[copy.deepcopy(hc.net)]]
    p_wind_opt = p_wind_min
    p_loss_opt = p_loss_min

    p_wind_next = p_wind_max
    p_loss_next = p_loss_max

    step = (p_loss_max - p_loss_min) / (n - 1)
    steps = np.linspace(p_loss_min + step, p_loss_max - step, (n - 2))

    # create output writer
    if kwargs.get('output_dir', False):
        output_dir = kwargs.get('output_dir')
        output_path = output_dir
        ow = create_output_writer(hc.net, range(100), output_path)
        ow.init_all(hc.net)
    else:
        ow = None

    step_change_tolerance = 1e-4
    step_change_counter = 0
    prev_step = None

    while step > 0.005:
        tqdm.write(str(steps))

        p_w = [p_wind_opt, p_wind_next]
        p_l = [p_loss_opt, p_loss_next]

        p_wind_eps, p_loss_eps, nets_eps, p_wind_opt, p_loss_opt, p_wind_next, p_loss_next, net_opt = epsilon_constraint(
            hc, steps, p_w, p_l, mode='opt', ow=ow)

        nets.append(nets_eps)
        p_wind.append(p_wind_eps)
        p_loss.append(p_loss_eps)

        if p_wind_opt is None:
            p_wind_opt = p_w[-1]
            p_loss_opt = p_l[-1]
            p_loss_next = p_loss_opt + step * 10
            p_wind_next = np.NaN
            net_opt = nets[-1][-1]
            # p_wind_opt = p_wind_eps[-1]
            # p_loss_opt = p_loss_eps[-1]
            # net_opt = nets_eps[-1]

        # fig, ax = plt.subplots()
        # for i in range(len(p_wind)):
        #     ax.plot(p_loss[i], p_wind[i], 'o-.')

        # Calculate the new steps
        step = abs(p_loss_opt - p_loss_next) / (n - 1)
        steps = np.linspace(p_loss_opt - step, p_loss_next + step, n - 2)

        # Check if step size change is within tolerance
        if prev_step is not None and abs(step - prev_step) <= step_change_tolerance:
            step_change_counter += 1
            # Break loop if step size change is within tolerance for 3 consecutive iterations
            if step_change_counter >= 3:
                break
        else:
            step_change_counter = 0

        # Update previous step size
        prev_step = step

    total_elements = sum([len(p) for p in p_wind])

    if kwargs.get('output_dir', False):
        ow.dump(hc.net)

    p_wind.append(p_wind_max)
    p_loss.append(p_loss_max)

    results = pd.DataFrame({'p_wind': p_wind, 'p_loss': p_loss})

    net = net_opt if net_opt else nets[-1][-1]

    return p_wind, p_loss, nets, results, p_wind_opt, p_loss_opt, net_opt


def epsilon_constraint(hc, steps, p_w, p_l, mode=None, ow=None):
    p_wind = [p_w[0]]
    p_loss = [p_l[0]]
    nets = []

    p_wind_opt = None
    p_loss_opt = None
    net_opt = None

    if ow:
        offset_timestep = ow.time_step + 1 if ow.time_step else 0

    for i, eps in enumerate(tqdm(steps)):
        hc.model.eps = eps
        hc.solve(solver='mindtpy', mip_solver='gurobi')

        p_wind.append(pe.value(sum(hc.model.pG[w] for w in hc.model.WIND)))
        p_loss.append(pe.value(sum(hc.model.pLfrom[l] + hc.model.pLto[l] for l in hc.model.L)))
        nets.append(copy.deepcopy(hc.net))

        if ow:
            pf_converged = hc.results.solver.termination_condition == (
                    pe.TerminationCondition.optimal or pe.TerminationCondition.feasible)

            ow.time_step = i + 1 + offset_timestep
            ow.save_results(hc.net, i + offset_timestep, pf_converged=pf_converged,
                            ctrl_converged=pe.check_optimal_termination(hc.results))

        if not mode:
            continue
        # calculate difference between last two values
        if len(p_wind) > 1:
            p_wind_opt, p_loss_opt, p_wind_next, p_loss_next = check_slope(p_wind, p_loss)

            if p_wind_opt:
                p_wind.pop(0)
                p_loss.pop(0)
                if len(p_wind) > 1:
                    net_opt = nets[-2]
                else:
                    net_opt = nets[-1]
                return np.array(p_wind), np.array(
                    p_loss), nets, p_wind_opt, p_loss_opt, p_wind_next, p_loss_next, net_opt

    p_wind.append(p_w[-1])
    p_loss.append(p_l[-1])

    if mode:
        p_wind_opt, p_loss_opt, p_wind_next, p_loss_next = check_slope(p_wind, p_loss)

        p_wind.pop(0)
        p_loss.pop(0)
        p_wind.pop()
        p_loss.pop()
        if p_wind_opt:
            net_opt = nets[-1]

        return np.array(p_wind), np.array(p_loss), nets, p_wind_opt, p_loss_opt, p_wind_next, p_loss_next, net_opt

    return np.array(p_wind), np.array(p_loss), nets, p_wind_opt, p_loss_opt


def check_slope(p_wind, p_loss):
    p_wind_opt, p_loss_opt, p_wind_next, p_loss_next = None, None, None, None

    p_wind_diff = p_wind[-1] - p_wind[-2]
    p_loss_diff = p_loss[-1] - p_loss[-2]
    sign = np.sign(p_wind_diff)
    if sign * (p_wind_diff / p_loss_diff) <= sign * 1.001:
        p_wind_opt = p_wind[-2]
        p_loss_opt = p_loss[-2]

        if p_wind_diff <= 0:
            p_wind_next = p_wind[-1]
            p_loss_next = p_loss[-1]
        else:
            p_wind_next = p_wind[-3]
            p_loss_next = p_loss[-3]

    return p_wind_opt, p_loss_opt, p_wind_next, p_loss_next


if __name__ == '__main__':
    with open('C:/Users/f.lohse/PycharmProjects/potpourri/potpourri/data/simbench_hv_grid_with_potential_pkl.pkl',
              'rb') as f:
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

    # # plot pareto front
    # fig, ax = plt.subplots()
    # ax.plot(p_loss, p_wind, 'o-.')
    # ax.set_xlabel('$P_{losses}$ [MW]')
    # ax.set_ylabel('$P_{wind}$ [MW]')
