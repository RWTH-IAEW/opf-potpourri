import copy
import math
import time

import numpy as np
import pandapower as pp
import matplotlib.pyplot as plt
import pyomo.environ as pe
import simbench as sb
import random
from potpourri.models.class_based.HC_ACOPF import HC_ACOPF
from potpourri.models.class_based.DC import DC
from potpourri.models.class_based.AC import AC

from potpourri.models.class_based.pyo_to_net import pyo_sol_to_net_res


def create_testnet():
    net = pp.create_empty_network()

    pp.create_buses(net, 3, name=["b1", "b2", "b3"], vn_kv=110, max_vm_pu=1.118, min_vm_pu=0.9)
    pp.create_lines(net, [0, 0, 1], [1, 2, 2], name=["l1", "l2", "l3"], length_km=[3, 3, 3],
                    std_type='305-AL1/39-ST1A 110.0')
    pp.create_sgens(net, [0, 1, 2], p_mw=[10, 10, 10])
    pp.create_ext_grid(net, 0)
    pp.create_loads(net, [0, 1, 2], name=["load1", "load2", "load3"], p_mw=[10, 10, 10])
    # pp.createsgen(net, 0, p_mw=1, controllable=True)
    # pp.create_sgen(net, 0, p_mw=1, controllable=True, max_p_mw=10)
    # pp.create_load(net, 1, p_mw=1, controllable=True)
    # pp.create_load(net, 2, p_mw=2, controllable=True, max_p_mw=5)

    pp.create_sgens(net, net.bus.index, p_mw=0, wind_hc=True)

    # pp.create_shunt(net, 2, p_mw=3, q_mvar=5)
    return net


def init_pGW_qGW(hc_model, pGW=1, qGW=0):
    for w in hc_model.WIND:
        hc_model.pG[w] = pGW / hc_model.baseMVA
        hc_model.qG[w] = qGW / hc_model.baseMVA


def vary_SLmax(hc):
    res = {'p_wind_mw': {}, 'q_wind_mvar': {}, 'v_bus_pu': {}, 'va_bus_degree': {}, 'line_loading_%': {},
           'p_line_mw': {}, 'q_line_mvar': {}, 'p_slack_mw': {}, 'q_slack_mvar': {}}
    obj = []
    for i in range(10, 110, 10):
        init_pGW_qGW(hc.model, pGW=10, qGW=0)
        #        hc = HC_ACOPF(net)
        hc.set_SLmax(i, 'percent')
        print("max loading: " + str(i))
        hc.solve()

        pyo_sol_to_net_res(hc.net, hc.model)
        store_values(hc, res)
        obj.append(pe.value(hc.model.OBJ))

    return res


def diff_SLmax(SLmax, hc=None, SWmax=None, peGmax=None):
    if hc:
        hc_act = hc
        init_pGW_qGW(hc_act.model, pGW=10, qGW=0)
        if SWmax:
            for w in hc_act.model.WIND:
                hc_act.model.SWmax[w] = SWmax
    else:
        if SWmax:
            hc_act = HC_ACOPF(net, SWmax=SWmax)
        else:
            hc_act = HC_ACOPF(net)
    if peGmax:
        hc_act.model.peGmax[0] = peGmax
    hc_act.set_SLmax(SLmax, 'MW')
    return hc_act


def diff_peGmax(peGmax, hc=None, SWmax=None, SLmax=None):
    if hc:
        hc_act = hc
        init_pGW_qGW(hc_act.model, pGW=10, qGW=0)
        if SWmax:
            for w in hc_act.model.WIND:
                hc_act.model.SWmax[w] = SWmax
        hc_act.model.peGmax[0] = peGmax
    else:
        if SWmax:
            hc_act = HC_ACOPF(net, SWmax=SWmax, peGmax=peGmax)
        else:
            hc_act = HC_ACOPF(net, peGmax=peGmax)
    if SLmax:
        hc_act.set_SLmax(SLmax, 'MW')
    return hc_act


def diff_SWmax(SWmax, hc=None, SLmax=None, peGmax=None):
    if hc:
        hc_act = hc
        init_pGW_qGW(hc_act.model, pGW=10, qGW=0)
        for w in hc_act.model.WIND:
            hc_act.model.SWmax[w] = SWmax
    else:
        hc_act = HC_ACOPF(net, SWmax=SWmax)

    if SLmax:
        hc_act.set_SLmax(SLmax)
    if peGmax:
        hc_act.model.peGmax[0] = peGmax

    return hc_act


def vary_slack_and_SLmax(values=None, key=None, hc=None, SL_percent=[100], SWmax=None, SLmax=None, peGmax=None,
                         save_net: bool = False,
                         folder=None, solver='ipopt', mip_solver=None):
    results = []
    if not values and not key:
        values = [0]
        key = 'SL_%'

    for val in values:
        if key == 'SL_%':
            val = []
        res = {'p_wind_mw': {}, 'q_wind_mvar': {}, 'v_bus_pu': {}, 'va_bus_degree': {}, 'line_loading_%': {},
               'p_line_mw': {}, 'q_line_mvar': {}, 'p_slack_mw': {}, 'q_slack_mvar': {}}
        obj = []
        x = []
        model = []
        limits = {'peGmax': [], 'SWmax': [], 'SLmax': []}

        for i in SL_percent:
            if key == 'peGmax':
                hc_act = diff_peGmax(val, hc, SWmax, SLmax)
            elif key == 'SLmax':
                hc_act = diff_SLmax(val, hc, SWmax, peGmax)
            elif key == 'SWmax':
                hc_act = diff_SWmax(val, hc, SLmax, peGmax)
            else:
                if hc:
                    hc_act = hc
                else:
                    hc_act = HC_ACOPF(net)
                if key == 'SL_%':
                    val.append(i)

            hc_act.set_SLmax(i, 'percent')

            print("max loading: " + str(i))
            if mip_solver:
                hc_act.solve(solver=solver, mip_solver=mip_solver)
            else:
                hc_act.solve(solver=solver)

            pyo_sol_to_net_res(hc_act.net, hc_act.model)
            store_values(hc_act, res)
            obj.append(pe.value(hc_act.model.OBJ))
            x.append(i)
            model.append(hc_act)
            baseMVA = hc_act.model.baseMVA.value
            limits['peGmax'].append(hc_act.model.peGmax[0].value * baseMVA)
            limits['SWmax'].append(hc_act.model.SWmax[3].value * baseMVA)
            limits['SLmax'].append(hc_act.model.SLmax[0].value * baseMVA)

            if save_net:
                pp.to_excel(hc_act.net, folder + key + '_' + str(val) + '_SWmax_' + str(SWmax) + '_SLmax_' + str(
                    SLmax) + '_peGmax_' + str(peGmax) + '.xlsx')

        results.append(
            {'varied': {'key': key, 'val': val}, 'measurements': res, 'obj': obj, 'x': x, 'model': model,
             'limits': limits})

    return results


def store_values(hc, res):
    for w in hc.model.WIND:
        res['p_wind_mw'].setdefault(w, []).append(hc.net.res_sgen.p_mw[w])
        res['q_wind_mvar'].setdefault(w, []).append(hc.net.res_sgen.q_mvar[w])

    for b in hc.model.B:
        res['v_bus_pu'].setdefault(b, []).append(hc.net.res_bus.vm_pu[b])
        res['va_bus_degree'].setdefault(b, []).append(hc.net.res_bus.va_degree[b])

    for l in hc.model.L:
        res['line_loading_%'].setdefault(l, []).append(hc.net.res_line.loading_percent[l])
        res['p_line_mw'].setdefault(l, []).append(hc.net.res_line.pl_mw[l])
        res['q_line_mvar'].setdefault(l, []).append(hc.net.res_line.ql_mvar[l])

    for g in hc.model.eG:
        res['p_slack_mw'].setdefault(g, []).append(hc.net.res_ext_grid.p_mw[g])
        res['q_slack_mvar'].setdefault(g, []).append(hc.net.res_ext_grid.q_mvar[g])

    return res


def plt_diff_SLmax(dicts, keys_to_plt=None, save: bool = False, folder='../hc_test_results/figures/'):
    marker = ['o', '*', '+', '.', ',']
    # x = np.linspace(100, 10, 10)
    res = dicts[0]
    x = res['x']

    if not keys_to_plt:
        keys_to_plt = list(res['measurements'].keys())

    plt.style.use('rwth-word')
    for (key, entries) in res['measurements'].items():
        if key in keys_to_plt:
            for (b, values) in entries.items():
                plt.plot(x, values, label=str(b) + ' - ' + res['varied']['key'] + ' ' + str(res['varied']['val']),
                         marker=marker[0], markersize=10)
                for i in range(len(dicts) - 1):
                    res_act = dicts[i + 1]
                    val = res_act['measurements'][key][b]
                    plt.plot(x, val,
                             label=str(b) + ' - ' + res_act['varied']['key'] + ' ' + str(res_act['varied']['val']),
                             marker=marker[i + 1], markersize=10)
            plt.title(str(key))
            # plt.suptitle('p_max slack = ' + str(hc.model.peGmax[0].value) + " MW")
            plt.xlabel("max line loading [%]")
            plt.grid()
            if save:
                plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
                plt.savefig(folder + str(key) + '_' + res['varied']['key'] + '.png', bbox_inches="tight")
            else:
                plt.legend()
            plt.show()


def plt_pq_wind_from_hc(hc):
    plt.style.use('rwth-word')

    x = np.arange(0, 1.1, 0.1)
    fig, ax = plt.subplots()

    y = math.tan(math.acos(0.95)) * x

    ax.plot(x, y, x, -y)
    ax.fill_between(x, y, -y, alpha=0.2)

    for w in hc.model.WIND:
        ax.plot(hc.model.pG[w].value / hc.model.SWmax[w].value, hc.model.qG[w].value / hc.model.SWmax[w].value,
                label=str(w), marker='o')
    plt.xlabel('P_wind / S_wind_max')
    plt.ylabel('Q_wind / S_wind_max')
    plt.legend()
    plt.show()


def plt_pq_wind_from_res(res, save: bool = False, folder=None):
    plt.style.use('rwth-word')

    markers = ['.', '*', '+', 'o', ',']
    clrs = ['#00549F', '#000000', '#E30066', '#FFED00', '#006165',
            '#0098A1', '#57AB27', '#BDCD00', '#F6A800', '#CC071E',
            '#A11035', '#612158', '#7A6FAC']

    fig, ax = plt.subplots()

    # --- PQ region from power factor ---
    # x = np.arange(0, 1.1, 0.1)
    # y = math.tan(math.acos(0.95)) * x
    #
    # ax.plot(x, y, x, -y, color=clrs[0])
    # ax.fill_between(x, y, -y, alpha=0.2)

    # --- PQ according to VDE Variante ---
    ax.vlines(1., -0.23, 0.48, color=clrs[0])  # Pmax
    ax.vlines(0.1, -0.1, 0.1, color=clrs[0])  # Pmin
    ax.plot([0.1, 0.2, 1], [-0.1, -0.23, -0.23], color=clrs[0])  # Qmin(P)
    ax.plot([0.1, 0.2, 1], [0.1, 0.48, 0.48], color=clrs[0])  # Qmax(P)

    ax.fill_between([0.1, 0.2, 1], [-0.1, -0.23, -0.23], [0.1, 0.48, 0.48], color=clrs[0], alpha=0.2)

    for key, dict in enumerate(res):
        for bus in dict['measurements']['p_wind_mw']:
            for i in range(len(dict['measurements']['p_wind_mw'][bus])):
                p = dict['measurements']['p_wind_mw'][bus][i]
                q = dict['measurements']['q_wind_mvar'][bus][i]
                SWmax = dict['limits']['SWmax'][i]

                if dict['varied']['key'] == 'SL_%':
                    value = dict['varied']['val'][i]
                    color = clrs[i + 1]
                else:
                    value = dict['varied']['val']
                    color = clrs[key + 1]

                ax.plot(p / p, q / p, label=str(bus) + ' - ' + dict['varied']['key'] + ' ' + str(value),
                        marker='.', color=color, markersize=10)

    plt.xlabel('P_wind / P_wind_inst')
    plt.ylabel('Q_wind / P_wind_inst')
    if save:
        plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
        plt.savefig(folder + 'pq_' + dict['varied']['key'] + '.png', bbox_inches="tight")
    else:
        plt.legend()

    plt.show()


# def plt_diff_SLmax(dicts):
#     marker = [',', '+', '.', 'o', '*']
#     x = np.linspace(100, 10, 10)
#     res = dicts[0]
#     for (key, entries) in res.items():
#         for (b, values) in entries.items():
#             plt.plot(x, values, label="bus " + str(b), marker=marker[0])
#             for i in range(len(dicts) - 1):
#                 res_act = dicts[i + 1]
#                 val = res_act[key][b]
#                 plt.plot(x, val, label=str(b) + ' ' + str(i+1))
#         plt.title(str(key))
#         # plt.suptitle('p_max slack = ' + str(hc.model.peGmax[0].value) + " MW")
#         plt.xlabel("max line loading [%]")
#         plt.grid()
#         plt.legend()
#         plt.show()


def plt_results(res, x, xlabel, suptitle):
    for (meas, entries) in res.items():
        for (b, values) in entries.items():
            plt.plot(x, values, label=str(b), marker='.')
        plt.title(str(meas))
        plt.suptitle(suptitle)
        plt.xlabel(xlabel)
        plt.grid()
        plt.legend()
        plt.show()


def test_pf(net, mode, save: bool = False, folder='../hc_test_results/'):
    pp.clear_result_tables(net)
    if mode == 'ac':
        pp.runpp(net)
        pf = AC(net)
    elif mode == 'dc':
        pp.rundcpp(net)
        pf = DC(net)

    pf.solve()
    pyo_sol_to_net_res(pf.net, pf.model)

    if save:
        pp.to_excel(pf.net, folder + mode + '_pyo.xlsx')
        pp.to_excel(net, folder + mode + '.xlsx')

    return net, pf, pp.nets_equal(net, pf.net, True)


def comp_hc_powerflow_to_pp(hc):
    pyo_sol_to_net_res(hc.net, hc.model)

    net = copy.deepcopy(hc.net)
    pp.runpp(net)

    return net, hc.net, pp.nets_equal(net, hc.net, True)


def plot_wind_hc_results(nets):
    # marker_trace = pp.plotting.create_weighted_marker_trace(net, "sgen", elm_ids=net.sgen.index[net.sgen.wind_hc],
    #                                                         column_to_plot="p_mw", marker_scaling=2, color="green")
    # pp.plotting.simple_plotly(net, bus_size=0.25, additional_traces=[marker_trace])

    # get maximum wind generation from all wind generator from all nets
    max_wind_gen = max([net.sgen[net.sgen.wind_hc].p_mw.max() for net in nets])

    for net in nets:
        sgen_wind = net.sgen[net.sgen.wind_hc]
        net.res_bus.loc[sgen_wind.bus, 'p_wind_mw'] = sgen_wind['p_mw'].values

        # create marker trace with marker size scaled according to wind generation
        marker_trace = pp.plotting.create_weighted_marker_trace(net, "res_bus",
                                                                elm_ids=net.res_bus.index[net.res_bus.p_wind_mw > 0],
                                                                column_to_plot="p_wind_mw", marker_scaling=2,
                                                                trace_name='P_wind')

        # create bus trace with color map for wind generation
        bt = pp.plotting.create_bus_trace(net, buses=net.res_bus.index[net.res_bus.p_wind_mw > 0], cmap='ylorrd',
                                          colormap_column='p_wind_mw', cbar_title='Wind generation [MW]',
                                          cmax=max_wind_gen)

        # add marker size information to bus trace
        bt[0]['marker']['size'] = marker_trace['marker']['size']
        bt[0]['marker']['sizemode'] = marker_trace['marker']['sizemode']
        bt[0]['text'] = marker_trace['text']

        pp.plotting.simple_plotly(net, bus_size=3, bus_color='black', additional_traces=bt, showlegend=False)


def get_limiting_constraints(model, tolerance=1e-4):
    lim_constr = []
    violated_constr = []
    for c in model.component_objects(pe.Constraint, active=True):
        for index in c:
            if c[index].equality:
                continue
            if c[index].slack() < tolerance:
                # don't consider wind constraints when y = 0
                if ('W_max_constraint' in c[index].name) | ('W_min_constraint' in c[index].name) | (
                        'U_max_constraint' in c[index].name) | ('U_min_constraint' in c[index].name):
                    if model.y[index].value == 0.:
                        continue
                lim_constr.append(c[index])
            # negative slack means constraint is violated
            if c[index].slack() < -tolerance:
                violated_constr.append(c[index])
    return lim_constr, violated_constr


def hc_nlp_swmin(net, SWmin=10):
    hc = HC_ACOPF(net, SWmax=10000, SWmin=0, peGmax=1000000)

    hc.solve()
    initial_solve = hc.results.solver.Time

    hc.add_OPF()
    hc.fix_vars('y', 1.)

    total_solve_time = 0
    start = time.perf_counter()

    hc.solve(to_net=False)
    total_solve_time += hc.results.solver.Time

    hc.change_vals('SWmin', SWmin / hc.baseMVA)

    for w in hc.model.WIND:
        if (hc.model.pG[w].value ** 2 + hc.model.qG[w].value ** 2) < (SWmin / hc.baseMVA) ** 2:
            hc.model.y[w].fix(0.)

    hc.solve(to_net=False)
    total_solve_time += hc.results.solver.Time

    end = time.perf_counter()
    total_time = end - start

    pyo_sol_to_net_res(hc.net, hc.model)

    return hc, {'total_time': total_time, 'total_solve_time': total_solve_time, 'initial_solve': initial_solve}


def hc_nlp_swmin_steps(net, SWmin=10, stepsize=1):
    hc = HC_ACOPF(net, SWmax=10000, SWmin=0, peGmax=1000000)

    hc.solve()
    initial_solve = hc.results.solver.Time

    hc.add_OPF()
    hc.fix_vars('y', 1.)

    total_solve_time = 0
    start = time.perf_counter()

    hc.solve(to_net=False)
    total_solve_time += hc.results.solver.Time

    for i in range(stepsize, SWmin + 1, stepsize):
        print(i / hc.baseMVA)
        hc.change_vals('SWmin', i / hc.baseMVA)

        for w in hc.model.WIND:
            if (hc.model.pG[w].value ** 2 + hc.model.qG[w].value ** 2) < (i / hc.baseMVA) ** 2:
                hc.model.y[w].fix(0.)

        hc.solve(to_net=False)
        total_solve_time += hc.results.solver.Time

    end = time.perf_counter()
    total_time = end - start
    pyo_sol_to_net_res(hc.net, hc.model)

    return hc, {'total_time': total_time, 'total_solve_time': total_solve_time, 'initial_solve': initial_solve}


def compare_hc_nlp_minlp_nlpstep(net):
    print('nlp')
    hc_nlp, time_nlp = hc_nlp_swmin_steps(net, 10, 10)

    print('minlp')
    hc_minlp = HC_ACOPF(net, SWmax=10000, SWmin=10, peGmax=1000000)
    hc_minlp.solve()
    initial_solve = hc_minlp.results.solver.Time
    hc_minlp.add_OPF()
    hc_minlp.solve(solver='mindtpy')
    total_solve_time = hc_minlp.results.solver[0]['Wallclock time']
    time_minlp = {'total_time': total_solve_time, 'total_solve_time': total_solve_time, 'initial_solve': initial_solve}

    print('nlp step')
    hc_nlp_step, time_nlp_step = hc_nlp_swmin_steps(net, SWmin=10, stepsize=1)
    pyo_sol_to_net_res(hc_nlp_step.net, hc_nlp_step.model)

    plot_wind_hc_results([hc_nlp.net, hc_minlp.net, hc_nlp_step.net])

    sl_loss_nlp = get_total_line_loss_mvar(hc_nlp.net)
    sl_loss_minlp = get_total_line_loss_mvar(hc_minlp.net)
    sl_loss_nlp_step = get_total_line_loss_mvar(hc_nlp_step.net)

    lim_c_nlp, viol_c_nlp = get_limiting_constraints(hc_nlp.model)
    lim_c_minlp, viol_c_minlp = get_limiting_constraints(hc_minlp.model)
    lim_c_nlp_step, viol_c_nlp_step = get_limiting_constraints(hc_nlp_step.model)

    return hc_nlp, hc_minlp, hc_nlp_step, {
        'time': {'time_nlp': time_nlp, 'time_minlp': time_minlp, 'time_nlp_step': time_nlp_step},
        'sl_loss': {'sl_loss_nlp': sl_loss_nlp, 'sl_loss_minlp': sl_loss_minlp,
                    'sl_loss_nlp_step': sl_loss_nlp_step},
        'lim_c': {'lim_c_nlp': lim_c_nlp, 'lim_c_minlp': lim_c_minlp, 'lim_c_minlp_step': lim_c_nlp_step},
        'viol_c': {'viol_c_nlp': viol_c_nlp, 'viol_c_minlp': viol_c_minlp, 'viol_c_minlp_step': viol_c_nlp_step}}


def get_total_line_loss_mvar(net):
    return np.sqrt(net.res_line.pl_mw.sum() ** 2 + net.res_line.ql_mvar.sum() ** 2)

def get_total_trafo_loss_mvar(net):
    return np.sqrt(net.res_trafo.pl_mw.sum() ** 2 + net.res_trafo.ql_mvar.sum() ** 2)


def vary_pg_pd(net):
    res_nets = []
    res_obj = []

    for i in range(10):
        net.sgen.scaling[net.sgen.wind_hc == False] = random.uniform(0.1, 1)
        hc = HC_ACOPF(net, SWmax=10000, SWmin=10, peGmax=1000000)
        hc.add_OPF()
        hc.solve(solver='mindtpy')
        pyo_sol_to_net_res(hc.net, hc.model)

        res_nets.append(hc.net)
        res_obj.append(hc.model.OBJ())

    for sgen in net.sgen.index:
        net.sgen.p_mw[sgen] = max([res.sgen.p_mw[sgen] for res in res_nets])
        net.sgen.q_mvar[sgen] = max([res.sgen.q_mvar[sgen] for res in res_nets])


def max_wind_min_loss(net):
    hc = HC_ACOPF(net, SWmax=10000, SWmin=0, peGmax=1000000)
    hc.add_OPF()
    hc.add_loss_obj()

    obj = []
    p_w = []
    p_line = []

    for i in range(11):
        hc.model.eps = 1 - i / 10
        hc.solve(solver='mindtpy')

        if pe.check_optimal_termination(hc.results):
            obj.append(hc.model.OBJ_with_loss)
            p_w.append(pe.value(sum(hc.model.pG[w] for w in hc.model.WIND)))
            p_line.append(pe.value(sum(hc.model.pLfrom[l] + hc.model.pLto[l] for l in hc.model.L)))
        else:
            obj.append(None)
            p_w.append(None)
            p_line.append(None)

    return obj, p_w, p_line


def plot_pq_gridcodes():
    plt.style.use('rwth-word')

    clrs = ['#00549F', '#000000', '#E30066', '#FFED00', '#006165',
            '#0098A1', '#57AB27', '#BDCD00', '#F6A800', '#CC071E',
            '#A11035', '#612158', '#7A6FAC']

    fig, ax = plt.subplots()

    # Variant 1
    ax.vlines(1., -0.23, 0.48, color=clrs[8])  # Pmax
    ax.vlines(0.1, -0.1, 0.1, color=clrs[8])  # Pmin
    ax.plot([0.1, 0.2, 1], [-0.1, -0.23, -0.23], color=clrs[8], label='Variante 1')  # Qmin(P)
    ax.plot([0.1, 0.2, 1], [0.1, 0.48, 0.48], color=clrs[8])  # Qmax(P)

    # Variant 2
    ax.vlines(1., -0.33, 0.41, color=clrs[7])  # Pmax
    ax.vlines(0.1, -0.1, 0.1, color=clrs[7])  # Pmin
    ax.plot([0.1, 0.2, 1], [-0.1, -0.33, -0.33], color=clrs[7], label='Variante 2')  # Qmin(P)
    ax.plot([0.1, 0.2, 1], [0.1, 0.41, 0.41], color=clrs[7])  # Qmax(P)

    # Variant 3
    ax.vlines(1., -0.41, 0.33, color=clrs[0])  # Pmax
    ax.vlines(0.1, -0.1, 0.1, color=clrs[0])  # Pmin
    ax.plot([0.1, 0.2, 1], [-0.1, -0.41, -0.41], color=clrs[0], label='Variante 3')  # Qmin(P)
    ax.plot([0.1, 0.2, 1], [0.1, 0.33, 0.33], color=clrs[0])  # Qmax(P)

    plt.xlabel('$P_{mom} / P_{inst}$')
    plt.ylabel('$Q / P_{inst}$')
    plt.yticks(np.arange(-0.5, 0.6, 0.1))
    plt.ylim(-0.5, 0.5)
    plt.legend(loc='center left', bbox_to_anchor=(0.6, 0.5))
    plt.show()


def plot_qu_gridcodes():
    plt.style.use('rwth-word')

    clrs = ['#00549F', '#000000', '#E30066', '#FFED00', '#006165',
            '#0098A1', '#57AB27', '#BDCD00', '#F6A800', '#CC071E',
            '#A11035', '#612158', '#7A6FAC']

    fig, ax = plt.subplots()

    # Variant 1
    ax.plot([96, 120, 127], [0.48, 0.48, -0.23], color=clrs[8], label='Variante 1')  # Qmin(P)
    ax.plot([96, 103, 127], [0.48, -0.23, -0.23], color=clrs[8])  # Qmax(P)

    # Variant 2
    ax.plot([96, 120, 127], [0.41, 0.41, -0.33], color=clrs[7], label='Variante 2')  # Qmin(P)
    ax.plot([96, 103, 127], [0.41, -0.33, -0.33], color=clrs[7])  # Qmax(P)

    # Variant 3
    ax.plot([96, 120, 127], [0.33, 0.33, -0.41], color=clrs[0], label='Variante 3')  # Qmin(P)
    ax.plot([96, 103, 127], [0.33, -0.41, -0.41], color=clrs[0])  # Qmax(P)

    plt.xlabel('$U_{b} [kV]$')
    plt.ylabel('$Q / P_{inst}$')
    plt.yticks(np.arange(-0.5, 0.6, 0.1))
    plt.ylim(-0.5, 0.5)
    plt.legend(loc='center')
    plt.show()


def plot_qu_res(nets, labels=None):
    # values for variant 1
    m_qu = (0.48 + 0.23) / (96 - 103) * 110
    qu_min = -m_qu * 96 / 110 + 0.48
    qu_max = -m_qu * 120 / 110 + 0.48
    ymax = m_qu * 1.1 + qu_max
    ymin = m_qu * 0.9 + qu_min

    plt.style.use('rwth-word')
    clrs = ['#00549F', '#E30066', '#BDCD00','#000000', '#FFED00', '#006165',
            '#0098A1', '#57AB27',  '#F6A800', '#CC071E',
            '#A11035', '#612158', '#7A6FAC']
    marker = ['o', '*', '+', '.', ',']

    fig, ax = plt.subplots()
    ax.fill_between([0.9, 103 / 110, 120 / 110, 1.1], [0.48, 0.48, 0.48, ymax], [ymin, -0.23, -0.23, -0.23],
                    color=clrs[8], alpha=0.2)

    ax.plot([0.9, 120 / 110, 1.1, 1.1], [0.48, 0.48, ymax, -0.23], color=clrs[8])  # Qmax(P)
    ax.plot([0.9, 0.9, 103 / 110, 1.1], [0.48, ymin, -0.23, -0.23], color=clrs[8])  # Qmin(P)

    plt.xticks(np.arange(0.9, 1.1, 0.05))
    plt.xlabel('$U_{b} [p.u.]$')
    plt.ylabel('$Q_{w} / P_{w}$')

    for i, net in enumerate(nets):
        line = ax.scatter(net.res_bus.vm_pu[net.sgen.bus[net.sgen.wind_hc]],
                          net.res_sgen.q_mvar[net.sgen.wind_hc].values / net.res_sgen.p_mw[net.sgen.wind_hc].values,
                          color=clrs[i], marker=marker[i])
        if labels:
            line.set_label(labels[i])

    if labels:
        plt.legend()


if __name__ == '__main__':
    net = create_testnet()
    net = pp.networks.simple_four_bus_system()

    net = sb.get_simbench_net("1-HV-mixed--0-no_sw")
    # hc_nlp, hc_minlp, hc_nlp_step, results = compare_hc_nlp_minlp_nlpstep(net)

    # start_0 = time.perf_counter()
    # times_nlp = []
    # times_minlp = []
    # times_nlp_step = []
    #
    # line_loss_nlp = []
    # line_loss_minlp = []
    # line_loss_nlp_step = []
    #
    # obj_nlp = []
    # obj_minlp = []
    # obj_nlp_step = []
    #
    # for i in range(20):
    #     hc_nlp, time_nlp = hc_nlp_swmin(net, SWmin=10)
    #     times_nlp.append(time_nlp)
    #     line_loss_nlp.append(get_total_line_loss_mvar(hc_nlp.net))
    #     obj_nlp.append(pe.value(hc_nlp.model.OBJ))
    #
    #     hc_minlp = HC_ACOPF(net, SWmax=10000, SWmin=10, peGmax=1000000)
    #     hc_minlp.solve()
    #     hc_minlp.add_OPF()
    #     hc_minlp.solve(solver='mindtpy')
    #     times_minlp.append(hc_minlp.results.solver[0]['Wallclock time'])
    #     line_loss_minlp.append(get_total_line_loss_mvar(hc_minlp.net))
    #     obj_minlp.append(pe.value(hc_minlp.model.OBJ))
    #
    #     hc_nlp_step, time_nlp_step = hc_nlp_swmin_steps(net, SWmin=10, stepsize=1)
    #     times_nlp_step.append(time_nlp_step)
    #     line_loss_nlp_step.append(get_total_line_loss_mvar(hc_nlp_step.net))
    #     obj_nlp_step.append(pe.value(hc_nlp_step.model.OBJ))
    # end_0 = time.perf_counter()
    #
    # get all entries of obj_nlp > 1136.815
    # obj_nlp = np.array(obj_nlp)
    # obj_nlp_step = np.array(obj_nlp_step)
    # obj_minlp = np.array(obj_minlp)
    #
    # # create array with all total_solve_times from times_nlp
    # total_solve_time_nlp = np.array([time['total_solve_time'] for time in times_nlp])