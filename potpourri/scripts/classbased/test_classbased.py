import copy
import math
import time

import pyomo.environ as pe
import simbench as sb
import random
from potpourri.models.class_based.HC_ACOPF import HC_ACOPF
from potpourri.models.class_based.DC import DC
from potpourri.models.class_based.AC import AC

from potpourri.models.class_based.pyo_to_net import pyo_sol_to_net_res

from potpourri.scripts.classbased.plot_functions import *


def create_testnet():
    net = pp.create_empty_network()

    pp.create_buses(net, 3, name=["b1", "b2", "b3"], vn_kv=110, max_vm_pu=1.118, min_vm_pu=0.9)
    pp.create_lines(net, [0, 0, 1], [1, 2, 2], name=["l1", "l2", "l3"], length_km=[3, 3, 3],
                    std_type='305-AL1/39-ST1A 110.0')
    pp.create_sgens(net, [0, 1, 2], p_mw=[10, 10, 10])
    pp.create_ext_grid(net, 0)
    pp.create_loads(net, [0, 1, 2], name=["load1", "load2", "load3"], p_mw=[1, 17, 10])
    # pp.createsgen(net, 0, p_mw=1, controllable=True)
    # pp.create_sgen(net, 0, p_mw=1, controllable=True, max_p_mw=10)
    # pp.create_load(net, 1, p_mw=1, controllable=True)
    # pp.create_load(net, 2, p_mw=2, controllable=True, max_p_mw=5)

    # pp.create_sgens(net, net.bus.index, p_mw=0, wind_hc=True)

    # pp.create_shunt(net, 2, p_mw=3, q_mvar=5)
    return net


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

    return net, pp.nets_equal(net, hc.net, True)


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
                    if model.pG[index].value + model.qG[index].value <= tolerance:
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


def vary_load(net):
    hcs = []
    res_obj = []

    for i in range(10):
        net.load.scaling = np.random.rand(len(net.load))
        hc = HC_ACOPF(net, SWmax=10000, SWmin=10, peGmax=1000000)
        hc.add_OPF()
        hc.solve(solver='mindtpy', mip_solver='gurobi')

        hcs.append(copy.deepcopy(hc))
        res_obj.append(hc.model.OBJ())


def diff_trafos(net):
    std_types_220 = pp.find_std_type_by_parameter(net, data={"vn_lv_kv": 110., "vn_hv_kv": 220.}, element='trafo')
    std_types_380 = pp.find_std_type_by_parameter(net, data={"vn_lv_kv": 110., "vn_hv_kv": 380.}, element='trafo')

    for i in len(std_types_220):
        net.trafo.std_type[net.trafo.vn_hv_kv == 220.] = std_types_220[i]
        net.trafo.std_type[net.trafo.vn_hv_kv == 380.] = std_types_380[i]


if __name__ == '__main__':
    # net = create_testnet()
    # net = pp.networks.simple_four_bus_system()
    #
    net = sb.get_simbench_net("1-HV-mixed--0-no_sw")

    # obj = []
    # hcs = []
    #
    # for i in range(10):
    #     net.trafo.parallel = i+1
    #     hc = HC_ACOPF(net, SWmin=10)
    #     hc.solve()
    #     hc.add_OPF()
    #     hc.solve(solver='mindtpy', mip_solver='gurobi')
    #     obj.append(pe.value(hc.model.OBJ))
    #     hcs.append(copy.deepcopy(hc))

    # plot hc.model.pG[w] where hc.model.y[w] == 1. for all w in hc.model.WINd and all hcs in a bar plot



    # net = pp.networks.simple_mv_open_ring_net()
    # net.switch.closed = True
    # net.trafo.shift_degree = 0.
    # pp.create_bus(net, 20.)
    # pp.create_switch(net, 2, 3, 'b')

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
