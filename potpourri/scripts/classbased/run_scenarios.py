from potpourri.models.class_based.HC_ACOPF import HC_ACOPF
from potpourri.scripts.classbased.init_pyo_from_pp_res import init_pyo_from_dcpp

import pandapower as pp
import pickle
import pyomo.environ as pe
from tqdm import tqdm
from pandapower.timeseries import OutputWriter
import numpy as np


def add_scenario_load_sgen_to_net(net, busses):
    net.load.in_service = False
    net.sgen.in_service = False

    load_ind = pp.create_loads(net, busses, p_mw=0, name='load_scenario')
    sgen_ind = pp.create_sgens(net, busses, p_mw=0, name='sgen_scenario')

    return load_ind, sgen_ind


def create_output_writer(net, time_steps, output_dir):
    ow = OutputWriter(net, time_steps, output_path=output_dir, output_file_type=".xlsx", log_variables=list())

    ow.log_variable('res_bus', 'vm_pu')
    ow.log_variable('res_bus', 'va_degree')
    ow.log_variable('res_sgen', 'q_mvar')
    ow.log_variable('res_line', 'loading_percent')
    ow.log_variable('res_line', 'pl_mw')
    ow.log_variable('res_line', 'pl_mw', eval_function=np.sum, eval_name='sum_p_line')
    ow.log_variable('res_line', 'ql_mvar')
    ow.log_variable('res_trafo', 'loading_percent')
    ow.log_variable('res_trafo', 'pl_mw')
    ow.log_variable('res_trafo', 'ql_mvar')
    ow.log_variable('res_ext_grid', 'p_mw')
    ow.log_variable('res_ext_grid', 'q_mvar')

    wind_index = net.sgen.index[net.sgen.wind_hc]
    ow.log_variable('res_sgen', 'p_mw', index=wind_index, eval_function=np.sum, eval_name='sum_p_wind')
    ow.log_variable('res_sgen', ['p_mw', 'y_wind'], eval_function=get_p_wind(ind=len(net.sgen.index)),
                    eval_name='wind_p')
    ow.log_variable('res_sgen', 'q_mvar')

    return ow


def get_p_wind(ind):
    def eval_func(results, n_columns=ind):
        p_wind = results[:, 0] * results[:, 1]
        return p_wind

    return eval_func


def run_scenarios(net, scenarios, output_dir, SWmin=10, peGmax=10000, cases=None):
    busses = scenarios['load_sgen_busses']
    if cases:
        pass
    else:
        cases = range(len(scenarios['load']))
    n_scenarios = len(scenarios['load'][0][0])

    load_ind, sgen_ind = add_scenario_load_sgen_to_net(net, busses)

    obj = {}

    for case in cases:
        obj[case] = []

        hc = HC_ACOPF(net, SWmin=SWmin, peGmax=peGmax)

        sc = 0
        for i, d in enumerate(load_ind):
            hc.model.pD[d] = scenarios['load'][case][i, sc]
            hc.model.qD[d] = scenarios['load_q'][case][i, sc]

        for i, g in enumerate(sgen_ind):
            hc.model.pG[g] = scenarios['sgen'][case][i, sc]
            hc.model.qG[g] = scenarios['sgen_q'][case][i, sc]

        if any(hc.net.trafo.shift_degree):
            init_pyo_from_dcpp(hc.net, hc.model)

        hc.solve()
        hc.add_OPF()
        # hc.add_loss_obj()
        # hc.model.eps = 0.9

        hc.model.OBJ.deactivate()

        def obj_wind_loss(model):
            return sum(model.pG[w] for w in model.WIND) - sum(
                model.pLfrom[l] + model.pLto[l] for l in model.L)

        hc.model.OBJ_wind_loss = pe.Objective(rule=obj_wind_loss, sense=pe.maximize)

        output_path = output_dir + '_' + str(case)
        ow = create_output_writer(hc.net, range(n_scenarios), output_path)
        ow.init_all(hc.net)

        hc.solve(solver='mindtpy', mip_solver='gurobi')

        # use solver.status == 'ok' instead?
        pf_converged = hc.results.solver.termination_condition == (
                pe.TerminationCondition.optimal or pe.TerminationCondition.feasible)

        ow.time_step = sc
        ow.save_results(hc.net, sc, pf_converged=pf_converged, ctrl_converged=pf_converged)

        obj[case].append(pe.value(hc.model.OBJ))

        for sc in tqdm(range(1, n_scenarios)):
            for i, d in enumerate(load_ind):
                hc.model.pD[d] = scenarios['load'][case][i, sc]
                hc.model.qD[d] = scenarios['load_q'][case][i, sc]

            for i, g in enumerate(sgen_ind):
                hc.model.pG[g] = scenarios['sgen'][case][i, sc]
                hc.model.qG[g] = scenarios['sgen_q'][case][i, sc]

            hc.solve(solver='mindtpy', mip_solver='gurobi')

            # use solver.status == 'ok' instead?
            pf_converged = hc.results.solver.termination_condition == (
                    pe.TerminationCondition.optimal or pe.TerminationCondition.feasible)

            ow.time_step = sc
            ow.save_results(hc.net, sc, pf_converged=pf_converged,
                            ctrl_converged=pe.check_optimal_termination(hc.results))

            obj[case].append(pe.value(hc.model.OBJ))

    return obj


if __name__ == '__main__':
    # with open("../../data/scenarios/1-HV-mixed--0-no_sw_scenario_types.pkl", "rb") as f:
    #     net = pickle.load(f)
    #
    # with open("../../data/scenarios/1-HV-mixed_loadcases_assumption_scenarios_with_q.pkl", "rb") as f:
    #     scenarios = pickle.load(f)
    #
    # results_dir = '../../results/scenarios/1-HV-mixed_loadcases_distribution_scenarios/'

    grid = '1-HV-mixed--0-no_sw'
    # grid = 'hv_grid'
    type = "assumption"

    # net_name = 'simbench_hv_grid_with_potential_pkl.pkl'
    # with open('C:\\Users\\f.lohse\PycharmProjects\potpourri\potpourri\data\\' + net_name,
    #           'rb') as f:
    #     net = pickle.load(f)
        # net.bus.windpot_p_mw = net.bus.windpot_p_mw.where(net.bus.windpot_p_mw <= 200, 200)

    net_name = grid + "_scenario_types.pkl"
    with open("../../data/scenarios/" + grid + "_scenario_types.pkl", "rb") as f:
        net = pickle.load(f)

    with open("../../data/scenarios/" + grid + "_loadcases_" + type + "_scenarios_with_q.pkl", "rb") as f:
        scenarios = pickle.load(f)

    # results_dir = '../../results/test_scenarios/' + grid + '_loadcases_' + type + '_scenarios/'

    case = 'lW'
    factors = net.loadcases.loc[case]
    net.load.p_mw *= factors['pload']
    net.load.q_mvar *= factors['qload']
    net.sgen.scaling[net.sgen.type == 'Wind'] = factors['Wind_p']
    net.sgen.scaling[net.sgen.type == 'PV'] = factors['PV_p']
    net.sgen.scaling[(net.sgen.type != 'Wind') & (net.sgen.type != 'Solar')] = factors['RES_p']
    net.ext_grid.vm_pu = factors['Slack_vm']

    results_dir = 'C:\\Users\\f.lohse\\PycharmProjects\\potpourri\\potpourri\\results\\scenarios_multiobj\\' + net_name + '_' + type + '\\'

    obj = run_scenarios(net, scenarios, results_dir + 'case', peGmax=10000, cases=[3, 4, 0, 1, 2], SWmin=10)

    # for i in range(len(scenarios['load'][case][0])):
    #     net.load.p_mw[net.load.name == 'load_scenario'] = scenarios['load'][case][:, i]
    #     net.load.q_mvar[net.load.name == 'load_scenario'] = scenarios['load_q'][case][:, i]
    #
    #     net.sgen.p_mw[net.sgen.name == 'sgen_scenario'] = scenarios['sgen'][case][:, i]
    #     net.sgen.q_mvar[net.sgen.name == 'sgen_scenario'] = scenarios['sgen_q'][case][:, i]
