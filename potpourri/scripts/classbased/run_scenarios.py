import copy
import datetime
import os

import matplotlib.pyplot as plt
import pandas as pd
from pandapower.timeseries import OutputWriter

from potpourri.models.class_based.HC_ACOPF import HC_ACOPF
from potpourri.scripts.classbased.init_pyo_from_pp_res import init_pyo_from_dcpp

import pandapower as pp
import pickle
import pyomo.environ as pe
import numpy as np
from tqdm import tqdm
import random

from potpourri.models.class_based.pyo_to_net import pyo_sol_to_net_res
from potpourri.scripts.classbased.plot_functions import plot_wind_hc_results


def add_scenario_load_sgen_to_net(net, busses):
    net.load.in_service = False
    net.sgen.in_service = False

    load_ind = pp.create_loads(net, busses, p_mw=0, name='load_scenario')
    sgen_ind = pp.create_sgens(net, busses, p_mw=0, name='sgen_scenario')

    return load_ind, sgen_ind


def create_output_writer(net, time_steps, output_dir):
    ow = OutputWriter(net, time_steps, output_path=output_dir, output_file_type=".xlsx", log_variables=list())

    ow.log_variable('res_bus', 'vm_pu')
    ow.log_variable('res_sgen', 'q_mvar')
    ow.log_variable('res_line', 'loading_percent')
    ow.log_variable('res_line', 'pl_mw')
    ow.log_variable('res_line', 'ql_mvar')
    ow.log_variable('res_trafo', 'loading_percent')
    ow.log_variable('res_ext_grid', 'p_mw')
    ow.log_variable('res_ext_grid', 'q_mvar')

    wind_index = net.sgen.index
    ow.log_variable('res_sgen', 'p_mw')
    ow.log_variable('res_sgen', ['p_mw', 'y_wind'], eval_function=get_p_wind(ind=len(wind_index)), eval_name='wind_p')
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


        # w = hc.model.WIND.data()
        # rand_w = random.sample(w, round(len(w)/4))
        # print(rand_w)
        # for w in rand_w:
        #     # hc.model.SWmin[w] = 20
        #     # hc.model.SWmax[w] = 150
        #     hc.model.y[w].fix(0.)

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


def scenario_supply(scenarios):
    n_scenarios = len(scenarios['load'][0][0])
    n_cases = len(scenarios['load'])

    scenarios['supply'] = [ [] for _ in range(n_cases) ]  # Initialize a list for each case
    scenarios['supply_sum'] = [ [] for _ in range(n_cases) ]  # Initialize a list for each case

    for case in range(n_cases):
        p_supply_case = []  # Initialize a list to store p_supply values for the current case
        for sc in range(n_scenarios):
            p_load = scenarios['load'][case][:, sc]
            p_sgen = scenarios['sgen'][case][:, sc]

            p_supply = p_sgen - p_load
            p_supply_case.append(p_supply)  # Append p_supply to the list for the current case

        scenarios['supply'][case] = np.array(p_supply_case)
        scenarios['supply_sum'][case] = np.sum(p_supply_case, axis=1)  # Sum over columns and store in 'supply_sum'

def plot_scenario_net(net, results_dir, n_cases):
    for case in range(n_cases):
        p_mw_file = os.path.join(results_dir, 'case_' + str(case), 'res_sgen', 'p_mw.xlsx')
        p_mw = pd.read_excel(p_mw_file, index_col=0)
        q_mvar_file = os.path.join(results_dir, 'case_' + str(case), 'res_sgen', 'q_mvar.xlsx')
        q_mvar = pd.read_excel(q_mvar_file, index_col=0)

        vm_pu_file = os.path.join(results_dir, 'case_' + str(case), 'res_bus', 'vm_pu.xlsx')
        vm_pu = pd.read_excel(vm_pu_file, index_col=0)

        ext_grid_p_mw_file = os.path.join(results_dir, 'case_' + str(case), 'res_ext_grid', 'p_mw.xlsx')
        ext_grid_p_mw = pd.read_excel(ext_grid_p_mw_file, index_col=0)
        ext_grid_q_mvar_file = os.path.join(results_dir, 'case_' + str(case), 'res_ext_grid', 'q_mvar.xlsx')
        ext_grid_q_mvar = pd.read_excel(ext_grid_q_mvar_file, index_col=0)

        line_loading_percent_file = os.path.join(results_dir, 'case_' + str(case), 'res_line', 'loading_percent.xlsx')
        line_loading_percent = pd.read_excel(line_loading_percent_file, index_col=0)
        line_pl_mw_file = os.path.join(results_dir, 'case_' + str(case), 'res_line', 'pl_mw.xlsx')
        line_pl_mw = pd.read_excel(line_pl_mw_file, index_col=0)
        line_ql_mvar_file = os.path.join(results_dir, 'case_' + str(case), 'res_line', 'ql_mvar.xlsx')
        line_ql_mvar = pd.read_excel(line_ql_mvar_file, index_col=0)

        trafo_loading_percent_file = os.path.join(results_dir, 'case_' + str(case), 'res_trafo', 'loading_percent.xlsx')
        trafo_loading_percent = pd.read_excel(trafo_loading_percent_file, index_col=0)

        for sc in range(len(p_mw)):
            net.res_sgen.p_mw = p_mw.loc[sc]
            net.res_sgen.q_mvar = q_mvar.loc[sc]

            net.res_bus.vm_pu = vm_pu.loc[sc]

            net.res_ext_grid.p_mw = ext_grid_p_mw.loc[sc]
            net.res_ext_grid.q_mvar = ext_grid_q_mvar.loc[sc]

            net.res_line.loading_percent = line_loading_percent.loc[sc]
            net.res_line.pl_mw = line_pl_mw.loc[sc]
            net.res_line.ql_mvar = line_ql_mvar.loc[sc]

            net.res_trafo.loading_percent = trafo_loading_percent.loc[sc]

            try:
                plot_wind_hc_results([net])
            except:
                print('Error plotting results of case ' + str(case) + ' scenario ' + str(sc))



def box_plot_scenarios(results_dir, n_cases, cases=None):
    p_wind = []

    # for case in range(n_cases):
    #     p_mw_file = os.path.join(results_dir, 'case_' + str(case), 'res_sgen', 'p_mw.xlsx')
    #     p_mw = pd.read_excel(p_mw_file, index_col=0)
    #     p_w = p_mw['wind_p']
    #     p_wind.append(p_w)

    for case in cases:
        p_mw_y_file = os.path.join(results_dir, 'case_' + str(case), 'res_sgen', '[\'p_mw\', \'y_wind\'].xlsx')
        p_mw_y = pd.read_excel(p_mw_y_file, index_col=0)
        p_w = p_mw_y.sum(axis=1)
        p_wind.append(p_w)


    fig, ax = plt.subplots()
    ax.boxplot(p_wind[0:3])

    ax.set_title("Windintegrationspotenzial")
    ax.set_ylabel('$P_{Wind}$ [MW]')
    ax.set_xlabel('Netznutzungsfall')
    ax.set_xticks([1, 2, 3], ['hL', 'hW', 'hPV'])

    fig, ax = plt.subplots()
    ax.boxplot(p_wind[3:5])

    ax.set_title("Windintegrationspotenzial")
    ax.set_ylabel('$P_{Wind}$ [MW]')
    ax.set_xlabel('Netznutzungsfall')
    ax.set_xticks([1, 2], ['lW', 'lPV'])

    return p_wind


def std_dev_scenarios(dir):
    with open(dir, "rb") as f:
        scenarios = pickle.load(f)

    scenarios['std_dev_load'] = {}
    scenarios['std_dev_sgen'] = {}
    scenarios['std_dev_supply'] = {}

    for case in range(len(scenarios['load'])):
        std_dev_load = np.std(scenarios['load'][case], axis=1)
        std_dev_sgen = np.std(scenarios['sgen'][case], axis=1)
        std_dev_supply = np.std(scenarios['supply'][case], axis=0)
        scenarios['std_dev_load'][case] = std_dev_load
        scenarios['std_dev_sgen'][case] = std_dev_sgen
        scenarios['std_dev_supply'][case] = std_dev_supply

    fig, ax = plt.subplots()
    ax.bar(range(len(scenarios['std_dev_sgen'][3])), scenarios['std_dev_sgen'][3])
    ax.set_xlabel("Knoten")
    ax.set_ylabel("Standardabweichung [MW]")

    return scenarios


def std_dev_p_wind_scenarios(dir, case):
    p_mw_file = os.path.join(results_dir, 'case_' + str(case), 'res_sgen', 'p_mw.xlsx')
    p_mw = pd.read_excel(p_mw_file, index_col=0)


def create_weighted_p_wind(probabilities, p_wind):
    min = np.min(probabilities)
    weights = probabilities / min * 3
    weights = weights.round()

    p_wind_weighted = []
    for p_w in p_wind:
        p_w_weighted = np.array([])
        for i, w in enumerate(weights):
            p_w_weighted = np.concatenate((p_w_weighted, np.repeat(p_w[i], w.astype(int))))

        p_wind_weighted.append(p_w_weighted)

    fig, ax = plt.subplots()
    ax.boxplot(p_wind_weighted)


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

    with open('C:\\Users\\f.lohse\PycharmProjects\potpourri\potpourri\data\simbench_hv_grid_with_potential_pkl.pkl',
              'rb') as f:
        net = pickle.load(f)
        net.bus.windpot_p_mw = net.bus.windpot_p_mw.where(net.bus.windpot_p_mw <= 200, 200)

    with open("../../data/scenarios/" + grid + "_scenario_types.pkl", "rb") as f:
        net = pickle.load(f)

    with open("../../data/scenarios/" + grid + "_loadcases_" + type + "_scenarios_with_q.pkl", "rb") as f:
        scenarios = pickle.load(f)

    results_dir = '../../results/scenarios/' + grid + '_loadcases_' + type + '_scenarios/'

    obj = run_scenarios(net, scenarios, results_dir + 'case', peGmax=10000, cases=[3], SWmin=10)

    # for i in range(len(scenarios['load'][0][0])):
    #     net.load.p_mw[net.load.name == 'load_scenario'] = scenarios['load'][0][:, i]
    #     net.load.q_mvar[net.load.name == 'load_scenario'] = scenarios['load_q'][0][:, i]
    #
    #     net.sgen.p_mw[net.sgen.name == 'sgen_scenario'] = scenarios['sgen'][0][:, i]
    #     net.sgen.q_mvar[net.sgen.name == 'sgen_scenario'] = scenarios['sgen_q'][0][:, i]
