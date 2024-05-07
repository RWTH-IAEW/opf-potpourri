import copy

from potpourri.models.HC_ACOPF import HC_ACOPF
from potpourri.models.pyo_to_net import pyo_sol_to_net_res
from potpourri.models.init_pyo_from_pp_res import init_pyo_from_dcpp
from scripts.plot_functions import set_plt_config

import pandapower as pp
import pickle
import pyomo.environ as pe
from tqdm import tqdm
from pandapower.timeseries import OutputWriter
import numpy as np
import matplotlib.pyplot as plt

from rwth_colors import colors

config = set_plt_config()


def add_scenario_load_sgen_to_net(net, busses):
    net.load.in_service = False
    net.sgen.in_service = False

    load_ind = pp.create_loads(net, busses, p_mw=0, name='load_scenario')
    sgen_ind = pp.create_sgens(net, busses, p_mw=0, name='sgen_scenario')

    return load_ind, sgen_ind


def write_scenario_to_model(model, load_ind, sgen_ind, scenarios, case, sc):
    for i, d in enumerate(load_ind):
        model.pD[d] = scenarios['load'][case][i, sc]
        model.qD[d] = scenarios['load_q'][case][i, sc]

    for i, g in enumerate(sgen_ind):
        model.psG[g] = scenarios['sgen'][case][i, sc]
        model.qsG[g] = scenarios['sgen_q'][case][i, sc]


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


def run_scenarios(net, scenarios, output_dir, SWmin=10, cases=None, solver='gurobi', tap='fixed'):
    busses = scenarios['load_sgen_busses']
    if cases:
        pass
    else:
        cases = range(len(scenarios['load']))
    n_scenarios = len(scenarios['load'][0][0])

    load_ind, sgen_ind = add_scenario_load_sgen_to_net(net, busses)

    hc_results = {}

    for case in cases:
        hc_results[case] = []

        hc = HC_ACOPF(net)

        sc = 0
        write_scenario_to_model(hc.model, load_ind, sgen_ind, scenarios, case, sc)

        if any(hc.net.trafo.shift_degree):
            init_pyo_from_dcpp(hc.net, hc.model)
        hc.solve()

        hc.add_OPF(SWmin=SWmin)
        if tap == 'discrete':
            hc.add_tap_changer_discrete()
        elif tap == 'linear':
            hc.add_tap_changer_linear()
        for w in hc.model.WIND_HC:
            if hc.model.pWmax[w].value == 0:
                hc.model.y[w].fix(0)
                hc.model.psG[w].fix(0.)
                hc.model.qsG[w].fix(0.)

        output_path = output_dir + 'case_' + str(case)
        # ow = create_output_writer(hc.net, range(n_scenarios), output_path)
        # ow.init_all(hc.net)

        hc.solve(solver='mindtpy', mip_solver=solver)
        if not pe.check_optimal_termination(hc.results):
            print('No optimal solution found for case ' + str(case) + ' and scenario ' + str(sc))
            pyo_sol_to_net_res(hc.net, hc.model)

        pickle.dump(hc.net, open(output_path + '\\sc_' + str(sc) + '.pkl', 'wb'))

        # use solver.status == 'ok' instead?
        pf_converged = hc.results.solver.termination_condition == (
                pe.TerminationCondition.optimal or pe.TerminationCondition.feasible)

        hc.net.res_sgen.y_wind.fillna(0, inplace=True)
        # ow.time_step = sc
        # ow.save_results(hc.net, sc, pf_converged=True, ctrl_converged=True)

        hc_results[case].append(pe.value(sum(hc.model.psG[w] for w in hc.model.WIND_HC)))

        for sc in tqdm(range(1, n_scenarios)):
            write_scenario_to_model(hc.model, load_ind, sgen_ind, scenarios, case, sc)

            hc.solve(solver='mindtpy', mip_solver=solver)
            if not pe.check_optimal_termination(hc.results):
                print('No optimal solution found for case ' + str(case) + ' and scenario ' + str(sc))
                pyo_sol_to_net_res(hc.net, hc.model)

            # use solver.status == 'ok' instead?
            pf_converged = hc.results.solver.termination_condition == (
                    pe.TerminationCondition.optimal or pe.TerminationCondition.feasible)

            hc.net.res_sgen.y_wind.fillna(0, inplace=True)
            # ow.time_step = sc
            # ow.save_results(hc.net, sc, pf_converged=True,
            #                 ctrl_converged=True)

            pickle.dump(hc.net, open(output_path + '\\sc_' + str(sc) + '.pkl', 'wb'))

            hc_results[case].append(pe.value(sum(hc.model.psG[w] for w in hc.model.WIND_HC)))

    return hc_results


def load_results(results_dir, scenarios):
    nets = []
    for sc in range(len(scenarios['load'][0][0])):
        with open(results_dir + '\\sc_' + str(sc) + '.pkl', "rb") as f:
            nets.append(pickle.load(f))

    return nets


def plot_results(nets, dir):
    buses = [net.sgen.bus[net.res_sgen.y_wind == 1] for net in nets]
    unique_buses = sorted(list(set(bus for sublist in buses for bus in sublist)))

    # Load all values of net.res_sgen.p_mw where net.res_sgen.y_wind == 1 for all net in nets
    p_mw_values = [net.res_sgen.p_mw[net.sgen.bus.isin(unique_buses) & net.sgen.wind_hc == True] for net in nets]

    # Transpose the list of lists
    p_wind_buses = list(map(list, zip(*p_mw_values)))

    # Create a new figure
    fig, ax = plt.subplots()

    # Create a boxplot for each bus
    ax.boxplot(p_wind_buses)
    ax.grid()
    ax.set_xticks(range(1, len(unique_buses) + 1), unique_buses)
    ax.set_ylabel('P [MW]')
    ax.set_xlabel('Knoten')
    fig.set_size_inches((config['textbreite'] * 0.8, 0.5 * config['textbreite']))
    fig.savefig(dir + '\\boxplot.pdf', format='pdf', bbox_inches='tight')

    # Get the max and min value for each index
    max_values_per_bus = [max(values) for values in p_wind_buses]
    min_values_per_bus = [min(values) for values in p_wind_buses]

    net_min_max = copy.deepcopy(nets[0])
    net_min_max.bus['p_wind_max'] = 0.
    net_min_max.bus['p_wind_min'] = 0.
    net_min_max.bus.loc[unique_buses, 'p_wind_max'] = max_values_per_bus
    net_min_max.bus.loc[unique_buses, 'p_wind_min'] = min_values_per_bus



    # net_max = copy.deepcopy(nets[0])
    # net_max.res_sgen.p_mw[net_max.sgen.bus.isin(unique_buses) & net_max.sgen.wind_hc == True] = max_values_per_bus
    # net_min = copy.deepcopy(nets[0])
    # net_min.res_sgen.p_mw[net_min.sgen.bus.isin(unique_buses) & net_min.sgen.wind_hc == True] = min_values_per_bus



def plot_min_max_hc_bus(net, dir=None):
    # create marker trace with marker size scaled according to wind generation
    p_wind_max_trace = pp.plotting.create_weighted_marker_trace(net, "bus",
                                                            elm_ids=net.bus.index[net.bus.p_wind_max > 0],
                                                            column_to_plot="p_wind_max",
                                                            trace_name='P_wind_max', color=colors['blue'])

    p_wind_min_trace = pp.plotting.create_weighted_marker_trace(net, "bus",
                                                            elm_ids=net.bus.index[net.bus.p_wind_min > 0],
                                                            column_to_plot="p_wind_min",
                                                            trace_name='P_wind_min', color=colors['red'])

    fig = pp.plotting.simple_plotly(net, bus_size=3, bus_color='black', additional_traces=[p_wind_max_trace, p_wind_min_trace], showlegend=False,
                                    auto_open=False, ext_grid_color=colors['green'])
    fig.update_layout(
        {
            "paper_bgcolor": "rgba(0, 0, 0, 0)",
            "plot_bgcolor": "rgba(0, 0, 0, 0)",
        }
    )
    fig.show(renderer="browser")

    x = 600
    fig.update_layout(width=x, height=x * 3 / 5)

    if dir:
        fig.write_image(dir + 'min_max_hc_bus.pdf')

    return fig


if __name__ == '__main__':
    # with open("../../data/scenarios/1-HV-mixed--0-no_sw_scenario_types.pkl", "rb") as f:
    #     net = pickle.load(f)
    #
    # with open("../../data/scenarios/1-HV-mixed_loadcases_assumption_scenarios_with_q.pkl", "rb") as f:
    #     scenarios = pickle.load(f)
    #
    # results_dir = '../../results/scenarios/1-HV-mixed_loadcases_distribution_scenarios/'

    # grid = '1-HV-mixed--0-no_sw'
    grid = 'sb_hv_grid'
    # grid = 'hv_grid'
    type = "distribution"

    # net_name = 'simbench_hv_grid_with_potential.pkl'
    # with open('C:\\Users\\f.lohse\PycharmProjects\potpourri\potpourri\data\\' + net_name,
    #           'rb') as f:
    #     net = pickle.load(f)
    #
    # net_name = grid + "_scenario_types.pkl"
    # with open("../../data/scenarios/" + grid + "_scenario_types.pkl", "rb") as f:
    #     net = pickle.load(f)
    # results_dir = '../../results/test_scenarios/' + grid + '_loadcases_' + type + '_scenarios/'

    with open("../../data/scenarios/" + grid + "_loadcases_" + type + "_scenarios_with_q.pkl", "rb") as f:
        scenarios = pickle.load(f)

    net_name = 'sb_hv_grid_urban'
    net_name = 'sb_hv_grid_with_potential_3MW_230m'
    with open('C:\\Users\\f.lohse\\PycharmProjects\\potpourri\\potpourri\\data\\windpot\\' + net_name + '.pkl',
              'rb') as f:
        net = pickle.load(f)

    case = 'lW'
    factors = net.loadcases.loc[case]
    net.ext_grid.vm_pu = factors['Slack_vm']

    results_dir = 'C:\\Users\\f.lohse\\PycharmProjects\\potpourri\\potpourri\\results\\scenarios\\' + net_name + '\\' + type + '\\'

    cases = [3]
    # obj = run_scenarios(net, scenarios, results_dir, cases=cases, SWmin=10)

    dir = results_dir + 'case_' + str(cases[0]) + '\\'

    # #check if all results are correct
    # for case in cases:
    #     hc_diff = np.diff(np.array(obj))
    #     sc_not_opt = np.where(hc_diff == 0)
    #     sc_not_opt = sc_not_opt[0]+1
    #
    #     output_path = results_dir + 'case_' + str(case) + '_glpk\\'
    #
    #     if sc_not_opt.size > 0:
    #         load_ind = net.load.index[net.load.name == 'load_scenario']
    #         sgen_ind = net.sgen.index[net.sgen.name == 'sgen_scenario']
    #
    #         hc = HC_ACOPF(net)
    #         hc.solve()
    #         hc.add_OPF()
    #
    #         for sc in sc_not_opt[0]+1:
    #             for i, d in enumerate(load_ind):
    #                 hc.model.pD[d] = scenarios['load'][case][i, sc]
    #                 hc.model.qD[d] = scenarios['load_q'][case][i, sc]
    #
    #             for i, g in enumerate(sgen_ind):
    #                 hc.model.psG[g] = scenarios['sgen'][case][i, sc]
    #                 hc.model.qsG[g] = scenarios['sgen_q'][case][i, sc]
    #
    #             hc.solve(solver='mindtpy', mip_solver='glpk')
    #             hc.add_OPF()
    #
    #             pickle.dump(hc.net, open(output_path + 'sc_' + str(sc) + '.pkl', 'wb'))


def plot_supply_task_sc(nets, scenarios, case):
    loads = scenarios['load'][case]
    load_sum_sc = []
    for sc in range(len(loads[0])):
        load_sum_sc.append(-sum(loads[:, sc]))

    hc_sum = []
    for net in nets:
        hc_sum.append(sum(net.res_sgen.p_mw[net.res_sgen.y_wind == 1]))

    # Get the indices that would sort load_sum_sc
    sorted_indices = np.argsort(load_sum_sc)[::-1]  # [::-1] is used to sort in descending order

    # Sort load_sum_sc and hc_sum using the sorted indices
    load_sum_sc = np.array(load_sum_sc)[sorted_indices]
    hc_sum = np.array(hc_sum)[sorted_indices]

    fig, ax = plt.subplots()
    ax.bar(range(len(load_sum_sc)), load_sum_sc, label='Last', color=colors['red'])
    ax.bar(range(len(hc_sum)), hc_sum, label='Netzintegrationspotenzial', color=colors['blue'])
    ax.plot(hc_sum+load_sum_sc, color='black', label='Residuum')

    ax.set_xlabel('Szenario')
    ax.set_ylabel('P [MW]')

    ax.grid()
    # ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.legend(bbox_to_anchor=(0.5, -0.15), loc='upper center')
    fig.set_size_inches((config['textbreite'] , 0.6 * config['textbreite']))