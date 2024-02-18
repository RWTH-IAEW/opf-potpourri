import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import pickle

from potpourri.scripts.classbased.plot_functions import plot_wind_hc_results


def scenario_supply(scenarios):
    n_scenarios = len(scenarios['load'][0][0])
    n_cases = len(scenarios['load'])

    scenarios['supply'] = [[] for _ in range(n_cases)]  # Initialize a list for each case
    scenarios['supply_sum'] = [[] for _ in range(n_cases)]  # Initialize a list for each case

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
    p_mw_file = os.path.join(dir, 'case_' + str(case), 'res_sgen', 'p_mw.xlsx')
    p_mw = pd.read_excel(p_mw_file, index_col=0)


def create_weighted_p_wind(probabilities, p_wind):
    min_prob = np.min(probabilities)
    weights = probabilities / min_prob * 3
    weights = weights.round()

    p_wind_weighted = []
    for p_w in p_wind:
        p_w_weighted = np.array([])
        for i, w in enumerate(weights):
            p_w_weighted = np.concatenate((p_w_weighted, np.repeat(p_w[i], w.astype(int))))

        p_wind_weighted.append(p_w_weighted)

    fig, ax = plt.subplots()
    ax.boxplot(p_wind_weighted)
