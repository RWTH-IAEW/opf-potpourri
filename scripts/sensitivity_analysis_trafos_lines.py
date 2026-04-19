import numpy as np
import pickle
from tqdm import tqdm

from src.potpourri.models.HC_ACOPF import HC_ACOPF
from src.potpourri.plotting.plot_functions import set_plt_config
import matplotlib.pyplot as plt


def plot_cmap_loadings(hc_results, step=5, x_ticks_max=120, y_ticks_max=120, decimals=0, unit='percent', ticks_step=2):
    config = set_plt_config()

    # Create the figure and axis
    fig, ax = plt.subplots()

    # Plot the heatmap
    cax = ax.imshow(hc_results, cmap='viridis', aspect='equal', origin='lower')

    # Set up the ticks
    yticks_len = len(hc_results)
    yticks_min = y_ticks_max - yticks_len * step
    xticks_len = len(hc_results[0])
    xticks_min = x_ticks_max - xticks_len * step

    yticks = np.around(np.linspace(y_ticks_max, yticks_min + step, yticks_len), decimals)
    xticks = np.around(np.linspace(xticks_min + step, x_ticks_max, xticks_len), decimals)

    if unit == 'percent':
        ax.set_xlabel("Max. Transformator Auslastung [%]")
        ax.set_ylabel("Max. Leitungsauslastung [%]")
        xticks = xticks.astype(int)
        yticks = yticks.astype(int)
    else:
        ax.set_xlabel("Max. Transformator Auslastung [p.u.]")
        ax.set_ylabel("Max. Leitungsauslastung [p.u.]")

    # Customize the ticks
    ax.set_xticks(np.arange(0, xticks_len, ticks_step))
    ax.set_xticklabels(xticks[::ticks_step])
    ax.set_yticks(np.arange(0, yticks_len, ticks_step))
    ax.set_yticklabels(yticks[::ticks_step])

    # Add colorbar
    cbar = fig.colorbar(cax, ax=ax)
    cbar.set_label('Netzintegrationspotenzial [MW]')

    # Adjust figure size
    fig.set_size_inches((config['textbreite'] * 0.9, 0.5 * config['textbreite']))

    return ax

def calc_hc_for_loading(hc, loading_steps):
    results_list = []
    for line_loading in tqdm(loading_steps):
        # set maximum line loadings
        for l in hc.model.L:
            hc.model.SLmax[l] = hc.line_data.loc[l, 'SLmax_data'] * line_loading

        results_trafo = []
        for trafo_loading in loading_steps:
            print('Line loading: ', line_loading, 'Trafo loading: ', trafo_loading)
            # set maximum transformer loadings
            for t in hc.model.TRANSF:
                hc.model.SLmaxT[t] = hc.trafo_data.loc[t, 'SLmaxT_data'] * trafo_loading

            # hc.solve(solver='mindtpy', mip_solver='gurobi', load_solutions=False)
            hc.solve()
            # save results
            # pickle.dump(hc.net, open(
            #     results_dir + 'hc_loading_line_' + str(line_loading) + '_trafo_' + str(trafo_loading) + '.pkl', 'wb'))

            # save hc result if optimal solution was found
            if hc.results.solver.termination_condition != 'optimal':
                hc_result = 0
            else:
                hc_result = sum(hc.model.psG[w]() for w in hc.model.WIND_HC)
            results_trafo.append(hc_result)

        results_list.append(results_trafo)

    return results_list

if __name__ == "__main__":
    results_dir = '../potpourri/results/node_potential/sensi'

    net_input_dir = '../potpourri/data/windpot/sb_hv_grid_with_potential_3MW_230m.pkl'
    net = pickle.load(open(net_input_dir, 'rb'))

    # apply loadcase 'lW' to net
    case = 'lW'
    factors = net.loadcases.loc[case]
    net.load.p_mw *= factors['pload']
    net.load.q_mvar *= factors['qload']
    net.sgen.scaling[net.sgen.type == 'Wind'] = factors['Wind_p']
    net.sgen.scaling[net.sgen.type == 'PV'] = factors['PV_p']
    net.sgen.scaling[(net.sgen.type != 'Wind') & (net.sgen.type != 'Solar')] = factors['RES_p']
    net.ext_grid.vm_pu = factors['Slack_vm']

    # make sure max loadings are at 100%
    net.line.max_loading_percent = 100.0
    net.trafo.max_loading_percent = 100.0

    # add wind control variant to existing wind generators
    net.sgen['controllable'] = False
    net.sgen['controllable'][net.sgen.type == 'Wind'] = True
    net.sgen['p_inst_mw'] = net.sgen['p_mw']
    net.sgen['var_q'] = None
    net.sgen['var_q'][net.sgen.type == 'Wind'] = 1

    # deactivate existing generators
    # net.sgen.in_service = False

    hc = HC_ACOPF(net)
    hc.solve()
    hc.add_OPF()
    hc.fix_vars('y', 1)
    # only Q variable for existing wind generators
    for w in hc.model.WINDc:
        hc.model.psG[w].fix()
    hc.add_tap_changer_linear()

    min_loading = 0.7
    max_loading = 2
    step = 0.05
    num_steps = int((max_loading-min_loading)/step +1)
    loading_steps = np.around(np.linspace(2, 0.7, num_steps), 2)
    hc_results = calc_hc_for_loading(hc, loading_steps)
    results_list = [lst[::-1] for lst in hc_results if any(lst)]

    ax = plot_cmap_loadings(results_list, step=int(step*100), decimals=2, ticks_max=max_loading*100)
    # ax.figure.savefig(results_dir + 'sensi_loading_3_230.svg', format='svg')

    res_trafos = [[] for _ in range(len(results_list[0]))]
    for list in results_list:
        for i, res in enumerate(list):
            res_trafos[i].append(res)
    res_trafos = [lst[::-1] for lst in res_trafos]
    res_trafos.reverse()
