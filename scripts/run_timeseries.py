import copy

import pandas as pd

from potpourri.models.ACOPF_base import ACOPF
import simbench as sb
import pyomo.environ as pe
from pandapower import timeseries as ts
import pickle
from pandapower.control.run_control import get_controller_order

from potpourri.models.HC_ACOPF import HC_ACOPF
from potpourri.models.init_pyo_from_pp_res import init_pyo_from_dcpp
from scripts.plot_functions import *


def create_run_pyomo_opf(net, **kwargs):
    opf = ACOPF(net)
    opf.add_OPF()
    opf.add_voltage_deviation_objective()

    for w in opf.model.WIND:
        opf.model.psG[w].fix()

        if opf.model.PsG[w] == 0:
            opf.model.QW_pos_constraint[w].deactivate()
            opf.model.QW_neg_constraint[w].deactivate()
            opf.model.PsG_Constraint[w].deactivate()
            opf.model.QsGmax[w] = 0.
            opf.model.QsGmin[w] = -0.1 * opf.model.PsG_inst[w]

    opf.solve()
    # plot_pq_res([opf.net])
    # plot_qu_res([opf.net])
    opf.net['converged'] = pe.check_optimal_termination(opf.results)
    # pp.plotting.pf_res_plotly(opf.net)
    net = opf.net
    print(pe.value(opf.model.obj_v_deviation))


def run_timeseries_pyomo_opf(net, **kwargs):
    opf = kwargs['opf']

    p_load = net.load.p_mw.values / opf.baseMVA
    q_load = net.load.q_mvar.values / opf.baseMVA
    p_sgen = net.sgen.p_mw.values / opf.baseMVA

    for d in opf.model.D:
        opf.model.pD[d].fix(p_load[d])
        opf.model.qD[d].fix(q_load[d])

    for g in opf.model.sG:
        opf.model.psG[g].fix(p_sgen[g])

    for w in opf.model.WINDc:
        if opf.model.psG[w].value < 0.1 * opf.model.PsG_inst[w]:
            opf.model.QW_pos_constraint[w].deactivate()
            opf.model.QW_neg_constraint[w].deactivate()
            opf.model.PsG_Constraint[w].deactivate()
            opf.model.QsGmax[w] = 0.
            opf.model.QsGmin[w] = -0.1 * opf.model.PsG_inst[w]

        else:
            opf.model.QW_pos_constraint[w].activate()
            opf.model.QW_neg_constraint[w].activate()
            opf.model.PsG_Constraint[w].activate()
            opf.model.QsGmax[w] = opf.static_generation_data['max_q'][w]
            opf.model.QsGmin[w] = opf.static_generation_data['min_q'][w]

    opf.solve()
    # net = opf.net
    # plot_pq_res([opf.net])
    # plot_qu_res([opf.net])
    if not pe.check_optimal_termination(opf.results):
        init_pyo_from_dcpp(opf.net, opf.model)
        opf.solve()
    net['converged'] = pe.check_optimal_termination(opf.results)
    # pp.plotting.pf_res_plotly(opf.net)

    # print(pe.value(opf.model.obj_v_deviation))
    return net


def create_output_writer(net, time_steps, output_dir=None):
    ow = ts.OutputWriter(net, time_steps, output_path=output_dir, output_file_type=".xlsx")
    ow.log_variable('res_sgen', 'p_mw')
    ow.log_variable('res_sgen', 'q_mvar')

    ow.log_variable('res_load', 'p_mw')
    ow.log_variable('res_load', 'q_mvar')

    ow.log_variable('res_bus', 'vm_pu')
    ow.log_variable('res_bus', 'va_degree')

    ow.log_variable('res_line', 'loading_percent')
    ow.log_variable('res_line', 'i_ka')

    ow.log_variable('res_trafo', 'loading_percent')

    ow.log_variable('res_ext_grid', 'p_mw')
    ow.log_variable('res_ext_grid', 'q_mvar')
    return ow


def add_wind_profile_to_net(net):
    s = [0, 1, 2, 3, 6, 7, 10, 11, 12, 13, 16, 17, 24, 25, 26, 27, 30, 31, 32, 33, 34, 35, 44, 45, 46, 47, 48, 49, 50,
         51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 70, 71, 72, 73, 74, 75, 78, 79, 88, 89, 94, 95, 104, 105,
         106, 107, 108, 109, 118, 119, 130, 131, 134, 135, 137, 139, 140, 141, 142, 152, 153, 154, 155, 156, 157, 160,
         161, 171, 172, 173, 174, 175, 176, 184, 185, 186, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210,
         211, 212, 214, 215, 216, 217, 218, 219, 222, 224, 227, 228, 229, 230, 233, 234, 241, 242, 243, 244, 245, 246,
         247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259, 260, 262, 269, 270, 272, 278, 279, 283, 286,
         287, 288, 289, 292, 293, 294, 295, 304, 305]

    n = [28, 29, 36, 37, 38, 39, 64, 65, 66, 67, 68, 69, 82, 83, 86, 87, 98, 99, 100, 101, 102, 103, 110, 111, 112, 113,
         122, 123, 124, 125, 126, 127, 128, 129, 138, 143, 144, 145, 146, 147, 148, 149, 150, 151, 164, 187, 188, 190,
         191, 192, 193, 194, 195, 196, 197, 213, 223, 235, 236, 237, 239, 240, 261, 263, 264, 265, 266, 267, 268, 273,
         274, 275, 280, 281, 282, 284, 285]

    o = [4, 5, 8, 9, 14, 15, 18, 19, 20, 21, 22, 23, 40, 41, 42, 43, 76, 77, 80, 81, 84, 85, 90, 91, 92, 93, 96, 97,
         114, 115, 116, 117, 120, 121, 132, 133, 136, 158, 159, 162, 163, 165, 166, 167, 168, 169, 170, 177, 178, 179,
         180, 181, 182, 183, 189, 198, 220, 221, 225, 226, 231, 232, 238, 271, 276, 277, 290, 291, 296, 297, 298, 299,
         300, 301, 302, 303]

    buses_s = net.bus.index[net.bus.index.isin(s)]
    buses_n = net.bus.index[net.bus.index.isin(n)]
    buses_o = net.bus.index[net.bus.index.isin(o)]

    net.bus['wind_profile'] = None
    net.bus['wind_profile'][buses_s] = 'WP10'
    net.bus['wind_profile'][buses_n] = 'WP7'
    net.bus['wind_profile'][buses_o] = 'WP4'

    net.sgen.profile[net.sgen.wind_hc] = net.bus['wind_profile'][net.sgen.bus[net.sgen.wind_hc]].values

    return net


if __name__ == "__main__":
    input_hc_net_dir = '../potpourri/results/node_potential/sb_hv_grid_with_potential_3MW_230m_var_1_hc_net.pkl'
    with open (input_hc_net_dir, 'rb') as f:
        net_hc = pickle.load(f)
    net = pickle.load(
        open('../potpourri/data/windpot/sb_hv_grid_with_potential_3MW_230m.pkl', 'rb'))
    pp.runpp(net)

    absolute_values = sb.get_absolute_values(net, False)  # get absolute values for loadcases
    for values in absolute_values.values():
        values.reset_index(inplace=True, drop=True)

    sb.apply_const_controllers(net, absolute_values)
    level_list, controller_order = get_controller_order(net, net.controller)
    ts.run_time_series.control_time_step(controller_order, 5)  # write values from loadcase 'lW' to net

    # calculate hosting capacity
    hc = HC_ACOPF(net)
    hc.solve()
    hc.add_OPF(SWmin=10)
    hc.solve(solver='mindtpy', mip_solver='gurobi')

    # reset net values to base case
    ts.run_time_series.control_time_step(controller_order, 0)

    net_wind = copy.deepcopy(net)

    # create wind generators in original net
    # wind_hc_index = hc.net.sgen.index[hc.net.res_sgen.y_wind == 1]
    # pp.create_sgens(net_wind, hc.net.sgen.bus[wind_hc_index], p_mw=hc.net.sgen.p_mw[wind_hc_index], var_q=0,
    #                 type='Wind', wind_hc=True)
    wind_hc_index = net_hc.sgen.index[net_hc.res_sgen.y_wind == 1]
    pp.create_sgens(net_wind, net_hc.sgen.bus[wind_hc_index], p_mw=net_hc.sgen.p_mw[wind_hc_index], var_q=0,
                    type='Wind', wind_hc=True)

    net_wind.sgen['wind_hc'].fillna(False, inplace=True)
    net_wind = add_wind_profile_to_net(net_wind)

    pp.drop_controllers_at_buses(net_wind, net_wind.bus.index)
    # net_wind.sgen['type'][net_wind.sgen.wind_hc] = 'Wind'

    net_wind.sgen['controllable'] = False
    net_wind.sgen['controllable'][net_wind.sgen.type == 'Wind'] = True
    # net_wind.sgen['in_service'][~net_wind.res_sgen.y_wind.astype(bool)] = False
    net_wind.sgen.wind_hc.fillna(False, inplace=True)
    net_wind.sgen['p_inst_mw'] = net_wind.sgen['p_mw']
    net_wind.sgen['var_q'][net_wind.sgen.type == 'Wind'] = 1
    net_wind.sgen['var_q'][net_wind.sgen.wind_hc] = 0

    pp.runpp(net_wind)

    absolute_values = sb.get_absolute_values(net_wind, True)  # get absolute values for loadcases
    for values in absolute_values.values():
        values.reset_index(inplace=True, drop=True)

    # ts.run_timeseries(net_wind, time_steps, run=create_run_pyomo_opf)
    opf = ACOPF(net_wind)
    opf.add_OPF()
    opf.add_voltage_deviation_objective()
    opf.add_tap_changer_linear()

    n_steps = len(values)
    time_steps = range(0, 35136)
    sb.apply_const_controllers(opf.net, absolute_values)

    output_dir = '../potpourri/results/timeseries/sb_hv_grid_with_potential_3MW_230m'
    ow = create_output_writer(opf.net, time_steps, output_dir=output_dir)

    ts.run_timeseries(opf.net, time_steps, run=run_timeseries_pyomo_opf, opf=opf, continue_on_divergence=True)


def evaluate_results(path):
    line_loading = pd.read_excel(path + '/res_line/loading_percent.xlsx', index_col=0)
    max_line_loading = [max(loadings) for loadings in line_loading.values]

    trafo_loading = pd.read_excel(path + '/res_trafo/loading_percent.xlsx', index_col=0)
    max_trafo_loading = [max(loadings) for loadings in trafo_loading.values]

    sgen_p_mw = pd.read_excel(path + '/res_sgen/p_mw.xlsx', index_col=0)
    sum_p_mw = [sum(p_mw) for p_mw in sgen_p_mw.values]
    load_p_mw = pd.read_excel(path + '/res_load/p_mw.xlsx', index_col=0)
    sum_load_p_mw = [sum(p_mw) for p_mw in load_p_mw.values]



    fig, ax = plt.subplots()
    ax.plot(max_trafo_loading[0:672], label='max trafo loading')
    ax.plot(sum_p_mw[0:672], label='sum p_mw')
    ax.plot(sum_load_p_mw[0:672], label='sum load p_mw')
    ax.legend()