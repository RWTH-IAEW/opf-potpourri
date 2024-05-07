from potpourri.models.HC_ACOPF import HC_ACOPF
from scripts.prepare_net import *

import pyomo.environ as pe
import pickle
import matplotlib.pyplot as plt

from scripts.plot_functions import set_plt_config


def add_wind_generation_box_constraints(hc):
    def q_p_box_constraint_lower(model, w):
        return model.qsG[w] >= - model.psG[w]

    def q_p_box_constraint_upper(model, w):
        return model.qsG[w] <= model.psG[w]

    hc.model.QP_hc_box_constraint_lower = pe.Constraint(hc.model.WIND_HC, rule=q_p_box_constraint_lower)
    hc.model.QP_hc_box_constraint_upper = pe.Constraint(hc.model.WIND_HC, rule=q_p_box_constraint_upper)

    for w in hc.model.WINDc:
        hc.model.QsGmin[w] = -hc.model.psG[w]()
        hc.model.QsGmax[w] = hc.model.psG[w]()

    return hc


def deactivate_wind_box_constraints(hc):
    try:
        hc.model.QP_hc_box_constraint_lower.deactivate()
        hc.model.QP_hc_box_constraint_upper.deactivate()
    except AttributeError:
        pass

    for w in hc.model.WINDc:
        hc.model.QsG_Constraint[w].deactivate()


def add_wind_cos_phi_control(hc, cos_phi):
    def cos_phi_control(model, w):
        return model.qsG[w] == model.psG[w] * pe.tan(pe.acos(cos_phi))

    hc.model.QP_hc_cos_phi_control = pe.Constraint(hc.model.WIND_HC, rule=cos_phi_control)

    for w in hc.model.WINDc:
        q_set_point = hc.model.psG[w]() * pe.tan(pe.acos(cos_phi))
        hc.model.qsG[w].fix(q_set_point)

    return hc

def deactivate_wind_cos_phi_control(hc):
    try:
        hc.model.QP_hc_cos_phi_control.deactivate()
    except AttributeError:
        pass

    for w in hc.model.WINDc:
        hc.model.qsG[w].unfix()

    return hc


def load_result_nets(dir, names):
    nets = []
    hc_results = []
    for name in names:
        with open(dir + name, 'rb') as f:
            net = pickle.load(f)
        nets.append(net)
        hc_results.append(sum(net.sgen.p_mw[net.sgen.wind_hc]))

    return nets, hc_results


def plot_hc_results(hc_results, path = None):
    config = set_plt_config()
    fig, ax = plt.subplots()
    ax.bar(['-P < Q < P','Q(U)', 'cos($\phi$) = 0,95', 'Q = 0'], hc_results)
    ax.grid()
    ax.set_ylabel('Netzintegrationspotenzial [MW]')
    fig.set_size_inches((config['textbreite'] * 0.7, 0.5 * config['textbreite']))
    plt.show()



if __name__ == "__main__":

    net_input_dir = 'potpourri/data/windpot/sb_hv_grid_with_potential_3MW_230m.pkl'
    with open(net_input_dir, 'rb') as f:
        net = pickle.load(f)

    results_dir = 'potpourri/results'

    # apply loadcase 'lW' to net
    case = 'lW'
    net = apply_loadcase_to_sb_net(net, case)

    add_regulatory_q_control_to_wind(net, 1)

    hc_res = []

    hc = HC_ACOPF(net)
    hc.solve()
    hc.add_OPF(SWmin=10)
    hc.add_tap_changer_linear()

    for w in hc.model.WINDc:
        hc.model.psG[w].fix()
        hc.model.PsG_Constraint[w].deactivate()

    hc.solve(solver='mindtpy', mip_solver='gurobi')
    hc_res.append(sum(hc.model.psG[w]() for w in hc.model.WIND_HC))
    with open(results_dir + 'var_q_1', 'wb') as f:
        pickle.dump(hc.net, f)

    # deactivate regulatory q control constraints
    hc.model.QW_min_constraint.deactivate()
    hc.model.QW_max_constraint.deactivate()
    hc.model.QU_min_hc_constraint.deactivate()
    hc.model.QU_max_hc_constraint.deactivate()

    hc.model.QU_min_constraint.deactivate()
    hc.model.QU_max_constraint.deactivate()
    hc.model.QW_pos_constraint.deactivate()
    hc.model.QW_neg_constraint.deactivate()

    # calculate hc for box constraints
    add_wind_generation_box_constraints(hc)
    hc.solve(solver='mindtpy', mip_solver='gurobi')
    hc_res.append(sum(hc.model.psG[w]() for w in hc.model.WIND_HC))
    with open(results_dir + 'box_constraints_net', 'wb') as f:
        pickle.dump(hc.net, f)
    deactivate_wind_box_constraints(hc)

    # calculate hc for fixed cos(phi) q control
    add_wind_cos_phi_control(hc, 0.95)
    hc.solve(solver='mindtpy', mip_solver='gurobi')
    hc_res.append(sum(hc.model.psG[w]() for w in hc.model.WIND_HC))
    with open(results_dir + 'cos_phi_095_net', 'wb') as f:
        pickle.dump(hc.net, f)
    deactivate_wind_cos_phi_control(hc)

    # calculate hc for 0 q
    for w in hc.model.WINDc:
        hc.model.qsG[w].fix(0.)

    for w in hc.model.WIND_HC:
        hc.model.qsG[w].fix(0.)

    hc.solve(solver='mindtpy', mip_solver='gurobi')
    hc_res.append(sum(hc.model.psG[w]() for w in hc.model.WIND_HC))
    with open(results_dir + 'q_0_net', 'wb') as f:
        pickle.dump(hc.net, f)


    names = ['box_constraints_net', 'var_q_1', 'cos_phi_095_net', 'q_0_net']
    #
    # for net in nets:
    #     net.sgen.p_inst_mw[net.res_sgen.y_wind == 1] = net.res_sgen.p_mw[net.res_sgen.y_wind == 1]
    #     plot_qu_res([net])



