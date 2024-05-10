import pandapower as pp
import math
def init_pyo_from_dcpp(net, model):
    pp.rundcpp(net)

    bus_lookup = net._pd2ppc_lookups["bus"]

    for b in net.bus.index:
        model.delta[bus_lookup[b]] = net.res_bus.va_degree[b] * math.pi / 180
        if 'AC' in model.name:
            model.v[bus_lookup[b]] = net.res_bus.vm_pu[b]

    for l in model.L:
        model.pLfrom[l] = net.res_line.p_from_mw[l] / model.baseMVA.value
        model.pLto[l] = net.res_line.p_to_mw[l] / model.baseMVA.value
        # if 'AC' in model.name:
        #     model.qLfrom[l] = net.res_line.q_from_mvar[l] / model.baseMVA.value
        #     model.qLto[l] = net.res_line.q_to_mvar[l] / model.baseMVA.value

    for m_ind, net_ind in enumerate(net.trafo.index):
        model.pThv[m_ind] = net.res_trafo.p_hv_mw[net_ind] / model.baseMVA.value
        model.pTlv[m_ind] = net.res_trafo.p_lv_mw[net_ind] / model.baseMVA.value
        # if 'AC' in model.name:
        #     model.qThv[t] = net.res_trafo.q_hv_mvar[t] / model.baseMVA.value
        #     model.qTlv[t] = net.res_trafo.q_lv_mvar[t] / model.baseMVA.value

    for m_ind, net_ind in enumerate(net.ext_grid.index):
        model.pG[m_ind] = net.res_ext_grid.p_mw[net_ind] / model.baseMVA.value
        # if 'AC' in model.name:
        #     model.qeG[e] = net.res_ext_grid.q_mvar[e] / model.baseMVA.value




