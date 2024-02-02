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

    for t in model.TRANSF:
        model.pThv[t] = net.res_trafo.p_hv_mw[t] / model.baseMVA.value
        model.pTlv[t] = net.res_trafo.p_lv_mw[t] / model.baseMVA.value
        # if 'AC' in model.name:
        #     model.qThv[t] = net.res_trafo.q_hv_mvar[t] / model.baseMVA.value
        #     model.qTlv[t] = net.res_trafo.q_lv_mvar[t] / model.baseMVA.value

    for e in model.eG:
        model.peG[e] = net.res_ext_grid.p_mw[e] / model.baseMVA.value
        # if 'AC' in model.name:
        #     model.qeG[e] = net.res_ext_grid.q_mvar[e] / model.baseMVA.value




