import math
import pandapower as pp


def init_pyo_from_pp_res(net, model):
    """Warm-start Pyomo variables from a pandapower DC power flow result.

    Runs pp.rundcpp() on the network and copies bus angles, line flows,
    transformer flows, and generator injections into the Pyomo model variables
    as initial values. For AC models, voltage magnitudes are also initialised.

    Args:
        net: pandapower network (will be solved in-place by rundcpp).
        model: Pyomo ConcreteModel with variables delta, pLfrom, pLto,
               pThv, pTlv, pG, and optionally v.
    """
    pp.rundcpp(net)

    bus_lookup = net._pd2ppc_lookups["bus"]

    # bus angles and voltage magnitudes
    for b in net.bus.index:
        model.delta[bus_lookup[b]] = net.res_bus.va_degree[b] * math.pi / 180
        if hasattr(model, 'v'):
            model.v[bus_lookup[b]] = net.res_bus.vm_pu[b]

    # line flows — Pyomo L index matches net.line.index
    for l in model.L:
        model.pLfrom[l] = net.res_line.p_from_mw[l] / model.baseMVA.value
        model.pLto[l] = net.res_line.p_to_mw[l] / model.baseMVA.value

    # transformer flows — trafo_data uses 0-based position indices
    for pos, net_ind in enumerate(net.trafo.index):
        if pos in model.TRANSF:
            model.pThv[pos] = net.res_trafo.p_hv_mw[net_ind] / model.baseMVA.value
            model.pTlv[pos] = net.res_trafo.p_lv_mw[net_ind] / model.baseMVA.value

    # generator injections — G index follows _gen_order contiguous ranges
    for element, (f, t) in net._gen_order.items():
        res_table = net['res_' + element]
        for i, net_ind in enumerate(net[element].index):
            pyo_ind = f + i
            if pyo_ind in model.G:
                model.pG[pyo_ind] = res_table.p_mw.iloc[i] / model.baseMVA.value
