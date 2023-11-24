import numpy as np


def pyo_sol_to_net_res(net, model):
    if 'DC' in model.name:
        net.res_bus.vm_pu = [1.] * len(net.bus.index)
    else:
        net.res_bus.vm_pu = model.v.get_values()

        # bus reactive power
        for b in model.B:
            net.res_bus.loc[b, 'q_mvar'] = -(sum(model.qLfrom[l].value for l in model.L if model.A[l, 1] == b) +
                                             sum(model.qLto[l].value for l in model.L if model.A[l, 2] == b) +
                                             sum(model.qLfromT[l].value for l in model.TRANSF if model.AT[l, 1] == b) +
                                             sum(model.qLtoT[l].value for l in model.TRANSF if model.AT[l, 2] == b))

        # reactive power on lines
        net.res_line.q_from_mvar = model.qLfrom.get_values()
        net.res_line.q_from_mvar *= model.baseMVA.value
        net.res_line.q_to_mvar = model.qLto.get_values()
        net.res_line.q_to_mvar *= model.baseMVA.value
        net.res_line.ql_mvar = net.res_line.q_from_mvar + net.res_line.q_to_mvar

        # external grid
        net.res_ext_grid.q_mvar = model.qeG.get_values()
        net.res_ext_grid.q_mvar *= model.baseMVA.value
        # load
        net.res_load.q_mvar = model.qD.get_values()
        net.res_load.q_mvar *= model.baseMVA.value
        # generator
        net.res_sgen.q_mvar = model.qG.get_values()
        net.res_sgen.q_mvar *= model.baseMVA.value
        # shunt
        for s in model.SHUNT:
            net.res_shunt.loc[s, 'q_mvar'] = -model.BB[s] * net.res_bus.vm_pu[net.shunt.bus].values ** 2

    net.res_bus.va_degree = model.delta.get_values()
    net.res_bus.va_degree *= 180 / np.pi

    # --- buses ---
    for b in model.B:
        net.res_bus.loc[b, 'p_mw'] = -(sum(model.pLfrom[l].value for l in model.L if model.A[l, 1] == b) +
                                       sum(model.pLto[l].value for l in model.L if model.A[l, 2] == b) +
                                       sum(model.pLfromT[l].value for l in model.TRANSF if model.AT[l, 1] == b) +
                                       sum(model.pLtoT[l].value for l in model.TRANSF if model.AT[l, 2] == b))

    # --- lines ---
    # voltage on lines
    net.res_line.vm_from_pu = net.res_bus.vm_pu[net.line.from_bus].values
    net.res_line.vm_to_pu = net.res_bus.vm_pu[net.line.to_bus].values

    # real power on lines
    net.res_line.p_from_mw = model.pLfrom.get_values()
    net.res_line.p_from_mw *= model.baseMVA.value
    net.res_line.p_to_mw = model.pLto.get_values()
    net.res_line.p_to_mw *= model.baseMVA.value
    net.res_line.pl_mw = net.res_line.p_from_mw + net.res_line.p_to_mw

    # voltage angle
    net.res_line.va_from_degree = net.res_bus.va_degree[net.line.from_bus].values
    net.res_line.va_to_degree = net.res_bus.va_degree[net.line.to_bus].values

    # current
    net.res_line.i_from_ka = np.sqrt(net.res_line.p_from_mw ** 2 + net.res_line.q_from_mvar.fillna(0) ** 2) / (
            net.res_line.vm_from_pu * np.sqrt(3) * net.bus.vn_kv[net.line.from_bus].values)
    net.res_line.i_to_ka = np.sqrt(net.res_line.p_to_mw ** 2 + net.res_line.q_to_mvar.fillna(0) ** 2) / (
            net.res_line.vm_to_pu * np.sqrt(3) * net.bus.vn_kv[net.line.to_bus].values)
    net.res_line.i_ka = net.res_line[['i_from_ka', 'i_to_ka']].max(axis=1)
    net.res_line.loading_percent = net.res_line.i_ka * 100 / (net.line.max_i_ka * net.line.df * net.line.parallel)

    # --- external grid ---
    net.res_ext_grid.p_mw = model.peG.get_values()
    net.res_ext_grid.p_mw *= model.baseMVA.value

    # --- load ---
    net.res_load.p_mw = model.pD.get_values()
    net.res_load.p_mw *= model.baseMVA.value

    # --- sgen ---
    net.res_sgen.p_mw = model.pG.get_values()
    net.res_sgen.p_mw *= model.baseMVA.value

    # --- shunt ---
    net.res_shunt.vm_pu = net.res_bus.vm_pu[net.shunt.bus].values
    for s in model.SHUNT:
        net.res_shunt.loc[s, 'p_mw'] = model.GB[s] * net.res_bus.vm_pu[net.shunt.bus].values ** 2
