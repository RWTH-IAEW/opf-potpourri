import numpy as np
import pandapower as pp
import pandas as pd

def pyo_sol_to_net_res(net, model):
    pp.clear_result_tables(net)

    if 'DC' in model.name:
        net.res_bus.vm_pu = pd.Series([1.] * len(net.bus.index), net.bus.index)
        # use values from net definition as voltage not calculated in DCLF calculation, to get same result tables as pandapower
        net.res_bus.vm_pu[net.gen.bus] = net.gen.vm_pu
        net.res_bus.vm_pu[net.ext_grid.bus] = net.ext_grid.vm_pu

    else:
        net.res_bus.vm_pu = model.v.get_values()

        # bus reactive power
        for b in model.B:
            net.res_bus.loc[b, 'q_mvar'] = -(sum(model.qLfrom[l].value for l in model.L if model.A[l, 1] == b) +
                                             sum(model.qLto[l].value for l in model.L if model.A[l, 2] == b) +
                                             sum(model.qThv[l].value for l in model.TRANSF if model.AT[l, 1] == b) +
                                             sum(model.qTlv[l].value for l in model.TRANSF if
                                                 model.AT[l, 2] == b)) * model.baseMVA.value

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

        # static generator
        for g in model.sG:
            net.res_sgen.loc[g, 'q_mvar'] = model.qG[g].value * model.baseMVA.value

        # generator
        ind_offset = len(model.sG)
        for g in model.gG:
            net.res_gen.loc[g, 'q_mvar'] = model.qG[g].value * model.baseMVA.value

        # transformer
        net.res_trafo.q_hv_mvar = model.qThv.get_values()
        net.res_trafo.q_hv_mvar *= model.baseMVA.value
        net.res_trafo.q_lv_mvar = model.qTlv.get_values()
        net.res_trafo.q_lv_mvar *= model.baseMVA.value
        net.res_trafo.ql_mvar = net.res_trafo.q_hv_mvar + net.res_trafo.q_lv_mvar

        # shunt
        net.res_shunt.vm_pu = pd.Series(net.res_bus.vm_pu[net.shunt.bus].values, net.shunt.index)
        for s in model.SHUNT:
            net.res_shunt.loc[s, 'q_mvar'] = -model.BB[s] * net.res_bus.vm_pu[
                net.shunt['bus'][s]] ** 2 * model.baseMVA.value

    # --- buses ---
    for b in model.B:
        net.res_bus.loc[b, 'p_mw'] = - (sum(model.pLfrom[l].value for l in model.L if model.A[l, 1] == b) +
                                        sum(model.pLto[l].value for l in model.L if model.A[l, 2] == b) +
                                        sum(model.pThv[l].value for l in model.TRANSF if model.AT[l, 1] == b) +
                                        sum(model.pTlv[l].value for l in model.TRANSF if
                                            model.AT[l, 2] == b)) * model.baseMVA.value

    net.res_bus.va_degree = model.delta.get_values()
    net.res_bus.va_degree *= 180 / np.pi

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

    # --- transformer ---
    # real power
    net.res_trafo.p_hv_mw = model.pThv.get_values()
    net.res_trafo.p_hv_mw *= model.baseMVA.value
    net.res_trafo.p_lv_mw = model.pTlv.get_values()
    net.res_trafo.p_lv_mw *= model.baseMVA.value
    net.res_trafo.pl_mw = net.res_trafo.p_hv_mw + net.res_trafo.p_lv_mw

    # voltage
    net.res_trafo.vm_hv_pu = net.res_bus.vm_pu[net.trafo.hv_bus].values
    net.res_trafo.vm_lv_pu = net.res_bus.vm_pu[net.trafo.lv_bus].values

    # current
    net.res_trafo.i_hv_ka = np.sqrt(net.res_trafo.p_hv_mw ** 2 + net.res_trafo.q_hv_mvar.fillna(0) ** 2) / (
            net.res_trafo.vm_hv_pu * np.sqrt(3) * net.bus.vn_kv[net.trafo.hv_bus].values)
    net.res_trafo.i_lv_ka = np.sqrt(net.res_trafo.p_lv_mw ** 2 + net.res_trafo.q_lv_mvar.fillna(0) ** 2) / (
            net.res_trafo.vm_lv_pu * np.sqrt(3) * net.bus.vn_kv[net.trafo.lv_bus].values)

    vns = np.vstack([net.trafo.vn_hv_kv.values, net.trafo.vn_lv_kv.values]).T
    lds_trafo = net.res_trafo[['i_hv_ka', 'i_lv_ka']] * vns * np.sqrt(3) / net.trafo.sn_mva.values[:, np.newaxis] * 100.
    with np.errstate(invalid='ignore'):
        ld_trafo = np.max(lds_trafo, axis=1)
    net.res_trafo.loading_percent = ld_trafo / net.trafo.parallel.values / net.trafo.df.values

    # voltage angle
    net.res_trafo.va_hv_degree = net.res_bus.va_degree[net.trafo.hv_bus].values
    net.res_trafo.va_lv_degree = net.res_bus.va_degree[net.trafo.lv_bus].values

    # --- external grid ---
    net.res_ext_grid.p_mw = model.peG.get_values()
    net.res_ext_grid.p_mw *= model.baseMVA.value

    # --- load ---
    net.res_load.p_mw = model.pD.get_values()
    net.res_load.p_mw *= model.baseMVA.value

    # --- sgen ---
    for g in model.sG:
        net.res_sgen.loc[g, 'p_mw'] = model.pG[g].value * model.baseMVA.value

    # --- gen ---
    ind_offset = len(model.sG)
    for g in model.gG:
        net.res_gen.loc[g, 'p_mw'] = model.pG[g].value * model.baseMVA.value

    net.res_gen.va_degree = net.res_bus.va_degree[net.gen.bus].values
    net.res_gen.vm_pu = net.res_bus.vm_pu[net.gen.bus].values

    # --- shunt ---
    for s in model.SHUNT:
        net.res_shunt.loc[s, 'p_mw'] = model.GB[s] * net.res_bus.vm_pu[net.shunt['bus'][s]] ** 2 * model.baseMVA.value

    net.res_line.q_from_mvar.fillna(0., inplace=True)
    net.res_line.q_to_mvar.fillna(0., inplace=True)
    net.res_line.ql_mvar.fillna(0., inplace=True)
    net.res_trafo.fillna(0., inplace=True)
