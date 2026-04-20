"""Post-processing for multi-period models: reads Pyomo solution at time t
into net.res_* DataFrames."""

import numpy as np
import pandapower as pp
import pandas as pd


def pyo_sol_to_net_res(net, model, t):
    """Write Pyomo solution for time step t to net.res_* DataFrames.

    Args:
        net: pandapower network whose res_* tables will be populated.
        model: solved Pyomo ConcreteModel.
        t: time step index to extract from multi-period variables.
    """
    if "HC" in model.name:
        for w in model.WIND_HC:
            net.sgen.p_mw[w] = (
                model.psG[w].value * model.baseMVA.value * model.y[w].value
            )
            net.sgen.q_mvar[w] = (
                model.qsG[w].value * model.baseMVA.value * model.y[w].value
            )

    pp.clear_result_tables(net)

    _bus_voltage_results_to_net(net, model, t)
    _line_results_to_net(net, model, t)
    _generation_results_to_net(net, model, t)
    _sgen_results_to_net(net, model, t)
    _load_results_to_net(net, model, t)
    _trafo_results_to_net(net, model, t)
    _shunt_results_to_net(net, model, t)
    bus_pq = _get_bus_power_results(net)
    net.res_bus.p_mw = bus_pq[:, 0]
    if "AC" in model.name:
        net.res_bus.q_mvar = bus_pq[:, 1]


def _bus_voltage_results_to_net(net, model, t):
    bus_lookup = net._pd2ppc_lookups["bus"]
    bus_idx = bus_lookup[net.bus.index.values]

    if "DC" in model.name:
        net.res_bus.vm_pu = pd.Series(
            [1.0] * len(net.bus.index), net.bus.index
        )
        # use values from net definition as voltage not calculated in DCLF
        # calculation, to get same result tables as pandapower
        net.res_bus.vm_pu[net.gen.bus] = net.gen.vm_pu
        net.res_bus.vm_pu[net.ext_grid.bus] = net.ext_grid.vm_pu

    else:
        v = model.v.get_values()
        v_res = [v[b, t] for b in bus_idx]
        net.res_bus.vm_pu = pd.Series(v_res, index=net.bus.index)

    va = model.delta.get_values()
    va_res = [va[b, t] for b in bus_idx]
    net.res_bus.va_degree = va_res
    net.res_bus.va_degree *= 180 / np.pi


def _get_bus_power_results(net):
    bus_pq = np.zeros(shape=(len(net["bus"].index), 2), dtype=np.float64)
    elements = ["load", "sgen", "gen", "ext_grid", "shunt"]

    b = []
    p = []
    q = []

    for element in elements:
        res_ = "res_" + element
        b_el = net[element]["bus"].values
        p_el = net[res_]["p_mw"].values
        q_el = net[res_]["q_mvar"].values

        if element.endswith("gen") or element.endswith("ext_grid"):
            p = np.hstack([p, -p_el])
            q = np.hstack([q, -q_el])
        else:
            p = np.hstack([p, p_el])
            q = np.hstack([q, q_el])
        b = np.hstack([b, b_el])

    dfgr = pd.DataFrame(np.vstack([b, p, q]).T, columns=["bus", "p", "q"])
    dfgr = dfgr.groupby("bus").sum(min_count=1)
    bus_idx = dfgr.index.to_numpy().astype(int)

    maxBus = max(net["bus"].index.values)
    bus_lookup_aranged = -np.ones(maxBus + 1, dtype=np.int64)
    bus_lookup_aranged[net["bus"].index.values] = np.arange(
        len(net["bus"].index.values)
    )

    b_i = bus_lookup_aranged[bus_idx]

    # assign p and q values to bus_pq according to dfgr['bus']
    bus_pq[b_i, 0] = dfgr["p"].values
    bus_pq[b_i, 1] = dfgr["q"].values

    return bus_pq


def _line_results_to_net(net, model, t):
    # --- lines ---
    # voltage on lines
    net.res_line.vm_from_pu = pd.Series(
        net.res_bus.vm_pu[net.line.from_bus].values, index=net.line.index
    )
    net.res_line.vm_to_pu = net.res_bus.vm_pu[net.line.to_bus].values

    # voltage angle
    net.res_line.va_from_degree = net.res_bus.va_degree[
        net.line.from_bus
    ].values
    net.res_line.va_to_degree = net.res_bus.va_degree[net.line.to_bus].values

    # real power on lines
    net.res_line.p_from_mw = pd.Series(
        [model.pLfrom[l, t].value for l in net.line.index],
        index=net.line.index,
    )
    net.res_line.p_from_mw *= model.baseMVA.value
    net.res_line.p_to_mw = pd.Series(
        [model.pLto[l, t].value for l in net.line.index], index=net.line.index
    )
    net.res_line.p_to_mw *= model.baseMVA.value
    net.res_line.pl_mw = net.res_line.p_from_mw + net.res_line.p_to_mw

    if "AC" in model.name:
        # reactive power on lines
        net.res_line.q_from_mvar = pd.Series(
            [model.qLfrom[l, t].value for l in net.line.index],
            index=net.line.index,
        )
        net.res_line.q_from_mvar *= model.baseMVA.value
        net.res_line.q_to_mvar = pd.Series(
            [model.qLto[l, t].value for l in net.line.index],
            index=net.line.index,
        )
        net.res_line.q_to_mvar *= model.baseMVA.value
        net.res_line.ql_mvar = (
            net.res_line.q_from_mvar + net.res_line.q_to_mvar
        )

    else:
        net.res_line.q_from_mvar.fillna(0.0, inplace=True)
        net.res_line.q_to_mvar.fillna(0.0, inplace=True)
        net.res_line.ql_mvar.fillna(0.0, inplace=True)

    # current
    net.res_line.i_from_ka = np.sqrt(
        net.res_line.p_from_mw**2 + net.res_line.q_from_mvar.fillna(0) ** 2
    ) / (
        net.res_line.vm_from_pu
        * np.sqrt(3)
        * net.bus.vn_kv[net.line.from_bus].values
    )
    net.res_line.i_to_ka = np.sqrt(
        net.res_line.p_to_mw**2 + net.res_line.q_to_mvar.fillna(0) ** 2
    ) / (
        net.res_line.vm_to_pu
        * np.sqrt(3)
        * net.bus.vn_kv[net.line.to_bus].values
    )
    net.res_line.i_ka = net.res_line[["i_from_ka", "i_to_ka"]].max(axis=1)
    net.res_line.loading_percent = (
        net.res_line.i_ka
        * 100
        / (net.line.max_i_ka * net.line.df * net.line.parallel)
    )

    net.res_line.fillna(0, inplace=True)


def _ext_grid_results_to_net(net, model, t):
    # --- external grid ---
    net.res_ext_grid.p_mw = pd.Series(
        model.peG[:, t].get_values(), index=net.ext_grid.index
    )
    net.res_ext_grid.p_mw *= model.baseMVA.value
    if "AC" in model.name:
        # external grid
        net.res_ext_grid.q_mvar = model.qeG[:, t].get_values()
        net.res_ext_grid.q_mvar *= model.baseMVA.value


def _generation_results_to_net(net, model, t):
    pg = model.pG.get_values()
    for gen, ord in net._gen_order.items():
        net["res_" + gen].p_mw = [pg[i, t] for i in range(ord[0], ord[1])]
        net["res_" + gen].p_mw *= model.baseMVA.value

        if "AC" in model.name:
            qg = model.qG.get_values()
            net["res_" + gen].q_mvar = [
                qg[i, t] for i in range(ord[0], ord[1])
            ]
            net["res_" + gen].q_mvar *= model.baseMVA.value

        net["res_" + gen].set_index(net[gen].index, inplace=True)

    net.res_gen.va_degree = net.res_bus.va_degree[net.gen.bus].values
    net.res_gen.vm_pu = net.res_bus.vm_pu[net.gen.bus].values


def _load_results_to_net(net, model, t):
    net.res_load = pd.DataFrame(
        columns=["p_mw", "q_mvar"], index=net.load.index
    )
    # --- load ---
    net.res_load.p_mw = pd.Series(
        [model.pD[d, t].value for d in net.load.index], index=net.load.index
    )
    net.res_load.p_mw *= model.baseMVA.value
    if "AC" in model.name:
        # load
        net.res_load.q_mvar = pd.Series(
            [model.qD[d, t].value for d in net.load.index],
            index=net.load.index,
        )
        net.res_load.q_mvar *= model.baseMVA.value
        net.res_load.q_mvar.fillna(
            net.load.q_mvar * net.load.scaling * net.load.in_service,
            inplace=True,
        )

    net.res_load.set_index(net.load.index, inplace=True)
    net.res_load.p_mw.fillna(
        net.load.p_mw * net.load.scaling * net.load.in_service, inplace=True
    )


def _sgen_results_to_net(net, model, t):
    # --- sgen ---
    net.res_sgen = pd.DataFrame(
        columns=["p_mw", "q_mvar"], index=net.sgen.index
    )
    for g in model.sG:
        net.res_sgen.iloc[g]["p_mw"] = (
            model.psG[g, t].value * model.baseMVA.value
        )

        if "AC" in model.name:
            net.res_sgen.iloc[g]["q_mvar"] = (
                model.qsG[g, t].value * model.baseMVA.value
            )

    if "HC" in model.name:
        y = model.y.get_values()
        net.res_sgen["y_wind"] = None
        for w in model.WIND_HC:
            net.res_sgen["y_wind"].iloc[w] = y[w]

    net.res_sgen.set_index(net.sgen.index, inplace=True)
    net.res_sgen.p_mw.fillna(
        net.sgen.p_mw * net.sgen.scaling * net.sgen.in_service, inplace=True
    )
    if "AC" in model.name:
        net.res_sgen.q_mvar.fillna(
            net.sgen.q_mvar * net.sgen.scaling * net.sgen.in_service,
            inplace=True,
        )


def _gen_results_to_net(net, model, t):
    # --- gen ---
    net.res_gen = pd.DataFrame(
        columns=["p_mw", "q_mvar", "va_degree", "vm_pu"], index=net.gen.index
    )
    for g in model.gG:
        # net.res_gen.loc[g, 'p_mw'] = model.pG[g].value * model.baseMVA.value

        if "AC" in model.name:
            net.res_gen.loc[g, "q_mvar"] = (
                model.qG[g, t].value * model.baseMVA.value
            )

    net.res_gen.va_degree = net.res_bus.va_degree[net.gen.bus].values
    net.res_gen.vm_pu = net.res_bus.vm_pu[net.gen.bus].values


def _trafo_results_to_net(net, model, t):
    net.res_trafo.p_hv_mw = pd.Series(
        [model.pThv[i, t].value for i in net.trafo.index],
        index=net.trafo.index,
    )
    net.res_trafo.p_hv_mw *= model.baseMVA.value
    net.res_trafo.p_lv_mw = pd.Series(
        [model.pTlv[i, t].value for i in net.trafo.index],
        index=net.trafo.index,
    )
    net.res_trafo.p_lv_mw *= model.baseMVA.value
    net.res_trafo.pl_mw = net.res_trafo.p_hv_mw + net.res_trafo.p_lv_mw

    # voltage
    net.res_trafo.vm_hv_pu = net.res_bus.vm_pu[net.trafo.hv_bus].values
    net.res_trafo.vm_lv_pu = net.res_bus.vm_pu[net.trafo.lv_bus].values

    if "AC" in model.name:
        # transformer
        net.res_trafo.q_hv_mvar = pd.Series(
            [model.qThv[i, t].value for i in net.trafo.index],
            index=net.trafo.index,
        )
        net.res_trafo.q_hv_mvar *= model.baseMVA.value
        net.res_trafo.q_lv_mvar = pd.Series(
            [model.qTlv[i, t].value for i in net.trafo.index],
            index=net.trafo.index,
        )
        net.res_trafo.q_lv_mvar *= model.baseMVA.value
        net.res_trafo.ql_mvar = (
            net.res_trafo.q_hv_mvar + net.res_trafo.q_lv_mvar
        )
    else:
        net.res_trafo.fillna(0.0, inplace=True)

    # current
    net.res_trafo.i_hv_ka = np.sqrt(
        net.res_trafo.p_hv_mw**2 + net.res_trafo.q_hv_mvar.fillna(0) ** 2
    ) / (
        net.res_trafo.vm_hv_pu
        * np.sqrt(3)
        * net.bus.vn_kv[net.trafo.hv_bus].values
    )
    net.res_trafo.i_lv_ka = np.sqrt(
        net.res_trafo.p_lv_mw**2 + net.res_trafo.q_lv_mvar.fillna(0) ** 2
    ) / (
        net.res_trafo.vm_lv_pu
        * np.sqrt(3)
        * net.bus.vn_kv[net.trafo.lv_bus].values
    )

    vns = np.vstack([net.trafo.vn_hv_kv.values, net.trafo.vn_lv_kv.values]).T
    lds_trafo = (
        net.res_trafo[["i_hv_ka", "i_lv_ka"]]
        * vns
        * np.sqrt(3)
        / net.trafo.sn_mva.values[:, np.newaxis]
        * 100.0
    )
    with np.errstate(invalid="ignore"):
        ld_trafo = np.max(lds_trafo, axis=1)
    net.res_trafo.loading_percent = (
        ld_trafo / net.trafo.parallel.values / net.trafo.df.values
    )

    # voltage angle
    net.res_trafo.va_hv_degree = net.res_bus.va_degree[net.trafo.hv_bus].values
    net.res_trafo.va_lv_degree = net.res_bus.va_degree[net.trafo.lv_bus].values

    if hasattr(model, "Tap_pos"):
        net.res_trafo["tap"] = pd.Series(
            [model.Tap[i, t].value for i in net.trafo.index],
            index=net.trafo.index,
        )
        net.res_trafo["tap_pos"] = pd.Series(
            [model.Tap_pos[i, t].value for i in net.trafo.index],
            index=net.trafo.index,
        )
    elif hasattr(model, "Tap_linear_constr"):
        net.res_trafo["tap"] = pd.Series(
            [model.Tap[i, t].value for i in net.trafo.index],
            index=net.trafo.index,
        )

    net.res_trafo.set_index(net.trafo.index, inplace=True)


def _shunt_results_to_net(net, model, t):
    for s in model.SHUNT:
        net.res_shunt.loc[s, "p_mw"] = (
            model.GB[s, t]
            * net.res_bus.vm_pu[net.shunt["bus"][s]] ** 2
            * model.baseMVA.value
        )

    if "AC" in model.name:
        net.res_shunt.vm_pu = pd.Series(
            net.res_bus.vm_pu[net.shunt.bus].values, net.shunt.index
        )
        for s in model.SHUNT:
            net.res_shunt.loc[s, "q_mvar"] = (
                -model.BB[s, t]
                * net.res_bus.vm_pu[net.shunt["bus"][s]] ** 2
                * model.baseMVA.value
            )
