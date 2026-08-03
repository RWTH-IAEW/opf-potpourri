"""Post-processing: reads Pyomo solution variables and writes results to
net.res_* DataFrames."""

import numpy as np
import pandas as pd

# pandapower 3.5 dropped clear_result_tables from the top-level namespace;
# pandapower.toolbox carries it on every version we support.
from pandapower.toolbox import clear_result_tables


def _is_ac(model):
    """Return True if model has voltage-magnitude variable v (AC model)."""
    return hasattr(model, "v")


def _is_hc(model):
    """Return True if model has binary placement variable y (HC model)."""
    return hasattr(model, "y")


def pyo_sol_to_net_res(net, model):
    """Write Pyomo solution to net.res_* DataFrames.

    Args:
        net: pandapower network whose res_* tables will be populated.
        model: solved Pyomo ConcreteModel.
    """
    if _is_hc(model):
        for w in model.WIND_HC:
            net.sgen.p_mw[w] = (
                model.psG[w].value * model.baseMVA.value * model.y[w].value
            )
            net.sgen.q_mvar[w] = (
                model.qsG[w].value * model.baseMVA.value * model.y[w].value
            )

    clear_result_tables(net)

    _bus_voltage_results_to_net(net, model)
    _line_results_to_net(net, model)
    _generation_results_to_net(net, model)
    _sgen_results_to_net(net, model)
    _load_results_to_net(net, model)
    _trafo_results_to_net(net, model)
    _shunt_results_to_net(net, model)
    bus_pq = _get_bus_power_results(net)
    net.res_bus.p_mw = bus_pq[:, 0]
    if _is_ac(model):
        net.res_bus.q_mvar = bus_pq[:, 1]


def _bus_voltage_results_to_net(net, model):
    bus_lookup = net._pd2ppc_lookups["bus"]
    bus_idx = bus_lookup[net.bus.index.values]

    if not _is_ac(model):
        net.res_bus.vm_pu = pd.Series(
            [1.0] * len(net.bus.index), net.bus.index
        )
        net.res_bus.vm_pu[net.gen.bus] = net.gen.vm_pu
        net.res_bus.vm_pu[net.ext_grid.bus] = net.ext_grid.vm_pu
    else:
        v = model.v.get_values()
        v_res = [v[b] for b in bus_idx]
        net.res_bus.vm_pu = pd.Series(v_res, index=net.bus.index)

    va = model.delta.get_values()
    va_res = [va[b] for b in bus_idx]
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

    bus_pq[b_i, 0] = dfgr["p"].values
    bus_pq[b_i, 1] = dfgr["q"].values

    return bus_pq


def _line_results_to_net(net, model):
    # model.L includes both net.line rows (indices < n_line) and any
    # net.impedance rows (indices >= n_line). pyo_to_net writes line flows
    # back to net.res_line for the native lines only; impedance results are
    # written to net.res_impedance.
    n_line = len(net.line.index)

    net.res_line.vm_from_pu = pd.Series(
        net.res_bus.vm_pu[net.line.from_bus].values, index=net.line.index
    )
    net.res_line.vm_to_pu = net.res_bus.vm_pu[net.line.to_bus].values

    net.res_line.va_from_degree = net.res_bus.va_degree[
        net.line.from_bus
    ].values
    net.res_line.va_to_degree = net.res_bus.va_degree[net.line.to_bus].values

    pL_from = model.pLfrom.get_values()
    pL_to = model.pLto.get_values()
    base = model.baseMVA.value

    # Index by net.line.index so out-of-service lines (absent from model.L)
    # get filled with 0.
    def _series_for_lines(d):
        return pd.Series(
            [d.get(int(idx), 0.0) * base for idx in net.line.index],
            index=net.line.index,
        )

    net.res_line.p_from_mw = _series_for_lines(pL_from)
    net.res_line.p_to_mw = _series_for_lines(pL_to)
    net.res_line.pl_mw = net.res_line.p_from_mw + net.res_line.p_to_mw

    if _is_ac(model):
        qL_from = model.qLfrom.get_values()
        qL_to = model.qLto.get_values()
        net.res_line.q_from_mvar = _series_for_lines(qL_from)
        net.res_line.q_to_mvar = _series_for_lines(qL_to)
        net.res_line.ql_mvar = (
            net.res_line.q_from_mvar + net.res_line.q_to_mvar
        )
    else:
        # Assignment rather than a chained `.fillna(..., inplace=True)`, which
        # mutates the temporary the column access returns. pandas warns today
        # and makes it a silent no-op in 3.0, leaving NaN for the current and
        # loading calculations below.
        for column in ("q_from_mvar", "q_to_mvar", "ql_mvar"):
            net.res_line[column] = net.res_line[column].fillna(0.0)

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

    # Write impedance branch flows back to net.res_impedance if present.
    imp = net.get("impedance")
    if imp is not None and not imp.empty:
        # impedance synthetic indices in model.L start at n_line
        def _imp_series(d):
            return [
                d.get(int(n_line + i), 0.0) * base for i in range(len(imp))
            ]

        rows_p_from = _imp_series(pL_from)
        rows_p_to = _imp_series(pL_to)
        if _is_ac(model):
            qL_from = model.qLfrom.get_values()
            qL_to = model.qLto.get_values()
            rows_q_from = _imp_series(qL_from)
            rows_q_to = _imp_series(qL_to)
        else:
            rows_q_from = [0.0] * len(imp)
            rows_q_to = [0.0] * len(imp)
        net.res_impedance = pd.DataFrame(
            {
                "p_from_mw": rows_p_from,
                "q_from_mvar": rows_q_from,
                "p_to_mw": rows_p_to,
                "q_to_mvar": rows_q_to,
            },
            index=imp.index,
        )


def _generation_results_to_net(net, model):
    pg = model.pG.get_values()
    qg = model.qG.get_values() if _is_ac(model) else None
    base = model.baseMVA.value
    G_keys = list(model.G)

    for gen, ord in net._gen_order.items():
        # Use the actual element indices that fall in this [f, t) slice,
        # rather than assuming pg/qg dicts have contiguous integer keys.
        keys = [k for k in G_keys if ord[0] <= k < ord[1]]
        keys.sort()
        net["res_" + gen].p_mw = [pg[k] * base for k in keys]
        if qg is not None:
            net["res_" + gen].q_mvar = [qg[k] * base for k in keys]

        net["res_" + gen].set_index(net[gen].index, inplace=True)

    net.res_gen.va_degree = net.res_bus.va_degree[net.gen.bus].values
    net.res_gen.vm_pu = net.res_bus.vm_pu[net.gen.bus].values


def _load_results_to_net(net, model):
    net.res_load = pd.DataFrame(
        columns=["p_mw", "q_mvar"], index=net.load.index, dtype=float
    )
    net.res_load.p_mw = model.pD.get_values()
    net.res_load.p_mw *= model.baseMVA.value
    if _is_ac(model):
        net.res_load.q_mvar = model.qD.get_values()
        net.res_load.q_mvar *= model.baseMVA.value
        net.res_load["q_mvar"] = net.res_load["q_mvar"].fillna(
            net.load["q_mvar"] * net.load["scaling"] * net.load["in_service"]
        )

    net.res_load.set_index(net.load.index, inplace=True)
    net.res_load["p_mw"] = net.res_load["p_mw"].fillna(
        net.load["p_mw"] * net.load["scaling"] * net.load["in_service"]
    )


def _sgen_results_to_net(net, model):
    net.res_sgen = pd.DataFrame(
        columns=["p_mw", "q_mvar"], index=net.sgen.index, dtype=float
    )
    for g in model.sG:
        net.res_sgen.loc[g, "p_mw"] = model.psG[g].value * model.baseMVA.value

        if _is_ac(model):
            net.res_sgen.loc[g, "q_mvar"] = (
                model.qsG[g].value * model.baseMVA.value
            )

    if _is_hc(model):
        y = model.y.get_values()
        net.res_sgen["y_wind"] = None
        for w in model.WIND_HC:
            net.res_sgen.loc[w, "y_wind"] = y[w]

    net.res_sgen.set_index(net.sgen.index, inplace=True)
    sgen_p_fallback = (
        net.sgen["p_mw"].astype(float)
        * net.sgen["scaling"].astype(float)
        * net.sgen["in_service"].astype(float)
    )
    net.res_sgen["p_mw"] = (
        net.res_sgen["p_mw"].astype(float).fillna(sgen_p_fallback)
    )

    if _is_ac(model):
        sgen_q_fallback = (
            net.sgen["q_mvar"].astype(float)
            * net.sgen["scaling"].astype(float)
            * net.sgen["in_service"].astype(float)
        )
        net.res_sgen["q_mvar"] = (
            net.res_sgen["q_mvar"].astype(float).fillna(sgen_q_fallback)
        )


def _trafo_results_to_net(net, model):
    net.res_trafo.p_hv_mw = model.pThv.get_values()
    net.res_trafo.p_hv_mw *= model.baseMVA.value
    net.res_trafo.p_lv_mw = model.pTlv.get_values()
    net.res_trafo.p_lv_mw *= model.baseMVA.value
    net.res_trafo.pl_mw = net.res_trafo.p_hv_mw + net.res_trafo.p_lv_mw

    net.res_trafo.vm_hv_pu = net.res_bus.vm_pu[net.trafo.hv_bus].values
    net.res_trafo.vm_lv_pu = net.res_bus.vm_pu[net.trafo.lv_bus].values

    if _is_ac(model):
        net.res_trafo.q_hv_mvar = model.qThv.get_values()
        net.res_trafo.q_hv_mvar *= model.baseMVA.value
        net.res_trafo.q_lv_mvar = model.qTlv.get_values()
        net.res_trafo.q_lv_mvar *= model.baseMVA.value
        net.res_trafo.ql_mvar = (
            net.res_trafo.q_hv_mvar + net.res_trafo.q_lv_mvar
        )
    else:
        net.res_trafo.fillna(0.0, inplace=True)

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

    net.res_trafo.va_hv_degree = net.res_bus.va_degree[net.trafo.hv_bus].values
    net.res_trafo.va_lv_degree = net.res_bus.va_degree[net.trafo.lv_bus].values

    if hasattr(model, "Tap_pos"):
        net.res_trafo["tap"] = model.Tap.get_values()
        net.res_trafo["tap_pos"] = model.Tap_pos.get_values()
    elif hasattr(model, "Tap_linear_constr"):
        net.res_trafo["tap"] = model.Tap.get_values()

    net.res_trafo.set_index(net.trafo.index, inplace=True)


def _shunt_results_to_net(net, model):
    for s in model.SHUNT:
        net.res_shunt.loc[s, "p_mw"] = (
            model.GB[s]
            * net.res_bus.vm_pu[net.shunt["bus"][s]] ** 2
            * model.baseMVA.value
        )

    if _is_ac(model):
        net.res_shunt.vm_pu = pd.Series(
            net.res_bus.vm_pu[net.shunt.bus].values, net.shunt.index
        )
        for s in model.SHUNT:
            net.res_shunt.loc[s, "q_mvar"] = (
                -model.BB[s]
                * net.res_bus.vm_pu[net.shunt["bus"][s]] ** 2
                * model.baseMVA.value
            )
