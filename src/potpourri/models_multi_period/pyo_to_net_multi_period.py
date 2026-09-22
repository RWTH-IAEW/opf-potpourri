# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Post-processing for multi-period models.

Reads Pyomo solution at time t into net.res_* DataFrames.
"""

import numpy as np
import pandas as pd

# pandapower 3.5 dropped clear_result_tables from the top-level namespace;
# pandapower.toolbox carries it on every version we support.
from pandapower.toolbox import clear_result_tables


def pyo_sol_to_net_res(net, model, t):
    """Write Pyomo solution for time step t to net.res_* DataFrames.

    Args:
        net: pandapower network whose res_* tables will be populated.
        model: solved Pyomo ConcreteModel.
        t: time step index to extract from multi-period variables.
    """
    if "HC" in model.name:
        # Dispatch is per (candidate, step) and the selection is per
        # candidate, so this writes the requested step's dispatch of the
        # units the model chose to build. The installed rating itself is
        # sqrt(SW2[w]) and has no time index -- a plant is built once.
        for w in model.WIND_HC:
            selected = model.y[w].value
            net.sgen.p_mw[w] = (
                model.psG[w, t].value * model.baseMVA.value * selected
            )
            net.sgen.q_mvar[w] = (
                model.qsG[w, t].value * model.baseMVA.value * selected
            )

    clear_result_tables(net)

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
    """Write bus voltages into `net.res_bus`.

    Magnitudes in p.u. and angles in degrees (the model holds radians). On a DC
    model there are no magnitudes to read, so the generator and external-grid
    set points are written and every other bus is left at 1.0 p.u.

    Args:
        net: The network to write results into.
        model: The solved Pyomo model to read.
        t: Time step to write. `net.res_*` has no time dimension, so one step
            is written at a time.

    Returns:
        None. The result tables are filled in place.
    """
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
    """Aggregate element powers onto their buses.

    Sums load, sgen, gen, ext_grid and shunt results per bus, so `res_bus.p_mw`
    reports the net injection there.

    Args:
        net: The network to write results into.

    Returns:
        None. The result tables are filled in place.
    """
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
    """Write line flows and loading into `net.res_line`.

    Fills both ends, the loss (`p_from + p_to`, which is what is left in the
    branch) and the loading percentage. Rows that the model carries as
    synthetic lines for `net.impedance` are written to `net.res_impedance`
    instead.

    Args:
        net: The network to write results into.
        model: The solved Pyomo model to read.
        t: Time step to write. `net.res_*` has no time dimension, so one step
            is written at a time.

    Returns:
        None. The result tables are filled in place.
    """
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
        # Assignment rather than a chained `.fillna(..., inplace=True)`, which
        # mutates the temporary the column access returns. pandas warns today
        # and makes it a silent no-op in 3.0, leaving NaN for the current and
        # loading calculations below.
        for column in ("q_from_mvar", "q_to_mvar", "ql_mvar"):
            net.res_line[column] = net.res_line[column].fillna(0.0)

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
    """Write the external-grid dispatch into `net.res_ext_grid`.

    Args:
        net: The network to write results into.
        model: The solved Pyomo model to read.
        t: Time step to write. `net.res_*` has no time dimension, so one step
            is written at a time.

    Returns:
        None. The result tables are filled in place.
    """
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
    """Write the external-grid dispatch into `net.res_ext_grid`.

    Args:
        net: The network to write results into.
        model: The solved Pyomo model to read.
        t: Time step to write. `net.res_*` has no time dimension, so one step
            is written at a time.

    Returns:
        None. The result tables are filled in place.
    """
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
    """Write the load dispatch into `net.res_load`.

    A controllable load may differ from its set point; a fixed one reproduces
    it.

    Args:
        net: The network to write results into.
        model: The solved Pyomo model to read.
        t: Time step to write. `net.res_*` has no time dimension, so one step
            is written at a time.

    Returns:
        None. The result tables are filled in place.
    """
    net.res_load = pd.DataFrame(
        columns=["p_mw", "q_mvar"], index=net.load.index, dtype=float
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
        net.res_load["q_mvar"] = net.res_load["q_mvar"].fillna(
            net.load.q_mvar * net.load.scaling * net.load.in_service
        )

    net.res_load.set_index(net.load.index, inplace=True)
    net.res_load["p_mw"] = net.res_load["p_mw"].fillna(
        net.load.p_mw * net.load.scaling * net.load.in_service
    )


def _sgen_results_to_net(net, model, t):
    """Write the static-generator dispatch into `net.res_sgen`.

    Generator sign convention: positive is injection.

    Args:
        net: The network to write results into.
        model: The solved Pyomo model to read.
        t: Time step to write. `net.res_*` has no time dimension, so one step
            is written at a time.

    Returns:
        None. The result tables are filled in place.
    """
    # --- sgen ---
    net.res_sgen = pd.DataFrame(
        columns=["p_mw", "q_mvar"], index=net.sgen.index, dtype=float
    )
    # Collect first, then assign whole columns. `res_sgen.iloc[g]["p_mw"] = x`
    # is chained assignment: it writes into the Series that iloc returns, which
    # is a copy whenever pandas decides to make one. It happens to land today
    # and will stop landing under copy-on-write in pandas 3.0 — silently, with
    # the optimised dispatch replaced by the profile via the fillna below, so a
    # curtailed sgen would report its uncurtailed output.
    positions = list(model.sG)
    p_column = net.res_sgen.columns.get_loc("p_mw")
    for g in positions:
        net.res_sgen.iloc[g, p_column] = (
            model.psG[g, t].value * model.baseMVA.value
        )

    if "AC" in model.name:
        q_column = net.res_sgen.columns.get_loc("q_mvar")
        for g in positions:
            net.res_sgen.iloc[g, q_column] = (
                model.qsG[g, t].value * model.baseMVA.value
            )

    if "HC" in model.name:
        y = model.y.get_values()
        net.res_sgen["y_wind"] = None
        y_column = net.res_sgen.columns.get_loc("y_wind")
        for w in model.WIND_HC:
            net.res_sgen.iloc[w, y_column] = y[w]

    net.res_sgen.set_index(net.sgen.index, inplace=True)
    net.res_sgen["p_mw"] = net.res_sgen["p_mw"].fillna(
        net.sgen.p_mw * net.sgen.scaling * net.sgen.in_service
    )
    if "AC" in model.name:
        net.res_sgen["q_mvar"] = net.res_sgen["q_mvar"].fillna(
            net.sgen.q_mvar * net.sgen.scaling * net.sgen.in_service
        )


def _gen_results_to_net(net, model, t):
    """Write the generator dispatch into `net.res_gen`.

    Args:
        net: The network to write results into.
        model: The solved Pyomo model to read.
        t: Time step to write. `net.res_*` has no time dimension, so one step
            is written at a time.

    Returns:
        None. The result tables are filled in place.
    """
    # --- gen ---
    net.res_gen = pd.DataFrame(
        columns=["p_mw", "q_mvar", "va_degree", "vm_pu"],
        index=net.gen.index,
        dtype=float,
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
    """Write transformer flows and loading into `net.res_trafo`.

    Both windings, plus the loading percentage against the nameplate rating.

    Args:
        net: The network to write results into.
        model: The solved Pyomo model to read.
        t: Time step to write. `net.res_*` has no time dimension, so one step
            is written at a time.

    Returns:
        None. The result tables are filled in place.
    """
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
    """Write shunt consumption into `net.res_shunt`.

    Computed from the solved bus voltage and the shunt's admittance, so it
    follows $v^2$ rather than the set point.

    Args:
        net: The network to write results into.
        model: The solved Pyomo model to read.
        t: Time step to write. `net.res_*` has no time dimension, so one step
            is written at a time.

    Returns:
        None. The result tables are filled in place.
    """
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
