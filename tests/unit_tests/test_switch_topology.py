"""Regression tests for switch handling and branch endpoint numbering.

Two defects were found while reproducing the linear AC-OPF study on the
SimBench MV networks (see ``docs/research/campaign_20260912/AUDIT.md``):

1. ``Basemodel`` built the line / transformer endpoint maps from
   pandapower bus indices while ``model.B`` and the load / generator maps
   use ppc bus numbers. Whenever pandapower's bus lookup is not the
   identity (orphan buses are moved to the end of the ppc numbering,
   open switches add auxiliary buses) lines attached to the wrong buses.
2. ``preprocess_grid`` merged closed bus-bus switches but left the
   merged-away bus behind as an orphan and did not move the switches
   sitting on it, which is precisely what made the lookup non-identity
   on ``1-MV-rural--0-sw``.
"""

import copy

import numpy as np
import pandapower as pp
import pyomo.environ as pyo
import pytest

from potpourri.models.ACOPF_base import ACOPF
from potpourri.models.basemodel import preprocess_grid


def _radial_four_bus_with_tie():
    """Ext grid at bus 0; bus 1 joined to bus 0 by a closed bus-bus switch;
    a normally-open tie line 1-3 whose switch sits on bus 1; bus 3 also
    supplied from bus 2."""
    net = pp.create_empty_network(sn_mva=1.0)
    b = [
        pp.create_bus(net, 20.0, max_vm_pu=1.1, min_vm_pu=0.9)
        for _ in range(4)
    ]
    pp.create_ext_grid(net, b[0], vm_pu=1.02)
    pp.create_line(net, b[0], b[2], 1.0, "NAYY 4x150 SE")
    pp.create_line(net, b[1], b[3], 1.0, "NAYY 4x150 SE")  # tie, open
    pp.create_line(net, b[2], b[3], 1.0, "NAYY 4x150 SE")
    pp.create_load(net, b[2], 0.3, 0.1)
    pp.create_load(net, b[3], 0.2, 0.05)
    pp.create_sgen(
        net,
        b[3],
        0.1,
        0.0,
        sn_mva=0.12,
        controllable=True,
        min_p_mw=0.0,
        max_p_mw=0.1,
        min_q_mvar=-0.05,
        max_q_mvar=0.05,
    )
    pp.create_switch(net, b[0], b[1], et="b", closed=True)
    pp.create_switch(net, b[1], 1, et="l", closed=False)
    return net


def _orphan_bus_net():
    """Bus 1 is in service but connected to nothing: pandapower numbers it
    last in the ppc, so pandapower and ppc bus numbers differ for buses
    2 and 3."""
    net = pp.create_empty_network(sn_mva=1.0)
    b = [
        pp.create_bus(net, 20.0, max_vm_pu=1.1, min_vm_pu=0.9)
        for _ in range(4)
    ]
    pp.create_ext_grid(net, b[0], vm_pu=1.02)
    pp.create_line(net, b[0], b[2], 1.0, "NAYY 4x150 SE")
    pp.create_line(net, b[0], b[3], 1.0, "NAYY 4x150 SE")
    pp.create_load(net, b[2], 0.3, 0.1)
    pp.create_load(net, b[3], 0.2, 0.05)
    return net


def test_preprocess_grid_drops_merged_bus_and_moves_its_switches():
    net = _radial_four_bus_with_tie()
    pre = preprocess_grid(copy.deepcopy(net))
    # merged-away bus gone, no bus-bus switch left, tie switch follows merge
    assert len(pre.bus) == len(net.bus) - 1
    assert not (pre.switch["et"] == "b").any()
    tie = pre.switch[pre.switch["et"] == "l"].iloc[0]
    line = pre.line.loc[int(tie["element"])]
    assert int(tie["bus"]) in (int(line["from_bus"]), int(line["to_bus"]))
    assert not bool(tie["closed"])
    # every element still references an existing bus
    for name, cols in (
        ("line", ("from_bus", "to_bus")),
        ("load", ("bus",)),
        ("sgen", ("bus",)),
        ("ext_grid", ("bus",)),
        ("switch", ("bus",)),
    ):
        for col in cols:
            assert set(pre[name][col].astype(int)) <= set(pre.bus.index)
    # topology unchanged: the tie stays open, power flow identical
    pp.runpp(net)
    pp.runpp(pre)
    assert abs(pre.res_line.p_from_mw.iloc[1]) < 1e-9
    np.testing.assert_allclose(
        sorted(pre.res_bus.vm_pu.values),
        sorted(net.res_bus.vm_pu.drop(1).values),
        atol=1e-9,
    )


def test_preprocess_grid_keeps_tightest_voltage_band_of_merged_pair():
    net = _radial_four_bus_with_tie()
    net.bus.at[1, "max_vm_pu"] = 1.05
    net.bus.at[1, "min_vm_pu"] = 0.95
    pre = preprocess_grid(copy.deepcopy(net))
    assert pre.bus.at[0, "max_vm_pu"] == pytest.approx(1.05)
    assert pre.bus.at[0, "min_vm_pu"] == pytest.approx(0.95)


def test_preprocess_grid_resolves_chained_bus_bus_switches():
    net = _orphan_bus_net()
    pp.create_switch(net, 0, 1, et="b", closed=True)
    pp.create_switch(net, 1, 2, et="b", closed=True)  # chain 0-1-2
    pre = preprocess_grid(copy.deepcopy(net))
    assert len(pre.bus) == 2
    assert set(pre.load.bus.astype(int)) <= set(pre.bus.index)
    pp.runpp(pre)
    assert pre.converged


def _assert_endpoints_match_ppc(opf):
    br = opf.net._ppc["branch"]
    for pos, l in enumerate(opf.line_data.index):
        if l not in opf.model.L:
            continue
        row = opf._line_rows_ppc[pos]
        assert opf.model.A[l, 1] == int(br[row, 0].real)
        assert opf.model.A[l, 2] == int(br[row, 1].real)
    for pos, t in enumerate(opf.trafo_data.index):
        if t not in opf.model.TRANSF:
            continue
        row = opf._trafo_rows_ppc[pos]
        assert opf.model.AT[t, 1] == int(br[row, 0].real)
        assert opf.model.AT[t, 2] == int(br[row, 1].real)
    for l in opf.model.L:
        assert opf.model.A[l, 1] in opf.model.B
        assert opf.model.A[l, 2] in opf.model.B


def test_branch_endpoints_follow_ppc_numbering_with_orphan_bus():
    net = _orphan_bus_net()
    opf = ACOPF(net)
    lookup = opf.net._pd2ppc_lookups["bus"][opf.net.bus.index.values]
    # the point of the test: pandapower's numbering differs from net.bus
    assert not np.array_equal(lookup, np.arange(len(opf.net.bus)))
    _assert_endpoints_match_ppc(opf)
    # every load sits on a bus that has a line attached in the model
    line_buses = {opf.model.A[l, k] for l in opf.model.L for k in (1, 2)}
    for b, _d in opf.model.Dbs:
        assert b in line_buses


def test_open_line_switch_is_respected_by_the_model():
    net = _radial_four_bus_with_tie()
    opf = ACOPF(net)
    _assert_endpoints_match_ppc(opf)
    # pandapower routes the open tie to an auxiliary bus: exactly one
    # endpoint of the tie is outside the pandapower-mapped buses
    tie = 1
    ends = {opf.model.A[tie, 1], opf.model.A[tie, 2]}
    assert len(ends - set(opf.model.Bpd)) == 1


def test_simple_four_bus_impedance_endpoints_still_pp_indices_when_identity():
    """On grids where pandapower and ppc numbering coincide the endpoint
    map is unchanged (guards the impedance-branch inclusion tests)."""
    net = pp.networks.simple_four_bus_system()
    b_new = pp.create_bus(net, net.bus.vn_kv.iloc[-1])
    pp.create_impedance(net, net.bus.index[-2], b_new, 0.01, 0.02, 1.0)
    pp.create_load(net, b_new, 0.01, 0.0)
    opf = ACOPF(net)
    n_line = len(net.line)
    assert opf.bus_line_dict[(n_line, 1)] == net.bus.index[-2]
    assert opf.bus_line_dict[(n_line, 2)] == b_new
    _assert_endpoints_match_ppc(opf)


@pytest.mark.slow
def test_acopf_on_merged_switch_network_matches_power_flow_at_fixed_dispatch():
    """With sgen dispatch pinned, the AC-OPF equations must reproduce the
    pandapower power flow of the switch-processed network."""
    net = _radial_four_bus_with_tie()
    net.sgen["controllable"] = False
    pp.runpp(net)
    opf = ACOPF(net)
    opf.add_OPF(thermal_limit="mva", free_slack_vm=False)
    opf.model.obj = pyo.Objective(expr=0.0)
    res = opf.solve(solver="ipopt", to_net=True, print_solver_output=False)
    assert pyo.check_optimal_termination(res)
    vm_pf = net.res_bus.vm_pu.drop(1).sort_index().values
    vm_opf = opf.net.res_bus.vm_pu.sort_index().values
    np.testing.assert_allclose(vm_opf, vm_pf, atol=2e-5)
