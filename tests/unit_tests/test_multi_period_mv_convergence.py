# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""The multi-period AC OPF on a SimBench MV network.

Regression tests for the warm start that left pandapower's auxiliary ppc
buses at the flat start. Those buses carry the open end of a line with an
open switch; they have no ``net.bus`` row, so a seed read through
``net.res_bus`` never reached them. On ``1-MV-rural--0-no_sw`` six of them
sit at about -148.7 deg behind the 110/20 kV Dyn5 transformers, and the
flat-start angle error alone put a residual of the order of 1e3 p.u. into
the branch-flow equations of short cables, from which IPOPT reported a
locally infeasible point on a feasible model.
"""

import copy
import warnings

import numpy as np
import pandapower as pp
import pyomo.environ as pyo
import pytest
import simbench as sb

from potpourri.models_multi_period.ACOPF_multi_period import (
    ACOPF_multi_period,
)
from potpourri.models_multi_period.init_pyo_from_pp_res_multi_period import (
    _seed_bus_state,
)

warnings.filterwarnings("ignore")

# The sunniest noon of the SimBench year (day 206, 12:00): highest PV feed-in,
# which is where the ±5 % band binds.
T0 = 19824

KVL_FAMILIES = (
    "KVL_real_from",
    "KVL_real_to",
    "KVL_reactive_from",
    "KVL_reactive_to",
    "KVL_real_fromTransf",
    "KVL_real_toTransf",
    "KVL_reactive_fromTransf",
    "KVL_reactive_toTransf",
)


@pytest.fixture(scope="module")
def mv_net():
    net = sb.get_simbench_net("1-MV-rural--0-no_sw")
    net.bus["min_vm_pu"] = 0.95
    net.bus["max_vm_pu"] = 1.05
    net.sgen["controllable"] = True
    net.sgen["min_p_mw"] = 0.0
    return net


def _model(mv_net, steps):
    mp = ACOPF_multi_period(copy.deepcopy(mv_net), fromT=T0, toT=T0 + steps)
    mp.add_OPF()
    mp.add_voltage_deviation_objective()
    return mp


def _auxiliary_buses(mp):
    """The ppc buses of the model that no pandapower bus maps onto."""
    reachable = {int(mp.bus_lookup[b]) for b in mp.net.bus.index}
    return sorted(set(mp.model.B) - reachable)


def _worst_residual(model, families):
    worst = 0.0
    for name in families:
        con = model.component(name)
        for idx in con:
            c = con[idx]
            body = pyo.value(c.body)
            if c.has_ub():
                worst = max(worst, body - pyo.value(c.upper))
            if c.has_lb():
                worst = max(worst, pyo.value(c.lower) - body)
    return worst


def test_the_mv_rural_network_has_auxiliary_buses(mv_net):
    """The premise: six open line switches, six ppc buses without a row."""
    mp = _model(mv_net, steps=1)
    aux = _auxiliary_buses(mp)
    assert len(aux) == 6
    open_line_switches = (mp.net.switch.et == "l") & ~mp.net.switch.closed
    assert int(open_line_switches.sum()) == 6
    # each auxiliary bus terminates exactly one in-service line
    lines_at_aux = [
        l
        for l in mp.model.L
        if mp.model.A[l, 1] in aux or mp.model.A[l, 2] in aux
    ]
    assert len(lines_at_aux) == 6


def test_warm_start_seeds_the_auxiliary_buses(mv_net):
    mp = _model(mv_net, steps=1)
    assert mp.warm_start_from_pf() == 1
    m = mp.model
    for b in _auxiliary_buses(mp):
        # about -148.7 deg, i.e. nowhere near the flat start
        assert abs(pyo.value(m.delta[b, T0])) > 2.0
        assert pyo.value(m.v[b, T0]) != 1.0
    # and the seed satisfies the branch-flow equations everywhere
    assert _worst_residual(m, KVL_FAMILIES) < 1e-4
    assert _worst_residual(m, ("KCL_real", "KCL_reactive")) < 1e-4


def test_seed_falls_back_to_res_bus_when_the_numbering_differs(mv_net):
    """A lookup that disagrees with the power flow must not mis-seed."""
    mp = _model(mv_net, steps=1)
    scratch = copy.deepcopy(mp.net)
    pp.runpp(scratch, voltage_depend_loads=False)
    wrong_lookup = np.array(mp.bus_lookup, copy=True)
    first, second = scratch.bus.index[:2]
    wrong_lookup[first], wrong_lookup[second] = (
        wrong_lookup[second],
        wrong_lookup[first],
    )
    seen = {}

    def record(component, index, value):
        seen[(component, index)] = value

    _seed_bus_state(mp.model, scratch, wrong_lookup, T0, record)
    # every pandapower bus seeded through the (wrong) lookup, no aux bus
    assert len(seen) == 2 * len(scratch.bus)
    for b in _auxiliary_buses(mp):
        assert ("v", (b, T0)) not in seen


def test_mv_rural_two_steps_converge_with_the_default_warm_start(mv_net):
    mp = _model(mv_net, steps=2)
    result = mp.solve(solver="ipopt", print_solver_output=False, to_net=False)
    assert pyo.check_optimal_termination(result)
    m = mp.model
    for b in m.Bpd:
        for t in m.T:
            assert 0.95 - 1e-5 <= pyo.value(m.v[b, t]) <= 1.05 + 1e-5


def test_cold_start_begins_at_the_base_power_flow(mv_net):
    """Without the seed ``v``/``delta`` hold the base power flow, aux buses too."""
    mp = _model(mv_net, steps=2)
    m = mp.model
    for b in _auxiliary_buses(mp):
        for t in m.T:
            assert abs(pyo.value(m.delta[b, t])) > 2.0
            assert pyo.value(m.v[b, t]) == pytest.approx(
                mp.bus_data.v_m[b], abs=1e-12
            )


def test_mv_rural_two_steps_converge_from_the_cold_start(mv_net):
    """``warm_start=False`` on a fresh model must converge too."""
    mp = _model(mv_net, steps=2)
    result = mp.solve(
        solver="ipopt",
        print_solver_output=False,
        to_net=False,
        warm_start=False,
    )
    assert pyo.check_optimal_termination(result)
