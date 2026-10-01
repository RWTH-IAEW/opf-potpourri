# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Opt-in switched-shunt control (`potpourri.models.shunt_control`).

A pandapower shunt is a constant admittance, so its reactive power follows
the voltage squared; the controllable model must keep that. The tests
build a two-bus feeder with a four-step capacitor bank and check, against
`pp.runpp` at every step, that the model with a *fixed* controllable step
reproduces bus voltages and shunt power — including a bank rated at a
different voltage than its bus, where pandapower applies a
`(V_bus / vn_kv)²` factor. Then an under-voltage case must switch the bank
in, the default model must stay untouched, and a multi-period schedule
must respect its movement limits.
"""

import copy

import numpy as np
import pandapower as pp
import pyomo.environ as pyo
import pytest
import simbench as sb

from potpourri.models.ACOPF_base import ACOPF
from potpourri.models.solver_guard import free_integer_variables
from potpourri.models.DCOPF import DCOPF
from potpourri.models.shunt_control import shunt_eligibility
from potpourri.models_multi_period.ACOPF_multi_period import (
    ACOPF_multi_period,
)


# Wide voltage bands keep every bank step admissible in the round trip; the
# grid-code Q(U) envelope warns about such a band, which does not matter
# for networks without a Q-controlled static generator.
pytestmark = pytest.mark.filterwarnings(
    "ignore::potpourri.technologies.q_control.EnvelopeRangeWarning"
)


def _minlp_available():
    try:
        return bool(pyo.SolverFactory("gurobi_direct_minlp").available())
    except Exception:  # noqa: BLE001 - any failure means "not available"
        return False


requires_minlp = pytest.mark.skipif(
    not _minlp_available(), reason="gurobi_direct_minlp not available"
)


def shunt_net(step=1, *, vn_kv=None, q_mvar=-0.5, p_mw=0.01, max_step=4):
    """Two 20 kV buses, a 10 km line, a load and a capacitor bank."""
    net = pp.create_empty_network(sn_mva=10.0)
    b0 = pp.create_bus(net, 20.0)
    b1 = pp.create_bus(net, 20.0)
    pp.create_ext_grid(net, b0, vm_pu=1.0)
    pp.create_line_from_parameters(
        net,
        b0,
        b1,
        length_km=10.0,
        r_ohm_per_km=0.3,
        x_ohm_per_km=0.35,
        c_nf_per_km=10.0,
        max_i_ka=0.4,
    )
    pp.create_load(net, b1, 4.0, 2.0)
    pp.create_shunt(
        net,
        b1,
        q_mvar=q_mvar,
        p_mw=p_mw,
        step=step,
        max_step=max_step,
        vn_kv=vn_kv,
    )
    net.bus["max_vm_pu"] = 1.5
    net.bus["min_vm_pu"] = 0.5
    return net


# ── eligibility and defaults ─────────────────────────────────────────────────


def test_eligibility_rules():
    assert shunt_eligibility(shunt_net()).at[0, "eligible"]
    net = shunt_net()
    net.shunt["step_dependency_table"] = True
    assert "step_dependency_table" in shunt_eligibility(net).at[0, "reason"]
    net = shunt_net(q_mvar=0.0, p_mw=0.0)
    assert "both zero" in shunt_eligibility(net).at[0, "reason"]
    net = shunt_net(step=5)
    assert "outside" in shunt_eligibility(net).at[0, "reason"]
    net = shunt_net()
    net.shunt["in_service"] = False
    assert "out of service" in shunt_eligibility(net).at[0, "reason"]


def test_default_model_has_no_controllable_shunt():
    opf = ACOPF(shunt_net())
    opf.add_OPF()
    assert not hasattr(opf.model, "SHUNT_CTRL")
    assert not hasattr(opf.model, "shunt_step")
    assert free_integer_variables(opf.model) == []


def test_enable_shunt_control_validates_and_is_refused_on_dc():
    opf = ACOPF(shunt_net())
    opf.add_OPF()
    with pytest.raises(ValueError, match="not a shunt"):
        opf.enable_shunt_control(shunts=[3])
    opf.enable_shunt_control()
    with pytest.raises(RuntimeError, match="already"):
        opf.enable_shunt_control()
    assert hasattr(opf.model, "shunt_step_movement_def")
    dc = DCOPF(shunt_net())
    dc.add_OPF()
    with pytest.raises(NotImplementedError):
        dc.enable_shunt_control()


def test_enable_shunt_control_with_no_shunt_explains(four_bus_net):
    opf = ACOPF(four_bus_net)
    opf.add_OPF()
    with pytest.raises(ValueError, match="No shunt is eligible"):
        opf.enable_shunt_control()


# ── pandapower round trip ────────────────────────────────────────────────────


@pytest.mark.integration
@pytest.mark.parametrize("vn_kv", [None, 21.0])
@pytest.mark.parametrize("step", [0, 1, 3, 4])
def test_fixed_step_reproduces_pandapower(vn_kv, step):
    ref = shunt_net(step, vn_kv=vn_kv)
    pp.runpp(ref)
    opf = ACOPF(shunt_net(1, vn_kv=vn_kv))
    opf.add_OPF(free_slack_vm=False)
    opf.enable_shunt_control(mode="discrete")
    opf.model.shunt_step[0].fix(step)
    opf.add_voltage_deviation_objective()
    assert pyo.check_optimal_termination(opf.solve(solver="ipopt"))
    assert opf.net.res_bus.vm_pu.values == pytest.approx(
        ref.res_bus.vm_pu.values, abs=1e-7
    )
    for column in ("p_mw", "q_mvar"):
        assert opf.net.res_shunt.at[0, column] == pytest.approx(
            ref.res_shunt.at[0, column], abs=1e-7
        ), column
    assert opf.net.res_shunt.at[0, "step"] == step
    # the bank's reactive power follows v², not the set point
    expected_q = (
        -0.5
        * step
        * (20.0 / (vn_kv or 20.0)) ** 2
        * ref.res_bus.vm_pu.at[1] ** 2
    )
    assert opf.net.res_shunt.at[0, "q_mvar"] == pytest.approx(
        expected_q, abs=1e-7
    )


# ── optimisation behaviour ───────────────────────────────────────────────────


V_MIN_UNDERVOLTAGE = 0.96


def _undervoltage_net():
    """Slack pinned at 1.0: the load bus sits at 0.950 with the bank off.

    Each step lifts it by about 0.0043 p.u. (0.954, 0.958, 0.963, 0.967),
    so a 0.96 p.u. limit is infeasible without the bank and needs at least
    three steps.
    """
    net = shunt_net(0)
    net.bus["min_vm_pu"] = V_MIN_UNDERVOLTAGE
    net.bus["max_vm_pu"] = 1.05
    return net


@pytest.mark.integration
def test_bank_is_switched_in_to_hold_the_voltage():
    fixed = ACOPF(_undervoltage_net())
    fixed.add_OPF(free_slack_vm=False)
    fixed.add_voltage_deviation_objective()
    assert not pyo.check_optimal_termination(fixed.solve(solver="ipopt"))

    opf = ACOPF(_undervoltage_net())
    opf.add_OPF(free_slack_vm=False)
    opf.enable_shunt_control(mode="discrete")
    opf.add_voltage_deviation_objective()
    results = opf.solve_shunt_round_and_fix(solver="ipopt")
    assert pyo.check_optimal_termination(results)
    step = opf.shunt_schedule()[0]
    assert 3 <= step <= 4
    assert opf.net.res_bus.vm_pu.at[1] >= V_MIN_UNDERVOLTAGE - 1e-6
    assert opf.shunt_operations()[0] == step
    check = opf.apply_shunt_steps(copy.deepcopy(_undervoltage_net()))
    assert check.shunt.at[0, "step"] == step
    pp.runpp(check)
    assert check.res_bus.vm_pu.at[1] == pytest.approx(
        opf.net.res_bus.vm_pu.at[1], abs=1e-6
    )


@pytest.mark.integration
@requires_minlp
def test_global_minlp_switches_the_bank():
    opf = ACOPF(_undervoltage_net())
    opf.add_OPF(free_slack_vm=False)
    opf.enable_shunt_control(mode="discrete")
    opf.add_voltage_deviation_objective()
    results = opf.solve(solver="gurobi_direct_minlp", time_limit=120)
    assert pyo.check_optimal_termination(results)
    assert 3 <= opf.shunt_schedule()[0] <= 4
    assert opf.net.res_bus.vm_pu.at[1] >= V_MIN_UNDERVOLTAGE - 1e-6


@pytest.mark.integration
def test_priced_switching_keeps_the_bank_off():
    net = shunt_net(0)
    net.bus["min_vm_pu"] = 0.9
    free = ACOPF(net)
    free.add_OPF(free_slack_vm=False)
    free.enable_shunt_control(mode="continuous")
    free.add_voltage_deviation_objective()
    assert pyo.check_optimal_termination(free.solve(solver="ipopt"))
    assert free.shunt_schedule()[0] > 0.05
    priced = ACOPF(net)
    priced.add_OPF(free_slack_vm=False)
    priced.enable_shunt_control(mode="continuous")
    priced.add_voltage_deviation_objective()
    priced.penalize_shunt_switching(cost=1.0)
    assert pyo.check_optimal_termination(priced.solve(solver="ipopt"))
    assert priced.shunt_schedule()[0] == pytest.approx(0.0, abs=1e-5)


# ── multi-period ─────────────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def lv_shunt_net():
    """SimBench LV rural1 with a small capacitor bank at its last bus."""
    net = sb.get_simbench_net("1-LV-rural1--0-sw")
    bus = int(net.load.bus.iloc[-1])
    pp.create_shunt(net, bus, q_mvar=-0.02, p_mw=0.0, step=0, max_step=3)
    net.bus["max_vm_pu"] = 1.05
    net.bus["min_vm_pu"] = 0.97
    return net


def test_multi_period_components_are_time_indexed(lv_shunt_net):
    mp = ACOPF_multi_period(lv_shunt_net, toT=3)
    mp.add_OPF()
    mp.enable_shunt_control(
        mode="discrete", max_change_per_step=1, max_operations=2
    )
    m = mp.model
    assert list(m.SHUNT_CTRL) == [0]
    assert len(m.shunt_step) == 3
    assert len(m.shunt_step_change_limit) == 3
    assert len(free_integer_variables(m)) == 3


@pytest.mark.integration
def test_multi_period_rounded_schedule_respects_the_limits(lv_shunt_net):
    mp = ACOPF_multi_period(lv_shunt_net, toT=3)
    mp.add_OPF()
    mp.enable_shunt_control(
        mode="discrete", max_change_per_step=1, max_operations=2
    )
    mp.add_voltage_deviation_objective()
    results = mp.solve_shunt_round_and_fix(
        solver="ipopt", print_solver_output=False
    )
    assert pyo.check_optimal_termination(results)
    schedule = mp.shunt_schedule()
    steps = [0] + schedule[0].tolist()
    assert all(0 <= s <= 3 for s in steps)
    assert all(abs(b - a) <= 1 for a, b in zip(steps[:-1], steps[1:]))
    assert mp.shunt_operations()[0] <= 2
    assert mp.net.res_shunt.at[0, "step"] == schedule[0].iloc[-1]
    assert np.isfinite(mp.net.res_shunt.at[0, "q_mvar"])
