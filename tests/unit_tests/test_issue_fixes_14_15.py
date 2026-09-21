# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Regression tests for GitLab issues #14 and #15.

* #14 — Battery, Heatpump and PV were not coupled to the power balance
* #15 — sgen power parameters were in MW while demand parameters were per-unit

Both are properties of the *formulation*, so most of these tests assert on
the constraint bodies and on a power-flow reference rather than on a solve.
"""

from __future__ import annotations

import copy
import warnings

import numpy as np
import pandapower as pp
import pyomo.environ as pyo
import pytest
import simbench as sb
from pyomo.core.expr import identify_variables

from potpourri.models_multi_period.ACOPF_multi_period import (
    ACOPF_multi_period,
)
from potpourri.models_multi_period.DCOPF_multi_period import (
    DCOPF_multi_period,
)
from potpourri.technologies.battery import Battery_multi_period
from potpourri.technologies.heat_pump import Heatpump_multi_period
from potpourri.technologies.pv import PV_multi_period

warnings.filterwarnings("ignore")

# Daylight, so the sgen and PV profiles are non-zero and a scaling error in
# either one shows up in the balance.
DAY_FROM, DAY_TO = 13860, 13872
PEAK = 13868

# Per-battery dispatch for the direction tests. Large enough that the voltage
# shift clears solver tolerance, small enough that 14 units in parallel do not
# drive the snapshot out of its (deliberately wide) voltage band.
DISPATCH_PU = 0.004


@pytest.fixture(scope="module")
def lv_net():
    return sb.get_simbench_net("1-LV-rural1--0-sw")


def _prepared(net, sn_mva=None):
    net = copy.deepcopy(net)
    if sn_mva is not None:
        net.sn_mva = sn_mva
    net.bus["max_vm_pu"] = 1.10
    net.bus["min_vm_pu"] = 0.90
    net.sgen["controllable"] = True
    net.sgen["max_p_mw"] = net.sgen["p_mw"]
    net.sgen["min_p_mw"] = 0.0
    net.ext_grid["max_q_mvar"] = 1e3
    net.ext_grid["min_q_mvar"] = -1e3
    return net


def _vars_in(constraint_data):
    if constraint_data.body is None:
        return set()
    return {
        v.name.split("[")[0] for v in identify_variables(constraint_data.body)
    }


def _constraints_referencing(model, names):
    """Names of active constraints whose body uses any of *names*."""
    hits = set()
    for con in model.component_objects(pyo.Constraint, active=True):
        for idx in con:
            if _vars_in(con[idx]) & names:
                hits.add(con.name)
                break
    return hits


# ---------------------------------------------------------------------------
# #15 — one per-unit convention across the whole multi-period model
# ---------------------------------------------------------------------------


def test_issue15_sgen_params_are_per_unit(lv_net):
    """``PsG`` / ``sPGmax`` were taken straight from the MW profile.

    Every other power quantity in the model is per-unit, so on a network with
    ``sn_mva != 1`` the balance mixed scales.
    """
    base = 10.0
    opf = ACOPF_multi_period(
        _prepared(lv_net, sn_mva=base), toT=DAY_TO, fromT=DAY_FROM
    )
    opf.add_OPF()
    profile_mw = opf.net.profiles[("sgen", "p_mw")][0].loc[PEAK]
    assert profile_mw > 0, "picked a time step with no generation"

    assert pyo.value(opf.model.PsG[0, PEAK]) == pytest.approx(
        profile_mw / base
    )
    assert pyo.value(opf.model.sPGmax[0, PEAK]) == pytest.approx(
        profile_mw / base
    )


def test_issue15_sgen_min_p_is_per_unit(lv_net):
    net = _prepared(lv_net, sn_mva=10.0)
    net.sgen["min_p_mw"] = 0.005
    opf = ACOPF_multi_period(net, toT=DAY_TO, fromT=DAY_FROM)
    opf.add_OPF()
    for g in opf.model.sGc:
        assert pyo.value(opf.model.sPGmin[g, PEAK]) == pytest.approx(
            0.005 / 10.0
        )


def test_issue15_sgen_and_load_params_share_a_scale(lv_net):
    """The two sides of the real-power balance must use the same base."""
    base = 10.0
    opf = ACOPF_multi_period(
        _prepared(lv_net, sn_mva=base), toT=DAY_TO, fromT=DAY_FROM
    )
    opf.add_OPF()
    sgen_mw = opf.net.profiles[("sgen", "p_mw")][0].loc[PEAK]
    load_col = opf.net.profiles[("load", "p_mw")].columns[0]
    load_mw = opf.net.profiles[("load", "p_mw")][load_col].loc[PEAK]

    sgen_ratio = pyo.value(opf.model.PsG[0, PEAK]) / sgen_mw
    load_ratio = pyo.value(opf.model.PD[0, PEAK]) / load_mw
    assert sgen_ratio == pytest.approx(load_ratio)
    assert sgen_ratio == pytest.approx(1.0 / base)


def test_issue15_pv_profile_is_per_unit(lv_net):
    base = 10.0
    opf = ACOPF_multi_period(
        _prepared(lv_net, sn_mva=base), toT=DAY_TO, fromT=DAY_FROM
    )
    pv = PV_multi_period(opf.net, T=DAY_TO - DAY_FROM, penetration=50.0)
    raw_mw = opf.net.pv_load_profiles["PV5"].abs().max()
    assert pv.pv_load_profile.abs().max() == pytest.approx(raw_mw / base)


def test_issue15_pv_p_inst_is_not_scaled_twice(lv_net):
    """``pv_p_inst`` divided an already-per-unit profile by the base again."""
    base = 10.0
    opf = ACOPF_multi_period(
        _prepared(lv_net, sn_mva=base), toT=DAY_TO, fromT=DAY_FROM
    )
    pv = PV_multi_period(
        opf.net, T=DAY_TO - DAY_FROM, penetration=50.0, q_control="qp"
    )
    assert pv.pv_p_inst == pytest.approx(pv.pv_load_profile.abs().max())


def test_issue15_heatpump_qloss_is_per_unit(lv_net):
    """``heat_load`` scaled a per-unit factor onto an MW profile."""
    base = 10.0
    net = _prepared(lv_net, sn_mva=base)
    opf = ACOPF_multi_period(net, toT=DAY_TO, fromT=DAY_FROM)
    hp = Heatpump_multi_period(
        opf.net, T=DAY_TO - DAY_FROM, penetration=50.0, qloss_max=0.01
    )
    # The peak of the scaled profile is exactly qloss_max by construction.
    assert hp.heat_load.max() == pytest.approx(0.01)


def test_issue15_solution_does_not_depend_on_the_system_base(lv_net):
    """A pinned model must reproduce pandapower's power flow at any base.

    With sgens non-controllable and the slack pinned, the OPF has nothing to
    decide: its answer is the power flow. Before the fix, sn_mva = 10 drove
    this model to a locally infeasible point roughly 6e-2 p.u. away.
    """
    results = {}
    for base in (1.0, 10.0):
        ref = copy.deepcopy(lv_net)
        ref.sn_mva = base
        prof = sb.get_absolute_values(
            ref, profiles_instead_of_study_cases=True
        )
        ref.load["p_mw"] = prof[("load", "p_mw")].iloc[PEAK].values
        ref.load["q_mvar"] = prof[("load", "q_mvar")].iloc[PEAK].values
        ref.sgen["p_mw"] = prof[("sgen", "p_mw")].iloc[PEAK].values
        ref.sgen["q_mvar"] = 0.0  # pf=1 in the model
        pp.runpp(ref, voltage_depend_loads=False)

        net = copy.deepcopy(lv_net)
        net.sn_mva = base
        net.bus["max_vm_pu"] = 1.10
        net.bus["min_vm_pu"] = 0.90
        net.sgen["controllable"] = False
        net.ext_grid["max_q_mvar"] = 1e3
        net.ext_grid["min_q_mvar"] = -1e3
        opf = ACOPF_multi_period(net, toT=PEAK + 1, fromT=PEAK)
        opf.add_OPF(free_slack_vm=False)
        opf.add_voltage_deviation_objective()
        res = opf.solve(solver="ipopt", print_solver_output=False)
        assert pyo.check_optimal_termination(res), f"sn_mva={base}"

        opf.map_to_net(PEAK)
        got = opf.net.res_bus.vm_pu.to_numpy()
        expected = ref.res_bus.vm_pu.to_numpy()
        n = min(len(got), len(expected))
        assert np.max(np.abs(got[:n] - expected[:n])) < 1e-4, (
            f"sn_mva={base} does not reproduce the power flow"
        )
        results[base] = got[:n]

    # Same physical answer whatever the base.
    assert np.max(np.abs(results[1.0] - results[10.0])) < 1e-6


# ---------------------------------------------------------------------------
# #14 — devices reach the nodal balance
# ---------------------------------------------------------------------------


def _with_battery(net, toT=4, fromT=0, **kw):
    opf = ACOPF_multi_period(net, toT=toT, fromT=fromT)
    bat = Battery_multi_period(
        opf.net,
        T=toT - fromT,
        penetration=50.0,
        power_pu=0.01,
        capacity_pu_h=0.03,
        **kw,
    )
    bat.get_all(opf.model)
    opf.add_OPF()
    return opf, bat


def test_issue14_battery_reaches_both_balances(lv_net):
    """``BAT_P`` used to appear only in the battery's own constraints."""
    opf, _ = _with_battery(_prepared(lv_net))
    hits = _constraints_referencing(
        opf.model, {"BAT_Pchg", "BAT_Pdis", "BAT_Q"}
    )
    assert "KCL_real" in hits
    assert "KCL_reactive" in hits


def test_issue14_heatpump_reaches_the_real_balance(lv_net):
    opf = ACOPF_multi_period(_prepared(lv_net), toT=4)
    hp = Heatpump_multi_period(opf.net, T=4, penetration=50.0)
    hp.get_all(opf.model)
    opf.add_OPF()
    assert "KCL_real" in _constraints_referencing(opf.model, {"hp_p"})


def test_issue14_pv_reaches_both_balances(lv_net):
    opf = ACOPF_multi_period(_prepared(lv_net), toT=4)
    pv = PV_multi_period(opf.net, T=4, penetration=50.0, q_control="both")
    pv.get_all(opf.model)
    opf.add_OPF()
    hits = _constraints_referencing(opf.model, {"pPV", "qPV"})
    assert "KCL_real" in hits
    assert "KCL_reactive" in hits


def test_issue14_battery_couples_in_the_dc_model(lv_net):
    """The DC formulation had no flexibility hook at all."""
    opf = DCOPF_multi_period(_prepared(lv_net), toT=4)
    bat = Battery_multi_period(
        opf.net, T=4, penetration=50.0, power_pu=0.01, capacity_pu_h=0.03
    )
    bat.get_all(opf.model)
    opf.add_OPF()
    assert "KCL_def" in _constraints_referencing(
        opf.model, {"BAT_Pchg", "BAT_Pdis"}
    )
    # No reactive balance in DC, so no reactive variable either.
    assert not hasattr(opf.model, "BAT_Q")


def test_issue14_double_registration_is_refused(lv_net):
    """Coupling the same device twice would silently double its power."""
    opf, bat = _with_battery(_prepared(lv_net))
    with pytest.raises(RuntimeError, match=r"already coupled"):
        bat.couple_to_power_balance(opf.model)


def test_issue14_rebuild_is_a_noop_without_devices(lv_net):
    """A model with no device must be unchanged by the rebuild."""
    opf = ACOPF_multi_period(_prepared(lv_net), toT=3)
    before = {
        idx: str(opf.model.KCL_real[idx].body) for idx in opf.model.KCL_real
    }
    opf.rebuild_kcl()
    after = {
        idx: str(opf.model.KCL_real[idx].body) for idx in opf.model.KCL_real
    }
    assert before == after
    assert opf.KCL_flexibility(opf.model, list(opf.model.B)[0], 0) == 0


def test_issue14_device_bus_is_mapped_through_the_bus_lookup(lv_net):
    """Placement sets hold pandapower buses; the balance uses ppc numbering."""
    opf, bat = _with_battery(_prepared(lv_net))
    for device, pd_bus in list(opf.model.BAT_bus):
        ppc_bus = int(bat.bus_lookup[int(pd_bus)])
        # The battery's variables must appear in the balance at its ppc bus.
        used = _vars_in(opf.model.KCL_real[ppc_bus, 0])
        assert {"BAT_Pchg", "BAT_Pdis"} & used, (
            f"battery {device} missing from the balance at ppc bus {ppc_bus}"
        )


def _pinned_snapshot(lv_net, dispatch=None):
    """Solve a single step with every degree of freedom removed.

    Sgens are non-controllable and the slack magnitude is pinned, so the model
    has a unique solution and comparisons do not depend on which local optimum
    IPOPT happens to find — the AC OPF is nonconvex, and an objective
    comparison between two runs is not a reliable probe of coupling.

    Args:
        lv_net: The low-voltage pandapower network to solve.
        dispatch: ``None`` for no battery, else ``(p_chg, p_dis, q)`` bounds to
            hold the battery at. Bounds rather than ``fix()``, so the
            capability constraints keep at least one free variable and the NL
            writer does not see a constraint with no variables at all.

    Returns:
        ``(voltages_by_ppc_bus, device_ppc_buses)``. The bus list excludes the
        slack: ``free_slack_vm=False`` pins its magnitude, so it cannot move
        between the two runs and asserting that it does would fail by
        construction rather than catch anything.
    """
    net = copy.deepcopy(lv_net)
    # A wide band deliberately: the voltage limits are only here to make the
    # model well posed, and the property under test is the *direction* in which
    # a dispatch moves the profile. With a tight band, an aggregate dispatch
    # large enough to measure pushes the snapshot outside it and the solve
    # fails for reasons that have nothing to do with the coupling.
    net.bus["max_vm_pu"] = 1.20
    net.bus["min_vm_pu"] = 0.80
    net.sgen["controllable"] = False
    net.ext_grid["max_q_mvar"] = 1e3
    net.ext_grid["min_q_mvar"] = -1e3

    opf = ACOPF_multi_period(net, toT=PEAK + 1, fromT=PEAK)
    bat = None
    if dispatch is not None:
        bat = Battery_multi_period(
            opf.net,
            T=1,
            penetration=100.0,
            power_pu=DISPATCH_PU,
            capacity_pu_h=1.0,
            soc_min=0.0,
            terminal_soc=None,
        )
        bat.get_all(opf.model)

    opf.add_OPF(free_slack_vm=False)
    opf.add_voltage_deviation_objective()

    if dispatch is not None:
        p_chg, p_dis, q = dispatch
        for b in opf.model.BAT:
            for t in opf.model.T:
                opf.model.BAT_Pchg[b, t].setlb(p_chg)
                opf.model.BAT_Pchg[b, t].setub(p_chg)
                opf.model.BAT_Pdis[b, t].setlb(p_dis)
                opf.model.BAT_Pdis[b, t].setub(p_dis)
                opf.model.BAT_Q[b, t].setlb(q)
                opf.model.BAT_Q[b, t].setub(q)

    res = opf.solve(solver="ipopt", print_solver_output=False)
    assert pyo.check_optimal_termination(res), f"dispatch={dispatch}"

    voltages = {b: pyo.value(opf.model.v[b, PEAK]) for b in opf.model.B}
    slack = set(opf.model.b0)
    device_buses = sorted(
        {
            int(bat.bus_lookup[int(p)])
            for _, p in list(opf.model.BAT_bus)
            if int(bat.bus_lookup[int(p)]) not in slack
        }
        if bat is not None
        else []
    )
    return voltages, device_buses


@pytest.mark.integration
def test_issue14_charging_battery_lowers_voltage(lv_net):
    """The point of the issue: the dispatch must move the grid state.

    A charging battery is extra load, so every bus carrying one must sit lower
    than in the same snapshot without it. While the device was inert, the two
    solves returned an identical voltage profile.
    """
    base, _ = _pinned_snapshot(lv_net)
    charging, buses = _pinned_snapshot(
        lv_net, dispatch=(DISPATCH_PU, 0.0, 0.0)
    )

    assert buses
    for b in buses:
        assert charging[b] < base[b] - 1e-6, f"bus {b} did not drop"


@pytest.mark.integration
def test_issue14_discharging_battery_raises_voltage(lv_net):
    """The opposite sign: discharging injects, so voltages must rise."""
    base, _ = _pinned_snapshot(lv_net)
    discharging, buses = _pinned_snapshot(
        lv_net, dispatch=(0.0, DISPATCH_PU, 0.0)
    )

    assert buses
    for b in buses:
        assert discharging[b] > base[b] + 1e-6, f"bus {b} did not rise"


@pytest.mark.integration
def test_issue14_capacitive_reactive_power_raises_voltage(lv_net):
    """Pins the reactive sign convention.

    ``BAT_Q`` is generator convention, so a positive value is capacitive
    injection and must raise voltages. A sign error here would show up as
    reactive support depressing the profile.
    """
    base, _ = _pinned_snapshot(lv_net)
    capacitive, buses = _pinned_snapshot(
        lv_net, dispatch=(0.0, 0.0, DISPATCH_PU)
    )
    inductive, _ = _pinned_snapshot(lv_net, dispatch=(0.0, 0.0, -DISPATCH_PU))

    assert buses
    for b in buses:
        assert capacitive[b] > base[b] + 1e-6, f"bus {b} did not rise"
        assert inductive[b] < base[b] - 1e-6, f"bus {b} did not drop"


# ---------------------------------------------------------------------------
# #14 — BESS reactive capability
# ---------------------------------------------------------------------------


def test_issue14_battery_has_reactive_power_on_ac_models(lv_net):
    opf, _ = _with_battery(_prepared(lv_net))
    assert hasattr(opf.model, "BAT_Q")
    assert hasattr(opf.model, "bat_inverter_s2")


def test_issue14_s_circle_bounds_apparent_power(lv_net):
    s_inv = 0.02
    opf, _ = _with_battery(_prepared(lv_net), s_inv_pu=s_inv)
    m = opf.model
    b, t = list(m.BAT)[0], list(m.T)[0]
    # At the rating in both P and Q the circle must be violated.
    m.BAT_Pchg[b, t].set_value(s_inv)
    m.BAT_Pdis[b, t].set_value(0.0)
    m.BAT_Q[b, t].set_value(s_inv)
    assert pyo.value(m.bat_inverter_s2[b, t].body) > s_inv**2
    # Pure reactive at the rating is exactly on the circle (STATCOM mode).
    m.BAT_Pchg[b, t].set_value(0.0)
    assert pyo.value(m.bat_inverter_s2[b, t].body) == pytest.approx(s_inv**2)


def test_issue14_reactive_var_is_boxed_by_the_rating(lv_net):
    """Implied by the circle, but IPOPT needs it stated to stay inside."""
    s_inv = 0.02
    opf, _ = _with_battery(_prepared(lv_net), s_inv_pu=s_inv)
    for idx in opf.model.BAT_Q:
        assert opf.model.BAT_Q[idx].lb == pytest.approx(-s_inv)
        assert opf.model.BAT_Q[idx].ub == pytest.approx(s_inv)


def test_issue14_converter_may_not_be_smaller_than_the_power_rating(lv_net):
    with pytest.raises(ValueError, match=r"below power_pu"):
        _with_battery(_prepared(lv_net), s_inv_pu=0.005)


def test_issue14_cos_phi_limits_reactive_power(lv_net):
    opf, _ = _with_battery(_prepared(lv_net), s_inv_pu=0.02, cos_phi_min=0.9)
    m = opf.model
    assert hasattr(m, "bat_cos_phi_upper")
    assert hasattr(m, "bat_cos_phi_lower")
    tan_phi = pyo.value(m.BAT_tan_phi[list(m.BAT)[0]])
    assert tan_phi == pytest.approx(np.tan(np.arccos(0.9)))


@pytest.mark.parametrize("cos_phi_min", [0.0, -0.5, 1.5])
def test_issue14_invalid_cos_phi_raises(lv_net, cos_phi_min):
    with pytest.raises(ValueError, match=r"cos_phi_min must be in"):
        _with_battery(_prepared(lv_net), cos_phi_min=cos_phi_min)


def test_issue14_invalid_q_control_raises(lv_net):
    with pytest.raises(ValueError, match=r"q_control must be"):
        _with_battery(_prepared(lv_net), q_control="bogus")


def test_issue14_qp_area_is_keyed_on_the_discharge_leg(lv_net):
    """A storage unit acts as a generating unit while discharging."""
    opf, _ = _with_battery(_prepared(lv_net), s_inv_pu=0.02, q_control="qp")
    m = opf.model
    assert hasattr(m, "bat_QP_pos")
    b, t = list(m.BAT)[0], list(m.T)[0]
    used = _vars_in(m.bat_QP_pos[b, t, 0])
    assert "BAT_Pdis" in used
    assert "BAT_Pchg" not in used


def test_issue14_qu_area_uses_the_bus_voltage(lv_net):
    opf, bat = _with_battery(_prepared(lv_net), s_inv_pu=0.02, q_control="qu")
    m = opf.model
    assert hasattr(m, "bat_QU_max")
    b, t = list(m.BAT)[0], list(m.T)[0]
    # The constraint must reference the voltage at the battery's own ppc bus.
    bus_of = {d: int(bat.bus_lookup[int(p)]) for d, p in list(m.BAT_bus)}
    body = str(m.bat_QU_max[b, t, 0].body)
    assert f"v[{bus_of[b]},{t}]" in body


@pytest.mark.integration
def test_issue14_reactive_dispatch_respects_the_circle(lv_net):
    """Solve and check the realised operating point, not just the constraint."""
    s_inv = 0.02
    opf, _ = _with_battery(_prepared(lv_net), toT=6, s_inv_pu=s_inv)
    opf.add_voltage_deviation_objective()
    assert pyo.check_optimal_termination(
        opf.solve(solver="ipopt", print_solver_output=False)
    )
    m = opf.model
    for b in m.BAT:
        for t in m.T:
            s = (
                pyo.value(m.BAT_P[b, t]) ** 2 + pyo.value(m.BAT_Q[b, t]) ** 2
            ) ** 0.5
            assert s <= s_inv + 1e-6, f"battery {b} at t={t} exceeds S_inv"
