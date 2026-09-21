# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Regression tests for the formulation-audit fixes (D1–D14).

Each test isolates one bug path identified in the audit and pins the
expected post-fix behaviour. Slower integration-style tests are marked.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandapower as pp
import pyomo.environ as pyo
import pytest

from potpourri.benchmarks.pglib import PGLIB_ROOT
from potpourri.models.ACOPF_base import ACOPF
from potpourri.models.DCOPF import DCOPF
from potpourri.models.basemodel import preprocess_grid
from potpourri.models.cost_objective import add_poly_cost_objective


warnings.filterwarnings("ignore")


# Auto-skip PGLib integration tests when the benchmarks/pglib-opf git
# submodule hasn't been initialised (e.g. on a fresh `git clone` without
# `--recurse-submodules`). Marks every PGLib-loading test in this module.
requires_pglib_submodule = pytest.mark.skipif(
    not (PGLIB_ROOT / "pglib_opf_case5_pjm.m").is_file(),
    reason=(
        "benchmarks/pglib-opf submodule not initialised; "
        "run `git submodule update --init` to enable PGLib tests"
    ),
)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _toy_two_bus_net(slack_q_range=(-100, 100)):
    """A minimal 2-bus net with one ext_grid, one load, one line."""
    net = pp.create_empty_network(sn_mva=100.0)
    b0 = pp.create_bus(net, vn_kv=110.0, max_vm_pu=1.1, min_vm_pu=0.9)
    b1 = pp.create_bus(net, vn_kv=110.0, max_vm_pu=1.1, min_vm_pu=0.9)
    pp.create_ext_grid(
        net,
        bus=b0,
        vm_pu=1.0,
        max_p_mw=500.0,
        min_p_mw=-500.0,
        max_q_mvar=slack_q_range[1],
        min_q_mvar=slack_q_range[0],
    )
    pp.create_load(net, bus=b1, p_mw=50.0, q_mvar=10.0)
    pp.create_line_from_parameters(
        net,
        from_bus=b0,
        to_bus=b1,
        length_km=10.0,
        r_ohm_per_km=0.1,
        x_ohm_per_km=0.3,
        c_nf_per_km=10.0,
        max_i_ka=1.0,
    )
    return net


# ---------------------------------------------------------------------------
# D5 — shunt-only buses survive preprocess_grid
# ---------------------------------------------------------------------------


def test_d5_shunt_only_bus_survives_preprocessing():
    """A bus referenced only by a shunt must not be dropped."""
    net = _toy_two_bus_net()
    b_shunt = pp.create_bus(net, vn_kv=110.0, max_vm_pu=1.1, min_vm_pu=0.9)
    pp.create_shunt(net, bus=b_shunt, q_mvar=5.0)
    # Attach b_shunt to b1 via a line so the network is still connected for
    # pp.runpp; the test asserts that shunt's reference is kept consistent.
    pp.create_line_from_parameters(
        net,
        from_bus=1,
        to_bus=b_shunt,
        length_km=1.0,
        r_ohm_per_km=0.1,
        x_ohm_per_km=0.3,
        c_nf_per_km=10.0,
        max_i_ka=1.0,
    )

    processed = preprocess_grid(net)
    assert len(processed.bus) == 3
    assert (processed.shunt["bus"] >= 0).all()
    # Shunt's bus must still point to a valid bus row.
    assert all(b in processed.bus.index for b in processed.shunt["bus"])


# ---------------------------------------------------------------------------
# D6 — bus-bus switch merge updates trafo + shunt + storage refs
# ---------------------------------------------------------------------------


def test_d6_switch_merge_updates_all_bus_references():
    net = pp.create_empty_network(sn_mva=100.0)
    bA = pp.create_bus(net, vn_kv=20.0)
    bB = pp.create_bus(net, vn_kv=20.0)
    bC = pp.create_bus(net, vn_kv=20.0)
    pp.create_ext_grid(net, bus=bA, vm_pu=1.0)
    pp.create_load(net, bus=bC, p_mw=1.0, q_mvar=0.2)
    pp.create_line_from_parameters(
        net,
        from_bus=bA,
        to_bus=bC,
        length_km=1.0,
        r_ohm_per_km=0.1,
        x_ohm_per_km=0.3,
        c_nf_per_km=10.0,
        max_i_ka=1.0,
    )
    pp.create_shunt(net, bus=bB, q_mvar=2.0)
    # closed zero-impedance bus-bus switch between bB and bA
    pp.create_switch(net, bus=bA, element=bB, et="b", closed=True, z_ohm=0.0)

    processed = preprocess_grid(net)
    # bB should be merged into bA; shunt's bus must now equal bA's new id
    # (after create_continuous_bus_index, indices are 0..N-1).
    assert len(processed.shunt) == 1
    assert int(processed.shunt.at[0, "bus"]) in processed.bus.index.tolist()


# ---------------------------------------------------------------------------
# D3 — out-of-service gen doesn't break P-limit broadcasting
# ---------------------------------------------------------------------------


def test_d3_out_of_service_gen_does_not_break_p_limits():
    net = _toy_two_bus_net()
    pp.create_gen(
        net,
        bus=1,
        p_mw=10.0,
        vm_pu=1.0,
        max_p_mw=100.0,
        min_p_mw=0.0,
        max_q_mvar=50.0,
        min_q_mvar=-50.0,
        controllable=True,
    )
    pp.create_gen(
        net,
        bus=1,
        p_mw=0.0,
        vm_pu=1.0,
        max_p_mw=80.0,
        min_p_mw=0.0,
        max_q_mvar=40.0,
        min_q_mvar=-40.0,
        controllable=True,
        in_service=False,
    )
    acopf = ACOPF(net)
    acopf.add_OPF()
    # Both gens (g0 + ext_grid) end up in model.G; the OOS gen is filtered.
    # Just constructing add_OPF without exception is the regression test.
    assert len(acopf.model.G) == 2


# ---------------------------------------------------------------------------
# D9 — same in-service filter for Q limits
# ---------------------------------------------------------------------------


def test_d9_out_of_service_gen_q_limits_filtered():
    net = _toy_two_bus_net()
    pp.create_gen(
        net,
        bus=1,
        p_mw=10.0,
        vm_pu=1.0,
        max_p_mw=100.0,
        min_p_mw=0.0,
        max_q_mvar=50.0,
        min_q_mvar=-50.0,
        controllable=True,
    )
    pp.create_gen(
        net,
        bus=1,
        p_mw=0.0,
        vm_pu=1.0,
        max_p_mw=80.0,
        min_p_mw=0.0,
        max_q_mvar=40.0,
        min_q_mvar=-40.0,
        controllable=True,
        in_service=False,
    )
    acopf = ACOPF(net)
    # Construction enough; nothing should raise on Q-limit broadcasting.
    acopf.add_OPF()


# ---------------------------------------------------------------------------
# D4 — degenerate Pmin == Pmax pins the variable instead of adding a
# redundant range constraint.
# ---------------------------------------------------------------------------


def test_d4_degenerate_gen_range_pins_variable():
    """A generator with Pmin == Pmax gets no range constraint.

    A gen with Pmin == Pmax must not contribute a range constraint
    (would trip IPOPT TOO_FEW_DOF) and the variable must be pinned via
    tight bounds, not Var.fix() (which would eliminate the var from the
    NL file and again reduce degrees of freedom).
    """
    net = _toy_two_bus_net()
    pp.create_gen(
        net,
        bus=1,
        p_mw=0.0,
        vm_pu=1.0,
        max_p_mw=0.0,
        min_p_mw=0.0,
        max_q_mvar=10.0,
        min_q_mvar=-10.0,
        controllable=True,
    )
    acopf = ACOPF(net)
    acopf.add_OPF()
    g_sync = max(acopf.model.G)
    # Range constraint must be skipped for this gen.
    assert g_sync not in acopf.model.PG_Constraint
    # Variable must be tightly bounded.
    lb = acopf.model.pG[g_sync].lb
    ub = acopf.model.pG[g_sync].ub
    assert lb is not None and ub is not None
    assert abs(ub - lb) < 1e-6
    # And not "fixed" in the Pyomo sense (otherwise the NL writer would
    # eliminate it and we'd lose a degree of freedom downstream).
    assert not acopf.model.pG[g_sync].is_fixed()


# ---------------------------------------------------------------------------
# D1 — slack voltage magnitude is free in ACOPF by default
# ---------------------------------------------------------------------------


def test_d1_slack_voltage_magnitude_is_free_by_default():
    net = _toy_two_bus_net()
    acopf = ACOPF(net)
    acopf.add_OPF()
    for b in acopf.model.b0:
        assert not acopf.model.v[b].is_fixed(), (
            "Slack v should be free in default ACOPF (free_slack_vm=True)"
        )


def test_d1_slack_voltage_pinned_when_requested():
    net = _toy_two_bus_net()
    acopf = ACOPF(net)
    acopf.add_OPF(free_slack_vm=False)
    for b in acopf.model.b0:
        assert acopf.model.v[b].is_fixed()


# ---------------------------------------------------------------------------
# D2 — 110 kV pinning is opt-in
# ---------------------------------------------------------------------------


def test_d2_no_hardcoded_110_kv_pin_by_default():
    net = _toy_two_bus_net()  # 110 kV buses
    acopf = ACOPF(net)
    acopf.add_OPF()
    # By default no 110-kV pinning constraint should be active.
    assert len(list(acopf.model.Bfix)) == 0


def test_d2_hv_bus_pin_can_be_enabled_explicitly():
    net = _toy_two_bus_net()  # both buses at 110 kV
    acopf = ACOPF(net)
    acopf.add_OPF(fix_hv_buses=True, hv_bus_kv=110.0)
    assert len(list(acopf.model.Bfix)) == 2


# ---------------------------------------------------------------------------
# D7 — thermal_limit flag toggles constraint form
# ---------------------------------------------------------------------------


def test_d7_thermal_limit_default_is_current():
    net = _toy_two_bus_net()
    acopf = ACOPF(net)
    acopf.add_OPF()
    assert acopf.thermal_limit_mode == "current"
    # An MVA-limit constraint would have a constant RHS; a current-limit
    # constraint contains the voltage variable squared.
    c0 = acopf.model.line_lim_from[0]
    body_str = str(c0.body)
    assert "v[" in body_str


def test_d7_thermal_limit_mva_mode():
    net = _toy_two_bus_net()
    acopf = ACOPF(net)
    acopf.add_OPF(thermal_limit="mva")
    assert acopf.thermal_limit_mode == "mva"
    c0 = acopf.model.line_lim_from[0]
    body_str = str(c0.body)
    # MVA mode shouldn't reference voltage variables in the constraint body.
    assert "v[" not in body_str


def test_d7_invalid_thermal_limit_raises():
    net = _toy_two_bus_net()
    acopf = ACOPF(net)
    with pytest.raises(ValueError):
        acopf.add_OPF(thermal_limit="bogus")


# ---------------------------------------------------------------------------
# D8 — branch angle-difference constraints opt-in
# ---------------------------------------------------------------------------


def test_d8_angle_limits_off_by_default():
    net = _toy_two_bus_net()
    acopf = ACOPF(net)
    acopf.add_OPF()
    assert not hasattr(acopf.model, "line_angle_diff")


def test_d8_angle_limits_enabled_when_data_present():
    net = _toy_two_bus_net()
    # attach explicit angle bounds to the single line
    net.line["angmin_degree"] = -30.0
    net.line["angmax_degree"] = 30.0
    acopf = ACOPF(net)
    acopf.add_OPF(angle_limits=True)
    assert hasattr(acopf.model, "line_angle_diff")
    # The constraint should be a range (lo ≤ body ≤ hi) tied to delta diff.
    c = acopf.model.line_angle_diff[0]
    assert c.has_lb() and c.has_ub()


def test_d8_dcopf_angle_limits_optional():
    net = _toy_two_bus_net()
    net.line["angmin_degree"] = -30.0
    net.line["angmax_degree"] = 30.0
    dcopf = DCOPF(net)
    dcopf.add_OPF(angle_limits=True)
    assert hasattr(dcopf.model, "line_angle_diff")


# ---------------------------------------------------------------------------
# D10 — sgen min_p falls back to p_mw, not 0
# ---------------------------------------------------------------------------


def test_d10_sgen_min_p_explicit_value_respected():
    """An sgen with explicit min_p_mw must keep it (no silent 0-flooring)."""
    net = _toy_two_bus_net()
    sg = pp.create_sgen(
        net,
        bus=1,
        p_mw=12.0,
        max_p_mw=20.0,
        min_p_mw=5.0,
        controllable=True,
    )
    acopf = ACOPF(net)
    acopf.add_OPF()
    val = float(acopf.static_generation_data.loc[sg, "min_p"])
    assert val == pytest.approx(5.0 / 100.0)


def test_d10_sgen_min_p_falls_back_to_zero_when_missing():
    """A missing min_p_mw column must not break the bounds.

    When min_p_mw is missing entirely (e.g. converted PGLib networks
    that don't carry the column on sgen), keep the historical 0-floor
    convention so OPF doesn't accidentally pin sgens to their setpoint.
    """
    net = _toy_two_bus_net()
    pp.create_sgen(
        net,
        bus=1,
        p_mw=12.0,
        max_p_mw=20.0,
        controllable=True,
    )
    if "min_p_mw" in net.sgen.columns:
        net.sgen.drop(columns=["min_p_mw"], inplace=True)
    acopf = ACOPF(net)
    acopf.add_OPF()
    val = float(acopf.static_generation_data.loc[0, "min_p"])
    assert val == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# D13 — PWL costs raise a clear error
# ---------------------------------------------------------------------------


def test_d13_pwl_cost_raises():
    """pwl_cost rows must trigger a clear error, not silent miscosting."""
    net = _toy_two_bus_net()
    pp.create_pwl_cost(
        net, element=0, et="ext_grid", points=[[-100, 0, 0], [0, 100, 100]]
    )
    acopf = ACOPF(net)
    acopf.add_OPF()
    with pytest.raises(ValueError, match="Piecewise-linear"):
        add_poly_cost_objective(acopf)


# ---------------------------------------------------------------------------
# D14 — pyo_to_net works with non-contiguous gen indices.
# (Smoke test: solve a tiny AC OPF and write back.)
# ---------------------------------------------------------------------------


@pytest.mark.integration
def test_d14_pyo_to_net_smoke():
    net = _toy_two_bus_net()
    pp.create_poly_cost(net, element=0, et="ext_grid", cp1_eur_per_mw=20.0)
    acopf = ACOPF(net)
    acopf.add_OPF(thermal_limit="mva")
    add_poly_cost_objective(acopf)
    res = acopf.solve(solver="ipopt", print_solver_output=False)
    assert pyo.check_optimal_termination(res)
    # res_ext_grid should now have one row with a finite p_mw
    assert np.isfinite(acopf.net.res_ext_grid.p_mw.iloc[0])


# ---------------------------------------------------------------------------
# Impedance-branch support (case89_pegase / case118_ieee in PGLib have
# pandapower net.impedance entries; Basemodel used to silently drop them)
# ---------------------------------------------------------------------------


def test_impedance_branch_included_in_model_L():
    """An impedance row appears in model.L with its constraints.

    A net.impedance row must show up in model.L at index >= len(net.line)
    with KVL constraints attached. Previously Basemodel only iterated
    net.line + net.trafo.
    """
    net = _toy_two_bus_net()
    # Add a third bus that is only reachable via an impedance branch
    b2 = pp.create_bus(net, vn_kv=33.0, max_vm_pu=1.1, min_vm_pu=0.9)
    pp.create_impedance(
        net,
        from_bus=1,
        to_bus=b2,
        rft_pu=0.01,
        xft_pu=0.05,
        sn_mva=100.0,
    )
    pp.create_load(net, bus=b2, p_mw=5.0, q_mvar=1.0)

    acopf = ACOPF(net)
    acopf.add_OPF(thermal_limit="mva")
    n_line = len(net.line)
    L = list(acopf.model.L)
    # The impedance branch should appear in model.L at index n_line.
    assert n_line in L, f"Expected impedance row at L index {n_line}; got {L}"
    # KVL constraints must exist for the impedance row.
    assert n_line in acopf.model.KVL_real_from
    assert n_line in acopf.model.KVL_real_to
    # bus_line_dict must point to the impedance endpoints.
    assert acopf.bus_line_dict[(n_line, 1)] == 1
    assert acopf.bus_line_dict[(n_line, 2)] == b2


@pytest.mark.integration
@requires_pglib_submodule
def test_pglib_case118_ieee_no_longer_isolated_bus():
    """case118's impedance branch leaves no bus infeasible.

    case118 has an impedance branch (67↔115) that previously left bus 115
    KCL-infeasible. With impedance branches included, AC OPF must converge
    within 1 % of the published reference.
    """
    from potpourri.benchmarks import load_pglib_case, PGLIB_BASELINE_TYP

    net = load_pglib_case("case118_ieee")
    acopf = ACOPF(net)
    acopf.add_OPF(thermal_limit="mva", angle_limits=True)
    add_poly_cost_objective(acopf)
    res = acopf.solve(solver="ipopt", print_solver_output=False)
    assert pyo.check_optimal_termination(res)
    obj = pyo.value(acopf.model.obj_poly_cost)
    ref = PGLIB_BASELINE_TYP["pglib_opf_case118_ieee"]["ac"]
    assert abs(obj - ref) / ref < 0.01, f"obj={obj} ref={ref}"


# ---------------------------------------------------------------------------
# Integration: tiny PGLib case5_pjm should match the published reference.
# ---------------------------------------------------------------------------


@pytest.mark.integration
@requires_pglib_submodule
def test_pglib_case5_pjm_matches_reference():
    from potpourri.benchmarks import load_pglib_case, PGLIB_BASELINE_TYP

    net = load_pglib_case("case5_pjm")
    acopf = ACOPF(net)
    acopf.add_OPF(thermal_limit="mva", angle_limits=True)
    add_poly_cost_objective(acopf)
    res = acopf.solve(solver="ipopt", print_solver_output=False)
    assert pyo.check_optimal_termination(res)
    obj = pyo.value(acopf.model.obj_poly_cost)
    ref = PGLIB_BASELINE_TYP["pglib_opf_case5_pjm"]["ac"]
    assert abs(obj - ref) / ref < 1e-3, f"obj={obj} ref={ref}"
