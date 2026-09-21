# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Unit and integration tests for the DCOPF model."""

import numpy as np
import pytest
import pandapower as pp
import pyomo.environ as pyo

from potpourri.models.DCOPF import DCOPF


@pytest.fixture()
def four_bus():
    return pp.networks.simple_four_bus_system()


# ── model structure (no solver required) ─────────────────────────────────────


def test_dcopf_model_structure(four_bus):
    """DCOPF builds expected Pyomo components after add_OPF()."""
    dcopf = DCOPF(four_bus)
    dcopf.add_OPF()

    assert hasattr(dcopf.model, "B"), "missing bus set"
    assert hasattr(dcopf.model, "L"), "missing line set"
    assert hasattr(dcopf.model, "SLmax"), "missing SLmax param"
    assert hasattr(dcopf.model, "line_lim_from"), "missing line upper limit"
    assert hasattr(dcopf.model, "line_lim_to"), "missing line lower limit"
    assert hasattr(dcopf.model, "PG_Constraint"), "missing generator limit"
    assert len(list(dcopf.model.B)) > 0
    assert len(list(dcopf.model.L)) > 0


def test_dcopf_has_no_reactive_vars(four_bus):
    """DC model must not declare reactive-power variables (qG, qsG, qD)."""
    dcopf = DCOPF(four_bus)
    for name in ("qG", "qsG", "qD", "qLfrom", "qLto"):
        assert not hasattr(dcopf.model, name), (
            f"DC model should not have reactive variable '{name}'"
        )


def test_dcopf_line_limits_non_negative(four_bus):
    """Thermal line limits SLmax must be non-negative."""
    dcopf = DCOPF(four_bus)
    dcopf.add_OPF()
    for line in dcopf.model.L:
        assert pyo.value(dcopf.model.SLmax[line]) >= 0, (
            f"negative SLmax for line {line}"
        )


# ── integration: local solver required ───────────────────────────────────────


@pytest.mark.integration
def test_dcopf_solves_with_glpk(four_bus):
    """DCOPF finds an optimal solution using the local GLPK solver."""
    dcopf = DCOPF(four_bus)
    dcopf.add_OPF()
    dcopf.model.obj = pyo.Objective(
        expr=sum(dcopf.model.pG[g] for g in dcopf.model.G),
        sense=pyo.minimize,
    )
    dcopf.solve(solver="glpk", to_net=False)
    assert pyo.check_optimal_termination(dcopf.results)


@pytest.mark.integration
@pytest.mark.slow
def test_dcopf_solves_with_neos():
    """DCOPF finds an optimal solution using the NEOS CPLEX solver."""
    import os

    os.environ.setdefault("NEOS_EMAIL", "test@example.com")
    net = pp.networks.simple_four_bus_system()
    dcopf = DCOPF(net)
    dcopf.add_OPF()
    dcopf.model.obj = pyo.Objective(
        expr=sum(dcopf.model.pG[g] for g in dcopf.model.G),
        sense=pyo.minimize,
    )
    dcopf.solve(solver="neos", neos_opt="cplex", to_net=False)
    assert pyo.check_optimal_termination(dcopf.results)


def test_dc_conventions():
    """``-1/x`` (MATPOWER, default) versus ``-x/(r²+x²)`` (PowerModels)."""
    import pandapower as pp
    import pytest

    from potpourri.models.DCOPF import DCOPF

    net = pp.create_empty_network(sn_mva=100.0)
    b0 = pp.create_bus(net, 110.0)
    b1 = pp.create_bus(net, 110.0)
    pp.create_ext_grid(net, b0)
    pp.create_line_from_parameters(
        net,
        b0,
        b1,
        1.0,
        r_ohm_per_km=6.05,
        x_ohm_per_km=12.1,
        c_nf_per_km=0.0,
        max_i_ka=1.0,
    )
    pp.create_load(net, b1, 20.0, 5.0)
    zn = 110.0**2 / 100.0  # per-unit base impedance
    r, x = 6.05 / zn, 12.1 / zn
    default = DCOPF(net)
    powermodels = DCOPF(net, dc_convention="powermodels")
    assert default.model.BL[0] == pytest.approx(-1 / x)
    assert powermodels.model.BL[0] == pytest.approx(-x / (r**2 + x**2))
    assert powermodels.model.BL[0] != pytest.approx(default.model.BL[0])
    with pytest.raises(ValueError):
        DCOPF(net, dc_convention="lossless")


def test_powermodels_convention_drops_the_phase_shift():
    """PowerModels' DC flow is p = -b (θ_from − θ_to): no phase shift."""
    import pandapower as pp
    import pyomo.environ as pyo

    from potpourri.models.DCOPF import DCOPF

    net = pp.create_empty_network(sn_mva=100.0)
    hv = pp.create_bus(net, 110.0)
    lv = pp.create_bus(net, 20.0)
    pp.create_ext_grid(net, hv)
    pp.create_transformer_from_parameters(
        net,
        hv,
        lv,
        sn_mva=40.0,
        vn_hv_kv=110.0,
        vn_lv_kv=20.0,
        vkr_percent=0.3,
        vk_percent=10.0,
        pfe_kw=0.0,
        i0_percent=0.0,
        shift_degree=30.0,
    )
    pp.create_load(net, lv, 10.0, 2.0)
    default = DCOPF(net)
    powermodels = DCOPF(net, dc_convention="powermodels")
    # the transformer's angle-difference constraint carries the shift only
    # under the MATPOWER convention: with every angle variable at zero its
    # body reduces to the shift term (deltaLT − (δ_hv − δ_lv − shift) = shift)
    shift = np.deg2rad(30.0)
    assert pyo.value(default.model.shift[0]) == pytest.approx(shift)
    for model in (default.model, powermodels.model):
        for var in (model.delta, model.deltaLT):
            for idx in var:
                var[idx].set_value(0.0)
    assert pyo.value(default.model.phase_diff2[0].body) == pytest.approx(shift)
    assert pyo.value(powermodels.model.phase_diff2[0].body) == pytest.approx(
        0.0
    )


def _trafo_feeder(tap_side, tap_pos):
    """110/20 kV transformer with a tap changer, load on the LV bus."""
    net = pp.create_empty_network(sn_mva=100.0)
    hv = pp.create_bus(net, 110.0)
    lv = pp.create_bus(net, 20.0)
    pp.create_ext_grid(net, hv)
    pp.create_transformer_from_parameters(
        net,
        hv,
        lv,
        sn_mva=40.0,
        vn_hv_kv=110.0,
        vn_lv_kv=20.0,
        vk_percent=12.0,
        vkr_percent=0.5,
        pfe_kw=0.0,
        i0_percent=0.0,
        tap_side=tap_side,
        tap_neutral=0,
        tap_pos=tap_pos,
        tap_step_percent=2.0,
        tap_min=-9,
        tap_max=9,
        tap_changer_type="Ratio",
    )
    pp.create_load(net, lv, 10.0, 2.0)
    return net


def test_powermodels_convention_ignores_the_tap_position():
    """An LV-side tap refers the impedance to the tapped voltage.

    Pandapower refers the series impedance to the tapped LV voltage when
    the tap sits on the LV side, (vn_trafo_lv / vn_lv_kv)². PowerModels reads
    the untapped reactance, so the susceptance must not move with the tap on
    either side; the default convention keeps pandapower's value.
    """
    neutral = DCOPF(_trafo_feeder("lv", 0), dc_convention="powermodels")
    b_neutral = float(neutral.trafo_data["BLT_data"].iloc[0])
    for side in ("lv", "hv"):
        tapped = DCOPF(_trafo_feeder(side, 5), dc_convention="powermodels")
        assert float(tapped.trafo_data["BLT_data"].iloc[0]) == pytest.approx(
            b_neutral, rel=1e-9
        ), side
    z_pu = 0.12 * 100.0 / 40.0  # vk % of the trafo rating, at the system base
    r_pu = 0.005 * 100.0 / 40.0
    x_pu = (z_pu**2 - r_pu**2) ** 0.5
    assert b_neutral == pytest.approx(-x_pu / (r_pu**2 + x_pu**2), rel=1e-6)
    scaled = DCOPF(_trafo_feeder("lv", 5))  # matpower convention: ppc value
    assert float(scaled.trafo_data["BLT_data"].iloc[0]) == pytest.approx(
        -1 / (x_pu * 1.1**2), rel=1e-6
    )
