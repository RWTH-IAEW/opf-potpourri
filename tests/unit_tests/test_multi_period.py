# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Unit tests for multi-period model construction (no solver required)."""

import pandapower as pp
import pytest
import simbench as sb

from potpourri.models_multi_period.ACOPF_multi_period import ACOPF_multi_period
from potpourri.models_multi_period.DCOPF_multi_period import DCOPF_multi_period


@pytest.fixture(scope="module")
def lv_net():
    """SimBench LV rural network (loaded once per module)."""
    return sb.get_simbench_net("1-LV-rural1--0-sw")


# ── time index ────────────────────────────────────────────────────────────────


def test_time_set_length(lv_net):
    """Model.T must have exactly toT − fromT elements."""
    fromT, toT = 0, 8
    opf = ACOPF_multi_period(lv_net, toT=toT, fromT=fromT)
    assert len(list(opf.model.T)) == toT - fromT


def test_time_set_range(lv_net):
    """Model.T must span [fromT, toT)."""
    fromT, toT = 4, 12
    opf = ACOPF_multi_period(lv_net, toT=toT, fromT=fromT)
    t_list = sorted(opf.model.T)
    assert t_list[0] == fromT
    assert t_list[-1] == toT - 1


# ── model structure ───────────────────────────────────────────────────────────


def test_multi_period_has_bus_set(lv_net):
    """Multi-period model must have a bus set B."""
    opf = ACOPF_multi_period(lv_net, toT=4)
    assert hasattr(opf.model, "B")
    assert len(list(opf.model.B)) == len(lv_net.bus)


def test_multi_period_voltage_variable_indexed(lv_net):
    """Voltage variable v must be indexed by (bus, time)."""
    opf = ACOPF_multi_period(lv_net, toT=4)
    assert hasattr(opf.model, "v")
    # Spot-check: v[b, t] should exist for first bus and first time step
    b0 = list(opf.model.B)[0]
    t0 = list(opf.model.T)[0]
    assert (b0, t0) in opf.model.v


def test_multi_period_kcl_count(lv_net):
    """KCL_real must have |B| × |T| constraints."""
    fromT, toT = 0, 4
    opf = ACOPF_multi_period(lv_net, toT=toT, fromT=fromT)
    n_buses = len(list(opf.model.B))
    n_steps = toT - fromT
    assert len(list(opf.model.KCL_real)) == n_buses * n_steps


def test_multi_period_add_opf_does_not_raise(lv_net):
    """add_OPF() on a multi-period model must not raise."""
    opf = ACOPF_multi_period(lv_net, toT=4)
    opf.add_OPF()  # should not raise


def test_multi_period_voltage_deviation_objective(lv_net):
    """add_voltage_deviation_objective creates obj_v_deviation."""
    opf = ACOPF_multi_period(lv_net, toT=4)
    opf.add_OPF()
    opf.add_voltage_deviation_objective()
    assert hasattr(opf.model, "obj_v_deviation")


# ── battery device attachment ─────────────────────────────────────────────────


def test_battery_attaches_to_model(lv_net):
    """Battery_multi_period.get_all() adds BAT and BAT_SOC to the model."""
    from potpourri.technologies.battery import Battery_multi_period

    fromT, toT = 0, 4
    opf = ACOPF_multi_period(lv_net, toT=toT, fromT=fromT)
    battery = Battery_multi_period(opf.net, T=toT - fromT, scenario=0)
    battery.get_all(opf.model)

    assert hasattr(opf.model, "BAT"), (
        "model.BAT set missing after battery attach"
    )
    assert hasattr(opf.model, "BAT_SOC"), "model.BAT_SOC var missing"


# ── impedance-branch (mirror of single-period fix) ───────────────────────────


def _net_with_impedance():
    """SimBench LV net with one extra pp.create_impedance branch added."""
    net = sb.get_simbench_net("1-LV-rural1--0-sw")
    a = int(net.bus.index[0])
    b = int(net.bus.index[3])
    pp.create_impedance(
        net, from_bus=a, to_bus=b, rft_pu=0.01, xft_pu=0.05, sn_mva=1.0
    )
    return net


def test_mp_ac_includes_impedance_branch_in_L():
    """Multi-period AC OPF must put net.impedance rows into model.L at
    synthetic indices >= len(net.line) (mirror of the single-period fix)."""
    net = _net_with_impedance()
    n_line = len(net.line)
    opf = ACOPF_multi_period(net, toT=2)
    L = list(opf.model.L)
    assert n_line in L, f"impedance index {n_line} missing from model.L: {L}"
    assert opf._n_native_lines == n_line
    # bus_line_dict must point to the impedance branch's bus pair
    assert (n_line, 1) in opf.bus_line_dict
    assert (n_line, 2) in opf.bus_line_dict


def test_mp_basemodel_includes_impedance_branch_in_line_data():
    """Multi-period ``Basemodel_multi_period`` must register impedance rows
    in ``line_data`` and ``bus_line_dict``. This is the shared shim both
    ``AC_multi_period`` and ``DC_multi_period`` rely on.
    """
    from potpourri.models_multi_period.basemodel_multi_period import (
        Basemodel_multi_period,
    )

    net = _net_with_impedance()
    n_line = len(net.line)
    bm = Basemodel_multi_period(net, toT=2)
    assert bm._n_native_lines == n_line
    assert n_line in bm.line_data.index
    assert (n_line, 1) in bm.bus_line_dict
    assert (n_line, 2) in bm.bus_line_dict


def test_mp_dc_constructs_with_time_steps():
    """``DC_multi_period`` must construct cleanly in the time-variant path
    (where ``self.model.T`` is the Pyomo Set and all branch / shunt
    parameters use the single-period indexing they were declared with).
    Regression test for the prior `self.T` int-vs-Set confusion and the
    spurious time index on `BL`, `BLT`, `shift`, `GB`.
    """
    from potpourri.models_multi_period.DC_multi_period import DC_multi_period

    net = sb.get_simbench_net("1-LV-rural1--0-sw")
    dc = DC_multi_period(net, toT=4)
    assert sorted(list(dc.model.T)) == [0, 1, 2, 3]
    # |KCL_def| == |B| × |T|
    n_b = len(list(dc.model.B))
    n_t = len(list(dc.model.T))
    assert len(list(dc.model.KCL_def)) == n_b * n_t
    # |KVL_real_fromend| == |L| × |T|
    n_l = len(list(dc.model.L))
    assert len(list(dc.model.KVL_real_fromend)) == n_l * n_t
    # |phase_angle_diff1| == |L| × |T|
    assert len(list(dc.model.phase_angle_diff1)) == n_l * n_t


def test_mp_dc_includes_impedance_branch_in_L():
    """``DC_multi_period`` must include net.impedance in model.L (mirror of
    the single-period DC fix)."""
    from potpourri.models_multi_period.DC_multi_period import DC_multi_period

    net = _net_with_impedance()
    n_line = len(net.line)
    dc = DC_multi_period(net, toT=2)
    L = list(dc.model.L)
    assert n_line in L, f"impedance index {n_line} missing from model.L: {L}"


def test_mp_dcopf_add_opf_does_not_raise(lv_net):
    """``DCOPF_multi_period`` must construct AND ``add_OPF`` without raising.

    Regression for two pre-existing bugs (both fixed together):

    * ``DC_multi_period.__init__`` did not accept the ``pf`` keyword, so
      ``DCOPF_multi_period.__init__`` (which forwards 4 args) crashed before
      the model was even built.
    * ``Sgens_multi_period.get_opf_parameters`` referenced
      ``self.QsGmax_tuple``, which is created only by
      ``static_generation_reactive_power_limits`` (an AC-only step). Calling
      it from the DC OPF path therefore raised an ``AttributeError``. Q
      parameters are now in a separate ``get_acopf_parameters`` hook.
    """

    dcopf = DCOPF_multi_period(lv_net, toT=4)
    dcopf.add_OPF()
    n_l = len(list(dcopf.model.L))
    n_t = len(list(dcopf.model.T))
    assert len(list(dcopf.model.line_lim_upper)) == n_l * n_t
    assert len(list(dcopf.model.line_lim_lower)) == n_l * n_t


def test_mp_acopf_add_opf_does_not_raise(lv_net):
    """``ACOPF_multi_period`` must still construct AND ``add_OPF`` after the
    sgens Q-parameter split (regression check that we didn't break the AC
    path while fixing DC)."""
    opf = ACOPF_multi_period(lv_net, toT=4)
    opf.add_OPF()
    # AC must keep its reactive-power bounds
    assert hasattr(opf.model, "QsGmax")
    assert hasattr(opf.model, "QsGmin")


def test_mp_static_generation_may_be_negative(lv_net):
    """``psG`` must not carry a non-negative domain.

    The domain would override ``net.sgen.min_p_mw`` wherever that bound is
    negative, exactly as it did in the single-period model (see
    ``test_sgen_negative_limits``). A missing or ``NaN`` bound still means 0.
    """
    import copy

    import pyomo.environ as pyo

    net = copy.deepcopy(lv_net)
    net.sgen["min_p_mw"] = -0.5
    net.sgen["max_p_mw"] = net.sgen.p_mw
    net.sgen["controllable"] = True
    opf = DCOPF_multi_period(net, toT=2)
    opf.add_OPF()
    index = list(opf.model.psG)[0]
    assert opf.model.psG[index].domain is pyo.Reals
    assert pyo.value(opf.model.sPGmin[index]) < 0
