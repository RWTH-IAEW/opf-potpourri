"""Unit tests for multi-period model construction (no solver required)."""

import pytest
import simbench as sb

from potpourri.models_multi_period.ACOPF_multi_period import ACOPF_multi_period


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
