"""Unit and integration tests for the ACOPF model."""

import pytest
import pandapower as pp
import pyomo.environ as pyo

from potpourri.models.ACOPF_base import ACOPF


@pytest.fixture()
def four_bus():
    return pp.networks.simple_four_bus_system()


# ── model structure (no solver required) ─────────────────────────────────────


def test_acopf_model_structure(four_bus):
    """ACOPF builds expected Pyomo components after add_OPF()."""
    acopf = ACOPF(four_bus)
    acopf.add_OPF()

    assert hasattr(acopf.model, "B"), "missing bus set"
    assert hasattr(acopf.model, "L"), "missing line set"
    assert hasattr(acopf.model, "Vmax"), "missing Vmax param"
    assert hasattr(acopf.model, "Vmin"), "missing Vmin param"
    assert hasattr(acopf.model, "v_pyo"), "missing voltage bound constraint"
    assert hasattr(acopf.model, "line_lim_from"), "missing line thermal limit"
    assert len(list(acopf.model.B)) > 0
    assert len(list(acopf.model.L)) > 0


def test_acopf_voltage_objective_added(four_bus):
    """add_voltage_deviation_objective creates obj_v_deviation on model."""
    acopf = ACOPF(four_bus)
    acopf.add_OPF()
    acopf.add_voltage_deviation_objective()
    assert hasattr(acopf.model, "obj_v_deviation")


def test_acopf_reactive_objective_added(four_bus):
    """add_reactive_power_flow_objective creates obj_reactive on model."""
    acopf = ACOPF(four_bus)
    acopf.add_OPF()
    acopf.add_reactive_power_flow_objective()
    assert hasattr(acopf.model, "obj_reactive")


def test_acopf_voltage_bounds_active_after_add_opf(four_bus):
    """v_pyo constraint is active after add_OPF()."""
    acopf = ACOPF(four_bus)
    acopf.add_OPF()
    assert acopf.model.v_pyo.active


def test_acopf_line_limits_active_after_add_opf(four_bus):
    """line_lim_from and line_lim_to are active after add_OPF()."""
    acopf = ACOPF(four_bus)
    acopf.add_OPF()
    assert acopf.model.line_lim_from.active
    assert acopf.model.line_lim_to.active


# ── integration: local solver required ───────────────────────────────────────


@pytest.mark.integration
def test_acopf_solves_with_ipopt(four_bus):
    """ACOPF finds an optimal solution using the local IPOPT solver."""
    acopf = ACOPF(four_bus)
    acopf.add_OPF()
    acopf.add_voltage_deviation_objective()
    acopf.solve(solver="ipopt", print_solver_output=False)

    assert pyo.check_optimal_termination(acopf.results)
    assert acopf.net.res_bus is not None
    assert not acopf.net.res_bus.vm_pu.isnull().any()


@pytest.mark.integration
@pytest.mark.slow
def test_acopf_solves_with_neos():
    """ACOPF finds an optimal solution using the NEOS IPOPT solver."""
    import os

    os.environ.setdefault("NEOS_EMAIL", "test@example.com")
    net = pp.networks.simple_four_bus_system()
    acopf = ACOPF(net)
    acopf.add_OPF()
    acopf.add_voltage_deviation_objective()
    acopf.solve(solver="neos", neos_opt="ipopt", print_solver_output=False)

    assert pyo.check_optimal_termination(acopf.results)
    assert acopf.net.res_bus is not None
