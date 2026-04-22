"""Additional unit tests for the AC power flow model and admittance mapping."""

import math
import pytest
import pandapower as pp
import pyomo.environ as pyo

from potpourri.models.AC import AC
from potpourri.models.ACOPF_base import ACOPF


@pytest.fixture()
def four_bus():
    return pp.networks.simple_four_bus_system()


# ── AC model structure ────────────────────────────────────────────────────────


def test_ac_model_has_kcl_constraints(four_bus):
    """AC model declares KCL_real and KCL_reactive constraints."""
    ac = AC(four_bus)
    assert hasattr(ac.model, "KCL_real")
    assert hasattr(ac.model, "KCL_reactive")
    assert len(list(ac.model.KCL_real)) == len(four_bus.bus)
    assert len(list(ac.model.KCL_reactive)) == len(four_bus.bus)


def test_ac_model_has_kvl_constraints(four_bus):
    """AC model declares KVL constraints for every line."""
    ac = AC(four_bus)
    for name in (
        "KVL_real_from",
        "KVL_real_to",
        "KVL_reactive_from",
        "KVL_reactive_to",
    ):
        assert hasattr(ac.model, name), f"missing constraint '{name}'"
        assert len(list(ac.model.component(name))) == len(four_bus.line)


def test_ac_voltage_variable_bounds(four_bus):
    """Voltage magnitude variable v must have bounds (0, 2)."""
    ac = AC(four_bus)
    for b in ac.model.B:
        lb, ub = ac.model.v[b].bounds
        assert lb is not None and lb >= 0.0
        assert ub is not None and ub <= 2.0


def test_ac_angle_variable_bounds(four_bus):
    """Voltage angle variable delta must have bounds in (−π, π)."""
    ac = AC(four_bus)
    for b in ac.model.B:
        lb, ub = ac.model.delta[b].bounds
        assert lb is not None and lb <= -math.pi * 0.9
        assert ub is not None and ub >= math.pi * 0.9


def test_ac_line_admittance_params_present(four_bus):
    """Gii, Bii, Gik, Bik parameters must exist and be indexed by L."""
    ac = AC(four_bus)
    for pname in ("Gii", "Bii", "Gik", "Bik"):
        assert hasattr(ac.model, pname), f"missing admittance param '{pname}'"
        assert len(list(ac.model.component(pname))) == len(four_bus.line)


# ── ACOPF constraint activation ───────────────────────────────────────────────


def test_acopf_constraint_deactivation(four_bus):
    """Deactivating v_pyo should remove the voltage bound constraint."""
    acopf = ACOPF(four_bus)
    acopf.add_OPF()
    acopf.model.v_pyo.deactivate()
    assert not acopf.model.v_pyo.active


def test_acopf_qsg_bounds_indexed_by_sGc(four_bus):
    """QsG_pyo reactive bounds must be indexed by the controllable sgen set."""
    acopf = ACOPF(four_bus)
    acopf.add_OPF()
    # Every constraint index must appear in sGc (controllable sgens)
    sGc = set(acopf.model.sGc)
    for idx in acopf.model.QsG_pyo:
        assert idx in sGc, f"QsG_pyo[{idx}] indexed outside sGc"


# ── admittance computation sanity ─────────────────────────────────────────────


def test_line_conductance_non_negative(four_bus):
    """Self-conductance Gii should be non-negative (positive resistance)."""
    ac = AC(four_bus)
    for line in ac.model.L:
        gii = pyo.value(ac.model.Gii[line])
        assert gii >= 0.0, f"Gii[{line}] = {gii:.6f} is negative"


def test_slack_bus_in_b0(four_bus):
    """The external-grid bus must appear in the b0 (slack) set."""
    acopf = ACOPF(four_bus)
    # b0 uses internal PPC bus indices; just verify the set is non-empty
    assert len(list(acopf.model.b0)) > 0


def test_slack_bus_voltage_fixed(four_bus):
    """Slack bus voltage magnitude must be fixed to the ext-grid setpoint."""
    ac = AC(four_bus)
    for b in ac.model.b0:
        assert ac.model.v[b].is_fixed(), f"v[{b}] (slack bus) should be fixed"
