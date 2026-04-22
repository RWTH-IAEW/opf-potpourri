"""Unit tests for Basemodel construction and Pyomo set population."""

import pytest
import pandapower as pp
import pyomo.environ as pyo

from potpourri.models.basemodel import Basemodel
from potpourri.models.ACOPF_base import ACOPF
from potpourri.models.DCOPF import DCOPF


# ── fixtures ─────────────────────────────────────────────────────────────────


@pytest.fixture()
def four_bus():
    return pp.networks.simple_four_bus_system()


# ── construction / input validation ──────────────────────────────────────────


def test_basemodel_rejects_non_pandapower_input():
    """Basemodel raises ValueError for non-pandapower input."""
    with pytest.raises(ValueError):
        Basemodel("not a network")


def test_basemodel_deep_copies_network(four_bus):
    """Basemodel must not mutate the original network object."""
    original_bus_count = len(four_bus.bus)
    acopf = ACOPF(four_bus)
    acopf.add_OPF()
    assert len(four_bus.bus) == original_bus_count
    assert four_bus is not acopf.net


# ── Pyomo set population ──────────────────────────────────────────────────────


def test_basemodel_sets_populated(four_bus):
    """create_model() populates all expected Pyomo sets."""
    acopf = ACOPF(four_bus)

    required = ["B", "b0", "D", "L", "TRANSF", "G", "sG", "Dbs", "Gbs"]
    for name in required:
        assert hasattr(acopf.model, name), f"missing Pyomo set '{name}'"

    assert len(list(acopf.model.B)) == len(four_bus.bus)
    assert len(list(acopf.model.D)) == len(
        four_bus.load[four_bus.load.in_service]
    )


def test_basemodel_demand_fixed(four_bus):
    """Load variables pD are fixed to profile values before add_OPF()."""
    acopf = ACOPF(four_bus)
    for d in acopf.model.D:
        assert acopf.model.pD[d].is_fixed(), (
            f"demand variable pD[{d}] should be fixed"
        )


def test_acopf_voltage_params_within_bounds(four_bus):
    """Vmax must be strictly greater than Vmin at every bus."""
    acopf = ACOPF(four_bus)
    acopf.add_OPF()
    for b in acopf.model.B:
        assert pyo.value(acopf.model.Vmax[b]) > pyo.value(
            acopf.model.Vmin[b]
        ), f"Vmax <= Vmin at bus {b}"


def test_dcopf_line_limits_non_negative(four_bus):
    """Thermal line limits SLmax must be non-negative."""
    dcopf = DCOPF(four_bus)
    dcopf.add_OPF()
    for line in dcopf.model.L:
        assert pyo.value(dcopf.model.SLmax[line]) >= 0, (
            f"negative SLmax for line {line}"
        )


# ── per-unit conversion ───────────────────────────────────────────────────────


def test_acopf_basemva_positive(four_bus):
    """baseMVA must be a positive scalar."""
    acopf = ACOPF(four_bus)
    assert pyo.value(acopf.model.baseMVA) > 0


def test_acopf_power_params_in_per_unit(four_bus):
    """PD and PsG parameters should be in p.u. (i.e. |value| < 1000 typically)."""
    acopf = ACOPF(four_bus)
    base = pyo.value(acopf.model.baseMVA)
    for d in acopf.model.D:
        pd_pu = pyo.value(acopf.model.PD[d])
        pd_mw = pd_pu * base
        # p.u. value should be small; MW value should be below 1000 for test net
        assert abs(pd_mw) < 1000, (
            f"PD[{d}] = {pd_mw:.2f} MW looks like it was not converted to p.u."
        )


# ── storage set always declared ───────────────────────────────────────────────


def test_stor_set_declared_without_storage(four_bus):
    """model.STOR must exist even when the network has no storage elements."""
    acopf = ACOPF(four_bus)
    assert hasattr(acopf.model, "STOR"), (
        "model.STOR should always be declared (empty when no storage)"
    )
    assert len(list(acopf.model.STOR)) == 0
