"""Unit and integration tests for the DCOPF model."""

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
