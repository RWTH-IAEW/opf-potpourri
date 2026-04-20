import os

os.environ["NEOS_EMAIL"] = "test-email@test.com"
import pandapower as pp
import pyomo.environ as pyo

from potpourri.models.DCOPF import DCOPF


def test_dcopf_model_structure():
    """DCOPF builds the expected Pyomo components after add_OPF() without solving."""
    net = pp.networks.simple_four_bus_system()
    dcopf = DCOPF(net)
    dcopf.add_OPF()

    assert hasattr(dcopf.model, "B"), "missing bus set"
    assert hasattr(dcopf.model, "L"), "missing line set"
    assert hasattr(dcopf.model, "SLmax"), "missing line thermal limit param"
    assert hasattr(dcopf.model, "line_lim_from"), (
        "missing line upper limit constraint"
    )
    assert hasattr(dcopf.model, "line_lim_to"), (
        "missing line lower limit constraint"
    )
    assert hasattr(dcopf.model, "PG_Constraint"), (
        "missing generator power limit"
    )
    assert len(list(dcopf.model.B)) > 0
    assert len(list(dcopf.model.L)) > 0


def test_dcopf_solves_with_cplex():
    """DCOPF finds an optimal solution using the NEOS CPLEX solver."""
    net = pp.networks.simple_four_bus_system()
    dcopf = DCOPF(net)
    dcopf.add_OPF()

    dcopf.model.obj = pyo.Objective(
        expr=sum(dcopf.model.pG[g] for g in dcopf.model.G),
        sense=pyo.minimize,
    )
    dcopf.solve(solver="neos", neos_opt="cplex", to_net=False)

    assert pyo.check_optimal_termination(dcopf.results)


if __name__ == "__main__":
    test_dcopf_model_structure()
    test_dcopf_solves_with_cplex()
