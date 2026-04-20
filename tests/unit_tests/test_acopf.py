import os
os.environ['NEOS_EMAIL'] = "test-email@test.com"
import pandapower as pp
import pyomo.environ as pyo

from potpourri.models.ACOPF_base import ACOPF


def test_acopf_model_structure():
    """ACOPF builds the expected Pyomo components after add_OPF() without solving."""
    net = pp.networks.simple_four_bus_system()
    acopf = ACOPF(net)
    acopf.add_OPF()

    assert hasattr(acopf.model, 'B'), "missing bus set"
    assert hasattr(acopf.model, 'L'), "missing line set"
    assert hasattr(acopf.model, 'Vmax'), "missing voltage upper bound param"
    assert hasattr(acopf.model, 'Vmin'), "missing voltage lower bound param"
    assert hasattr(acopf.model, 'v_pyo'), "missing voltage bound constraint"
    assert hasattr(acopf.model, 'line_lim_from'), "missing line thermal limit"
    assert len(list(acopf.model.B)) > 0
    assert len(list(acopf.model.L)) > 0


def test_acopf_solves_with_ipopt():
    """ACOPF finds an optimal solution using the local IPOPT solver."""
    net = pp.networks.simple_four_bus_system()
    acopf = ACOPF(net)
    acopf.add_OPF()
    acopf.add_voltage_deviation_objective()
    acopf.solve(solver='neos', neos_opt='ipopt')

    assert pyo.check_optimal_termination(acopf.results)
    assert acopf.net.res_bus is not None
    assert not acopf.net.res_bus.vm_pu.isnull().any()


if __name__ == '__main__':
    test_acopf_model_structure()
    test_acopf_solves_with_ipopt()