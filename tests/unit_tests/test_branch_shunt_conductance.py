"""The AC branch model must include pandapower's BR_G (iron losses)."""

import numpy as np
import pandapower as pp
import pyomo.environ as pyo
import pytest

from potpourri.models.ACOPF_base import ACOPF
from potpourri.models.basemodel import branch_charging_admittance


def _trafo_feeder():
    net = pp.create_empty_network(sn_mva=1.0)
    hv = pp.create_bus(net, 20.0, max_vm_pu=1.1, min_vm_pu=0.9)
    lv = pp.create_bus(net, 0.4, max_vm_pu=1.1, min_vm_pu=0.9)
    lv2 = pp.create_bus(net, 0.4, max_vm_pu=1.1, min_vm_pu=0.9)
    pp.create_ext_grid(net, hv, vm_pu=1.02)
    pp.create_transformer(net, hv, lv, "0.4 MVA 20/0.4 kV")
    pp.create_line(net, lv, lv2, 0.3, "NAYY 4x150 SE")
    pp.create_load(net, lv2, 0.15, 0.03)
    pp.create_sgen(net, lv2, 0.05, 0.0, sn_mva=0.06, controllable=False)
    return net


def test_charging_admittance_reads_br_g_column():
    net = _trafo_feeder()
    pp.runpp(net)
    y = branch_charging_admittance(net._ppc["branch"])
    n_line = len(net.line)
    assert y[n_line].real > 0.0  # trafo iron losses present
    assert y[:n_line].real.max() == 0.0  # lines carry no conductance
    np.testing.assert_allclose(y.imag, net._ppc["branch"][:, 4].real)
    # a MATPOWER-style array without BR_G yields a purely imaginary result
    y2 = branch_charging_admittance(net._ppc["branch"][:, :13])
    assert np.all(y2.real == 0.0)


@pytest.mark.slow
def test_acopf_fixed_dispatch_reproduces_pandapower_trafo_losses():
    net = _trafo_feeder()
    pp.runpp(net, voltage_depend_loads=False)
    opf = ACOPF(net)
    opf.add_OPF(thermal_limit="mva", free_slack_vm=False)
    opf.model.obj = pyo.Objective(expr=0.0)
    res = opf.solve(solver="ipopt", to_net=True, print_solver_output=False)
    assert pyo.check_optimal_termination(res)
    assert opf.net.res_trafo.pl_mw.iloc[0] == pytest.approx(
        net.res_trafo.pl_mw.iloc[0], abs=1e-7
    )
    assert opf.net.res_ext_grid.p_mw.iloc[0] == pytest.approx(
        net.res_ext_grid.p_mw.iloc[0], abs=1e-7
    )
    np.testing.assert_allclose(
        opf.net.res_bus.vm_pu.values, net.res_bus.vm_pu.values, atol=1e-7
    )
