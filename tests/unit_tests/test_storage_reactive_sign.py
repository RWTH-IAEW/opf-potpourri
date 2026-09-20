"""Storage reactive power must follow pandapower's load convention in the
AC balance (positive q_mvar = consumption), like its active power.

Before the fix ``KCL_reactive`` added ``qSTOR`` on the generation side while
``KCL_real`` subtracted ``pSTOR``, so an AC-OPF state with non-zero storage Q
violated the exact reactive balance when realised with ``pp.runpp``
(0.41 Mvar on SimBench ``1-MV-rural--2-sw``).
"""

import numpy as np
import pandapower as pp
import pyomo.environ as pyo
import pytest

from potpourri.models.ACOPF_base import ACOPF


def _feeder_with_storage():
    net = pp.create_empty_network(sn_mva=1.0)
    hv = pp.create_bus(net, 20.0, max_vm_pu=1.1, min_vm_pu=0.9)
    lv = pp.create_bus(net, 0.4, max_vm_pu=1.1, min_vm_pu=0.9)
    lv2 = pp.create_bus(net, 0.4, max_vm_pu=1.1, min_vm_pu=0.9)
    pp.create_ext_grid(net, hv, vm_pu=1.02)
    pp.create_transformer(net, hv, lv, "0.4 MVA 20/0.4 kV")
    pp.create_line(net, lv, lv2, 0.4, "NAYY 4x150 SE")
    pp.create_load(net, lv2, 0.12, 0.03)
    # charging 30 kW and consuming 40 kvar (both positive = consumption)
    pp.create_storage(
        net,
        lv2,
        p_mw=0.03,
        q_mvar=0.04,
        max_e_mwh=0.2,
        soc_percent=50.0,
        sn_mva=0.06,
        efficiency_percent=95.0,
    )
    return net


def _fix_storage(model, p, q):
    for s in model.STOR:
        model.STOR_Pchg[s].fix(max(p, 0.0))
        model.STOR_Pdis[s].fix(max(-p, 0.0))
        model.qSTOR[s].fix(q)


@pytest.mark.slow
def test_acopf_with_pinned_storage_q_reproduces_pandapower():
    net = _feeder_with_storage()
    pp.runpp(net, voltage_depend_loads=False)
    opf = ACOPF(net)
    opf.add_OPF(thermal_limit="mva", free_slack_vm=False)
    _fix_storage(opf.model, 0.03, 0.04)
    opf.model.obj = pyo.Objective(expr=0.0)
    res = opf.solve(solver="ipopt", to_net=True, print_solver_output=False)
    assert pyo.check_optimal_termination(res)
    np.testing.assert_allclose(
        opf.net.res_bus.vm_pu.values, net.res_bus.vm_pu.values, atol=1e-6
    )
    assert opf.net.res_ext_grid.q_mvar.iloc[0] == pytest.approx(
        net.res_ext_grid.q_mvar.iloc[0], abs=1e-6
    )
