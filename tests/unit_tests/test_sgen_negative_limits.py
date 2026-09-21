# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""A static generator whose ``min_p_mw`` is negative must be allowed to consume.

``psG`` used to be declared over the non-negative reals, which overrode the
network's own lower bound without a word: a unit that can take power off the
grid was held at zero, and the OPF returned a dearer dispatch (PGLib
``case588_sdet`` has four such units and came out 1.6 % above the reference).
The default for a missing or ``NaN`` bound is still zero.
"""

import pandapower as pp
import pyomo.environ as pyo
import pytest

from potpourri.models.cost_objective import add_poly_cost_objective
from potpourri.models.DCOPF import DCOPF

pytestmark = pytest.mark.filterwarnings("ignore")


def _two_unit_net(sgen_min_p_mw):
    """Cheap external grid, expensive sgen: the optimum sends the sgen to its
    lower bound and buys the difference from the grid."""
    net = pp.create_empty_network(sn_mva=100.0)
    b0 = pp.create_bus(net, 110.0)
    b1 = pp.create_bus(net, 110.0)
    pp.create_ext_grid(net, b0, min_p_mw=-500.0, max_p_mw=500.0)
    pp.create_line(net, b0, b1, 10.0, "149-AL1/24-ST1A 110.0")
    sg = pp.create_sgen(
        net,
        b1,
        p_mw=20.0,
        min_p_mw=sgen_min_p_mw,
        max_p_mw=80.0,
        controllable=True,
    )
    pp.create_load(net, b1, 30.0, 0.0)
    pp.create_poly_cost(net, 0, "ext_grid", cp1_eur_per_mw=10.0)
    pp.create_poly_cost(net, sg, "sgen", cp1_eur_per_mw=50.0)
    return net, sg


def _solve(net):
    model = DCOPF(net, dc_convention="powermodels")
    model.add_OPF()
    add_poly_cost_objective(model, allow_quadratic=True)
    res = model.solve(solver="ipopt", print_solver_output=False, to_net=False)
    assert pyo.check_optimal_termination(res)
    return model


@pytest.mark.integration
def test_negative_sgen_bound_is_honoured():
    net, sg = _two_unit_net(-40.0)
    model = _solve(net)
    p_mw = pyo.value(model.model.psG[sg]) * 100.0
    assert p_mw == pytest.approx(-40.0, abs=1e-4)


@pytest.mark.integration
def test_default_lower_bound_is_still_zero():
    net, sg = _two_unit_net(float("nan"))
    model = _solve(net)
    p_mw = pyo.value(model.model.psG[sg]) * 100.0
    assert p_mw == pytest.approx(0.0, abs=1e-4)


def test_psg_domain_allows_negative_values():
    """Solver-free: the variable itself must not carry a non-negative domain."""
    net, _ = _two_unit_net(-40.0)
    model = DCOPF(net, dc_convention="powermodels")
    model.add_OPF()
    for index in model.model.psG:
        assert model.model.psG[index].domain is pyo.Reals
