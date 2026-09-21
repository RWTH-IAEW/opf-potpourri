# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Tests for the 'gurobi_direct_minlp' solver path.

Pyomo's ``gurobi_direct_minlp`` interface (Pyomo >= 6.10, gurobipy >= 12) maps
Pyomo's ``sin``/``cos``/``exp``/``log``/``sqrt`` onto ``gurobipy.nlfunc``, so
the polar-form AC power flow in :mod:`potpourri.models.AC` can be handed to
Gurobi directly instead of going through IPOPT or MindtPy. The older
``gurobi``, ``gurobi_direct`` and ``gurobi_persistent`` interfaces are limited
to expressions of degree 2 and reject the model.
"""

import pandapower as pp
import pytest
import pyomo.environ as pyo

from potpourri.models.ACOPF_base import ACOPF

pytestmark = pytest.mark.filterwarnings("ignore")


def _minlp_available():
    try:
        return pyo.SolverFactory("gurobi_direct_minlp").available()
    except Exception:
        return False


requires_minlp = pytest.mark.skipif(
    not _minlp_available(), reason="gurobi_direct_minlp not available"
)


def _four_bus_acopf():
    opf = ACOPF(pp.networks.simple_four_bus_system())
    opf.add_OPF()
    opf.add_voltage_deviation_objective()
    return opf


@pytest.mark.integration
@requires_minlp
def test_gurobi_minlp_solves_polar_acopf():
    """Gurobi accepts the sin/cos power flow and reaches optimality."""
    opf = _four_bus_acopf()
    res = opf.solve(solver="gurobi_direct_minlp")
    assert pyo.check_optimal_termination(res)


@pytest.mark.integration
@requires_minlp
def test_gurobi_minlp_option_names_are_translated():
    """max_iter/time_limit must map to Gurobi's own parameter names.

    Passing IPOPT's 'max_iter' straight through raises
    GurobiError("Unknown parameter 'max_iter'").
    """
    opf = _four_bus_acopf()
    res = opf.solve(solver="gurobi_direct_minlp", max_iter=1000, time_limit=30)
    assert pyo.check_optimal_termination(res)


@pytest.mark.integration
@requires_minlp
def test_gurobi_minlp_matches_ipopt_objective():
    """Global solution is at least as good as IPOPT's local one."""
    ipopt = pyo.SolverFactory("ipopt")
    if not ipopt.available():
        pytest.skip("IPOPT not available")

    opf_g = _four_bus_acopf()
    opf_g.solve(solver="gurobi_direct_minlp")
    obj_g = pyo.value(opf_g.model.obj_v_deviation)

    opf_i = _four_bus_acopf()
    opf_i.solve(solver="ipopt")
    obj_i = pyo.value(opf_i.model.obj_v_deviation)

    assert obj_g <= obj_i + 1e-6
