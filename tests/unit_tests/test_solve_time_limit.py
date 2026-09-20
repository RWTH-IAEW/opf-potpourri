"""``solve(time_limit=...)`` must reach the solver.

It was honoured by Gurobi and mindtpy only and silently dropped for IPOPT, so a
benchmark asking for a one-hour cap saw solves run for hours. Solver-free: the
SolverFactory is replaced by a stub that records the options it was given.
"""

import pandapower as pp
import pyomo.environ as pyo
import pytest

import potpourri.models.basemodel as bm
from potpourri.models.ACOPF_base import ACOPF

pytestmark = pytest.mark.filterwarnings("ignore")


class _Results:
    class solver:
        status = pyo.SolverStatus.ok
        termination_condition = pyo.TerminationCondition.optimal


class _Optimizer:
    def __init__(self, name):
        self.name = name
        self.options = {}

    def solve(self, model, **kwargs):
        return _Results()


@pytest.fixture
def created(monkeypatch):
    made = []

    def factory(name):
        opt = _Optimizer(name)
        made.append(opt)
        return opt

    monkeypatch.setattr(bm.pyo, "SolverFactory", factory)
    return made


def _opf():
    opf = ACOPF(pp.networks.simple_four_bus_system())
    opf.add_OPF()
    return opf


def test_time_limit_reaches_ipopt_as_max_wall_time(created):
    _opf().solve(solver="ipopt", time_limit=42, to_net=False)
    assert created[-1].options == {"max_wall_time": 42.0}


def test_ipopt_gets_no_limit_unless_asked(created):
    _opf().solve(solver="ipopt", to_net=False)
    assert created[-1].options == {}
    _opf().solve(solver="ipopt", max_iter=7, to_net=False)
    assert created[-1].options == {"max_iter": 7}


def test_gurobi_keeps_its_default_and_own_option_names(created):
    _opf().solve(solver="gurobi_direct_minlp", to_net=False)
    assert created[-1].options == {"TimeLimit": 600}
    _opf().solve(
        solver="gurobi_direct_minlp", time_limit=30, max_iter=5, to_net=False
    )
    assert created[-1].options == {"TimeLimit": 30, "IterationLimit": 5}
