"""Regression tests for GitLab issues #18 and #19.

* #18 — the shunt parameters ``GB`` and ``BB`` were declared with opposite index
  sets, so any network carrying a shunt raised ``KeyError`` on construction.
  Invisible because every SimBench network has zero shunts.
* #19 — chained ``inplace=True`` fills in the result mappers, which pandas turns
  into silent no-ops under copy-on-write, leaving ``NaN`` in ``net.res_*``.
"""

from __future__ import annotations

import ast
import copy
import logging
import warnings
from pathlib import Path

import pandapower as pp
import pyomo.environ as pyo
import pytest
import simbench as sb

from potpourri.models.ACOPF_base import ACOPF
from potpourri.models_multi_period.ACOPF_multi_period import (
    ACOPF_multi_period,
)
from potpourri.models_multi_period.DCOPF_multi_period import (
    DCOPF_multi_period,
)

warnings.filterwarnings("ignore")
logging.getLogger("pyomo.core").setLevel(logging.ERROR)

REPO_ROOT = Path(__file__).resolve().parents[2]
MAPPERS = (
    REPO_ROOT / "src" / "potpourri" / "models" / "pyo_to_net.py",
    REPO_ROOT
    / "src"
    / "potpourri"
    / "models_multi_period"
    / "pyo_to_net_multi_period.py",
)


@pytest.fixture(scope="module")
def lv_net():
    return sb.get_simbench_net("1-LV-rural1--0-sw")


def _with_shunt(lv_net):
    """The LV net plus one shunt — the case no shipped network covers."""
    net = copy.deepcopy(lv_net)
    pp.create_shunt(net, bus=int(net.bus.index[3]), q_mvar=-0.01, p_mw=0.001)
    net.bus["max_vm_pu"] = 1.10
    net.bus["min_vm_pu"] = 0.90
    net.sgen["controllable"] = True
    net.sgen["max_p_mw"] = net.sgen["p_mw"]
    net.sgen["min_p_mw"] = 0.0
    net.ext_grid["max_q_mvar"] = 500.0
    net.ext_grid["min_q_mvar"] = -500.0
    return net


# ---------------------------------------------------------------------------
# #18 — shunt parameter indexing
# ---------------------------------------------------------------------------


def test_issue18_no_shipped_network_has_a_shunt(lv_net):
    """Why this went unnoticed, asserted so the reason stays visible.

    ``for s in model.SHUNT`` iterates zero times on every network in the
    examples and the tests, so neither the bad balance term nor the bad mapper
    term was ever evaluated.
    """
    assert len(lv_net.shunt) == 0


@pytest.mark.parametrize(
    "cls",
    [ACOPF_multi_period, DCOPF_multi_period],
    ids=lambda c: c.__name__,
)
def test_issue18_model_builds_with_a_shunt(lv_net, cls):
    """Construction used to raise KeyError before returning a model."""
    opf = cls(_with_shunt(lv_net), toT=3)
    opf.add_OPF()
    assert list(opf.model.SHUNT) == [0]


def test_issue18_gb_and_bb_share_the_time_index(lv_net):
    """The two disagreed: GB over SHUNT, BB over SHUNT x T."""
    opf = ACOPF_multi_period(_with_shunt(lv_net), toT=3)
    steps = list(opf.model.T)
    expected = {(0, t) for t in steps}
    assert set(opf.model.GB.keys()) == expected
    assert set(opf.model.BB.keys()) == expected


def test_issue18_shunt_results_are_mapped(lv_net):
    """The mapper half, reachable only since #9 wired it in."""
    opf = ACOPF_multi_period(_with_shunt(lv_net), toT=3)
    opf.add_OPF()
    opf.add_voltage_deviation_objective()
    res = opf.solve(solver="ipopt", print_solver_output=False)
    assert pyo.check_optimal_termination(res)

    assert len(opf.net.res_shunt) == 1
    assert opf.net.res_shunt["p_mw"].notna().all()
    assert opf.net.res_shunt["q_mvar"].notna().all()


def test_issue18_shunt_free_networks_still_work(lv_net):
    """The fix must not disturb the case that always worked."""
    net = copy.deepcopy(lv_net)
    net.bus["max_vm_pu"] = 1.10
    net.bus["min_vm_pu"] = 0.90
    opf = ACOPF_multi_period(net, toT=3)
    opf.add_OPF()
    assert list(opf.model.SHUNT) == []


# ---------------------------------------------------------------------------
# #19 — chained inplace assignment in the result mappers
# ---------------------------------------------------------------------------


def _chained_inplace_calls(path):
    """Find ``<frame>.<column>.<method>(..., inplace=True)`` in a module.

    A chained call has an attribute access on an attribute access — the column
    is taken off the frame and then mutated, which operates on a temporary.
    ``frame.method(..., inplace=True)`` is fine and is not reported.
    """
    tree = ast.parse(path.read_text(), filename=str(path))
    found = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        if not any(
            kw.arg == "inplace"
            and isinstance(kw.value, ast.Constant)
            and kw.value.value is True
            for kw in node.keywords
        ):
            continue
        func = node.func
        # func is <something>.method; chained when <something> is itself an
        # attribute access whose own base is an attribute access or a call.
        if isinstance(func, ast.Attribute) and isinstance(
            func.value, ast.Attribute
        ):
            base = func.value.value
            if isinstance(base, (ast.Attribute, ast.Subscript)):
                found.append(f"{path.name}:{node.lineno}")
    return found


def test_issue19_mappers_have_no_chained_inplace_calls():
    """These become silent no-ops under copy-on-write in pandas 3.0."""
    offenders = []
    for path in MAPPERS:
        offenders.extend(_chained_inplace_calls(path))
    assert offenders == [], (
        "chained inplace mutation writes into a temporary; assign the column "
        f"instead: {offenders}"
    )


def test_issue19_pandas_pin_has_an_upper_bound():
    """An unbounded pin is how pandapower 3.5 broke things in #8."""
    assert '"pandas>=2.0,<3"' in (REPO_ROOT / "pyproject.toml").read_text()


@pytest.mark.integration
def test_issue19_multi_period_mapping_raises_no_future_warning(lv_net):
    """The whole mapper path must be clean, not merely warning-suppressed."""
    net = copy.deepcopy(lv_net)
    net.bus["max_vm_pu"] = 1.10
    net.bus["min_vm_pu"] = 0.90
    net.sgen["controllable"] = True
    net.sgen["max_p_mw"] = net.sgen["p_mw"]
    net.sgen["min_p_mw"] = 0.0
    net.ext_grid["max_q_mvar"] = 500.0
    net.ext_grid["min_q_mvar"] = -500.0

    opf = ACOPF_multi_period(net, toT=3)
    opf.add_OPF()
    opf.add_voltage_deviation_objective()
    with warnings.catch_warnings():
        warnings.simplefilter("error", FutureWarning)
        res = opf.solve(solver="ipopt", print_solver_output=False)
    assert pyo.check_optimal_termination(res)


@pytest.mark.integration
def test_issue19_single_period_mapping_raises_no_future_warning(lv_net):
    net = copy.deepcopy(lv_net)
    net.bus["max_vm_pu"] = 1.10
    net.bus["min_vm_pu"] = 0.90
    net.ext_grid["max_q_mvar"] = 500.0
    net.ext_grid["min_q_mvar"] = -500.0

    opf = ACOPF(net)
    opf.add_OPF()
    opf.add_voltage_deviation_objective()
    with warnings.catch_warnings():
        warnings.simplefilter("error", FutureWarning)
        res = opf.solve(solver="ipopt", print_solver_output=False)
    assert pyo.check_optimal_termination(res)


@pytest.mark.integration
def test_issue19_result_tables_have_no_leftover_nan(lv_net):
    """The point of the fills: no NaN should survive into net.res_*."""
    net = copy.deepcopy(lv_net)
    net.bus["max_vm_pu"] = 1.10
    net.bus["min_vm_pu"] = 0.90
    net.sgen["controllable"] = True
    net.sgen["max_p_mw"] = net.sgen["p_mw"]
    net.sgen["min_p_mw"] = 0.0
    net.ext_grid["max_q_mvar"] = 500.0
    net.ext_grid["min_q_mvar"] = -500.0

    opf = ACOPF_multi_period(net, toT=3)
    opf.add_OPF()
    opf.add_voltage_deviation_objective()
    assert pyo.check_optimal_termination(
        opf.solve(solver="ipopt", print_solver_output=False)
    )

    for table in ("res_bus", "res_line", "res_load", "res_sgen", "res_trafo"):
        frame = opf.net[table]
        assert not frame.isna().any().any(), f"{table} still holds NaN"


@pytest.mark.integration
def test_issue19_optimised_sgen_dispatch_reaches_res_sgen(lv_net):
    """Guards the chained ``res_sgen.iloc[g]["p_mw"] = ...`` assignment.

    If that write lands on a copy, the fill below replaces the curtailed value
    with the uncurtailed profile — so a curtailed sgen would report its full
    output and the error would look like a physical result.
    """
    net = copy.deepcopy(lv_net)
    net.bus["max_vm_pu"] = 1.10
    net.bus["min_vm_pu"] = 0.90
    net.sgen["controllable"] = True
    net.sgen["max_p_mw"] = net.sgen["p_mw"]
    net.sgen["min_p_mw"] = 0.0
    net.ext_grid["max_q_mvar"] = 500.0
    net.ext_grid["min_q_mvar"] = -500.0

    peak = 13868
    opf = ACOPF_multi_period(net, toT=peak + 1, fromT=peak)
    opf.add_OPF()
    opf.add_voltage_deviation_objective()
    # Force curtailment well below the profile, so the optimised value and the
    # profile fallback are far apart.
    for g in opf.model.sGc:
        opf.model.psG[g, peak].setub(0.001)
    assert pyo.check_optimal_termination(
        opf.solve(solver="ipopt", print_solver_output=False)
    )

    base = pyo.value(opf.model.baseMVA)
    for g in opf.model.sGc:
        expected = pyo.value(opf.model.psG[g, peak]) * base
        assert opf.net.res_sgen.p_mw.iloc[g] == pytest.approx(
            expected, abs=1e-9
        )
        profile = opf.net.profiles[("sgen", "p_mw")][g].loc[peak]
        assert opf.net.res_sgen.p_mw.iloc[g] < profile, (
            "res_sgen reports the profile, not the curtailed dispatch"
        )
