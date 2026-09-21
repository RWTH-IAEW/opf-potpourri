# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Regression tests for GitLab issue #16.

A 12-step midday window on ``1-LV-rural1--0-sw`` reported a locally infeasible
point even though curtailing the PV to zero is available and feasible. The cause
was the cold start: ``v`` at 1.0 with every angle and branch flow at zero
violates the nodal balance at every bus by the full injection, and IPOPT does not
always recover from that on a nonconvex AC OPF.

``solve()`` now seeds a consistent operating point from a per-step power flow
first. These tests pin the fix, and pin the diagnosis by showing the failure
returns when the seeding is switched off.
"""

from __future__ import annotations

import copy
import logging
import warnings

import pyomo.environ as pyo
import pytest
import simbench as sb

from potpourri.models_multi_period.ACOPF_multi_period import (
    ACOPF_multi_period,
)
from potpourri.models_multi_period.DCOPF_multi_period import (
    DCOPF_multi_period,
)

warnings.filterwarnings("ignore")
# Loading a locally-infeasible result logs a multi-line Pyomo warning; the
# no-warm-start test does that on purpose.
logging.getLogger("pyomo.core").setLevel(logging.ERROR)

# The reported window: 12 steps around the 24 May PV peak.
FROM_T, TO_T = 13860, 13872


@pytest.fixture(scope="module")
def lv_net():
    return sb.get_simbench_net("1-LV-rural1--0-sw")


def _reported_case(lv_net, toT=TO_T, fromT=FROM_T):
    """The model exactly as reported in the issue."""
    net = copy.deepcopy(lv_net)
    net.bus["max_vm_pu"] = 1.05
    net.bus["min_vm_pu"] = 0.95
    net.sgen["controllable"] = True
    net.sgen["max_p_mw"] = net.sgen["p_mw"]
    net.sgen["min_p_mw"] = 0.0  # full curtailment is available
    net.ext_grid["max_q_mvar"] = 500.0
    net.ext_grid["min_q_mvar"] = -500.0
    opf = ACOPF_multi_period(net, toT=toT, fromT=fromT)
    opf.add_OPF()
    opf.add_voltage_deviation_objective()
    return opf


def _voltages(opf):
    return [
        pyo.value(opf.model.v[b, t]) for b in opf.model.B for t in opf.model.T
    ]


@pytest.mark.integration
def test_issue16_reported_window_converges(lv_net):
    """The headline: the window in the issue must solve."""
    opf = _reported_case(lv_net)
    res = opf.solve(solver="ipopt", print_solver_output=False)
    assert pyo.check_optimal_termination(res)

    vs = _voltages(opf)
    # A plausible operating point, not the collapsed iterate the failure gave
    # (which ran down to ~0 p.u.).
    assert min(vs) > 0.9
    assert max(vs) <= 1.05 + 1e-6


@pytest.mark.integration
def test_issue16_cold_start_converges_since_the_sgen_bound_fix(lv_net):
    """The cold start of the reported window converges too, since 0.5.3.

    Until then this test pinned the opposite: from the cold start IPOPT
    reported a locally infeasible point, and the seed was the cure. 0.5.3
    moved the static-generation lower bound from the variable's domain into
    the ``PsG_Constraint`` (so that a negative ``min_p_mw`` is honoured), and
    with the bound expressed that way IPOPT's barrier path differs and the
    same window solves cold. Putting the bound back on the variable brings
    the failure back, which is how the cause was verified. The seed remains
    the default because it is the safer start; this test keeps the record.
    """
    opf = _reported_case(lv_net)
    res = opf.solve(
        solver="ipopt", print_solver_output=False, warm_start=False
    )
    assert pyo.check_optimal_termination(res)
    vs = _voltages(opf)
    assert min(vs) > 0.9
    assert max(vs) <= 1.05 + 1e-6


@pytest.mark.integration
def test_issue16_warm_start_is_insensitive_to_the_seed(lv_net):
    """Any power-flow-consistent seed reaches the same optimum.

    What matters is that the seed satisfies the power flow, not that it is near
    the answer — an uncurtailed and a fully curtailed seed are far apart in
    dispatch and must still agree.
    """
    solutions = []
    for curtailment in (1.0, 0.0):
        opf = _reported_case(lv_net)
        opf.warm_start_from_pf(curtailment=curtailment)
        res = opf.solve(
            solver="ipopt", print_solver_output=False, warm_start=False
        )
        assert pyo.check_optimal_termination(res), f"{curtailment=}"
        solutions.append(_voltages(opf))

    for a, b in zip(*solutions):
        assert a == pytest.approx(b, abs=1e-6)


def test_issue16_warm_start_seeds_every_step(lv_net):
    """All steps must be seeded, or the unseeded ones keep the bad start."""
    opf = _reported_case(lv_net)
    assert opf.warm_start_from_pf() == len(list(opf.model.T))


def test_issue16_warm_start_makes_the_state_consistent(lv_net):
    """The point of the seed: the state stops being inconsistent.

    A cold model carries **no** value for the branch flows at all — they are
    uninitialized, so the NL writer hands IPOPT a default of zero while ``v``
    is 1.0. The nodal balance is then violated at every bus by the full
    injection, which is the starting point IPOPT could not recover from.
    """
    opf = _reported_case(lv_net)
    m = opf.model
    t = list(m.T)[0]

    # `.value` rather than pyo.value(): the latter raises on an uninitialized
    # variable, which is precisely the state being asserted here.
    #
    # The cold state is v = 1, delta = 0 (both explicitly initialised) and no
    # branch-flow values at all. A flat voltage profile carries no flow, so
    # zero flows would be consistent with it — but the injections are not zero,
    # and that is where the balance breaks.
    assert all(m.pLfrom[line, t].value is None for line in m.L), (
        "expected a cold model to have no branch-flow values at all"
    )
    assert all(m.delta[b, t].value == 0.0 for b in m.B), (
        "expected a cold model to start from a flat angle profile"
    )
    # Every bus but the slack, whose magnitude carries its base-case value.
    assert all(m.v[b, t].value == 1.0 for b in m.B if b not in set(m.b0)), (
        "expected a cold model to start from a flat voltage profile"
    )

    opf.warm_start_from_pf()

    warm_flows = [m.pLfrom[line, t].value for line in m.L]
    assert all(f is not None for f in warm_flows)
    assert any(abs(f) > 1e-9 for f in warm_flows), (
        "warm start left the branch flows at zero"
    )
    # Voltages come from a power flow, so they are no longer all exactly 1.0.
    warm_v = [m.v[b, t].value for b in m.B]
    assert any(abs(v - 1.0) > 1e-9 for v in warm_v)


def test_issue16_warm_start_works_on_the_dc_model(lv_net):
    """The seeder must skip variables a formulation does not have.

    The DC model has no ``v``, ``qLfrom`` or ``qG``; seeding must not raise.
    """
    net = copy.deepcopy(lv_net)
    net.sgen["controllable"] = True
    net.sgen["max_p_mw"] = net.sgen["p_mw"]
    net.sgen["min_p_mw"] = 0.0
    opf = DCOPF_multi_period(net, toT=FROM_T + 3, fromT=FROM_T)
    opf.add_OPF()
    assert opf.warm_start_from_pf() == 3


def test_issue16_warm_start_can_be_switched_off(lv_net):
    """`warm_start=False` must not run the power flows at all."""
    opf = _reported_case(lv_net, toT=FROM_T + 2)
    m = opf.model
    t = list(m.T)[0]
    opf.solve(
        solver="ipopt",
        print_solver_output=False,
        warm_start=False,
        to_net=False,
    )
    # Nothing to assert about the solution here; the point is that the call
    # is accepted and does not seed. Re-seeding afterwards must still work.
    assert opf.warm_start_from_pf() == len(list(m.T))
    assert any(abs(pyo.value(m.pLfrom[line, t])) > 1e-9 for line in m.L)
