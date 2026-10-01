# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Opt-in on-load tap changer control (`potpourri.models.oltc`).

Three things are pinned down here.

1. **Nothing changes unless asked.** A model built without `enable_oltc`
   has no control components, no free integer variable, and transformer
   equations that evaluate to the MATPOWER/PowerModels form with the
   pandapower ratio on the HV side.
2. **pandapower is the oracle.** For HV- and LV-side tap changers, a
   transformer whose rated voltages differ from the bus voltages, a 150°
   vector-group shift and iron/magnetising losses, the model built at the
   neutral position and moved to another one must reproduce `pp.runpp`
   at that position — residuals of the transformer equations at
   pandapower's solution first (no solver), then a full IPOPT solve
   compared field by field.
3. **The decisions are right.** An under-voltage case moves the tap the
   right way on each side, a PV case curtails less with the tap free, a
   priced tap stays put, limits hold, the continuous relaxation bounds the
   discrete solution, and a discrete schedule is reproducible with
   `pp.runpp` after `apply_tap_positions`. The multi-period model adds
   movement and operation limits.

IPOPT suffices for everything except the tests marked with
`requires_minlp`, which use Gurobi's global MINLP and skip without it.
"""

import copy
import math

import numpy as np
import pandapower as pp
import pyomo.environ as pyo
import pytest
import simbench as sb

from potpourri.models.AC import AC
from potpourri.models.ACOPF_base import ACOPF
from potpourri.models.basemodel import free_integer_variables
from potpourri.models.DCOPF import DCOPF
from potpourri.models.oltc import oltc_eligibility
from potpourri.models_multi_period.ACOPF_multi_period import (
    ACOPF_multi_period,
)
from potpourri.models_multi_period.DCOPF_multi_period import (
    DCOPF_multi_period,
)

# The feeders below use a deliberately wide voltage band so that every tap
# position of the round trip is admissible; the grid-code Q(U) envelope
# warns about such a band, which is irrelevant for networks without a
# Q-controlled static generator.
pytestmark = [
    pytest.mark.filterwarnings("ignore::DeprecationWarning"),
    pytest.mark.filterwarnings(
        "ignore::potpourri.technologies.q_control.EnvelopeRangeWarning"
    ),
]


def _minlp_available():
    try:
        return bool(pyo.SolverFactory("gurobi_direct_minlp").available())
    except Exception:  # noqa: BLE001 - any failure means "not available"
        return False


requires_minlp = pytest.mark.skipif(
    not _minlp_available(), reason="gurobi_direct_minlp not available"
)

# Default feeder: a 40 MVA 110/21 kV unit on 110/20 kV buses (nominal
# mismatch r0 = 20/21), Dyn5 shift, iron and magnetising losses, and a
# ±9 × 2 % tap changer. Every one of those details has bitten a tap model
# somewhere.
TAP_STEP_PERCENT = 2.0


def feeder(
    tap_side="hv",
    tap_pos=0,
    *,
    load_mw=20.0,
    load_mvar=5.0,
    vn_lv_kv=21.0,
    shift_degree=150.0,
    pfe_kw=30.0,
    i0_percent=0.1,
    vm_pu=1.02,
    tap_changer_type="Ratio",
    v_band=(0.5, 1.5),
    **trafo_kwargs,
):
    net = pp.create_empty_network(sn_mva=100.0)
    hv = pp.create_bus(net, 110.0)
    lv = pp.create_bus(net, 20.0)
    pp.create_ext_grid(net, hv, vm_pu=vm_pu)
    pp.create_transformer_from_parameters(
        net,
        hv,
        lv,
        sn_mva=40.0,
        vn_hv_kv=110.0,
        vn_lv_kv=vn_lv_kv,
        vk_percent=12.0,
        vkr_percent=0.5,
        pfe_kw=pfe_kw,
        i0_percent=i0_percent,
        shift_degree=shift_degree,
        tap_side=tap_side,
        tap_neutral=0,
        tap_pos=tap_pos,
        tap_step_percent=TAP_STEP_PERCENT,
        tap_min=-9,
        tap_max=9,
        tap_changer_type=tap_changer_type,
        **trafo_kwargs,
    )
    pp.create_load(net, lv, load_mw, load_mvar)
    net.bus["min_vm_pu"] = v_band[0]
    net.bus["max_vm_pu"] = v_band[1]
    return net


def _controlled(net, mode="discrete", **opf_kwargs):
    """ACOPF with the single transformer under OLTC control."""
    opf = ACOPF(net)
    opf.add_OPF(**opf_kwargs)
    opf.enable_oltc(transformers=[0], mode=mode)
    return opf


def _equality_residuals(model, names):
    """Largest |body − rhs| over the indices of the named constraints."""
    worst = 0.0
    for name in names:
        component = model.component(name)
        if component is None:
            continue
        for idx in component:
            con = component[idx]
            worst = max(worst, abs(pyo.value(con.body) - pyo.value(con.upper)))
    return worst


def _load_pandapower_state(opf, ref, tap_pos):
    """Copy pandapower's solution at `tap_pos` into the OPF variables."""
    model, lookup, base = opf.model, opf.bus_lookup, opf.baseMVA
    for b in ref.bus.index:
        model.v[lookup[b]].set_value(ref.res_bus.vm_pu.at[b])
        model.delta[lookup[b]].set_value(
            math.radians(ref.res_bus.va_degree.at[b])
        )
    model.pThv[0].set_value(ref.res_trafo.p_hv_mw.at[0] / base)
    model.pTlv[0].set_value(ref.res_trafo.p_lv_mw.at[0] / base)
    model.qThv[0].set_value(ref.res_trafo.q_hv_mvar.at[0] / base)
    model.qTlv[0].set_value(ref.res_trafo.q_lv_mvar.at[0] / base)
    data = opf.oltc_setup.data[0]
    model.trafo_tap_position[0].set_value(tap_pos)
    model.trafo_tap_factor[0].set_value(data.factor(tap_pos))
    if data.side == "hv":
        model.Tap[0].set_value(data.ratio_nominal * data.factor(tap_pos))
    else:
        model.Tap_lv[0].set_value(data.factor(tap_pos) / data.factor_base)


# ── eligibility ──────────────────────────────────────────────────────────────


@pytest.mark.parametrize("tap_side", ["hv", "lv"])
def test_eligibility_accepts_a_ratio_changer_on_either_side(tap_side):
    report = oltc_eligibility(feeder(tap_side))
    assert bool(report.at[0, "eligible"])
    assert report.at[0, "reason"] == ""
    assert report.at[0, "n_positions"] == 19


def test_eligibility_rejects_a_transformer_without_tap_changer_type():
    """A None type means pandapower ignores tap_pos; so do we."""
    report = oltc_eligibility(feeder(tap_changer_type=None))
    assert not report.at[0, "eligible"]
    assert "Ratio" in report.at[0, "reason"]


@pytest.mark.parametrize("kind", ["Ideal", "Symmetrical", "Tabular"])
def test_eligibility_rejects_other_tap_changer_types(kind):
    kwargs = {"tap_step_degree": 1.0} if kind != "Tabular" else {}
    report = oltc_eligibility(feeder(tap_changer_type=kind, **kwargs))
    assert not report.at[0, "eligible"]
    assert kind in report.at[0, "reason"]


def test_eligibility_rejects_a_cross_regulator():
    report = oltc_eligibility(feeder(tap_step_degree=30.0))
    assert not report.at[0, "eligible"]
    assert "tap_step_degree" in report.at[0, "reason"]


def test_eligibility_rejects_tabular_characteristics_and_second_changers():
    net = feeder()
    net.trafo["tap_dependency_table"] = True
    assert not oltc_eligibility(net).at[0, "eligible"]
    net = feeder()
    net.trafo["tap2_pos"] = 0.0
    report = oltc_eligibility(net)
    assert not report.at[0, "eligible"]
    assert "second tap changer" in report.at[0, "reason"]


def test_eligibility_rejects_incomplete_or_inconsistent_tap_data():
    net = feeder()
    net.trafo["tap_min"] = np.nan
    assert "tap_min" in oltc_eligibility(net).at[0, "reason"]
    net = feeder(tap_pos=12)
    assert "outside" in oltc_eligibility(net).at[0, "reason"]
    net = feeder()
    net.trafo["tap_step_percent"] = 0.0
    assert "zero" in oltc_eligibility(net).at[0, "reason"]
    net = feeder()
    net.trafo["in_service"] = False
    assert "out of service" in oltc_eligibility(net).at[0, "reason"]


def test_eligibility_reads_pre_3_0_tap_phase_shifter_data():
    """A pandapower 2.x table has no type column but a phase-shifter flag."""
    net = feeder()
    net.trafo = net.trafo.drop(columns=["tap_changer_type"])
    net.trafo["tap_phase_shifter"] = False
    assert oltc_eligibility(net).at[0, "eligible"]
    net.trafo["tap_phase_shifter"] = True
    report = oltc_eligibility(net)
    assert not report.at[0, "eligible"]
    assert "Ideal" in report.at[0, "reason"]


def test_simbench_transformers_need_the_type_set_first(lv_rural_net):
    """SimBench delivers tap data but no tap_changer_type."""
    assert not oltc_eligibility(lv_rural_net)["eligible"].any()
    net = copy.deepcopy(lv_rural_net)
    net.trafo["tap_changer_type"] = "Ratio"
    assert oltc_eligibility(net)["eligible"].all()


def test_enable_oltc_rejects_an_ineligible_explicit_selection():
    opf = ACOPF(feeder(tap_changer_type=None))
    opf.add_OPF()
    with pytest.raises(ValueError, match="Ratio"):
        opf.enable_oltc(transformers=[0])
    with pytest.raises(ValueError, match="not a transformer"):
        opf.enable_oltc(transformers=[7])


def test_enable_oltc_with_nothing_eligible_explains_why():
    opf = ACOPF(feeder(tap_changer_type=None))
    opf.add_OPF()
    with pytest.raises(ValueError, match="No transformer is eligible"):
        opf.enable_oltc()


def test_enable_oltc_rejects_bad_mode_and_a_second_call():
    opf = ACOPF(feeder())
    opf.add_OPF()
    with pytest.raises(ValueError, match="mode"):
        opf.enable_oltc(mode="integer")
    opf.enable_oltc()
    with pytest.raises(RuntimeError, match="already"):
        opf.enable_oltc()


def test_enable_oltc_is_refused_on_dc_models():
    dc = DCOPF(feeder())
    dc.add_OPF()
    with pytest.raises(NotImplementedError):
        dc.enable_oltc()


# ── default behaviour unchanged ──────────────────────────────────────────────


@pytest.mark.parametrize("tap_side", ["hv", "lv"])
def test_default_model_has_no_controllable_tap(tap_side):
    opf = ACOPF(feeder(tap_side, 3))
    opf.add_OPF()
    model = opf.model
    assert not hasattr(model, "TRANSF_OLTC")
    assert not hasattr(model, "trafo_tap_position")
    assert model.Tap[0].fixed and model.Tap_lv[0].fixed
    assert pyo.value(model.Tap_lv[0]) == 1.0
    assert free_integer_variables(model) == []
    # the HV-side ratio is pandapower's ppc ratio, including the nominal
    # mismatch (20/21) and the tap factor on the tapped side
    n = 1 + 3 * TAP_STEP_PERCENT / 100
    r0 = (110.0 / 21.0) / (110.0 / 20.0)
    expected = r0 * n if tap_side == "hv" else r0 / n
    assert pyo.value(model.Tap[0]) == pytest.approx(expected, rel=1e-12)


@pytest.mark.parametrize("tap_side", ["hv", "lv"])
def test_default_transformer_equations_are_the_matpower_form(tap_side):
    """At any voltages and angles the four rules equal the textbook form.

    With `Tap_lv` fixed at 1 the equations must be MATPOWER's: tap on the
    from side only, i.e. `1/tau²` on the HV self term, `1/tau` on the
    mutual terms, nothing on the LV self term.
    """
    ac = AC(feeder(tap_side, -4))
    m = ac.model
    rng = np.random.default_rng(1)
    v = {b: 0.9 + 0.2 * rng.random() for b in m.B}
    d = {b: 0.3 * rng.random() for b in m.B}
    for b in m.B:
        m.v[b].set_value(v[b])
        m.delta[b].set_value(d[b])
    i, k = m.AT[0, 1], m.AT[0, 2]
    tau = pyo.value(m.Tap[0])
    phi = pyo.value(m.shift[0])
    Gii, Bii, Gik, Bik = (
        pyo.value(m.GiiT[0]),
        pyo.value(m.BiiT[0]),
        pyo.value(m.GikT[0]),
        pyo.value(m.BikT[0]),
    )
    th = d[i] - d[k] - phi
    p_hv = Gii / tau**2 * v[i] ** 2 + v[i] * v[k] / tau * (
        Gik * math.cos(th) + Bik * math.sin(th)
    )
    q_hv = -Bii / tau**2 * v[i] ** 2 + v[i] * v[k] / tau * (
        Gik * math.sin(th) - Bik * math.cos(th)
    )
    p_lv = Gii * v[k] ** 2 + v[i] * v[k] / tau * (
        Gik * math.cos(-th) + Bik * math.sin(-th)
    )
    q_lv = -Bii * v[k] ** 2 + v[i] * v[k] / tau * (
        Gik * math.sin(-th) - Bik * math.cos(-th)
    )
    m.pThv[0].set_value(p_hv)
    m.qThv[0].set_value(q_hv)
    m.pTlv[0].set_value(p_lv)
    m.qTlv[0].set_value(q_lv)
    names = (
        "KVL_real_fromTransf",
        "KVL_real_toTransf",
        "KVL_reactive_fromTransf",
        "KVL_reactive_toTransf",
    )
    assert _equality_residuals(m, names) < 1e-12


# ── pandapower round trip ────────────────────────────────────────────────────

TRANSFORMER_RULES = (
    "KVL_real_fromTransf",
    "KVL_real_toTransf",
    "KVL_reactive_fromTransf",
    "KVL_reactive_toTransf",
    "trafo_tap_factor_def",
    "trafo_tap_ratio_hv_def",
    "trafo_tap_ratio_lv_def",
)


@pytest.mark.parametrize("tap_side", ["hv", "lv"])
@pytest.mark.parametrize("tap_pos", [-9, -4, 0, 5, 9])
def test_controlled_equations_vanish_at_pandapowers_solution(
    tap_side, tap_pos
):
    """Model built at neutral, moved to `tap_pos`: residuals are zero.

    No solver: the transformer equations with the OLTC machinery enabled
    are evaluated at pandapower's own solution for that position.
    """
    ref = feeder(tap_side, tap_pos)
    pp.runpp(ref)
    opf = _controlled(feeder(tap_side, 0))
    _load_pandapower_state(opf, ref, tap_pos)
    assert _equality_residuals(opf.model, TRANSFORMER_RULES) < 1e-10


@pytest.mark.integration
@pytest.mark.parametrize("tap_side", ["hv", "lv"])
@pytest.mark.parametrize("tap_pos", [-9, 0, 5, 9])
def test_fixed_position_solve_reproduces_pandapower(tap_side, tap_pos):
    """Fix the controlled position and solve: every result matches runpp."""
    ref = feeder(tap_side, tap_pos)
    pp.runpp(ref)
    opf = _controlled(feeder(tap_side, 0), free_slack_vm=False)
    opf.model.trafo_tap_position[0].fix(tap_pos)
    opf.add_voltage_deviation_objective()
    results = opf.solve(solver="ipopt")
    assert pyo.check_optimal_termination(results)
    res = opf.net
    assert res.res_bus.vm_pu.values == pytest.approx(
        ref.res_bus.vm_pu.values, abs=1e-6
    )
    assert res.res_bus.va_degree.values == pytest.approx(
        ref.res_bus.va_degree.values, abs=1e-5
    )
    for column in ("p_hv_mw", "q_hv_mvar", "p_lv_mw", "q_lv_mvar", "pl_mw"):
        assert res.res_trafo.at[0, column] == pytest.approx(
            ref.res_trafo.at[0, column], abs=1e-5
        ), column
    assert res.res_trafo.at[0, "loading_percent"] == pytest.approx(
        ref.res_trafo.at[0, "loading_percent"], abs=1e-4
    )
    assert res.res_trafo.at[0, "tap_pos"] == tap_pos
    assert res.res_trafo.at[0, "tap_factor"] == pytest.approx(
        1 + tap_pos * TAP_STEP_PERCENT / 100
    )


# ── optimisation behaviour ───────────────────────────────────────────────────


def _undervoltage_net(tap_side):
    """Heavy load, slack at 1.0 p.u.: the LV bus sits below 0.98."""
    return feeder(
        tap_side,
        0,
        load_mw=30.0,
        vn_lv_kv=20.0,
        shift_degree=0.0,
        vm_pu=1.0,
        v_band=(0.98, 1.06),
    )


@pytest.mark.integration
@pytest.mark.parametrize("tap_side", ["hv", "lv"])
def test_undervoltage_is_cured_by_moving_the_tap_the_right_way(tap_side):
    fixed = ACOPF(_undervoltage_net(tap_side))
    fixed.add_OPF(free_slack_vm=False)
    fixed.add_voltage_deviation_objective()
    assert not pyo.check_optimal_termination(fixed.solve(solver="ipopt"))

    opf = _controlled(
        _undervoltage_net(tap_side), mode="continuous", free_slack_vm=False
    )
    opf.add_voltage_deviation_objective()
    assert pyo.check_optimal_termination(opf.solve(solver="ipopt"))
    k = opf.tap_schedule()[0]
    # fewer HV turns (hv side) or more LV turns (lv side) raise the LV voltage
    assert k < 0 if tap_side == "hv" else k > 0
    assert opf.net.res_bus.vm_pu.at[1] >= 0.98 - 1e-6
    assert -9 <= k <= 9


@pytest.mark.integration
@pytest.mark.parametrize("tap_side", ["hv", "lv"])
def test_round_and_fix_gives_an_integer_position_and_is_reproducible(
    tap_side,
):
    opf = _controlled(_undervoltage_net(tap_side), free_slack_vm=False)
    opf.add_voltage_deviation_objective()
    with pytest.raises(ValueError, match="integer variable"):
        opf.solve(solver="ipopt")
    results = opf.solve_oltc_round_and_fix(solver="ipopt")
    assert pyo.check_optimal_termination(results)
    k = opf.tap_schedule()[0]
    assert isinstance(k, (int, np.integer))
    assert k == (-1 if tap_side == "hv" else 1)
    assert opf.tap_operations()[0] == 1
    relaxed = opf.rounding_info["relaxed"]["trafo_tap_position"][0]
    assert abs(relaxed - k) <= 0.5 + 1e-9
    # the relaxation is a bound on the rounded solution
    assert (
        opf.rounding_info["objective_relaxed"]
        <= pyo.value(opf.model.obj_v_deviation) + 1e-9
    )
    assert opf.model.trafo_tap_position[0].fixed

    # pandapower reproduces the rounded solution once the tap is applied
    check = opf.apply_tap_positions(copy.deepcopy(_undervoltage_net(tap_side)))
    assert check.trafo.at[0, "tap_pos"] == k
    pp.runpp(check)
    assert check.res_bus.vm_pu.at[1] == pytest.approx(
        opf.net.res_bus.vm_pu.at[1], abs=1e-6
    )


@pytest.mark.integration
@requires_minlp
@pytest.mark.parametrize("tap_side", ["hv", "lv"])
def test_global_minlp_agrees_with_the_rounded_solution(tap_side):
    opf = _controlled(_undervoltage_net(tap_side), free_slack_vm=False)
    opf.add_voltage_deviation_objective()
    results = opf.solve(solver="gurobi_direct_minlp", time_limit=120)
    assert pyo.check_optimal_termination(results)
    assert opf.tap_schedule()[0] == (-1 if tap_side == "hv" else 1)
    assert opf.net.res_bus.vm_pu.at[1] >= 0.98 - 1e-6


@pytest.mark.integration
def test_relax_integrality_lets_ipopt_solve_the_relaxation():
    opf = _controlled(_undervoltage_net("hv"), free_slack_vm=False)
    opf.add_voltage_deviation_objective()
    results = opf.solve(solver="ipopt", relax_integrality=True)
    assert pyo.check_optimal_termination(results)
    k = pyo.value(opf.model.trafo_tap_position[0])
    assert abs(k - round(k)) > 1e-3  # fractional: it was a relaxation


PV_MW = 10.0


def _pv_net():
    """PV at the end of a 12 km MV line: at neutral tap the cap curtails it.

    Active power through the transformer alone hardly lifts the LV voltage
    (its impedance is almost purely reactive), so the rise that matters
    happens along a resistive overhead line: 10 MW give 1.057 p.u. at the
    PV bus with the tap at neutral, 1.020 at tap +2.
    """
    net = feeder(
        "hv",
        0,
        load_mw=2.0,
        load_mvar=0.5,
        vn_lv_kv=20.0,
        shift_degree=0.0,
        vm_pu=1.0,
        v_band=(0.95, 1.04),
    )
    pv_bus = pp.create_bus(net, 20.0, min_vm_pu=0.95, max_vm_pu=1.04)
    pp.create_line_from_parameters(
        net,
        1,
        pv_bus,
        length_km=12.0,
        r_ohm_per_km=0.3,
        x_ohm_per_km=0.35,
        c_nf_per_km=10.0,
        max_i_ka=0.5,
    )
    pp.create_sgen(
        net,
        pv_bus,
        p_mw=PV_MW,
        q_mvar=0.0,
        controllable=True,
        max_p_mw=PV_MW,
        min_p_mw=0.0,
        max_q_mvar=0.0,
        min_q_mvar=0.0,
    )
    return net


def _maximise_pv(opf):
    """Objective: as much PV infeed as possible."""

    @opf.model.Objective(sense=pyo.maximize)
    def obj_pv(model):
        """Total static-generator infeed.

        Args:
            model: The Pyomo model.

        Returns:
            A Pyomo expression to maximise.
        """
        return sum(model.psG[g] for g in model.sG)


@pytest.mark.integration
def test_overvoltage_tap_permits_more_pv():
    fixed = ACOPF(_pv_net())
    fixed.add_OPF(free_slack_vm=False)
    _maximise_pv(fixed)
    assert pyo.check_optimal_termination(fixed.solve(solver="ipopt"))
    p_fixed = fixed.net.res_sgen.p_mw.at[0]
    assert p_fixed < PV_MW - 1e-3  # the cap binds without the tap
    assert fixed.net.res_bus.vm_pu.at[2] == pytest.approx(1.04, abs=1e-5)

    opf = _controlled(_pv_net(), mode="continuous", free_slack_vm=False)
    _maximise_pv(opf)
    assert pyo.check_optimal_termination(opf.solve(solver="ipopt"))
    assert opf.tap_schedule()[0] > 0  # more HV turns lower the LV voltage
    assert opf.net.res_sgen.p_mw.at[0] > p_fixed + 0.5
    assert opf.net.res_bus.vm_pu.max() <= 1.04 + 1e-6
    assert opf.net.res_bus.vm_pu.min() >= 0.95 - 1e-6


@pytest.mark.integration
def test_priced_tap_movement_keeps_the_initial_position():
    """A tiny voltage gain does not pay for a tap operation once priced."""
    net = feeder(
        "hv",
        0,
        load_mw=8.0,
        vn_lv_kv=20.0,
        shift_degree=0.0,
        vm_pu=1.0,
        v_band=(0.9, 1.1),
    )
    free = _controlled(net, mode="continuous", free_slack_vm=False)
    free.add_voltage_deviation_objective()
    assert pyo.check_optimal_termination(free.solve(solver="ipopt"))
    assert abs(free.tap_schedule()[0]) > 0.05

    priced = _controlled(net, mode="continuous", free_slack_vm=False)
    priced.add_voltage_deviation_objective()
    priced.penalize_tap_movement(cost=1.0)
    assert pyo.check_optimal_termination(priced.solve(solver="ipopt"))
    assert priced.tap_schedule()[0] == pytest.approx(0.0, abs=1e-5)
    assert priced.tap_operations()[0] == pytest.approx(0.0, abs=1e-5)


def test_penalize_tap_movement_needs_an_objective():
    opf = _controlled(feeder())
    with pytest.raises(ValueError, match="objective"):
        opf.penalize_tap_movement(1.0)


@pytest.mark.integration
def test_tap_limits_hold_when_the_optimum_wants_more():
    """Slack at 0.75 p.u.: reaching 1.0 would need tap −12.5, so −9 binds."""
    net = feeder(
        "hv", 0, load_mw=5.0, vn_lv_kv=20.0, shift_degree=0.0, vm_pu=0.75
    )
    opf = _controlled(net, mode="continuous", free_slack_vm=False)
    opf.add_voltage_deviation_objective()
    assert pyo.check_optimal_termination(opf.solve(solver="ipopt"))
    k = opf.tap_schedule()[0]
    assert k == pytest.approx(-9.0, abs=1e-5)  # interior-point slack only
    assert k >= -9.0 - 1e-6  # never beyond tap_min
    assert pyo.value(opf.model.Tap[0]) >= opf.model.Tap[0].lb - 1e-6


def test_result_columns_and_legacy_methods():
    opf = _controlled(feeder("lv", 2), mode="continuous")
    assert hasattr(opf.model, "TRANSF_OLTC_LV")
    assert not opf.model.Tap_lv[0].fixed and opf.model.Tap[0].fixed
    legacy = ACOPF(feeder())
    legacy.add_OPF()
    with pytest.warns(DeprecationWarning):
        legacy.add_tap_changer_linear()
    assert hasattr(legacy.model, "Tap_linear_constr")


def test_apply_tap_positions_warns_when_pandapower_would_ignore_them():
    opf = _controlled(feeder())
    target = feeder(tap_changer_type=None)
    with pytest.warns(UserWarning, match="tap_changer_type"):
        opf.apply_tap_positions(target)
    assert target.trafo.at[0, "tap_pos"] == 0


# ── multi-period ─────────────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def lv_oltc_net():
    """SimBench LV rural1 with its tap changer typed and a tight band."""
    net = sb.get_simbench_net("1-LV-rural1--0-sw")
    net.trafo["tap_changer_type"] = "Ratio"
    net.bus["max_vm_pu"] = 1.03
    net.bus["min_vm_pu"] = 0.97
    return net


def test_multi_period_components_are_time_indexed(lv_oltc_net):
    mp = ACOPF_multi_period(lv_oltc_net, toT=3)
    mp.add_OPF()
    mp.enable_oltc(mode="discrete", max_change_per_step=1, max_operations=2)
    m = mp.model
    assert list(m.TRANSF_OLTC) == [0]
    assert len(m.trafo_tap_position) == 3
    assert len(m.trafo_tap_movement_def) == 3
    assert len(m.trafo_tap_change_limit) == 3
    assert len(m.trafo_tap_operations_limit) == 1
    assert all(not m.Tap[0, t].fixed for t in m.T)
    assert all(m.Tap_lv[0, t].fixed for t in m.T)
    assert len(free_integer_variables(m)) == 3


def test_multi_period_dc_refuses_oltc(lv_oltc_net):
    mp = DCOPF_multi_period(lv_oltc_net, toT=2)
    mp.add_OPF()
    with pytest.raises(NotImplementedError):
        mp.enable_oltc()


@pytest.mark.integration
def test_multi_period_rounded_schedule_respects_the_limits(lv_oltc_net):
    mp = ACOPF_multi_period(lv_oltc_net, toT=3)
    mp.add_OPF()
    mp.enable_oltc(mode="discrete", max_change_per_step=1, max_operations=2)
    mp.add_voltage_deviation_objective()
    results = mp.solve_oltc_round_and_fix(
        solver="ipopt", print_solver_output=False
    )
    assert pyo.check_optimal_termination(results)
    schedule = mp.tap_schedule()
    assert schedule.shape == (3, 1)
    positions = [int(mp.net.trafo.tap_pos.iloc[0])] + schedule[0].tolist()
    assert all(-2 <= k <= 2 for k in positions)
    assert all(abs(b - a) <= 1 for a, b in zip(positions[:-1], positions[1:]))
    assert mp.tap_operations()[0] <= 2
    assert mp.net.res_trafo.at[0, "tap_pos"] == schedule[0].iloc[-1]


@pytest.mark.integration
def test_multi_period_switching_cost_reduces_operations(lv_oltc_net):
    free = ACOPF_multi_period(lv_oltc_net, toT=3)
    free.add_OPF()
    free.enable_oltc(mode="continuous")
    free.add_voltage_deviation_objective()
    assert pyo.check_optimal_termination(
        free.solve(solver="ipopt", print_solver_output=False)
    )
    priced = ACOPF_multi_period(lv_oltc_net, toT=3)
    priced.add_OPF()
    priced.enable_oltc(mode="continuous")
    priced.add_voltage_deviation_objective()
    priced.penalize_tap_movement(cost=10.0)
    assert pyo.check_optimal_termination(
        priced.solve(solver="ipopt", print_solver_output=False)
    )
    assert priced.tap_operations()[0] <= free.tap_operations()[0] + 1e-6
    assert priced.tap_operations()[0] == pytest.approx(0.0, abs=1e-4)


@pytest.mark.integration
def test_multi_period_relaxation_bounds_the_discrete_and_fixed_models(
    lv_oltc_net,
):
    objective = {}
    fixed = ACOPF_multi_period(lv_oltc_net, toT=3)
    fixed.add_OPF()
    fixed.add_voltage_deviation_objective()
    assert pyo.check_optimal_termination(
        fixed.solve(solver="ipopt", print_solver_output=False)
    )
    objective["fixed"] = pyo.value(fixed.model.obj_v_deviation)
    for mode in ("continuous", "discrete"):
        mp = ACOPF_multi_period(lv_oltc_net, toT=3)
        mp.add_OPF()
        mp.enable_oltc(mode=mode)
        mp.add_voltage_deviation_objective()
        if mode == "continuous":
            results = mp.solve(solver="ipopt", print_solver_output=False)
        else:
            results = mp.solve_oltc_round_and_fix(
                solver="ipopt", print_solver_output=False
            )
        assert pyo.check_optimal_termination(results)
        objective[mode] = pyo.value(mp.model.obj_v_deviation)
    assert objective["continuous"] <= objective["discrete"] + 1e-7
    assert objective["discrete"] <= objective["fixed"] + 1e-7
