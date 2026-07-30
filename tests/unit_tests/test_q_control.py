"""Unit tests for grid-code reactive-power control (no solver required).

Three groups:

* the :mod:`~potpourri.technologies.q_control` grid-code registry and the
  capability curves it computes;
* single-period activation — which ``add_OPF`` arguments actually create
  which Pyomo blocks.  A missing precondition is silent (the constraint is
  simply not added and no warning is raised), so these tests pin the
  preconditions documented in ``docs/user-guide/reactive-power-control.md``;
* multi-period activation, which is driven by ``net.sgen`` columns rather
  than keyword arguments.

Every test deep-copies the network fixture before annotating it: the
fixtures are session-scoped and adding columns such as ``var_q`` would
otherwise leak into unrelated tests.
"""

import copy
import warnings

import numpy as np
import pytest

from potpourri.models.ACOPF_base import ACOPF
from potpourri.models_multi_period.ACOPF_multi_period import (
    ACOPF_multi_period,
)
from potpourri.technologies import q_control as qc

MP_STEPS = 3  # keep multi-period construction cheap


# ── helpers ───────────────────────────────────────────────────────────────


def _annotate(
    base_net,
    *,
    var_q=0,
    sn_mva=False,
    cos_phi_min=None,
    pu_curtail=False,
    fixed_cos_phi=None,
    cos_phi_p_profile=False,
    drop_sn_mva=False,
):
    """Return a deep copy of ``base_net`` with the requested annotations.

    ``drop_sn_mva`` removes the column entirely, which is needed because
    SimBench networks ship an ``sn_mva`` column by default — simply not
    setting it does not test its absence.
    """
    net = copy.deepcopy(base_net)
    mask = (
        net.sgen["type"].astype(str).str.contains("PV", case=False, na=False)
    )
    if not mask.any():
        mask = net.sgen.index == net.sgen.index[0]

    net.sgen["controllable"] = True
    net.sgen["p_inst_mw"] = net.sgen["p_mw"].abs()
    net.sgen["max_p_mw"] = net.sgen["p_mw"].abs()
    net.sgen["min_p_mw"] = 0.0
    net.bus["max_vm_pu"] = 1.05
    net.bus["min_vm_pu"] = 0.95

    net.sgen["var_q"] = None
    if var_q is not None:
        net.sgen.loc[mask, "var_q"] = int(var_q)

    if drop_sn_mva:
        if "sn_mva" in net.sgen:
            net.sgen.drop(columns=["sn_mva"], inplace=True)
    elif sn_mva:
        net.sgen["sn_mva"] = net.sgen["p_mw"].abs() / 0.9

    if cos_phi_min is not None:
        net.sgen["cos_phi_min"] = float(cos_phi_min)
    elif "cos_phi_min" in net.sgen:
        net.sgen.drop(columns=["cos_phi_min"], inplace=True)

    if pu_curtail:
        net.sgen["pu_curtail"] = False
        net.sgen.loc[mask, "pu_curtail"] = True
    if fixed_cos_phi is not None:
        net.sgen["fixed_cos_phi"] = float("nan")
        net.sgen.loc[mask, "fixed_cos_phi"] = float(fixed_cos_phi)
    if cos_phi_p_profile:
        net.sgen["cos_phi_p_profile"] = False
        net.sgen.loc[mask, "cos_phi_p_profile"] = True

    return net


def _build_sp(net, **add_opf_kwargs):
    """Build a single-period ACOPF, suppressing the provisional warning."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", qc.ProvisionalGridCodeWarning)
        opf = ACOPF(net)
        opf.add_OPF(**add_opf_kwargs)
    return opf


def _build_mp(net, **add_opf_kwargs):
    """Build a multi-period ACOPF, suppressing the provisional warning."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", qc.ProvisionalGridCodeWarning)
        mp = ACOPF_multi_period(net, toT=MP_STEPS)
        mp.add_OPF(**add_opf_kwargs)
    return mp


# ── registry: capability curves ───────────────────────────────────────────


def test_compute_q_curves_columns():
    """All seven slope/intercept columns must be present."""
    curves = qc.compute_q_curves()
    assert set(curves.columns) == {
        "m_qv",
        "b_qv_min",
        "b_qv_max",
        "m_qp_max",
        "m_qp_min",
        "b_qp_max",
        "b_qp_min",
    }


def test_compute_q_curves_one_row_per_variant():
    """The frame is indexed by var_q variant."""
    curves = qc.compute_q_curves()
    assert len(curves) == qc.VDE_AR_N_4105.n_variants == 3


@pytest.mark.parametrize(
    ("variant", "q_max", "q_min"),
    [(0, 0.48, -0.23), (1, 0.41, -0.33), (2, 0.33, -0.41)],
)
def test_q_envelope_at_reference_point(variant, q_max, q_min):
    """At P = 0.2 Pn the envelope must hit the VDE-AR-N 4105 table values."""
    code = qc.VDE_AR_N_4105
    curves = qc.compute_q_curves(code)
    p = code.qp_p_low
    upper = curves.b_qp_max[variant] + curves.m_qp_max[variant] * p
    lower = curves.b_qp_min[variant] + curves.m_qp_min[variant] * p
    assert upper == pytest.approx(q_max, abs=1e-9)
    assert lower == pytest.approx(q_min, abs=1e-9)


@pytest.mark.parametrize("variant", [0, 1, 2])
def test_q_envelope_narrows_at_lower_breakpoint(variant):
    """At P = 0.1 Pn the envelope collapses to ±0.1 Pn for every variant."""
    code = qc.VDE_AR_N_4105
    curves = qc.compute_q_curves(code)
    p = code.qp_p_high
    upper = curves.b_qp_max[variant] + curves.m_qp_max[variant] * p
    lower = curves.b_qp_min[variant] + curves.m_qp_min[variant] * p
    assert upper == pytest.approx(code.qp_p_high, abs=1e-9)
    assert lower == pytest.approx(-code.qp_p_high, abs=1e-9)


def test_qp_bound_is_not_clipped_above_reference_point():
    """The Q(P) bound is one unclipped linear segment.

    It keeps widening beyond P = 0.2 Pn rather than holding at the table
    value, which is why Q(P) alone does not limit reactive power at high
    active output — the S² circle and cos(phi) cone do.  Documented in the
    user guide; pinned here so the behaviour cannot change silently.
    """
    curves = qc.compute_q_curves(qc.VDE_AR_N_4105)
    at_ref = curves.b_qp_max[0] + curves.m_qp_max[0] * 0.2
    at_full = curves.b_qp_max[0] + curves.m_qp_max[0] * 1.0
    assert at_full > at_ref
    assert at_full > 1.0  # far beyond any physical inverter rating


def test_q_curves_monotonic_in_variant():
    """Variant 0 is the widest capacitive envelope, variant 2 the narrowest."""
    curves = qc.compute_q_curves(qc.VDE_AR_N_4105)
    p = qc.VDE_AR_N_4105.qp_p_low
    upper = [curves.b_qp_max[v] + curves.m_qp_max[v] * p for v in curves.index]
    assert upper[0] > upper[1] > upper[2]


# ── registry: grid-code resolution ────────────────────────────────────────


def test_resolve_default_is_4105():
    """No selector means VDE-AR-N 4105."""
    assert qc.resolve_grid_code(None) is qc.VDE_AR_N_4105
    assert qc.DEFAULT_GRID_CODE is qc.VDE_AR_N_4105


def test_resolve_accepts_grid_code_instance():
    """A GridCode instance is returned unchanged."""
    assert qc.resolve_grid_code(qc.VDE_AR_N_4105) is qc.VDE_AR_N_4105


@pytest.mark.parametrize(
    "selector", ["4110", "VDE-AR-N 4110", "vde-ar-n-4110", "VDEARN4110"]
)
def test_resolve_accepts_aliases(selector):
    """Short name, full designation and punctuation variants all resolve."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", qc.ProvisionalGridCodeWarning)
        assert qc.resolve_grid_code(selector) is qc.VDE_AR_N_4110


def test_resolve_unknown_raises_value_error():
    """An unknown grid code is an error, not a silent fallback."""
    with pytest.raises(ValueError, match="Unknown grid code"):
        qc.resolve_grid_code("9999")


def test_registry_contains_both_codes():
    """Both shipped grid codes are registered under their short names."""
    assert set(qc.GRID_CODES) == {"4105", "4110"}


# ── registry: provisional-value contract ──────────────────────────────────


def test_provisional_code_warns():
    """Selecting a placeholder grid code must warn, not fail silently."""
    with pytest.warns(qc.ProvisionalGridCodeWarning):
        qc.resolve_grid_code("4110")


def test_provisional_warning_names_the_code():
    """The warning has to say which code is provisional and why."""
    with pytest.warns(qc.ProvisionalGridCodeWarning) as record:
        qc.resolve_grid_code("4110")
    message = str(record[0].message)
    assert "4110" in message
    assert "placeholder" in message.lower()
    assert "not" in message.lower()


def test_normative_code_does_not_warn():
    """VDE-AR-N 4105 carries real values and must stay quiet."""
    with warnings.catch_warnings():
        warnings.simplefilter("error", qc.ProvisionalGridCodeWarning)
        qc.resolve_grid_code("4105")  # must not raise


def test_4110_is_still_flagged_provisional():
    """Tripwire: VDE-AR-N 4110 currently holds placeholder values.

    When its normative medium-voltage parameters are entered and the
    ``provisional`` flag is cleared, this test fails on purpose — that is
    the reminder to update the user guide and the CHANGELOG, which both
    state that 4110 results are not compliant.
    """
    assert qc.VDE_AR_N_4110.provisional is True
    assert qc.VDE_AR_N_4105.provisional is False


def test_backcompat_constants_mirror_default_code():
    """The legacy module constants must stay aliases of the default code."""
    code = qc.DEFAULT_GRID_CODE
    assert np.array_equal(qc.VQU_V_POINTS, code.vqu_v_points)
    assert np.array_equal(qc.VQU_Q_MAX, code.vqu_q_max)
    assert qc.QP_P_HIGH == code.qp_p_high
    assert qc.QP_P_LOW == code.qp_p_low
    assert qc.VPU_V_CURTAIL == code.vpu_v_curtail
    assert qc.VPU_V_MAX == code.vpu_v_max
    assert qc.CPP_P_THRESHOLD_PU == code.cpp_p_threshold_pu


# ── single-period: Q(P) / Q(U) activation ─────────────────────────────────


@pytest.mark.parametrize(
    ("mode", "expect_qp", "expect_qu"),
    [
        (None, False, False),
        (False, False, False),
        ("qp", True, False),
        ("qu", False, True),
        ("both", True, True),
        (True, True, True),
    ],
)
def test_pv_q_control_modes(lv_rural_net, mode, expect_qp, expect_qu):
    """pv_q_control selects exactly the requested characteristic."""
    opf = _build_sp(_annotate(lv_rural_net), pv_q_control=mode)
    assert hasattr(opf.model, "PV_QP_pos") is expect_qp
    assert hasattr(opf.model, "PV_QP_neg") is expect_qp
    assert hasattr(opf.model, "PV_QU_max") is expect_qu
    assert hasattr(opf.model, "PV_QU_min") is expect_qu


def test_pv_q_control_needs_var_q(lv_rural_net):
    """Without var_q there is nothing to constrain."""
    net = _annotate(lv_rural_net, var_q=None)
    opf = _build_sp(net, pv_q_control="both")
    assert not hasattr(opf.model, "PV_QP_pos")


# ── single-period: inverter S² circle and cos(phi) cone ───────────────────


def test_inverter_s2_is_opt_in(lv_rural_net):
    """inverter_s2 defaults to False, so sn_mva alone changes nothing."""
    opf = _build_sp(_annotate(lv_rural_net, sn_mva=True))
    assert not hasattr(opf.model, "sgen_inverter_s2")


def test_inverter_s2_enabled(lv_rural_net):
    """inverter_s2=True plus sn_mva creates the apparent-power circle."""
    opf = _build_sp(_annotate(lv_rural_net, sn_mva=True), inverter_s2=True)
    assert hasattr(opf.model, "sgen_inverter_s2")


def test_inverter_s2_requires_sn_mva(lv_rural_net):
    """Without sn_mva there is no rating to bound against."""
    opf = _build_sp(
        _annotate(lv_rural_net, drop_sn_mva=True), inverter_s2=True
    )
    assert not hasattr(opf.model, "sgen_inverter_s2")


def test_cos_phi_cone_requires_inverter_s2(lv_rural_net):
    """The cone is nested inside the S² block.

    Setting net.sgen["cos_phi_min"] on its own adds no constraint and
    raises no warning — the failure mode this pins down.
    """
    net = _annotate(lv_rural_net, sn_mva=True, cos_phi_min=0.9)
    opf = _build_sp(net)  # inverter_s2 left at its default
    assert not hasattr(opf.model, "sgen_cos_phi_upper")
    assert not hasattr(opf.model, "sgen_cos_phi_lower")


def test_cos_phi_cone_from_column(lv_rural_net):
    """With inverter_s2=True the per-row column creates the cone."""
    net = _annotate(lv_rural_net, sn_mva=True, cos_phi_min=0.9)
    opf = _build_sp(net, inverter_s2=True)
    assert hasattr(opf.model, "sgen_cos_phi_upper")
    assert hasattr(opf.model, "sgen_cos_phi_lower")


def test_cos_phi_cone_from_scalar_kwarg(lv_rural_net):
    """A scalar cos_phi_min works when no column is present."""
    net = _annotate(lv_rural_net, sn_mva=True)
    opf = _build_sp(net, inverter_s2=True, cos_phi_min=0.9)
    assert hasattr(opf.model, "sgen_cos_phi_upper")


def test_cos_phi_cone_tan_phi_value(lv_rural_net):
    """tan_phi must equal tan(arccos(cos_phi_min))."""
    net = _annotate(lv_rural_net, sn_mva=True, cos_phi_min=0.9)
    opf = _build_sp(net, inverter_s2=True)
    expected = float(np.tan(np.arccos(0.9)))
    for g in opf.model.sGpf:
        assert opf.model.tan_phi[g] == pytest.approx(expected)


def test_cos_phi_cone_needs_sn_mva(lv_rural_net):
    """No sn_mva means no S² block, hence no cone either."""
    net = _annotate(lv_rural_net, drop_sn_mva=True, cos_phi_min=0.9)
    opf = _build_sp(net, inverter_s2=True)
    assert not hasattr(opf.model, "sgen_cos_phi_upper")


# ── single-period: inverter controller modes ──────────────────────────────


def test_pu_curtail_activation(lv_rural_net):
    """pu_curtail=True creates the P(U) curtailment constraint."""
    opf = _build_sp(_annotate(lv_rural_net), pu_curtail=True)
    assert hasattr(opf.model, "sgen_pu_curtail")


def test_pu_curtail_is_opt_in(lv_rural_net):
    """Without the flag no curtailment constraint appears."""
    opf = _build_sp(_annotate(lv_rural_net))
    assert not hasattr(opf.model, "sgen_pu_curtail")


def test_pu_curtail_thresholds_from_grid_code(lv_rural_net):
    """Absent per-row columns, thresholds come from the grid code."""
    opf = _build_sp(_annotate(lv_rural_net), pu_curtail=True)
    code = qc.VDE_AR_N_4105
    for g in opf.model.sGpu:
        assert opf.model.V_curtail[g] == pytest.approx(code.vpu_v_curtail)
        assert opf.model.V_max_curtail[g] == pytest.approx(code.vpu_v_max)


def test_fixed_cos_phi_scalar(lv_rural_net):
    """A scalar fixed_cos_phi creates the equality for every sgen."""
    opf = _build_sp(_annotate(lv_rural_net), fixed_cos_phi=0.95)
    assert hasattr(opf.model, "sgen_fixed_cos_phi")


def test_fixed_cos_phi_value(lv_rural_net):
    """fixed_tan_phi must equal tan(arccos(cos_phi))."""
    opf = _build_sp(_annotate(lv_rural_net), fixed_cos_phi=0.95)
    expected = float(np.tan(np.arccos(0.95)))
    for g in opf.model.sGfcf:
        assert opf.model.fixed_tan_phi[g] == pytest.approx(expected)


def test_cos_phi_p_profile_activation(lv_rural_net):
    """cos_phi_p_profile=True plus cos_phi_min creates the CPP equality."""
    net = _annotate(lv_rural_net, cos_phi_min=0.9)
    opf = _build_sp(net, cos_phi_p_profile=True)
    assert hasattr(opf.model, "sgen_cpp")


def test_cos_phi_p_profile_needs_cos_phi_min(lv_rural_net):
    """Without a power factor at full output there is no profile."""
    net = _annotate(lv_rural_net)  # cos_phi_min column dropped
    opf = _build_sp(net, cos_phi_p_profile=True)
    assert not hasattr(opf.model, "sgen_cpp")


def test_equality_modes_are_mutually_exclusive_in_practice(lv_rural_net):
    """Stacking both equality modes builds two constraints on the same Q.

    Fixed cos(phi) and the cos(phi)(P) profile are both equalities on
    qsG, so together they over-determine reactive power.  This test
    documents that the model does not guard against it — the user guide
    tells callers to keep the two on disjoint sgens.
    """
    net = _annotate(lv_rural_net, cos_phi_min=0.9)
    opf = _build_sp(net, fixed_cos_phi=0.95, cos_phi_p_profile=True)
    assert hasattr(opf.model, "sgen_fixed_cos_phi")
    assert hasattr(opf.model, "sgen_cpp")


# ── single-period: grid-code threading ────────────────────────────────────


def test_grid_code_resolved_onto_model(lv_rural_net):
    """add_OPF stores the resolved code so later steps can read it."""
    opf = _build_sp(_annotate(lv_rural_net), pv_q_control="both")
    assert opf._grid_code is qc.VDE_AR_N_4105


def test_grid_code_selection_reaches_single_period(lv_rural_net):
    """Selecting 4110 must actually resolve to 4110 in the model."""
    opf = _build_sp(
        _annotate(lv_rural_net), pv_q_control="both", grid_code="4110"
    )
    assert opf._grid_code is qc.VDE_AR_N_4110


def test_grid_code_curves_match_selected_code(lv_rural_net):
    """q_limit_parameter must come from the selected grid code."""
    opf = _build_sp(_annotate(lv_rural_net), pv_q_control="both")
    expected = qc.compute_q_curves(qc.VDE_AR_N_4105)
    assert np.allclose(opf.q_limit_parameter.values, expected.values)


def test_unknown_grid_code_raises_from_add_opf(lv_rural_net):
    """An unknown grid code fails at add_OPF rather than being ignored."""
    with pytest.raises(ValueError, match="Unknown grid code"):
        ACOPF(_annotate(lv_rural_net)).add_OPF(grid_code="nope")


def test_wind_var_q_populates_reactive_bounds(lv_rural_net):
    """var_q sgens get finite Q bounds from the capability table.

    Regression test: this path once referenced a table that had been
    removed during a refactor, and the whole suite still passed because
    nothing exercised it.
    """
    opf = _build_sp(_annotate(lv_rural_net), pv_q_control="both")
    data = opf.static_generation_data
    idx = data.index[data.var_q.notna()]
    assert len(idx) > 0
    assert np.isfinite(data["max_q"][idx].astype(float)).all()
    assert np.isfinite(data["min_q"][idx].astype(float)).all()
    assert (
        data["max_q"][idx].astype(float) > data["min_q"][idx].astype(float)
    ).all()


def test_wind_var_q_bounds_scale_with_capability_table(lv_rural_net):
    """max_q / min_q must be Pn times the grid-code Q/Pn entries."""
    opf = _build_sp(_annotate(lv_rural_net, var_q=0), pv_q_control="both")
    data = opf.static_generation_data
    idx = data.index[data.var_q.notna()]
    table = qc.VDE_AR_N_4105.vqu_q_max
    for g in idx:
        pn = float(data["p_inst"][g])
        assert float(data["max_q"][g]) == pytest.approx(table[0, 0] * pn)
        assert float(data["min_q"][g]) == pytest.approx(table[1, 0] * pn)


# ── multi-period: activation is column-driven ─────────────────────────────


def test_mp_q_control_auto_detected(lv_rural_net):
    """var_q alone is enough; the multi-period model needs no flag."""
    mp = _build_mp(_annotate(lv_rural_net))
    assert hasattr(mp.model, "sG_QP_pos")
    assert hasattr(mp.model, "sG_QP_neg")
    assert hasattr(mp.model, "sG_QU_max")
    assert hasattr(mp.model, "sG_QU_min")


def test_mp_no_var_q_no_q_control(lv_rural_net):
    """Without var_q the Q-control blocks stay absent."""
    mp = _build_mp(_annotate(lv_rural_net, var_q=None))
    assert not hasattr(mp.model, "sG_QP_pos")


def test_mp_inverter_s2_auto_detected(lv_rural_net):
    """Unlike the single-period model, sn_mva alone enables the circle."""
    mp = _build_mp(_annotate(lv_rural_net, sn_mva=True))
    assert hasattr(mp.model, "sgen_inverter_s2")


def test_mp_inverter_s2_absent_without_sn_mva(lv_rural_net):
    """No rating column, no circle."""
    mp = _build_mp(_annotate(lv_rural_net, drop_sn_mva=True))
    assert not hasattr(mp.model, "sgen_inverter_s2")


def test_mp_cone_auto_detected(lv_rural_net):
    """sn_mva plus cos_phi_min is sufficient — no inverter_s2 flag exists."""
    mp = _build_mp(_annotate(lv_rural_net, sn_mva=True, cos_phi_min=0.9))
    assert hasattr(mp.model, "sgen_cos_phi_upper")
    assert hasattr(mp.model, "sgen_cos_phi_lower")


def test_mp_cone_needs_cos_phi_min(lv_rural_net):
    """The circle can exist without the cone."""
    mp = _build_mp(_annotate(lv_rural_net, sn_mva=True))
    assert hasattr(mp.model, "sgen_inverter_s2")
    assert not hasattr(mp.model, "sgen_cos_phi_upper")


def test_mp_pu_curtail_auto_detected(lv_rural_net):
    """The pu_curtail column drives the P(U) constraint."""
    mp = _build_mp(_annotate(lv_rural_net, pu_curtail=True))
    assert hasattr(mp.model, "sgen_pu_curtail")


def test_mp_fixed_cos_phi_auto_detected(lv_rural_net):
    """The fixed_cos_phi column drives the equality."""
    mp = _build_mp(_annotate(lv_rural_net, fixed_cos_phi=0.95))
    assert hasattr(mp.model, "sgen_fixed_cos_phi")


def test_mp_cpp_auto_detected(lv_rural_net):
    """cos_phi_p_profile plus cos_phi_min drives the CPP equality."""
    mp = _build_mp(
        _annotate(lv_rural_net, cos_phi_p_profile=True, cos_phi_min=0.9)
    )
    assert hasattr(mp.model, "sgen_cpp")


def test_mp_time_indexed_constraints(lv_rural_net):
    """Multi-period Q-control constraints are indexed by (sgen, time)."""
    mp = _build_mp(_annotate(lv_rural_net))
    keys = list(mp.model.sG_QP_pos)
    assert keys, "expected at least one Q-controlled sgen"
    assert len(keys[0]) == 2
    times = {k[1] for k in keys}
    assert times == set(mp.model.T)


def test_mp_grid_code_selection_threaded(lv_rural_net):
    """grid_code reaches the multi-period sgen device module."""
    mp = _build_mp(_annotate(lv_rural_net), grid_code="4110")
    sgens = next(obj for obj in mp.flexibilities if hasattr(obj, "grid_code"))
    assert sgens.grid_code is qc.VDE_AR_N_4110


def test_mp_per_sgen_strategies_are_independent(lv_rural_net):
    """Different sgens may follow different controllers in one model.

    Mirrors scripts/grid_code_q_strategies.py: the strategies are assigned
    to disjoint sgens, so no unit gets two equalities on its reactive power.
    """
    net = copy.deepcopy(lv_rural_net)
    mask = (
        net.sgen["type"].astype(str).str.contains("PV", case=False, na=False)
    )
    pv = list(net.sgen.index[mask]) or list(net.sgen.index[:1])
    if len(pv) < 3:
        pytest.skip("needs at least three sgens to separate strategies")

    net.sgen["controllable"] = True
    net.sgen["p_inst_mw"] = net.sgen["p_mw"].abs()
    net.sgen["max_p_mw"] = net.sgen["p_mw"].abs()
    net.sgen["min_p_mw"] = 0.0
    net.bus["max_vm_pu"] = 1.05
    net.bus["min_vm_pu"] = 0.95

    net.sgen["var_q"] = None
    net.sgen["fixed_cos_phi"] = float("nan")
    net.sgen["cos_phi_p_profile"] = False
    net.sgen["cos_phi_min"] = float("nan")

    net.sgen.at[pv[0], "var_q"] = 0
    net.sgen.at[pv[1], "fixed_cos_phi"] = 0.95
    net.sgen.at[pv[2], "cos_phi_p_profile"] = True
    net.sgen.at[pv[2], "cos_phi_min"] = 0.9

    mp = _build_mp(net)
    assert hasattr(mp.model, "sG_QP_pos")
    assert hasattr(mp.model, "sgen_fixed_cos_phi")
    assert hasattr(mp.model, "sgen_cpp")
    # each equality applies to exactly the sgen it was assigned to
    assert {k[0] for k in mp.model.sgen_fixed_cos_phi} == {pv[1]}
    assert {k[0] for k in mp.model.sgen_cpp} == {pv[2]}
