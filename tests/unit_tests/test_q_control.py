# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

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
import pathlib
import warnings
from dataclasses import replace

import numpy as np
import pyomo.environ as pyo
import pytest

from potpourri.models.ACOPF_base import DEFAULT_PV_SGEN_TYPES, ACOPF
from potpourri.models_multi_period.ACOPF_multi_period import (
    ACOPF_multi_period,
)
from potpourri.technologies import q_control as qc
from potpourri.technologies import windpower

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
    assert len(curves) == qc.DEFAULT_GRID_CODE.n_variants == 3


@pytest.mark.parametrize(
    ("variant", "q_max", "q_min"),
    [
        (0, 0.484322, -0.227902),
        (1, 0.410775, -0.328684),
        (2, 0.328684, -0.410775),
    ],
)
def test_q_envelope_at_reference_point(variant, q_max, q_min):
    """At P = 0.2 Pn the envelope hits the VDE-AR-N 4120 table values."""
    code = qc.VDE_AR_N_4120
    lower, upper = code.pq_area.q_flexibility(code.qp_p_low, variant)
    assert upper == pytest.approx(q_max, abs=1e-9)
    assert lower == pytest.approx(q_min, abs=1e-9)


@pytest.mark.parametrize("variant", [0, 1, 2])
def test_q_envelope_narrows_at_lower_breakpoint(variant):
    """At P = 0.1 Pn the envelope collapses to ±0.1 Pn for every variant."""
    code = qc.VDE_AR_N_4120
    lower, upper = code.pq_area.q_flexibility(code.qp_p_high, variant)
    assert upper == pytest.approx(code.qp_p_high, abs=1e-9)
    assert lower == pytest.approx(-code.qp_p_high, abs=1e-9)


def test_qp_bound_saturates_above_the_reference_point():
    """The Q(P) bound holds its limit above P = 0.2 Pn instead of climbing.

    Before 0.4.1 this bound was a single unclipped line that reached
    +3.52 Pn at full output against a grid-code limit of +0.484 — the
    saturation shelf was missing.  Pinned so it cannot regress.
    """
    code = qc.VDE_AR_N_4120
    _, at_ref = code.pq_area.q_flexibility(0.2, 0)
    _, at_full = code.pq_area.q_flexibility(1.0, 0)
    assert at_full == pytest.approx(at_ref, abs=1e-12)
    assert at_full == pytest.approx(0.484322, abs=1e-9)


def test_qu_bound_saturates_at_nominal_voltage():
    """At v = 1.0 the Q(U) band is the plateau, not an extrapolated line.

    Before 0.4.1 the two bounds were unclipped parallel lines, giving
    [-0.940, +1.494] at nominal voltage — 3.2x the grid-code range.
    """
    lower, upper = qc.VDE_AR_N_4120.qv_area.q_flexibility(1.0, 0)
    assert lower == pytest.approx(-0.227902, abs=1e-9)
    assert upper == pytest.approx(0.484322, abs=1e-9)


def test_q_curves_monotonic_in_variant():
    """Variant 0 is the widest capacitive envelope, variant 2 the narrowest."""
    code = qc.VDE_AR_N_4120
    upper = [
        code.pq_area.q_flexibility(code.qp_p_low, v)[1]
        for v in range(code.n_variants)
    ]
    assert upper[0] > upper[1] > upper[2]


# ── registry: grid-code resolution ────────────────────────────────────────


def test_resolve_default_is_4120():
    """No selector means VDE-AR-N 4120.

    The pre-0.4.1 constants were labelled 4105 but carried the 110 kV
    breakpoints (96/103/120/127 kV) and the three 4120 variants, so 4120 is
    the code those defaults actually described.  Keeping it as the default
    also keeps existing nets working, whose var_q spans 0..2.
    """
    assert qc.resolve_grid_code(None) is qc.VDE_AR_N_4120
    assert qc.DEFAULT_GRID_CODE is qc.VDE_AR_N_4120


def test_resolve_accepts_grid_code_instance():
    """A GridCode instance is returned unchanged."""
    assert qc.resolve_grid_code(qc.VDE_AR_N_4105) is qc.VDE_AR_N_4105


@pytest.mark.parametrize(
    "selector", ["4110", "VDE-AR-N 4110", "vde-ar-n-4110", "VDEARN4110"]
)
def test_resolve_accepts_aliases(selector):
    """Short name, full designation and punctuation variants all resolve."""
    assert qc.resolve_grid_code(selector) is qc.VDE_AR_N_4110


def test_resolve_unknown_raises_value_error():
    """An unknown grid code is an error, not a silent fallback."""
    with pytest.raises(ValueError, match="Unknown grid code"):
        qc.resolve_grid_code("9999")


def test_registry_contains_every_shipped_code():
    """All three shipped grid codes are registered under their short names."""
    assert set(qc.GRID_CODES) == {"4105", "4110", "4120"}


@pytest.mark.parametrize(
    ("name", "n_variants", "level"),
    [
        ("4105", 2, "low voltage"),
        ("4110", 1, "medium voltage"),
        ("4120", 3, "high voltage"),
    ],
)
def test_variant_count_per_code(name, n_variants, level):
    """Each code exposes exactly the variants the standard defines."""
    code = qc.GRID_CODES[name]
    assert code.n_variants == n_variants
    assert code.voltage_level == level


# ── registry: provisional-value contract ──────────────────────────────────


def test_no_shipped_code_is_provisional():
    """Every shipped code now carries values sourced from pandapower.

    4110 was a placeholder copy of the (mislabelled) 4105 entry until
    0.4.1.  If a future code ships provisional again, flip this and restore
    a tripwire like the one this replaced.
    """
    assert not any(c.provisional for c in qc.GRID_CODES.values())


def test_no_shipped_code_warns_on_resolve():
    """Resolving any registered code must stay quiet."""
    with warnings.catch_warnings():
        warnings.simplefilter("error", qc.ProvisionalGridCodeWarning)
        for name in qc.GRID_CODES:
            qc.resolve_grid_code(name)  # must not raise


def test_provisional_mechanism_still_warns():
    """The provisional machinery works, even with nothing shipped using it."""
    stub = replace(
        qc.VDE_AR_N_4105,
        name="stub",
        provisional=True,
        provisional_note="stub carries placeholder values, not normative.",
    )
    with pytest.warns(qc.ProvisionalGridCodeWarning) as record:
        qc.resolve_grid_code(stub)
    assert "placeholder" in str(record[0].message).lower()


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


# ── single-period: sgen type selection ────────────────────────────────────


def _retype(base_net, new_type):
    """Deep copy with every sgen relabelled to ``new_type``.

    SimBench names medium-voltage PV ``PV_MV`` rather than ``PV``, so
    relabelling the LV fixture reproduces the MV naming without downloading
    an MV grid.
    """
    net = _annotate(base_net)
    net.sgen["type"] = new_type
    return net


def test_default_types_include_lv_and_mv_pv():
    """The default must cover both SimBench PV spellings."""
    assert "PV" in DEFAULT_PV_SGEN_TYPES
    assert "PV_MV" in DEFAULT_PV_SGEN_TYPES


@pytest.mark.parametrize("sgen_type", ["PV", "PV_MV"])
def test_default_types_reach_both_pv_spellings(lv_rural_net, sgen_type):
    """Q-control reaches PV at LV and at MV without extra arguments.

    Regression test: matching only ``type == "PV"`` silently reached zero
    sgens on every SimBench MV grid, building no constraints and raising no
    warning.
    """
    opf = _build_sp(_retype(lv_rural_net, sgen_type), pv_q_control="both")
    assert hasattr(opf.model, "PV_QP_pos")
    assert len(list(opf.model.PVc)) > 0


def test_unlisted_type_is_not_reached(lv_rural_net):
    """A category outside the list stays unconstrained."""
    opf = _build_sp(_retype(lv_rural_net, "Wind_MV"), pv_q_control="both")
    assert not hasattr(opf.model, "PV_QP_pos")


def test_sgen_types_widens_selection(lv_rural_net):
    """An explicit list can pull in further categories."""
    net = _retype(lv_rural_net, "lv_RES")
    opf = _build_sp(
        net, pv_q_control="both", sgen_types=("PV", "PV_MV", "lv_RES")
    )
    assert hasattr(opf.model, "PV_QP_pos")
    assert len(list(opf.model.PVc)) > 0


def test_sgen_types_narrows_selection(lv_rural_net):
    """An explicit list can also exclude the defaults."""
    opf = _build_sp(
        _retype(lv_rural_net, "PV"),
        pv_q_control="both",
        sgen_types=("PV_MV",),
    )
    assert not hasattr(opf.model, "PV_QP_pos")


def test_sgen_types_matching_is_exact(lv_rural_net):
    """Matching is exact, not substring: 'PV' must not catch 'PV_MV_extra'."""
    opf = _build_sp(
        _retype(lv_rural_net, "PV_MV_extra"),
        pv_q_control="both",
        sgen_types=("PV",),
    )
    assert not hasattr(opf.model, "PV_QP_pos")


def test_sgen_types_accepts_any_iterable(lv_rural_net):
    """A list works as well as a tuple."""
    opf = _build_sp(
        _retype(lv_rural_net, "PV_MV"),
        pv_q_control="both",
        sgen_types=["PV_MV"],
    )
    assert hasattr(opf.model, "PV_QP_pos")


# ── wind path: sgen type selection ────────────────────────────────────────


def test_default_wind_types_cover_every_simbench_spelling():
    """SimBench spells wind four ways across the voltage levels."""
    assert set(qc.DEFAULT_WIND_SGEN_TYPES) == {
        "Wind",
        "Wind_MV",
        "wind onshore",
        "wind offshore",
    }


@pytest.mark.parametrize(
    "sgen_type", ["Wind", "Wind_MV", "wind onshore", "wind offshore"]
)
def test_wind_path_reaches_every_spelling(lv_rural_net, sgen_type):
    """The wind Q path must reach wind at MV and EHV, not only HV.

    Regression test: matching only ``type == "Wind"`` reached nothing on any
    SimBench MV or EHV grid, so no wind Q-control was built and nothing
    warned.
    """
    opf = _build_sp(_retype(lv_rural_net, sgen_type))
    assert len(list(opf.model.WINDc)) > 0


def test_wind_path_ignores_unrelated_types(lv_rural_net):
    """A non-wind category stays out of the wind set."""
    opf = _build_sp(_retype(lv_rural_net, "Biomass_MV"))
    assert len(list(opf.model.WINDc)) == 0


def test_wind_sgen_types_is_overridable(lv_rural_net):
    """An explicit list replaces the default."""
    opf = _build_sp(
        _retype(lv_rural_net, "Wind_MV"), wind_sgen_types=("Wind",)
    )
    assert len(list(opf.model.WINDc)) == 0


def test_pv_and_wind_paths_are_disjoint_by_default(lv_rural_net):
    """With the defaults no sgen can be claimed by both paths."""
    assert not (set(DEFAULT_PV_SGEN_TYPES) & set(qc.DEFAULT_WIND_SGEN_TYPES))


def test_overlapping_selection_warns_and_defers_to_wind(lv_rural_net):
    """Listing a wind category in sgen_types must not double-constrain qsG.

    Both paths impose the same grid-code characteristic on the same
    variable, so an sgen in both sets would get duplicate constraints.  The
    wind path owns those types; PVc gives them up and the caller is told.
    """
    net = _retype(lv_rural_net, "Wind_MV")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        opf = ACOPF(net)
        opf.add_OPF(
            pv_q_control="both",
            sgen_types=("PV", "PV_MV", "pv", "Wind_MV"),
        )
    overlaps = [
        w for w in caught if issubclass(w.category, qc.SgenTypeOverlapWarning)
    ]
    assert overlaps, "expected an overlap warning"
    wind = set(opf.model.WINDc)
    pvc = set(opf.model.PVc) if hasattr(opf.model, "PVc") else set()
    assert wind, "wind path should still claim them"
    assert not (pvc & wind), "PVc must not overlap WINDc"


def test_no_overlap_warning_with_default_types(lv_rural_net):
    """The default configuration must stay silent."""
    with warnings.catch_warnings():
        warnings.simplefilter("error", qc.SgenTypeOverlapWarning)
        _build_sp(_retype(lv_rural_net, "Wind_MV"), pv_q_control="both")


# ── hosting-capacity path uses the registry ───────────────────────────────


def test_windpower_has_no_private_capability_table():
    """Windpower must read the registry, not its own copy of the table.

    It previously carried a private VDE-AR-N 4105 table and a third
    byte-identical copy of the Q-curve maths, so ``grid_code`` never reached
    the wind or hosting-capacity paths and the numbers could drift from the
    registry unnoticed.
    """
    src = pathlib.Path(windpower.__file__).read_text()
    for name in (
        "_VQU_V_POINTS",
        "_VQU_Q_MAX",
        "_QP_P_BREAK_HIGH",
        "_QP_HC_MAX",
        "_QP_HC_MIN",
    ):
        assert name not in src, f"{name} still defined in windpower"


def test_hc_q_bounds_are_the_widest_envelope():
    """The simplified HC check uses the widest band the code offers."""
    hc_max, hc_min = windpower._hc_q_bounds(qc.VDE_AR_N_4120)
    assert hc_max == pytest.approx(0.484322)
    assert hc_min == pytest.approx(-0.410775)


def test_hc_defaults_track_the_default_grid_code():
    """The public HC keyword defaults are derived, not literals.

    ``qp_max`` / ``qp_min`` on ``Windpower_multi_period.__init__`` were
    hard-coded 0.48 / -0.41; they now come from the default grid code, which
    supplies the same limits at the standard's own precision.
    """
    expected_max, expected_min = windpower._hc_q_bounds(qc.DEFAULT_GRID_CODE)
    assert windpower._DEFAULT_HC_Q_MAX == expected_max
    assert windpower._DEFAULT_HC_Q_MIN == expected_min
    # Same limits the literals approximated, to two decimals.
    assert windpower._DEFAULT_HC_Q_MAX == pytest.approx(0.48, abs=5e-3)
    assert windpower._DEFAULT_HC_Q_MIN == pytest.approx(-0.41, abs=5e-3)


def test_hc_slopes_are_pinned():
    """Regression: the Q(U) hosting-capacity slopes.

    These shifted in 0.4.1 only because the capability table moved from
    two-decimal literals to the standard's own values (0.48 -> 0.484322);
    the derivation is unchanged.
    """
    code = qc.DEFAULT_GRID_CODE
    x, y = code.vqu_v_points, code.vqu_q_max
    hc_max, _ = windpower._hc_q_bounds(code)
    last = y.shape[1] - 1

    m_qu_max = (hc_max + abs(y[1, 0])) / (x[0, 0] - x[1, 0])
    qu_max = -m_qu_max * x[1, 0] + hc_max
    m_qu_min = (abs(y[0, last]) + abs(y[1, last])) / (x[0, 0] - x[1, 0])
    qu_min = -m_qu_min * x[0, 0] + y[0, last]

    assert m_qu_max == pytest.approx(-3.2643600000000004, abs=1e-12)
    assert qu_max == pytest.approx(4.045442, abs=1e-12)
    assert m_qu_min == pytest.approx(-3.3891870833333337, abs=1e-12)
    assert qu_min == pytest.approx(3.2865200000000003, abs=1e-12)


def test_windpower_q_curves_come_from_the_registry():
    """compute_q_curves must report the envelope's own ramp segment."""
    code = qc.VDE_AR_N_4120
    curves = qc.compute_q_curves(code)
    for v in range(code.n_variants):
        ramp_m, ramp_b = min(
            code.qv_area.lower_pieces(v), key=lambda mb: mb[0]
        )
        assert curves.m_qv[v] == pytest.approx(ramp_m, abs=1e-12)
        assert curves.b_qv_min[v] == pytest.approx(ramp_b, abs=1e-12)


def test_hc_bounds_track_the_selected_grid_code():
    """A different code yields its own HC envelope, not the default's."""
    custom = qc.GridCode(
        name="test-tar",
        title="Test TAR",
        voltage_level="medium voltage",
        pq_area=qc.Envelope(
            x_points=np.array([0.1, 0.2, 1.0]),
            q_min=np.array([[-0.1, -0.10, -0.10], [-0.1, -0.25, -0.25]]),
            q_max=np.array([[0.1, 0.30, 0.30], [0.1, 0.20, 0.20]]),
        ),
        qv_area=qc.Envelope(
            x_points=qc.VDE_AR_N_4120.qv_area.x_points,
            q_min=np.array(
                [[0.30, -0.10, -0.10, -0.10], [0.20, -0.25, -0.25, -0.25]]
            ),
            q_max=np.array(
                [[0.30, 0.30, 0.30, -0.10], [0.20, 0.20, 0.20, -0.25]]
            ),
        ),
        vpu_v_curtail=1.06,
        vpu_v_max=1.10,
        cpp_p_threshold_pu=0.2,
    )
    hc_max, hc_min = windpower._hc_q_bounds(custom)
    assert hc_max == pytest.approx(0.30)
    assert hc_min == pytest.approx(-0.25)
    assert custom.n_variants == 2  # the derivation must not assume three


# ── bus numbering: ppc space vs pandapower space ──────────────────────────


def test_model_covers_every_ppc_bus(lv_rural_net):
    """model.B must span the whole ppc bus table, not a truncated slice.

    Regression test: the bus set was built from the first ``len(net.bus)``
    ppc rows.  Where pandapower's conversion inserts auxiliary buses for
    node-node switches the ppc table is longer, and that slice kept
    auxiliary buses while dropping real pandapower buses together with the
    branches attached to them.
    """
    opf = _build_sp(_annotate(lv_rural_net), pv_q_control="both")
    assert len(list(opf.model.B)) == len(opf.net._ppc["bus"])


def test_bpd_is_exactly_the_pandapower_backed_buses(lv_rural_net):
    """Bpd must equal the set of ppc buses that a pandapower bus maps onto."""
    opf = _build_sp(_annotate(lv_rural_net), pv_q_control="both")
    mapped = {int(b) for b in opf.pd_bus_to_ppc}
    assert set(opf.model.Bpd) == mapped
    assert set(opf.model.Bpd) <= set(opf.model.B)


def test_bpd_equals_b_without_auxiliary_buses(lv_rural_net):
    """On a grid with no auxiliary ppc buses nothing changes.

    This is the compatibility guarantee for the bus-numbering fix: whenever
    the ppc bus count equals ``len(net.bus)``, ``Bpd`` and ``B`` coincide and
    the voltage bounds cover every bus exactly as before.
    """
    opf = _build_sp(_annotate(lv_rural_net), pv_q_control="both")
    if len(opf.net._ppc["bus"]) == len(opf.net.bus):
        assert set(opf.model.Bpd) == set(opf.model.B)
    else:  # pragma: no cover - fixture is an auxiliary-free LV grid
        pytest.skip("fixture unexpectedly has auxiliary ppc buses")


def test_voltage_limits_are_indexed_by_ppc_bus(lv_rural_net):
    """get_v_limits must key on ppc bus numbers, matching bus_lookup.

    Positional arrays were correct only while pandapower and ppc numbering
    coincided; every consumer resolves through ``bus_lookup``.
    """
    opf = _build_sp(_annotate(lv_rural_net), pv_q_control="both")
    vmax, vmin = opf.v_limits
    assert set(vmax.index) == set(opf.model.Bpd)
    assert set(vmin.index) == set(opf.model.Bpd)
    # values must still be the ones from net.bus, resolved via the lookup
    got_max = vmax.loc[opf.pd_bus_to_ppc].to_numpy()
    assert np.allclose(got_max, opf.net.bus.max_vm_pu.values)


def test_voltage_bounds_applied_to_pandapower_backed_buses(lv_rural_net):
    """The bound constraint is indexed over Bpd, not B."""
    opf = _build_sp(_annotate(lv_rural_net), pv_q_control="both")
    assert set(opf.model.v_pyo) == set(opf.model.Bpd)


def test_results_write_back_resolves_every_bus(lv_rural_net):
    """Every bus_lookup target must exist in the voltage variable.

    This is the precondition ``pyo_to_net._bus_voltage_results_to_net``
    relies on; when it failed the write-back raised KeyError.
    """
    opf = _build_sp(_annotate(lv_rural_net), pv_q_control="both")
    v_keys = set(opf.model.v.keys())
    targets = {int(b) for b in opf.pd_bus_to_ppc}
    assert targets <= v_keys


def test_multi_period_defines_bpd(lv_rural_net):
    """The multi-period model exposes Bpd as well."""
    mp = _build_mp(_annotate(lv_rural_net))
    assert hasattr(mp.model, "Bpd")
    assert set(mp.model.Bpd) <= set(mp.model.B)


def test_multi_period_voltage_bounds_over_bpd(lv_rural_net):
    """Multi-period voltage bounds are indexed (Bpd, T)."""
    mp = _build_mp(_annotate(lv_rural_net))
    buses = {k[0] for k in mp.model.v_constraint}
    assert buses == set(mp.model.Bpd)


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
    assert opf._grid_code is qc.VDE_AR_N_4120


def test_grid_code_selection_reaches_single_period(lv_rural_net):
    """Selecting 4110 must actually resolve to 4110 in the model."""
    opf = _build_sp(
        _annotate(lv_rural_net), pv_q_control="both", grid_code="4110"
    )
    assert opf._grid_code is qc.VDE_AR_N_4110


def test_grid_code_curves_match_selected_code(lv_rural_net):
    """q_limit_parameter must come from the selected grid code."""
    opf = _build_sp(_annotate(lv_rural_net), pv_q_control="both")
    expected = qc.compute_q_curves(qc.VDE_AR_N_4120)
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
    table = qc.DEFAULT_GRID_CODE.vqu_q_max
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
    """Multi-period Q-control constraints are (sgen, time, piece)-indexed.

    The trailing piece index arrived in 0.4.1: a capability bound is a
    piecewise-linear envelope, so it needs one inequality per affine piece
    rather than one per element.
    """
    mp = _build_mp(_annotate(lv_rural_net))
    keys = list(mp.model.sG_QP_pos)
    assert keys, "expected at least one Q-controlled sgen"
    assert len(keys[0]) == 3
    times = {k[1] for k in keys}
    assert times == set(mp.model.T)
    pieces = {k[2] for k in keys}
    assert pieces == set(mp.model.sG_QP_PIECE)


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


# ── parity with pandapower's own capability areas ─────────────────────────
# Every area here should reproduce the matching class in
# pandapower.control.controller.DERController.  4105 and 4110 are shapely
# polygons upstream, so those skip without shapely; 4120 is plain branch
# logic and always runs.

PP_AREAS = "pandapower.control.controller.DERController"


def _pp_pair(code_name, variant):
    """Return pandapower's (pq_area, qv_area) for one variant of a code."""
    import importlib

    m = importlib.import_module(PP_AREAS)
    if code_name == "4105":
        return m.PQArea4105(variant + 1), m.QVArea4105(variant + 1)
    if code_name == "4110":
        return m.PQArea4110(), m.QVArea4110()
    q = [
        (-0.227902, 0.484322),
        (-0.328684, 0.410775),
        (-0.410775, 0.328684),
    ][variant]
    # potpourri encodes the 2015 active-power breakpoints (0.1 / 0.2).
    return m.PQArea4120(*q, version=2015), m.QVArea4120(*q)


@pytest.mark.parametrize(
    ("code_name", "variant"),
    [
        ("4105", 0),
        ("4105", 1),
        ("4110", 0),
        ("4120", 0),
        ("4120", 1),
        ("4120", 2),
    ],
)
def test_areas_match_pandapower(code_name, variant):
    """Our envelopes must equal pandapower's q_flexibility to machine eps.

    Sampled strictly inside each breakpoint span: the polygon classes return
    [0, 0] outside it, where numpy.interp saturates instead, and shapely's
    ``contains`` excludes the boundary itself.
    """
    if code_name in ("4105", "4110"):
        pytest.importorskip("shapely", reason="pandapower polygon areas")
    pq_pp, qv_pp = _pp_pair(code_name, variant)
    code = qc.GRID_CODES[code_name]

    lo, hi = code.pq_area.exact_range()
    p = np.linspace(lo + 1e-6, hi - 1e-6, 201)
    ours = np.column_stack(code.pq_area.q_flexibility(p, variant))
    assert np.allclose(pq_pp.q_flexibility(p), ours, atol=1e-9)

    lo, hi = code.qv_area.exact_range()
    v = np.linspace(lo + 1e-6, hi - 1e-6, 201)
    ours = np.column_stack(code.qv_area.q_flexibility(v, variant))
    assert np.allclose(
        qv_pp.q_flexibility(np.ones_like(v), v), ours, atol=1e-9
    )


def test_4110_is_no_longer_a_placeholder():
    """4110 carries its own values, not a copy of another code's.

    This replaces the tripwire that guarded the placeholder: the medium
    voltage rule now has a 0.05 p.u. active-power threshold that neither
    other code shares.
    """
    code = qc.VDE_AR_N_4110
    assert code.qp_p_high == pytest.approx(0.05)
    assert code.n_variants == 1
    assert not np.array_equal(
        code.pq_area.x_points, qc.VDE_AR_N_4120.pq_area.x_points
    )


# ── linear pieces: exactness and low-P feasibility ────────────────────────


@pytest.mark.parametrize("code_name", ["4105", "4110", "4120"])
def test_pieces_reproduce_the_envelope_inside_its_span(code_name):
    """Inside the breakpoints, min/max-of-pieces *is* the envelope."""
    code = qc.GRID_CODES[code_name]
    for area in (code.pq_area, code.qv_area):
        lo, hi = area.exact_range()
        for v in range(code.n_variants):
            for x in np.linspace(lo, hi, 51):
                env_lo, env_hi = area.q_flexibility(x, v)
                got_hi = min(m * x + b for m, b in area.upper_pieces(v))
                got_lo = max(m * x + b for m, b in area.lower_pieces(v))
                assert got_hi == pytest.approx(env_hi, abs=1e-12)
                assert got_lo == pytest.approx(env_lo, abs=1e-12)


@pytest.mark.parametrize("code_name", ["4105", "4110", "4120"])
def test_q_stays_feasible_down_to_zero_output(code_name):
    """Regression: the Q(P) bounds must not cross below the first breakpoint.

    Before 0.4.1 they did.  The lower bound extrapolated upward and the
    upper one downward, so for VDE-AR-N 4120 no Q at all satisfied both
    below P = 0.061 Pn — any curtailed sgen made the model infeasible.
    The pieces are now taken over the operating range, which replaces each
    bound by its hull where the exact area is non-convex.
    """
    code = qc.GRID_CODES[code_name]
    area = code.pq_area
    rng = qc.DEFAULT_P_RANGE_PU
    for v in range(code.n_variants):
        for p in np.linspace(*rng, 101):
            hi = min(m * p + b for m, b in area.upper_pieces(v, rng))
            lo = max(m * p + b for m, b in area.lower_pieces(v, rng))
            assert lo <= hi + 1e-12, f"{code_name} v{v}: empty band at P={p}"


@pytest.mark.parametrize("code_name", ["4105", "4110", "4120"])
def test_hull_only_ever_relaxes(code_name):
    """The hull may permit more than the grid code, never less.

    Erring toward permissive is the deliberate choice: the exact area is
    non-convex outside the breakpoints, and the strict alternative makes
    the model infeasible rather than conservative.
    """
    code = qc.GRID_CODES[code_name]
    area = code.pq_area
    rng = qc.DEFAULT_P_RANGE_PU
    for v in range(code.n_variants):
        for p in np.linspace(*rng, 101):
            env_lo, env_hi = area.q_flexibility(p, v)
            hi = min(m * p + b for m, b in area.upper_pieces(v, rng))
            lo = max(m * p + b for m, b in area.lower_pieces(v, rng))
            assert hi >= env_hi - 1e-12
            assert lo <= env_lo + 1e-12


def test_hull_is_exact_at_full_output():
    """Relaxation is confined below the reference point, not at rated P."""
    rng = qc.DEFAULT_P_RANGE_PU
    for name in ("4105", "4110", "4120"):
        code = qc.GRID_CODES[name]
        area = code.pq_area
        for v in range(code.n_variants):
            env_lo, env_hi = area.q_flexibility(1.0, v)
            hi = min(m + b for m, b in area.upper_pieces(v, rng))
            lo = max(m + b for m, b in area.lower_pieces(v, rng))
            assert hi == pytest.approx(env_hi, abs=1e-12)
            assert lo == pytest.approx(env_lo, abs=1e-12)


def test_range_outside_the_envelope_warns():
    """Leaving the exact span is flagged, not silent."""
    area = qc.VDE_AR_N_4105.qv_area  # spans 0.90-1.10
    with pytest.warns(qc.EnvelopeRangeWarning, match="exact range"):
        qc.warn_if_outside_exact_range(area, 0.85, 1.15)


def test_range_inside_the_envelope_is_quiet():
    """The usual 0.9/1.1 bus limits must not warn against a wider area."""
    area = qc.VDE_AR_N_4120.qv_area  # spans 0.87-1.15
    with warnings.catch_warnings():
        warnings.simplefilter("error", qc.EnvelopeRangeWarning)
        assert not qc.warn_if_outside_exact_range(area, 0.90, 1.10)


# ── dead-band Q(U) characteristic ─────────────────────────────────────────


def test_deadband_curve_is_zero_through_the_band():
    """Q is held at zero across the dead band and ramps outside it."""
    curve = qc.VDE_AR_N_4110.deadband_curve()
    assert curve.deadband(0) == (0.95, 1.05)
    for v in (0.95, 0.98, 1.00, 1.02, 1.05):
        assert float(curve.step(v, 0)) == pytest.approx(0.0, abs=1e-12)
    assert float(curve.step(0.90, 0)) == pytest.approx(0.484322, abs=1e-9)
    assert float(curve.step(1.10, 0)) == pytest.approx(-0.484322, abs=1e-9)


def test_deadband_curve_is_monotonic():
    """More voltage never means more capacitive reactive power."""
    curve = qc.VDE_AR_N_4120.deadband_curve()
    q = [float(curve.step(v, 0)) for v in np.linspace(0.85, 1.20, 71)]
    assert all(b <= a + 1e-12 for a, b in zip(q, q[1:]))


def test_deadband_width_is_configurable():
    """An operator parameterisation overrides the code's own plateau."""
    curve = qc.VDE_AR_N_4120.deadband_curve(deadband=(0.98, 1.02))
    assert curve.deadband(0) == (0.98, 1.02)
    assert float(curve.step(1.00, 0)) == pytest.approx(0.0, abs=1e-12)
    assert float(curve.step(0.96, 0)) > 0.0  # outside the narrower band


def test_deadband_outside_the_curve_range_is_rejected():
    """A dead band wider than the characteristic is an error."""
    with pytest.raises(ValueError, match="must lie inside"):
        qc.VDE_AR_N_4105.deadband_curve(deadband=(0.5, 1.5))


def test_deadband_curve_covers_every_variant():
    """One curve per var_q variant, matching the code."""
    for name in ("4105", "4110", "4120"):
        code = qc.GRID_CODES[name]
        assert code.deadband_curve().n_variants == code.n_variants


def test_padding_saturates_rather_than_extrapolating():
    """A curve narrower than the bus limits holds its end value."""
    curve = qc.VDE_AR_N_4105.deadband_curve()  # spans 0.90-1.10
    wide = curve.padded(0.80, 1.20)
    assert float(wide.step(0.80, 0)) == pytest.approx(
        float(curve.step(0.90, 0)), abs=1e-12
    )
    assert float(wide.step(1.20, 0)) == pytest.approx(
        float(curve.step(1.10, 0)), abs=1e-12
    )
    # Already-covering ranges are returned untouched.
    assert curve.padded(0.95, 1.05) is curve


@pytest.mark.parametrize(
    ("spec", "expected"),
    [
        (None, None),
        (False, None),
        (True, (0.95, 1.05)),
        ((0.98, 1.02), (0.98, 1.02)),
    ],
)
def test_resolve_qu_curve_accepts_each_spelling(spec, expected):
    """None/False keep the area; True and a pair build a curve."""
    curve = qc.resolve_qu_curve(spec, "4110")
    if expected is None:
        assert curve is None
    else:
        assert curve.deadband(0) == expected


def test_deadband_replaces_the_area_single_period(lv_rural_net):
    """With a dead band the Q(U) bounds give way to the characteristic."""
    opf = _build_sp(
        _annotate(lv_rural_net), pv_q_control="both", qu_deadband=(0.98, 1.02)
    )
    assert hasattr(opf.model, "pv_qu_db_pw")
    assert not hasattr(opf.model, "PV_QU_min")
    assert not hasattr(opf.model, "PV_QU_max")
    # Q(P) is unaffected: it does not depend on voltage.
    assert hasattr(opf.model, "PV_QP_pos")


def test_area_is_kept_without_a_deadband(lv_rural_net):
    """The default stays the convex area, solvable with IPOPT alone."""
    opf = _build_sp(_annotate(lv_rural_net), pv_q_control="both")
    assert hasattr(opf.model, "PV_QU_min")
    assert not hasattr(opf.model, "pv_qu_db_pw")


def test_deadband_introduces_integrality(lv_rural_net):
    """The dead band is genuinely non-convex, so it needs binaries."""
    import pyomo.environ as pyo

    opf = _build_sp(
        _annotate(lv_rural_net), pv_q_control="both", qu_deadband=True
    )
    binaries = [
        v for v in opf.model.component_data_objects(pyo.Var) if v.is_binary()
    ]
    assert binaries, "expected binary variables from the piecewise block"


def test_deadband_multi_period_is_time_indexed(lv_rural_net):
    """Every (sgen, time) pair gets its own point on the characteristic."""
    mp = _build_mp(_annotate(lv_rural_net), qu_deadband=(0.98, 1.02))
    assert hasattr(mp.model, "sG_qu_db_pw")
    assert not hasattr(mp.model, "sG_QU_min")
    keys = list(mp.model.sG_qu_db_IDX)
    assert keys and len(keys[0]) == 2
    assert {k[1] for k in keys} == set(mp.model.T)


# ── var_q validated against the selected code ─────────────────────────────


def test_var_q_beyond_the_codes_variants_is_rejected():
    """4110 defines one variant, so var_q=1 is an error, not a silent clamp."""
    with pytest.raises(ValueError, match="not valid for VDE-AR-N 4110"):
        qc.check_var_q([0, 1], "4110")


def test_var_q_within_range_is_accepted():
    """Valid indices pass, and NaN (no Q-control) is ignored."""
    qc.check_var_q([0, 1], "4105")
    qc.check_var_q([0, 1, 2, np.nan], "4120")
    qc.check_var_q([], "4110")


def test_var_q_error_names_the_code_and_the_limit():
    """The message has to say what to do about it."""
    with pytest.raises(ValueError) as err:
        qc.check_var_q([2], "4105", context="net.sgen.var_q")
    text = str(err.value)
    assert "net.sgen.var_q" in text
    assert "VDE-AR-N 4105" in text
    assert "0..1" in text


def test_var_q_is_validated_when_building_a_model(lv_rural_net):
    """The check fires from the model, not just the helper."""
    with pytest.raises(ValueError, match="not valid for VDE-AR-N 4110"):
        _build_sp(
            _annotate(lv_rural_net, var_q=2),
            pv_q_control="both",
            grid_code="4110",
        )


# ── regressions from the independent review ───────────────────────────────


@pytest.mark.parametrize("bad", [0.9, 1.7, 2.5, -0.5])
def test_fractional_var_q_is_rejected_not_rounded(bad):
    """A fractional var_q must fail, not silently pick a neighbour.

    ``int(0.9)`` is 0, so truncating would quietly select variant 0 and
    change dispatch with no error anywhere.
    """
    with pytest.raises(ValueError, match="whole number"):
        qc.check_var_q([bad], "4120")


def test_integral_float_var_q_is_accepted():
    """Pandas stores var_q as float64 whenever the column holds NaN."""
    qc.check_var_q([0.0, 1.0, 2.0], "4120")
    qc.check_var_q(np.array([0.0, 2.0]), "4120")


def test_fractional_variant_rejected_at_the_envelope():
    """The same guard on the lower-level accessor."""
    with pytest.raises(IndexError, match="whole number"):
        qc.VDE_AR_N_4120.qv_area.q_flexibility(1.0, 1.5)


def test_qu_pieces_never_demand_q_the_code_does_not(lv_rural_net):
    """Below the code's voltage span the pieces must relax, not extrapolate.

    VDE-AR-N 4105's QV area starts at 0.90 p.u.  Extrapolating its lower
    ramp to 0.85 requires Q >= +0.329 Pn, where the standard requires
    nothing at all — the model would forbid a dispatch the grid code
    permits.  Passing the bus voltage range switches to the hull.
    """
    area = qc.VDE_AR_N_4105.qv_area
    v_span = (0.85, 1.10)
    for v in np.linspace(*v_span, 51):
        env_lo, env_hi = area.q_flexibility(v, 0)
        lo = max(m * v + b for m, b in area.lower_pieces(0, v_span))
        hi = min(m * v + b for m, b in area.upper_pieces(0, v_span))
        assert lo <= env_lo + 1e-12, f"over-restrictive lower bound at v={v}"
        assert hi >= env_hi - 1e-12, f"over-restrictive upper bound at v={v}"
        assert lo <= hi + 1e-12


def test_qu_pieces_unchanged_when_limits_sit_inside_the_span():
    """Grids inside the code's own voltage range are unaffected by the hull."""
    for name in ("4105", "4110", "4120"):
        area = qc.GRID_CODES[name].qv_area
        span = area.exact_range()
        for v in range(qc.GRID_CODES[name].n_variants):
            assert area.lower_pieces(v, span) == area.lower_pieces(v)
            assert area.upper_pieces(v, span) == area.upper_pieces(v)


def test_bus_voltage_range_reads_the_network(lv_rural_net):
    """The span comes from net.bus, widened past the default when needed."""
    net = copy.deepcopy(lv_rural_net)
    net.bus["min_vm_pu"] = 0.85
    net.bus["max_vm_pu"] = 1.15
    assert qc.bus_voltage_range(net) == (0.85, 1.15)
    # Tighter-than-default limits must not narrow the span below the
    # default, or the pieces would stop covering the code's own range.
    net.bus["min_vm_pu"] = 0.97
    net.bus["max_vm_pu"] = 1.03
    assert qc.bus_voltage_range(net) == (0.9, 1.1)


def test_bus_voltage_range_falls_back_without_columns(lv_rural_net):
    """A network with no limits gets pandapower's own default band."""
    net = copy.deepcopy(lv_rural_net)
    net.bus.drop(
        columns=["min_vm_pu", "max_vm_pu"], errors="ignore", inplace=True
    )
    assert qc.bus_voltage_range(net) == qc.DEFAULT_V_RANGE_PU


def test_wide_voltage_limits_keep_the_model_feasible(lv_rural_net):
    """A grid allowing 0.85-1.15 must still build and stay satisfiable."""
    net = _annotate(lv_rural_net)
    net.bus["min_vm_pu"] = 0.85
    net.bus["max_vm_pu"] = 1.15
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", qc.EnvelopeRangeWarning)
        opf = _build_sp(net, pv_q_control="both", grid_code="4105")
    code = opf._grid_code
    v_span = qc.bus_voltage_range(net)
    for v in np.linspace(*v_span, 41):
        lo = max(m * v + b for m, b in code.qv_area.lower_pieces(0, v_span))
        hi = min(m * v + b for m, b in code.qv_area.upper_pieces(0, v_span))
        assert lo <= hi + 1e-12, f"empty Q band at v={v}"


def test_multi_period_stores_the_resolved_grid_code(lv_rural_net):
    """``_grid_code`` must mean the same on both models.

    The single-period model stored the resolved GridCode while the
    multi-period one kept the raw selector, so ``mp._grid_code.pq_area``
    raised AttributeError for anyone who read it.
    """
    mp = _build_mp(_annotate(lv_rural_net), grid_code="4110")
    assert mp._grid_code is qc.VDE_AR_N_4110
    assert mp._grid_code.pq_area is qc.VDE_AR_N_4110.pq_area
    opf = _build_sp(_annotate(lv_rural_net), grid_code="4110")
    assert type(mp._grid_code) is type(opf._grid_code)


# ── multi-period Q-control must not be a no-op ────────────────────────────


def test_multi_period_reactive_bounds_come_from_the_grid_code(lv_rural_net):
    """Regression: Q-controlled sgens must not be pinned to Q = 0.

    ``static_generation_reactive_power_limits`` derives QsGmax / QsGmin from
    the ``q_mvar`` profile, which SimBench ships as zero for PV.  Those
    bounds pinned ``qsG`` to zero, so every multi-period Q-control
    constraint — Q(P), Q(U), the inverter circle — was satisfied trivially
    and no reactive power was ever dispatched.  The model looked
    Q-controlled and did nothing.
    """
    net = _annotate(lv_rural_net)
    assert (net.sgen["q_mvar"].abs() < 1e-12).all(), (
        "fixture no longer has zero q_mvar; the regression it guards is gone"
    )
    mp = _build_mp(net)
    qcs = list(mp.model.sGqc)
    assert qcs, "expected Q-controlled sgens"
    for g in qcs:
        hi = float(pyo.value(mp.model.QsGmax[g, 0]))
        lo = float(pyo.value(mp.model.QsGmin[g, 0]))
        assert hi > 1e-9, f"sgen {g} cannot inject reactive power at all"
        assert lo < -1e-9, f"sgen {g} cannot absorb reactive power at all"


def test_multi_period_bounds_match_the_single_period_path(lv_rural_net):
    """The two paths must grant the same capability for the same net.

    They diverged: the single-period model overrode the profile-derived
    bounds from the capability table, the multi-period one did not.
    """
    net = _annotate(lv_rural_net)
    mp = _build_mp(net)
    sp = _build_sp(net, pv_q_control="both")
    sp_data = sp.static_generation_data
    for g in mp.model.sGqc:
        assert float(pyo.value(mp.model.QsGmax[g, 0])) == pytest.approx(
            float(sp_data["max_q"][g]), rel=1e-9
        )
        assert float(pyo.value(mp.model.QsGmin[g, 0])) == pytest.approx(
            float(sp_data["min_q"][g]), rel=1e-9
        )


def test_multi_period_bounds_scale_with_the_capability_table(lv_rural_net):
    """The bounds are Pn times the grid code's own Q/Pn entries."""
    net = _annotate(lv_rural_net)
    mp = _build_mp(net)
    table = mp._grid_code.vqu_q_max
    p_inst = net.sgen["p_inst_mw"].values / mp.baseMVA
    for g in mp.model.sGqc:
        v = int(net.sgen["var_q"][g])
        assert float(pyo.value(mp.model.QsGmax[g, 0])) == pytest.approx(
            table[0, v] * p_inst[g], rel=1e-9
        )
        assert float(pyo.value(mp.model.QsGmin[g, 0])) == pytest.approx(
            table[1, v] * p_inst[g], rel=1e-9
        )


def test_multi_period_bounds_track_the_selected_code(lv_rural_net):
    """A different grid code gives different reactive bounds."""
    net = _annotate(lv_rural_net)
    a = _build_mp(net, grid_code="4110")
    b = _build_mp(net, grid_code="4120")
    g = list(a.model.sGqc)[0]
    assert float(pyo.value(a.model.QsGmin[g, 0])) != pytest.approx(
        float(pyo.value(b.model.QsGmin[g, 0]))
    )


def test_sgens_without_var_q_keep_their_profile_bounds(lv_rural_net):
    """The override applies only to Q-controlled sgens."""
    net = _annotate(lv_rural_net, var_q=None)
    mp = _build_mp(net)
    assert not getattr(mp, "sgen_qc_indices", [])
    # Nothing was overridden, so the profile-derived bounds stand.
    for g in list(mp.model.sG)[:3]:
        assert float(pyo.value(mp.model.QsGmax[g, 0])) == pytest.approx(
            abs(float(net.sgen["q_mvar"][g])) / mp.baseMVA, abs=1e-12
        )


# ── regressions from the review of the shipped 0.4.1 ──────────────────────


def test_deadband_conflicting_with_qp_area_warns():
    """A characteristic the Q(P) area cannot satisfy must be flagged.

    Both are imposed on the same reactive power. For VDE-AR-N 4110 with its
    default dead band the curve assigns +0.484 Pn at 0.90 p.u., which the
    Q(P) area permits only at rated output — so a bus reaching that voltage
    at lower active power makes the model infeasible, and the solver says
    only "infeasible".
    """
    code = qc.VDE_AR_N_4110
    curve = code.deadband_curve()
    with pytest.warns(qc.QuCurveOutsidePqAreaWarning, match="disagree"):
        assert qc.warn_if_curve_leaves_pq_area(curve, code.pq_area)


def test_deadband_conflict_warning_names_the_operating_point():
    """The message has to say where it breaks and what to do."""
    code = qc.VDE_AR_N_4110
    with pytest.warns(qc.QuCurveOutsidePqAreaWarning) as rec:
        qc.warn_if_curve_leaves_pq_area(code.deadband_curve(), code.pq_area)
    text = str(rec[0].message)
    assert "0.9000" in text
    assert "+0.4843" in text
    assert "dead band" in text.lower()


def test_shallow_deadband_curve_does_not_warn():
    """A curve the Q(P) area can satisfy everywhere must stay quiet.

    Widening the dead band does not help: the curve still assigns its full
    reactive limit at its outermost breakpoints, and the Q(P) area grants
    that only at rated output.  What removes the conflict is a curve whose
    amplitude fits inside the area at zero active power.
    """
    code = qc.VDE_AR_N_4110
    at_zero = min(
        m * 0.0 + b
        for m, b in code.pq_area.upper_pieces(0, qc.DEFAULT_P_RANGE_PU)
    )
    curve = qc.deadband_qv_curve(code.qv_area.x_points, q_max=at_zero * 0.9)
    with warnings.catch_warnings():
        warnings.simplefilter("error", qc.QuCurveOutsidePqAreaWarning)
        assert not qc.warn_if_curve_leaves_pq_area(curve, code.pq_area)


def test_widening_the_deadband_does_not_remove_the_conflict():
    """Pinned so the misconception does not creep back into the docs."""
    code = qc.VDE_AR_N_4110
    v = code.qv_area.x_points
    curve = code.deadband_curve(
        deadband=(float(v[0]) + 1e-3, float(v[-1]) - 1e-3)
    )
    with pytest.warns(qc.QuCurveOutsidePqAreaWarning):
        assert qc.warn_if_curve_leaves_pq_area(curve, code.pq_area)


def test_deadband_conflict_checked_when_building_a_model(lv_rural_net):
    """The check fires from add_OPF, not just the helper."""
    with pytest.warns(qc.QuCurveOutsidePqAreaWarning):
        _build_sp(
            _annotate(lv_rural_net),
            pv_q_control="both",
            grid_code="4110",
            qu_deadband=True,
        )


def test_no_conflict_warning_without_a_deadband(lv_rural_net):
    """The area formulation cannot produce this conflict."""
    with warnings.catch_warnings():
        warnings.simplefilter("error", qc.QuCurveOutsidePqAreaWarning)
        _build_sp(
            _annotate(lv_rural_net), pv_q_control="both", grid_code="4110"
        )


def test_reactive_bound_override_refuses_to_run_out_of_order(lv_rural_net):
    """Calling the Q-control step first must fail loudly, not silently.

    The override mutates the profile-derived `QsGmax` / `QsGmin` dicts.  If
    they do not exist yet it used to return quietly, which would restore
    exactly the no-op this release fixes.
    """
    from potpourri.technologies.sgens import Sgens_multi_period

    net = _annotate(lv_rural_net)
    mp = _build_mp(net)
    sg = next(o for o in mp.flexibilities if isinstance(o, Sgens_multi_period))
    # Reproduce the out-of-order state: the profile-derived bounds the
    # override mutates have not been built.
    del sg.QsGmax_data_dict
    with pytest.raises(RuntimeError, match="must run before"):
        sg.static_generation_q_ctrl_data(net)


def test_reactive_bound_override_in_the_documented_order(lv_rural_net):
    """In the order the model uses, the override applies."""
    mp = _build_mp(_annotate(lv_rural_net))
    for g in mp.model.sGqc:
        assert float(pyo.value(mp.model.QsGmax[g, 0])) > 1e-9
