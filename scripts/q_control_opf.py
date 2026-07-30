"""VDE-AR-N 4105 reactive-power control — single- and multi-period example.

Demonstrates how to annotate a simbench network with grid-code Q-control data
(``var_q``, ``p_inst_mw``) and solve three variants:

* **Single-period** AC OPF — compares four Q-control modes side-by-side:
  uncontrolled, Q(P) only, Q(U) only, and Q(P)+Q(U) combined.  For each mode
  the objective value, voltage band, and per-sgen Q dispatch are shown; the
  Q(P) bounds are verified analytically.

* **PV inverter control modes** — three additional single-period modes that
  further shape the inverter (P, Q) operating region:

  * P(U) active-power curtailment (VDE-AR-N 4105 §8.5): reduces active output
    when bus voltage exceeds ``V_CURTAIL_PU``, reaching zero at
    ``V_MAX_CURTAIL_PU``.
  * Fixed cos(φ): locks the Q/P ratio at the prescribed power factor
    (equality constraint — inverter tracks the target exactly).
  * cos(φ)(P) profile: Q = 0 below a threshold *P*_t, then ramps with P up
    to Q_max = P_n · tan(arccos(cos_phi_min)) at full output.

* **Multi-period** AC OPF over a 24-step daytime horizon.  Q-control
  constraints are added automatically for all sgens with ``var_q`` set via
  :class:`~potpourri.models_multi_period.ACOPF_multi_period.ACOPF_multi_period`
  — no extra flag needed.

Grid-code background
--------------------
VDE-AR-N 4105 / BDEW requires LV/MV generators to provide reactive power
according to either a Q(P) characteristic (reactive power proportional to
active output) or a Q(U) droop (reactive power proportional to local voltage).
The ``var_q`` variant index (0–2) selects the Q/P envelope column:

    variant 0  Qmax =  0.48 Pn  (most permissive)
    variant 1  Qmax =  0.41 Pn
    variant 2  Qmax =  0.33 Pn  (strictest)

Network: 1-LV-rural1--0-sw (low-voltage rural feeder, 4 rooftop PV units).

Institut für Elektrische Anlagen und Netze, Digitalisierung und
Energiewirtschaft (IAEW)
(c) 2024, Steffen Kortmann
"""

import copy
import warnings

import pyomo.environ as pe
import simbench as sb

from potpourri.models.ACOPF_base import ACOPF
from potpourri.models_multi_period.ACOPF_multi_period import ACOPF_multi_period
from potpourri.technologies.q_control import compute_q_curves

warnings.filterwarnings("ignore")

# ── Configuration ────────────────────────────────────────────────────────────
SOLVER = "ipopt"
NET_NAME = "1-LV-rural1--0-sw"
SGEN_TYPE = "PV"  # annotate only sgens of this type with Q-control
VAR_Q = 0  # VDE-AR-N 4105 variant (0=Qmax 0.48 Pn, 1=0.41, 2=0.33)

# Single-period: profile index to use as the snapshot
# 13868 is a peak-PV midday step in the 35 136-step simbench year series.
PROFILE_IDX = 13868

# Multi-period: time window (15-min steps, 24 steps ≈ 6 h centred on solar noon)
# Q(P) constraints only make sense when P > 0.  Window 13855–13878 is the first
# 24-step run where every PV sgen exceeds its Q(P) crossover threshold (~3–6 kW).
FROM_T = 13855  # first step where every PV sgen exceeds the Q(P) crossover threshold
TO_T = 13879  # 24 steps = 6 hours through solar noon

# Section 3: PV inverter controller modes
COS_PHI_FIXED = 0.95  # fixed power factor for fixed-cos(φ) mode
COS_PHI_MIN_CPP = 0.90  # power factor at full output for cos(φ)(P) profile
V_CURTAIL_PU = 1.06  # P(U): voltage at which curtailment begins (p.u.)
V_MAX_CURTAIL_PU = 1.10  # P(U): voltage at which active output reaches zero


# ── Helpers ──────────────────────────────────────────────────────────────────
def annotate_sgens_with_q_control(net, sgen_type: str = "PV", var_q: int = 0):
    """Add ``var_q`` and ``p_inst_mw`` columns to sgens of the given type.

    Non-matching rows get ``None`` so that the underlying ``pd.notna`` guard
    in :meth:`Sgens_multi_period.static_generation_q_ctrl_data` skips them.
    The ``p_inst_mw`` column is set to the nominal ``p_mw`` value (the
    per-unit installed capacity used as Pn in the Q-control characteristic).
    """
    # Initialise as object-dtype column so non-annotated rows carry None,
    # not NaN — this keeps the `pd.notna` filter correct regardless of solver.
    net.sgen["var_q"] = None
    mask = net.sgen["type"] == sgen_type
    net.sgen.loc[mask, "var_q"] = int(var_q)
    net.sgen["p_inst_mw"] = net.sgen["p_mw"].abs()
    return net


def solve_sp(net_snapshot, mode, solver):
    """Build and solve a single-period ACOPF for one Q-control mode.

    Args:
        net_snapshot: pandapower network with snapshot values already applied.
        mode: ``pv_q_control`` argument — ``None``, ``"qp"``, ``"qu"``, or
            ``"both"``.
        solver: Pyomo solver name.

    Returns:
        dict with keys ``term``, ``obj``, ``v_min``, ``v_max``, ``q_dispatch``
        (list of ``(g, p_kw, q_kvar, qmax_kvar, qmin_kvar)`` per PV sgen).
    """
    ac = ACOPF(copy.deepcopy(net_snapshot))
    ac.add_OPF(pv_q_control=mode)
    ac.add_voltage_deviation_objective()
    res = ac.solve(solver=solver, print_solver_output=False)
    term = (
        res.solver.termination_condition.value
        if res is not None
        else "no_result"
    )
    if term != "optimal":
        return {"term": term}

    base = ac.model.baseMVA
    qc = compute_q_curves()
    v_vals = [pe.value(ac.model.v[b]) for b in ac.model.B]

    q_dispatch = []
    pvc = list(ac.model.PVc) if hasattr(ac.model, "PVc") else []
    # When mode is None, report all controllable sgens instead.
    sgen_list = pvc if pvc else list(ac.model.sGc)
    for g in sgen_list:
        p_pu = pe.value(ac.model.psG[g])
        q_pu = pe.value(ac.model.qsG[g])
        # Q(P) bounds (uses p_inst if available, else p_mw)
        if pvc:
            pi_pu = pe.value(ac.model.PV_p_inst[g])
            v_idx = int(pe.value(ac.model.PV_var_q[g]))
        else:
            pi_pu = p_pu  # fallback: no p_inst param in uncontrolled case
            v_idx = VAR_Q
        qmax_pu = qc.b_qp_max[v_idx] * pi_pu + qc.m_qp_max[v_idx] * p_pu
        qmin_pu = qc.b_qp_min[v_idx] * pi_pu + qc.m_qp_min[v_idx] * p_pu
        q_dispatch.append(
            (
                g,
                p_pu * base * 1e3,
                q_pu * base * 1e3,
                qmax_pu * base * 1e3,
                qmin_pu * base * 1e3,
            )
        )

    return {
        "term": term,
        "obj": pe.value(ac.model.obj_v_deviation),
        "v_min": min(v_vals),
        "v_max": max(v_vals),
        "q_dispatch": q_dispatch,
    }


def solve_sp_ctrl(net_snapshot, add_opf_kwargs, solver):
    """Build and solve a single-period ACOPF with arbitrary add_OPF kwargs.

    Args:
        net_snapshot: pandapower network with snapshot values applied.
        add_opf_kwargs: keyword arguments forwarded verbatim to ``add_OPF()``.
        solver: Pyomo solver name.

    Returns:
        dict with ``term``, ``obj``, ``v_min``, ``v_max``, ``p_q``
        (list of ``(g, p_kw, q_kvar)`` for every controllable sgen).
    """
    ac = ACOPF(copy.deepcopy(net_snapshot))
    ac.add_OPF(**add_opf_kwargs)
    ac.add_voltage_deviation_objective()
    res = ac.solve(solver=solver, print_solver_output=False)
    term = (
        res.solver.termination_condition.value
        if res is not None
        else "no_result"
    )
    if term != "optimal":
        return {"term": term}

    base = ac.model.baseMVA
    v_vals = [pe.value(ac.model.v[b]) for b in ac.model.B]
    p_q = [
        (
            g,
            pe.value(ac.model.psG[g]) * base * 1e3,
            pe.value(ac.model.qsG[g]) * base * 1e3,
        )
        for g in sorted(ac.model.sGc)
    ]
    return {
        "term": term,
        "obj": pe.value(ac.model.obj_v_deviation),
        "v_min": min(v_vals),
        "v_max": max(v_vals),
        "p_q": p_q,
    }


if __name__ == "__main__":
    # ── Load network ─────────────────────────────────────────────────────────
    net_base = sb.get_simbench_net(NET_NAME)
    profiles = sb.get_absolute_values(
        net_base, profiles_instead_of_study_cases=True
    )

    net_base.bus["max_vm_pu"] = 1.05
    net_base.bus["min_vm_pu"] = 0.95
    net_base.line["max_loading_percent"] = 80.0
    net_base.sgen["controllable"] = True
    net_base.ext_grid["max_q_mvar"] = 500.0
    net_base.ext_grid["min_q_mvar"] = -500.0

    n_pv = (net_base.sgen["type"] == SGEN_TYPE).sum()
    print(f"Network  : {NET_NAME}")
    print(f"PV sgens : {n_pv}  (type='{SGEN_TYPE}', var_q variant {VAR_Q})")
    print(f"Solver   : {SOLVER}\n")

    # ── 1. Single-period comparison ───────────────────────────────────────────
    print("=" * 62)
    print(f"1. SINGLE-PERIOD COMPARISON  (snapshot t={PROFILE_IDX})")
    print("=" * 62)

    net_sp = copy.deepcopy(net_base)
    net_sp.sgen["p_mw"] = profiles[("sgen", "p_mw")].iloc[PROFILE_IDX]
    net_sp.load["p_mw"] = profiles[("load", "p_mw")].iloc[PROFILE_IDX]
    net_sp.load["q_mvar"] = profiles[("load", "q_mvar")].iloc[PROFILE_IDX]
    annotate_sgens_with_q_control(net_sp, SGEN_TYPE, VAR_Q)

    qc = compute_q_curves()
    modes = [
        (None, "uncontrolled"),
        ("qp", "Q(P) only  "),
        ("qu", "Q(U) only  "),
        ("both", "Q(P)+Q(U)  "),
    ]

    results = {}
    for mode, label in modes:
        print(f"  Solving {label} ...", end=" ", flush=True)
        results[label] = solve_sp(net_sp, mode, SOLVER)
        print(results[label]["term"])

    # ── Summary table ─────────────────────────────────────────────────────────
    print()
    print(f"  {'Mode':<14}  {'Objective':>12}  {'vm_min':>8}  {'vm_max':>8}")
    print("  " + "-" * 48)
    for _, label in modes:
        r = results[label]
        if r["term"] != "optimal":
            print(f"  {label}  {'— ' + r['term']:>12}")
            continue
        print(
            f"  {label}  {r['obj']:>12.6f}  {r['v_min']:>8.4f}  {r['v_max']:>8.4f}"
        )

    # ── Per-sgen Q dispatch table ─────────────────────────────────────────────
    print()
    print(
        f"  Per-sgen reactive dispatch  "
        f"[Q(P) bounds from variant {VAR_Q}: "
        f"Qmax = {qc.b_qp_max[VAR_Q]:.2f}·Pn,"
        f"  Qmin = {qc.b_qp_min[VAR_Q]:.2f}·Pn]"
    )
    print(
        f"  {'g':>3}  {'P (kW)':>7}  "
        + "  ".join(f"{'Q ' + lbl.strip():>12}" for _, lbl in modes)
        + f"  {'Qmax_QP':>10}  {'Qmin_QP':>10}"
    )
    print("  " + "-" * (3 + 2 + 7 + 2 + 12 * len(modes) + 2 * len(modes) + 24))

    # Gather sgen indices from the uncontrolled result (all sGc present there)
    ref = results["uncontrolled"]
    if ref["term"] == "optimal":
        for row in ref["q_dispatch"]:
            g, p_kw = row[0], row[1]
            qmax_kvar, qmin_kvar = row[3], row[4]
            q_vals = []
            for _, label in modes:
                r = results[label]
                if r["term"] != "optimal":
                    q_vals.append(float("nan"))
                    continue
                match = next(
                    (x[2] for x in r["q_dispatch"] if x[0] == g), float("nan")
                )
                q_vals.append(match)
            print(
                f"  {g:>3}  {p_kw:>7.2f}  "
                + "  ".join(f"{q:>12.3f}" for q in q_vals)
                + f"  {qmax_kvar:>10.3f}  {qmin_kvar:>10.3f}"
            )

    # ── 2. PV inverter control modes ──────────────────────────────────────────
    print()
    print("=" * 62)
    print(f"2. PV INVERTER CONTROL MODES  (snapshot t={PROFILE_IDX})")
    print("=" * 62)
    print(
        f"  P(U) thresholds : V_curtail={V_CURTAIL_PU} p.u., "
        f"V_max={V_MAX_CURTAIL_PU} p.u."
    )
    print(f"  Fixed cos(φ)    : {COS_PHI_FIXED}")
    print(
        f"  cos(φ)(P) profile: cos_phi_min={COS_PHI_MIN_CPP}"
        f"  (Q=0 at P≤0.2·Pn, ramps to Q_max at full output)\n"
    )

    # net_sp already has p_inst_mw; add per-sgen controller columns.
    net_sp2 = copy.deepcopy(net_sp)
    pv_mask2 = net_sp2.sgen["type"] == SGEN_TYPE
    net_sp2.sgen.loc[pv_mask2, "cos_phi_min"] = COS_PHI_MIN_CPP
    net_sp2.sgen.loc[pv_mask2, "v_curtail_pu"] = V_CURTAIL_PU
    net_sp2.sgen.loc[pv_mask2, "v_max_curtail_pu"] = V_MAX_CURTAIL_PU

    ctrl_modes = [
        ("P(U) curtail   ", {"pu_curtail": True}),
        (f"Fixed cos(φ)={COS_PHI_FIXED}", {"fixed_cos_phi": COS_PHI_FIXED}),
        (f"cos(φ)(P) pf={COS_PHI_MIN_CPP} ", {"cos_phi_p_profile": True}),
    ]

    ctrl_results = {}
    for label, kwargs in ctrl_modes:
        print(f"  Solving {label} ...", end=" ", flush=True)
        ctrl_results[label] = solve_sp_ctrl(net_sp2, kwargs, SOLVER)
        print(ctrl_results[label]["term"])

    print()
    print(f"  {'Mode':<22}  {'Objective':>12}  {'vm_min':>8}  {'vm_max':>8}")
    print("  " + "-" * 56)
    for label, _ in ctrl_modes:
        r = ctrl_results[label]
        if r["term"] != "optimal":
            print(f"  {label:<22}  {'— ' + r['term']:>12}")
            continue
        print(
            f"  {label:<22}  {r['obj']:>12.6f}"
            f"  {r['v_min']:>8.4f}  {r['v_max']:>8.4f}"
        )

    print()
    print(
        f"  {'g':>3}  {'P (kW)':>7}  "
        + "  ".join(f"{'Q ' + l.strip()[:11]:>13}" for l, _ in ctrl_modes)
    )
    print("  " + "-" * (3 + 2 + 7 + 2 + 15 * len(ctrl_modes)))
    for g in sorted(net_sp2.sgen.index):
        if not net_sp2.sgen.at[g, "controllable"]:
            continue
        p_ref = None
        q_vals_row = []
        for label, _ in ctrl_modes:
            r = ctrl_results[label]
            if r["term"] != "optimal":
                q_vals_row.append(float("nan"))
                continue
            match = next((row for row in r["p_q"] if row[0] == g), None)
            if match is not None:
                if p_ref is None:
                    p_ref = match[1]
                q_vals_row.append(match[2])
            else:
                q_vals_row.append(float("nan"))
        p_kw_disp = p_ref if p_ref is not None else 0.0
        print(
            f"  {g:>3}  {p_kw_disp:>7.2f}  "
            + "  ".join(f"{q:>13.3f}" for q in q_vals_row)
        )

    # ── 3. Multi-period AC OPF ────────────────────────────────────────────────
    print()
    print("=" * 62)
    print(
        f"3. MULTI-PERIOD  (t={FROM_T}…{TO_T - 1}, "
        f"{TO_T - FROM_T} × 15 min = "
        f"{(TO_T - FROM_T) * 15 // 60} h)"
    )
    print("=" * 62)

    net_mp = copy.deepcopy(net_base)
    net_mp.sgen["max_p_mw"] = net_base.sgen["p_mw"]
    net_mp.sgen["min_p_mw"] = 0.0
    # Allow non-zero reactive dispatch: the simbench default q_mvar=0 would
    # pin QsGmin=QsGmax=0 via static_generation_reactive_power_bounds, making
    # the Q(P) upper bound (which goes negative at low P) infeasible.
    net_mp.sgen["q_mvar"] = net_base.sgen["p_mw"].abs() * 0.5

    # Setting var_q on net.sgen is enough — ACOPF_multi_period picks it up
    # automatically in _calc_opf_parameters via static_generation_q_ctrl_data.
    annotate_sgens_with_q_control(net_mp, SGEN_TYPE, VAR_Q)

    mpopf = ACOPF_multi_period(net_mp, toT=TO_T, fromT=FROM_T)
    mpopf.add_OPF()
    mpopf.add_voltage_deviation_objective()

    print("Solving ...", flush=True)
    res_mp = mpopf.solve(solver=SOLVER, print_solver_output=False)
    term_mp = (
        res_mp.solver.termination_condition.value
        if res_mp is not None
        else "no_result"
    )
    print(f"Termination : {term_mp}")

    if term_mp == "optimal":
        base = mpopf.model.baseMVA

        # Q-controlled sgens are listed in model.sGqc (subset of model.sGc).
        qctrl = (
            sorted(mpopf.model.sGqc) if hasattr(mpopf.model, "sGqc") else []
        )
        print(
            f"Q-controlled sgens (model.sGqc): {qctrl}"
            f"  [variant {VAR_Q}: Qmax = {qc.b_qp_max[VAR_Q]:.2f}·Pn]\n"
        )

        # Time-series dispatch for the first Q-controlled sgen
        if qctrl:
            g0 = qctrl[0]
            print(f"Time-series dispatch — sgen {g0}:")
            print(
                f"  {'t':>5}  {'P (kW)':>8}  {'Q (kvar)':>10}  "
                f"{'vm_max (p.u.)':>14}"
            )
            for t in sorted(mpopf.model.T):
                p_kw = pe.value(mpopf.model.psG[g0, t]) * base * 1e3
                q_kvar = pe.value(mpopf.model.qsG[g0, t]) * base * 1e3
                vm_max = max(
                    pe.value(mpopf.model.v[b, t]) for b in mpopf.model.B
                )
                print(
                    f"  {t:>5}  {p_kw:>8.2f}  {q_kvar:>10.2f}  {vm_max:>14.4f}"
                )

        v_max = max(
            pe.value(mpopf.model.v[b, t])
            for b in mpopf.model.B
            for t in mpopf.model.T
        )
        v_min = min(
            pe.value(mpopf.model.v[b, t])
            for b in mpopf.model.B
            for t in mpopf.model.T
        )
        obj_mp = pe.value(mpopf.model.obj_v_deviation)
        print(
            f"\nVoltage band over full horizon: [{v_min:.4f}, {v_max:.4f}] p.u."
        )
        print(f"Objective Σ(v−1)²: {obj_mp:.6f}")
        print(
            "\nKey takeaway: adding var_q + q_mvar to net.sgen is sufficient — "
            "ACOPF_multi_period automatically adds sG_QP_pos/neg and "
            "sG_QU_min/max constraints for all annotated sgens.  "
            "(q_mvar > 0 is required so that the default Q=0 box constraint "
            "does not conflict with the Q(P) upper bound at low P dispatch.)"
        )
