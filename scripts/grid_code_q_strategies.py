# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Grid-code selection and per-sgen Q-control strategies.

Two things this example shows:

* **Selecting a grid code.**  The technical connection rules live in a
  registry (``potpourri.technologies.q_control.GRID_CODES``) and are chosen
  model-wide with ``add_OPF(grid_code=...)``.  Section 1 solves the same
  snapshot under every registered code and reports the resulting Q(P)
  envelope and dispatch.

  VDE-AR-N 4110 is registered but **provisional**: it currently reuses the
  VDE-AR-N 4105 values as a placeholder, so selecting it emits a
  ``ProvisionalGridCodeWarning`` and its results are *not* 4110-compliant.
  This example surfaces that warning rather than silencing it.

* **Assigning a different strategy per sgen.**  Section 2 gives each PV unit
  its own controller: Q(P)/Q(U), fixed cos(phi), cos(phi)(P), and P(U)
  curtailment.  The sets are kept **disjoint** on purpose — fixed cos(phi)
  and cos(phi)(P) are both equality constraints on the same reactive power,
  so applying both to one sgen over-determines Q.

  Per-sgen assignment is a multi-period feature: there every strategy is
  driven by a ``net.sgen`` column.  In the single-period model ``pu_curtail``
  and ``cos_phi_p_profile`` are model-wide switches, so only
  ``fixed_cos_phi`` and ``var_q`` can vary per row.  Section 2 therefore uses
  the multi-period model.

Network: 1-LV-rural1--0-sw (low-voltage rural feeder, 4 rooftop PV units).

Institut für Elektrische Anlagen und Netze, Digitalisierung und
Energiewirtschaft (IAEW)
Author: Steffen Kortmann (2024)
"""

import copy
import warnings

import pyomo.environ as pyo
import simbench as sb

from potpourri.models.ACOPF_base import ACOPF
from potpourri.models_multi_period.ACOPF_multi_period import ACOPF_multi_period
from potpourri.technologies.q_control import (
    GRID_CODES,
    ProvisionalGridCodeWarning,
    compute_q_curves,
    resolve_grid_code,
)

# ── Configuration ────────────────────────────────────────────────────────────
SOLVER = "ipopt"
NET_NAME = "1-LV-rural1--0-sw"
SGEN_TYPE = "PV"  # annotate only sgens of this type
VAR_Q = 0  # capability variant (0 = widest capacitive envelope)

# Grid codes compared in section 1; None means "the default" (VDE-AR-N 4105).
GRID_CODES_TO_COMPARE = (None, "4105", "4110")

# Section 1 snapshot: peak-PV midday step of the simbench year series.
PROFILE_IDX = 13868

# Section 2 window (15-min steps). Q(P) only bites while P > 0, so this
# window is chosen to sit through solar noon.
FROM_T = 13855
TO_T = 13879

# Section 2 controller parameters
COS_PHI_FIXED = 0.95  # fixed cos(phi) target
COS_PHI_MIN_CPP = 0.90  # cos(phi) at full output for the cos(phi)(P) profile
V_CURTAIL_PU = 1.06  # P(U): curtailment onset
V_MAX_CURTAIL_PU = 1.10  # P(U): zero-output voltage

VOLTAGE_BAND = (0.95, 1.05)


# ── Helpers ──────────────────────────────────────────────────────────────────
def _pv_mask(net):
    """Boolean mask of the sgens this example controls."""
    return net.sgen["type"] == SGEN_TYPE


def _prepare(net):
    """Common annotation: installed capacity, bounds, voltage band."""
    net.sgen["p_inst_mw"] = net.sgen["p_mw"].abs()
    net.sgen["controllable"] = True
    net.bus["min_vm_pu"], net.bus["max_vm_pu"] = (
        VOLTAGE_BAND[0],
        VOLTAGE_BAND[1],
    )
    return net


def compare_grid_codes(
    base_net, grid_codes=GRID_CODES_TO_COMPARE, solver=SOLVER
):
    """Solve one snapshot under each grid code and report the envelope.

    Args:
        base_net: pandapower network to copy per run.
        grid_codes: iterable of grid-code selectors (``None`` = default).
        solver: Pyomo solver name.

    Returns:
        dict keyed by the selector, holding objective, voltage band and the
        Q(P) plateau of the selected code.
    """
    print("=" * 70)
    print("1. GRID-CODE SELECTION")
    print("=" * 70)
    print(f"registered grid codes: {sorted(GRID_CODES)}\n")

    results = {}
    for selector in grid_codes:
        net = _prepare(copy.deepcopy(base_net))
        net.sgen["var_q"] = None
        net.sgen.loc[_pv_mask(net), "var_q"] = VAR_Q

        # Capture the provisional warning instead of hiding it.
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            opf = ACOPF(net)
            opf.add_OPF(pv_q_control="both", grid_code=selector)
            provisional = [
                w
                for w in caught
                if issubclass(w.category, ProvisionalGridCodeWarning)
            ]

        code = resolve_grid_code_quietly(selector)
        opf.add_voltage_deviation_objective()
        opf.solve(solver=solver, print_solver_output=False)

        v = [pyo.value(opf.model.v[b]) for b in opf.model.B]
        curves = compute_q_curves(code)
        # Q(P) envelope evaluated at the capability reference point.
        q_at_ref = (
            curves.b_qp_max[VAR_Q] + curves.m_qp_max[VAR_Q] * code.qp_p_low
        )
        label = "default" if selector is None else str(selector)
        results[label] = {
            "code": code.title,
            "provisional": bool(provisional),
            "obj": pyo.value(opf.model.obj_v_deviation),
            "vmin": min(v),
            "vmax": max(v),
            "q_over_pn_at_ref": q_at_ref,
        }

        print(f"selector {label!r} -> {code.title} ({code.voltage_level})")
        print(f"  variants available     : {code.n_variants}")
        print(f"  Q/Pn at P = {code.qp_p_low:.2f} Pn : {q_at_ref:+.4f}")
        print(
            f"  objective sum (v-1)^2  : "
            f"{pyo.value(opf.model.obj_v_deviation):.6f}"
        )
        print(f"  voltage band           : [{min(v):.4f}, {max(v):.4f}] p.u.")
        if provisional:
            print(
                "  PROVISIONAL: placeholder values, results are NOT "
                f"{code.title}-compliant"
            )
        print()

    return results


def resolve_grid_code_quietly(selector):
    """Resolve a selector without re-emitting the provisional warning."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ProvisionalGridCodeWarning)
        return resolve_grid_code(selector)


def assign_per_sgen_strategies(base_net, solver=SOLVER):
    """Give each PV sgen a different Q-control strategy and solve.

    The strategy sets are disjoint: no sgen receives two equality
    constraints on its reactive power.

    Returns:
        dict with the assignment table and the solved dispatch.
    """
    print("=" * 70)
    print("2. ONE STRATEGY PER SGEN  (multi-period)")
    print("=" * 70)

    net = _prepare(copy.deepcopy(base_net))
    pv = list(net.sgen.index[_pv_mask(net)])
    if len(pv) < 4:
        print(
            f"note: only {len(pv)} {SGEN_TYPE} sgens; strategies are "
            "assigned round-robin over the ones available"
        )

    # Round-robin so the example still works on networks with fewer sgens.
    strategies = ["qp_qu", "fixed_cos_phi", "cos_phi_p", "pu_curtail"]
    assignment = {g: strategies[i % len(strategies)] for i, g in enumerate(pv)}

    # Every column starts empty, then each sgen gets exactly one strategy.
    net.sgen["var_q"] = None
    net.sgen["fixed_cos_phi"] = float("nan")
    net.sgen["cos_phi_p_profile"] = False
    net.sgen["cos_phi_min"] = float("nan")
    net.sgen["pu_curtail"] = False
    net.sgen["v_curtail_pu"] = V_CURTAIL_PU
    net.sgen["v_max_curtail_pu"] = V_MAX_CURTAIL_PU
    net.sgen["max_p_mw"] = net.sgen["p_mw"].abs()
    net.sgen["min_p_mw"] = 0.0

    for g, strategy in assignment.items():
        if strategy == "qp_qu":
            net.sgen.at[g, "var_q"] = VAR_Q
        elif strategy == "fixed_cos_phi":
            net.sgen.at[g, "fixed_cos_phi"] = COS_PHI_FIXED
        elif strategy == "cos_phi_p":
            net.sgen.at[g, "cos_phi_p_profile"] = True
            net.sgen.at[g, "cos_phi_min"] = COS_PHI_MIN_CPP
        elif strategy == "pu_curtail":
            # P(U) limits active power, so it may coexist with a Q rule.
            net.sgen.at[g, "pu_curtail"] = True
            net.sgen.at[g, "var_q"] = VAR_Q

    print("\nassignment:")
    for g, strategy in assignment.items():
        print(f"  sgen {g}: {strategy}")

    mpopf = ACOPF_multi_period(net, toT=TO_T, fromT=FROM_T)
    mpopf.add_OPF()
    mpopf.add_voltage_deviation_objective()

    built = [
        name
        for name in (
            "sG_QP_pos",
            "sG_QU_max",
            "sgen_fixed_cos_phi",
            "sgen_cpp",
            "sgen_pu_curtail",
        )
        if hasattr(mpopf.model, name)
    ]
    print("\nconstraint blocks built:", ", ".join(built) or "none")

    res = mpopf.solve(solver=SOLVER, print_solver_output=False, to_net=False)
    print("termination :", res.solver.termination_condition)

    base = mpopf.model.baseMVA
    t0 = min(mpopf.model.T)
    print(f"\ndispatch at first step (t={t0}):")
    print(
        f"  {'sgen':>5}  {'strategy':<15} {'P (kW)':>9} {'Q (kvar)':>10}"
        f" {'cos(phi)':>9}"
    )
    for g, strategy in assignment.items():
        p = pyo.value(mpopf.model.psG[g, t0]) * base * 1e3
        q = pyo.value(mpopf.model.qsG[g, t0]) * base * 1e3
        s = (p**2 + q**2) ** 0.5
        cos_phi = abs(p) / s if s > 1e-9 else float("nan")
        print(f"  {g:>5}  {strategy:<15} {p:>9.3f} {q:>10.3f} {cos_phi:>9.4f}")

    print(
        "\nReading the table: this feeder sits comfortably inside its voltage\n"
        "band, so the objective has little reason to spend reactive power and\n"
        "Q stays near zero for the bound-type strategies. The fixed-cos(phi)\n"
        "unit behaves differently on purpose: its equality ties Q rigidly to\n"
        "P, and since every bus is already above 1.0 p.u. the only way to\n"
        "reduce its reactive injection is to reduce active power, so it\n"
        "curtails towards zero. That coupling is the practical cost of a\n"
        "fixed power factor, and it is why the bound-type rules (Q(P), Q(U))\n"
        "are usually preferred when active yield matters."
    )

    v_all = [
        pyo.value(mpopf.model.v[b, t])
        for b in mpopf.model.B
        for t in mpopf.model.T
    ]
    print(
        f"\nvoltage band over horizon: "
        f"[{min(v_all):.4f}, {max(v_all):.4f}] p.u."
    )
    print(
        f"objective sum (v-1)^2    : "
        f"{pyo.value(mpopf.model.obj_v_deviation):.6f}"
    )

    return {"assignment": assignment, "built": built}


def main():
    """Run the analysis this script demonstrates.

    Configuration comes from the module-level constants above, not from the
    command line. Edit those, or import and call this function, to change what
    is run.

    Returns:
        None. Results are printed, and written to the paths named in the
        configuration block where the script produces files.
    """
    base_net = sb.get_simbench_net(NET_NAME)
    compare_grid_codes(base_net)
    print()
    assign_per_sgen_strategies(base_net)
    print("\n" + "=" * 70)
    print(
        "Takeaway: grid_code selects the capability envelope model-wide, "
        "while\nthe net.sgen columns select which controller each sgen "
        "follows. Keep the\ntwo equality strategies (fixed cos(phi), "
        "cos(phi)(P)) on disjoint sgens."
    )
    print("=" * 70)


if __name__ == "__main__":
    main()
