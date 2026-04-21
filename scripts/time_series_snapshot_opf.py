"""Time-series snapshot OPF demo for potpourri.

Demonstrates how potpourri can solve an independent AC-OPF for each of a
small set of representative time snapshots drawn from the simbench yearly
profile.  Eight snapshots are chosen to sample a range of seasonal and
diurnal operating conditions on the 0.4 kV LV rural1 network.

Workflow per snapshot
---------------------
1. Pull absolute load and PV generation values from the simbench
   15-min profile (96 slots per day, 35 040 slots per year).
2. Derive PV Q-capability from a fixed cos(φ) = 0.95 envelope.
3. Allow active-power curtailment: P ∈ [0, p_rated].
4. Build an ACOPF model with voltage bounds [0.95, 1.06] p.u.
5. Solve with IPOPT, minimising Σ_b (v[b] − 1)².
6. Record ext-grid P/Q, voltage band, max line loading, and PV dispatch.

Important: every snapshot is solved as a fully independent optimisation
problem.  There is no warm-starting, no inter-snapshot coupling, and no
temporal linkage (batteries, ramp constraints, storage state of charge).
This is the "snapshot OPF" paradigm — a lightweight way to test a handful
of representative operating points.  For coupled multi-period optimisation
see ACOPF_multi_period in src/potpourri/models_multi_period/.

Institut für Elektrische Anlagen und Netze, Digitalisierung und
Energiewirtschaft (IAEW)
(c) 2023, Steffen Kortmann
"""

import copy
import os
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import simbench as sb

from potpourri.models.ACOPF_base import ACOPF

warnings.filterwarnings("ignore")

NET_NAME = "1-LV-rural1--0-sw"
SOLVER = "ipopt"
RESULTS_DIR = "results"

PV_SCALE = 3.0  # multiply profile PV to create realistic voltage pressure
PF_MIN = 0.95  # PV inverter minimum power factor
VM_MAX = 1.06  # upper voltage bound [p.u.]
VM_MIN = 0.95  # lower voltage bound [p.u.]

SLOTS_PER_DAY = 96  # 15-min resolution: 4 slots/hour × 24 h


def _slot(day: int, hour: int) -> int:
    """Day-of-year (0-based) + hour (0–23) → 15-min slot index."""
    return day * SLOTS_PER_DAY + hour * 4


# ── snapshot catalogue ─────────────────────────────────────────────────────
# Eight time slots chosen to represent diverse seasonal / diurnal conditions.
# Slot indices are for a 365-day, 15-min dataset (35 040 steps).
SNAPSHOTS = [
    (_slot(0, 12), "Jan 1 noon"),  # winter midday
    (_slot(12, 18), "Jan 13 eve"),  # winter evening peak load
    (_slot(102, 12), "Apr 13 noon"),  # spring noon, moderate PV
    (_slot(177, 0), "Jun 26 mid"),  # summer midnight, near-zero PV
    (_slot(196, 12), "Jul 15 noon"),  # summer noon, high PV
    (_slot(268, 12), "Sep 25 noon"),  # autumn noon, fading PV
    (_slot(307, 18), "Nov 3 eve"),  # autumn evening, low PV + load
    (_slot(349, 18), "Dec 15 eve"),  # winter evening, high load
]


# ── per-snapshot network setup ─────────────────────────────────────────────


def _build_snapshot_net(base_net, profiles, slot_idx):
    """Return a deep-copied net with profile values and OPF bounds applied."""
    net = copy.deepcopy(base_net)

    net.load["p_mw"] = profiles[("load", "p_mw")].iloc[slot_idx]
    net.load["q_mvar"] = profiles[("load", "q_mvar")].iloc[slot_idx]
    net.sgen["p_mw"] = profiles[("sgen", "p_mw")].iloc[slot_idx] * PV_SCALE

    # PV Q-capability at minimum power factor
    tan_phi = np.sqrt(1 - PF_MIN**2) / PF_MIN
    net.sgen["max_q_mvar"] = tan_phi * net.sgen["p_mw"]
    net.sgen["min_q_mvar"] = -tan_phi * net.sgen["p_mw"]

    # P bounds: allow curtailment to zero
    net.sgen["max_p_mw"] = net.sgen["p_mw"]
    net.sgen["min_p_mw"] = 0.0
    net.sgen["controllable"] = True

    net.bus["max_vm_pu"] = VM_MAX
    net.bus["min_vm_pu"] = VM_MIN

    return net


# ── result extraction ──────────────────────────────────────────────────────


def _extract(slot_idx, label, ac, status):
    """Return a results dict from a solved (or failed) ACOPF instance."""
    if status != "optimal":
        return {
            "slot": slot_idx,
            "label": label,
            "status": status,
            "eg_p_mw": np.nan,
            "eg_q_mvar": np.nan,
            "vm_min_pu": np.nan,
            "vm_max_pu": np.nan,
            "max_loading_%": np.nan,
            "pv_mw": np.nan,
            "load_mw": np.nan,
        }

    net = ac.net
    slack_bus = net.ext_grid.bus.iloc[0]
    vm = net.res_bus.vm_pu.drop(index=slack_bus)

    return {
        "slot": slot_idx,
        "label": label,
        "status": status,
        "eg_p_mw": round(net.res_ext_grid.p_mw.sum(), 4),
        "eg_q_mvar": round(net.res_ext_grid.q_mvar.sum(), 4),
        "vm_min_pu": round(vm.min(), 4),
        "vm_max_pu": round(vm.max(), 4),
        "max_loading_%": round(net.res_line.loading_percent.max(), 2),
        "pv_mw": round(net.res_sgen.p_mw.sum(), 4),
        "load_mw": round(net.res_load.p_mw.sum(), 4),
    }


# ── plot ───────────────────────────────────────────────────────────────────


def _plot(df, save_path):
    """Four-panel chart: active import, PV dispatch, line loading, voltage band."""
    x = np.arange(len(df))
    tick_labels = df["label"].tolist()
    bar_color = "#00549F"

    fig, axes = plt.subplots(2, 2, figsize=(11, 7))
    fig.suptitle(
        f"AC-OPF time-series snapshots  ·  {NET_NAME}"
        f"  (PV × {PV_SCALE}, min voltage deviation, Vmax = {VM_MAX})",
        fontsize=11,
        fontweight="bold",
    )

    def _bar(ax, col, ylabel, title):
        vals = df[col]
        ax.bar(x, vals, color=bar_color, width=0.6, edgecolor="white")
        ax.set_xticks(x)
        ax.set_xticklabels(tick_labels, fontsize=8, rotation=15, ha="right")
        ax.set_ylabel(ylabel, fontsize=9)
        ax.set_title(title, fontsize=10, fontweight="bold")
        ax.grid(axis="y", linestyle="--", alpha=0.5)
        for i, v in enumerate(vals):
            ax.text(
                i,
                v + 0.001 * abs(v + 1e-9),
                f"{v:.3f}",
                ha="center",
                va="bottom",
                fontsize=7.5,
            )

    _bar(axes[0, 0], "eg_p_mw", "Ext-grid P [MW]", "Active import / export")
    _bar(
        axes[0, 1],
        "pv_mw",
        "PV generation [MW]",
        "PV dispatch (after curtailment)",
    )
    _bar(
        axes[1, 0], "max_loading_%", "Max line loading [%]", "Max line loading"
    )

    # voltage band panel
    ax = axes[1, 1]
    ax.fill_between(
        x,
        df["vm_min_pu"],
        df["vm_max_pu"],
        alpha=0.25,
        color=bar_color,
        label="V band",
    )
    ax.plot(x, df["vm_max_pu"], "^--", color=bar_color, ms=6, label="vm_max")
    ax.plot(x, df["vm_min_pu"], "v--", color=bar_color, ms=6, label="vm_min")
    ax.axhline(VM_MAX, color="red", lw=0.9, ls=":", label=f"Vmax = {VM_MAX}")
    ax.axhline(VM_MIN, color="red", lw=0.9, ls=":", label=f"Vmin = {VM_MIN}")
    ax.set_xticks(x)
    ax.set_xticklabels(tick_labels, fontsize=8, rotation=15, ha="right")
    ax.set_ylabel("Voltage [p.u.]", fontsize=9)
    ax.set_title("Bus voltage band", fontsize=10, fontweight="bold")
    ax.set_ylim(0.93, 1.09)
    ax.grid(axis="y", linestyle="--", alpha=0.5)
    ax.legend(fontsize=8, loc="lower right")

    fig.tight_layout()
    os.makedirs(RESULTS_DIR, exist_ok=True)
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    return fig


# ── main ───────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    # ── 1. load base network and full-year profiles ───────────────────────
    print(f"Network  : {NET_NAME}")
    print("Loading simbench profiles ...", end="  ", flush=True)
    base_net = sb.get_simbench_net(NET_NAME)
    profiles = sb.get_absolute_values(
        base_net, profiles_instead_of_study_cases=True
    )
    n_slots = len(profiles[("load", "p_mw")])
    print(f"OK  ({n_slots} slots × 15 min = {n_slots // SLOTS_PER_DAY} days)")
    print(
        f"Objective: min Σ(v−1)²  |  V ∈ [{VM_MIN}, {VM_MAX}] p.u."
        f"  |  PV × {PV_SCALE}  |  Solver = {SOLVER}\n"
    )

    # ── 2. solve each snapshot independently ──────────────────────────────
    rows = []

    for slot_idx, label in SNAPSHOTS:
        print(f"  [{slot_idx:5d}]  {label:<15s}", end="  ", flush=True)

        snap_net = _build_snapshot_net(base_net, profiles, slot_idx)
        pv_rated = snap_net.sgen["p_mw"].sum()
        load_mw = snap_net.load["p_mw"].sum()

        ac = ACOPF(snap_net)
        ac.add_OPF()
        ac.add_voltage_deviation_objective()
        result = ac.solve(solver=SOLVER, print_solver_output=False)

        term = (
            result.solver.termination_condition.value
            if result is not None
            else "no_result"
        )
        status_ok = term == "optimal"
        print(
            f"PV = {pv_rated:.4f} MW  load = {load_mw:.5f} MW  "
            f"→ {'OK' if status_ok else f'FAILED ({term})'}"
        )

        rows.append(_extract(slot_idx, label, ac, term))

    # ── 3. results table ──────────────────────────────────────────────────
    df = pd.DataFrame(rows)

    display_cols = [
        "label",
        "status",
        "eg_p_mw",
        "eg_q_mvar",
        "vm_min_pu",
        "vm_max_pu",
        "max_loading_%",
        "pv_mw",
        "load_mw",
    ]

    pd.set_option("display.float_format", "{:.4f}".format)
    pd.set_option("display.width", 130)
    pd.set_option("display.max_colwidth", 18)

    print("\n" + "─" * 105)
    print(
        "RESULTS — AC-OPF per snapshot  (each snapshot solved independently, no inter-snapshot coupling)"
    )
    print("─" * 105)
    print(df[display_cols].to_string(index=False))
    print("─" * 105)

    n_ok = (df["status"] == "optimal").sum()
    print(f"\n{n_ok}/{len(df)} snapshots converged.")

    # ── 4. plot (converged snapshots only) ────────────────────────────────
    df_ok = df[df["status"] == "optimal"].reset_index(drop=True)
    save_path = os.path.join(RESULTS_DIR, "time_series_snapshot_opf.png")
    _plot(df_ok, save_path)
    plt.close("all")
    print(f"Plot saved → {save_path}")

    # ── 5. learning summary ───────────────────────────────────────────────
    print("""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
WHAT THIS DEMO SHOWS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Each row in the table is a fully independent AC-OPF solve.  No
information flows between snapshots: no ramp constraints, no battery
state of charge, no warm-starting from the previous solution.

The seasonal pattern in a LV rural network is visible across the
eight snapshots:

  • Summer noon  — PV production is high, load is low, ext-grid
    exports surplus to the HV grid (eg_p_mw < 0).  The OPF must
    curtail PV or absorb reactive power to hold vm ≤ Vmax.

  • Winter evening — PV is near zero, residential load peaks,
    the ext-grid imports significant active power (eg_p_mw > 0),
    and vm_min approaches Vmin.  Max line loading is at its
    seasonal high.

  • Spring / autumn noon — intermediate operating point: moderate
    PV, moderate load, the OPF has room to regulate voltages
    with minimal curtailment.

  • Summer midnight — PV is zero, load is low, the network is
    almost idle.  Voltages stay near 1.0 pu with no intervention
    needed.

The "snapshot OPF" paradigm is useful for:
  – hosting-capacity studies (can the grid carry this generation
    level at this load?),
  – planning sensitivity scans across seasons, and
  – quick validation that operating constraints are met.

For coupled multi-period optimisation with storage, EVs, or ramp
constraints use ACOPF_multi_period in
src/potpourri/models_multi_period/.
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
""")
