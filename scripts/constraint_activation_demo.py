"""Constraint activation demo for potpourri.

Demonstrates how activating different groups of AC-OPF constraints
changes the optimal dispatch on the same network and scenario.

Five variants are compared — each activates a different subset of the
constraint groups that make up a full AC-OPF:

  Variant       Voltage  Line  Q-limits
  ─────────────────────────────────────
  Base (P only)    ✗       ✗       ✗
  + Voltage        ✓       ✗       ✗
  + Lines          ✗       ✓       ✗
  + Q limits       ✗       ✗       ✓
  All limits       ✓       ✓       ✓

Objective: minimise active power drawn from the external grid,
           i.e. maximise local PV generation use.  A small quadratic
           Q-regularisation (1e-6) is added to all variants to prevent
           degenerate reactive solutions in the unconstrained cases.

Scenario: 1-LV-rural1--0-sw, lW loadcase, PV × 3.
  Without any limits the optimizer keeps all PV at its rated output
  and exports the surplus to the HV grid.  Each constraint group
  restricts the feasible space and forces additional curtailment or
  reactive-power adjustments.

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
import pyomo.environ as pyo
import simbench as sb

from potpourri.models.ACOPF_base import ACOPF

warnings.filterwarnings("ignore")

NET_NAME = "1-LV-rural1--0-sw"
SOLVER = "ipopt"
RESULTS_DIR = "results"

PV_SCALE = 3.0  # multiply rated PV to create genuine over-generation
PF_MIN = 0.95  # PV inverter minimum power factor
VM_MAX = 1.06  # upper voltage limit [p.u.]
VM_MIN = 0.95  # lower voltage limit [p.u.]

# ── variant catalogue ──────────────────────────────────────────────────────
_VARIANTS = [
    {
        "label": "Base\n(P only)",
        "short": "P-only",
        "v": False,
        "lines": False,
        "q": False,
    },
    {
        "label": "+ Voltage\nlimits",
        "short": "+V",
        "v": True,
        "lines": False,
        "q": False,
    },
    {
        "label": "+ Line\nloading",
        "short": "+L",
        "v": False,
        "lines": True,
        "q": False,
    },
    {
        "label": "+ Q\nlimits",
        "short": "+Q",
        "v": False,
        "lines": False,
        "q": True,
    },
    {
        "label": "All limits\n(full OPF)",
        "short": "All",
        "v": True,
        "lines": True,
        "q": True,
    },
]


# ── scenario ───────────────────────────────────────────────────────────────


def _build_scenario():
    """Return a pandapower network set up for the high-PV demo scenario."""
    net = sb.get_simbench_net(NET_NAME)
    factors = net.loadcases.loc["lW"]
    net.load.p_mw *= factors["pload"]
    net.load.q_mvar *= factors["qload"]
    net.ext_grid.vm_pu = factors["Slack_vm"]

    net.sgen["p_mw"] = net.sgen["p_mw"] * PV_SCALE
    tan_phi = np.sqrt(1 - PF_MIN**2) / PF_MIN
    net.sgen["max_q_mvar"] = tan_phi * net.sgen["p_mw"]
    net.sgen["min_q_mvar"] = -tan_phi * net.sgen["p_mw"]
    net.sgen["max_p_mw"] = net.sgen["p_mw"]
    net.sgen["min_p_mw"] = 0.0
    net.sgen["controllable"] = True

    net.bus["max_vm_pu"] = VM_MAX
    net.bus["min_vm_pu"] = VM_MIN
    return net


# ── model construction ─────────────────────────────────────────────────────


def _build_acopf(net, *, activate_v, activate_lines, activate_q):
    """Build ACOPF, deactivate specified constraint groups, set objective."""
    ac = ACOPF(copy.deepcopy(net))
    ac.add_OPF()

    # selectively deactivate constraint groups
    if not activate_v:
        ac.model.v_pyo.deactivate()
        if len(list(ac.model.Bfix)) > 0:
            ac.model.v_fixed.deactivate()

    if not activate_lines:
        ac.model.line_lim_from.deactivate()
        ac.model.line_lim_to.deactivate()
        ac.model.transf_lim1.deactivate()
        ac.model.transf_lim2.deactivate()

    if not activate_q:
        ac.model.QsG_pyo.deactivate()
        ac.model.QG_pyo.deactivate()

    # objective: min active import + tiny Q regularisation to avoid
    # degenerate reactive solutions when Q constraints are absent
    ac.model.obj_import = pyo.Objective(
        expr=(
            sum(ac.model.pG[g] for g in ac.model.eG)
            + 1e-6 * sum(ac.model.qsG[g] ** 2 for g in ac.model.sG)
        ),
        sense=pyo.minimize,
    )
    return ac


# ── result extraction ──────────────────────────────────────────────────────


def _extract(variant, ac, status):
    """Return a results dict from a solved (or failed) ACOPF instance."""
    label = variant["label"].replace("\n", " ")
    short = variant["short"]
    active = (
        ("V" if variant["v"] else " ")
        + ("L" if variant["lines"] else " ")
        + ("Q" if variant["q"] else " ")
    )

    if status != "optimal":
        return {
            "variant": label,
            "short": short,
            "active": active,
            "status": status,
            **{
                k: np.nan
                for k in [
                    "eg_p_mw",
                    "eg_q_mvar",
                    "vm_min_pu",
                    "vm_max_pu",
                    "max_loading_%",
                    "pv_mw",
                    "pv_curtail_%",
                ]
            },
        }

    net = ac.net
    slack_bus = net.ext_grid.bus.iloc[0]
    vm = net.res_bus.vm_pu.drop(index=slack_bus)
    pv_rated = net.sgen["p_mw"].sum()
    pv_actual = net.res_sgen.p_mw.sum()
    curtail_pct = 100.0 * (1.0 - pv_actual / pv_rated) if pv_rated > 0 else 0.0

    return {
        "variant": label,
        "short": short,
        "active": active,
        "status": status,
        "eg_p_mw": round(net.res_ext_grid.p_mw.sum(), 4),
        "eg_q_mvar": round(net.res_ext_grid.q_mvar.sum(), 4),
        "vm_min_pu": round(vm.min(), 4),
        "vm_max_pu": round(vm.max(), 4),
        "max_loading_%": round(net.res_line.loading_percent.max(), 2),
        "pv_mw": round(pv_actual, 4),
        "pv_curtail_%": round(curtail_pct, 1),
    }


# ── plot ───────────────────────────────────────────────────────────────────


def _plot(df, save_path):
    """Four-panel chart showing the effect of each constraint group."""
    x = np.arange(len(df))
    labels = df["short"].tolist()

    # colour: highlight "All limits" in gold, rest in IAEW blue
    colors = ["#F6A800" if s == "All" else "#00549F" for s in df["short"]]

    fig, axes = plt.subplots(2, 2, figsize=(10, 7))
    fig.suptitle(
        f"Constraint activation  ·  {NET_NAME}"
        f"  (min active import, PV × {PV_SCALE})",
        fontsize=11,
        fontweight="bold",
    )

    def _bar(ax, col, ylabel, title, ref=None, ref_label=None):
        vals = df[col].fillna(0)
        ax.bar(x, vals, color=colors, width=0.55, edgecolor="white")
        if ref is not None:
            ax.axhline(
                ref,
                color="red",
                lw=1.0,
                ls="--",
                label=ref_label or f"limit = {ref}",
            )
            ax.legend(fontsize=8)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=9)
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
                fontsize=8,
            )

    _bar(axes[0, 0], "pv_curtail_%", "PV curtailment [%]", "PV curtailment")
    _bar(
        axes[0, 1],
        "vm_max_pu",
        "Voltage [p.u.]",
        "Max bus voltage",
        ref=VM_MAX,
        ref_label=f"Vmax = {VM_MAX}",
    )
    _bar(
        axes[1, 0],
        "max_loading_%",
        "Loading [%]",
        "Max line loading",
        ref=100.0,
        ref_label="100 % rating",
    )
    _bar(
        axes[1, 1],
        "eg_q_mvar",
        "Ext-grid Q [MVAr]",
        "Reactive import (ext-grid)",
    )

    fig.tight_layout()
    os.makedirs(RESULTS_DIR, exist_ok=True)
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    return fig


# ── main ───────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    net = _build_scenario()

    print(f"Network  : {NET_NAME}")
    print(
        f"Scenario : lW loadcase  |  PV × {PV_SCALE}  |  Vm ∈ [{VM_MIN}, {VM_MAX}]"
    )
    print(
        f"PV total : {net.sgen['p_mw'].sum():.4f} MW"
        f"  (rated {net.sgen['p_mw'].sum() / PV_SCALE:.4f} MW × {PV_SCALE})"
    )
    print(f"Load     : {net.load.p_mw.sum():.5f} MW\n")

    # ── solve five variants ────────────────────────────────────────────────
    rows = []
    for spec in _VARIANTS:
        label_short = spec["short"]
        flags = (
            f"V={'on' if spec['v'] else 'off'}"
            f"  L={'on' if spec['lines'] else 'off'}"
            f"  Q={'on' if spec['q'] else 'off'}"
        )
        print(f"  [{label_short:<8s}]  {flags}  ...", end="  ", flush=True)

        ac = _build_acopf(
            net,
            activate_v=spec["v"],
            activate_lines=spec["lines"],
            activate_q=spec["q"],
        )
        result = ac.solve(solver=SOLVER, print_solver_output=False)
        term = (
            result.solver.termination_condition.value
            if result is not None
            else "no_result"
        )
        print("OK" if term == "optimal" else f"FAILED ({term})")
        rows.append(_extract(spec, ac, term))

    # ── results table ──────────────────────────────────────────────────────
    df = pd.DataFrame(rows)

    pd.set_option("display.float_format", "{:.4f}".format)
    pd.set_option("display.width", 130)
    pd.set_option("display.max_colwidth", 24)

    display_cols = [
        "variant",
        "active",
        "eg_p_mw",
        "eg_q_mvar",
        "vm_min_pu",
        "vm_max_pu",
        "max_loading_%",
        "pv_curtail_%",
    ]

    print("\n" + "─" * 100)
    print(
        "RESULTS  (active flags: V = voltage bounds, "
        "L = line loading, Q = reactive limits)"
    )
    print("─" * 100)
    print(df[display_cols].to_string(index=False))
    print("─" * 100)

    # annotate violations
    print()
    for _, row in df[df["status"] == "optimal"].iterrows():
        flags = []
        if row["vm_max_pu"] > VM_MAX + 1e-4:
            flags.append(
                f"⚠  vm_max = {row['vm_max_pu']:.4f} > Vmax = {VM_MAX}"
            )
        if row["max_loading_%"] > 100.0 + 1e-2:
            flags.append(f"⚠  loading = {row['max_loading_%']:.1f}% > 100%")
        if flags:
            print(f"  {row['short']:8s}  " + "  |  ".join(flags))
    if all(
        (df["vm_max_pu"] <= VM_MAX + 1e-4)
        & (df["max_loading_%"] <= 100.0 + 0.01)
    ):
        print("  (no constraint violations in any variant)")

    # ── plot ───────────────────────────────────────────────────────────────
    df_ok = df[df["status"] == "optimal"].reset_index(drop=True)
    save_path = os.path.join(RESULTS_DIR, "constraint_activation_demo.png")
    _plot(df_ok, save_path)
    plt.close("all")
    print(f"\nPlot saved → {save_path}")

    # ── learning summary ───────────────────────────────────────────────────
    print("""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
WHAT THIS DEMO SHOWS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

The objective is the same in every variant: minimise active power
drawn from the external grid (= maximise local PV use).  Only the
active constraints differ.

Base (P only)
  No voltage, line, or Q constraints.  The optimizer keeps all PV
  at rated output and exports the surplus.  Bus voltages likely
  exceed the physical limit; reactive power is unconstrained and
  may be unrealistically large.  This is the physically infeasible
  "optimistic" solution.

+ Voltage limits
  Bus voltages must stay in [Vmin, Vmax].  The optimizer must
  curtail PV on buses that would otherwise exceed Vmax.  Reactive
  power injection from the PV inverters helps push down voltages,
  reducing the required active-power curtailment.  Line loading
  and Q remain unconstrained.

+ Line loading
  Line apparent-power ratings enforced (≤ 100 %).  In this low-load
  high-PV scenario the feeder current is dominated by the surplus
  export flow; limiting it forces additional curtailment.  Voltages
  remain unconstrained, so vm_max may still violate Vmax.

+ Q limits
  Reactive power of PV inverters and the external grid is bounded
  by the cos(φ) = 0.95 envelope.  Without voltage bounds, the
  optimizer uses reactive power to compensate voltage deviations.
  Restricting Q removes that lever and forces more active curtailment
  to keep operating conditions feasible.

All limits (full OPF)
  Every constraint active simultaneously.  This is the most
  conservative feasible dispatch: least PV, highest curtailment,
  but fully respecting all physical and regulatory limits.

Key insight: each constraint group "protects" one observable:
  voltage limits → vm_max
  line limits    → max loading
  Q limits       → eg_q_mvar
Adding all three produces the true constrained optimum.  Running
the OPF with only a subset gives results that look numerically
better but are physically infeasible.
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
""")
