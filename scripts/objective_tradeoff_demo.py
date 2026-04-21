"""AC-OPF objective function tradeoff demo for potpourri.

Shows how four different objectives steer the same AC-OPF to measurably
different operating points on the same 0.4 kV LV network under the same
high-PV / light-load scenario.

Objectives compared
-------------------
1. Voltage deviation   minimize  Σ_b (v[b] - 1)²
2. Reactive generation minimize  Σ_g qsG[g]²
3. Active import       minimize  Σ_g pG[g]  for g in external grids
4. Network losses      minimize  Σ_l (pLfrom[l] + pLto[l])
                                + Σ_t (pThv[t] + pTlv[t])

Scenario: 1-LV-rural1--0-sw, light-load / high-PV.
  PV is scaled to 5× rated to push voltages above 1.06 pu without OPF,
  forcing genuine curtailment decisions.  Voltage bounds are tightened to
  [0.95, 1.06] (tighter than the default 0.9–1.1).
  PV sgens are controllable: P ∈ [0, p_mw], Q ∈ [−cos(φ=0.95)·S, +…].
  Transformer tap is optimised continuously via add_tap_changer_linear().

The differences in the result table illustrate a key practical insight:
which objective you choose determines how the OPF allocates the "PV
curtailment budget" and reactive power across the network.

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
import pandapower as pp
import pyomo.environ as pyo
import simbench as sb

from potpourri.models.ACOPF_base import ACOPF

warnings.filterwarnings("ignore")

NET_NAME = "1-LV-rural1--0-sw"
SOLVER = "ipopt"
RESULTS_DIR = "results"

# ── scenario ─────────────────────────────────────────────────────────────────

PV_SCALE = 5.0  # multiply PV output to create voltage pressure
VM_MAX = 1.06  # tightened upper voltage bound (p.u.)
VM_MIN = 0.95  # lower voltage bound (p.u.)
PF_MIN = 0.95  # minimum power factor → Q capability of PV inverter


def _build_scenario():
    """Return a pandapower network configured for the high-PV demo scenario."""
    net = sb.get_simbench_net(NET_NAME)

    # light-load operating point
    factors = net.loadcases.loc["lW"]
    net.load.p_mw *= factors["pload"]
    net.load.q_mvar *= factors["qload"]
    net.ext_grid.vm_pu = factors["Slack_vm"]

    # boost PV to create real voltage pressure
    net.sgen["p_mw"] = net.sgen["p_mw"] * PV_SCALE

    # PV inverter Q capability at cos(φ) = PF_MIN
    tan_phi = np.sqrt(1 - PF_MIN**2) / PF_MIN
    net.sgen["max_q_mvar"] = tan_phi * net.sgen["p_mw"]
    net.sgen["min_q_mvar"] = -tan_phi * net.sgen["p_mw"]

    # P bounds: can curtail to zero
    net.sgen["max_p_mw"] = net.sgen["p_mw"]
    net.sgen["min_p_mw"] = 0.0
    net.sgen["controllable"] = True

    # tighten voltage bounds for more interesting optimization
    net.bus["max_vm_pu"] = VM_MAX
    net.bus["min_vm_pu"] = VM_MIN

    return net


# ── model construction ───────────────────────────────────────────────────────


def _build_acopf(net):
    """Construct ACOPF with OPF constraints and linear tap changer."""
    ac = ACOPF(copy.deepcopy(net))
    ac.add_OPF()
    ac.add_tap_changer_linear()
    return ac


def _add_import_objective(model):
    """Minimize total active power drawn from external grids."""
    model.obj_import = pyo.Objective(
        expr=sum(model.pG[g] for g in model.eG),
        sense=pyo.minimize,
    )


def _add_losses_objective(model):
    """Minimize total active power losses in lines and transformers."""
    model.obj_losses = pyo.Objective(
        expr=(
            sum(model.pLfrom[l] + model.pLto[l] for l in model.L)
            + sum(model.pThv[t] + model.pTlv[t] for t in model.TRANSF)
        ),
        sense=pyo.minimize,
    )


# ── result extraction ─────────────────────────────────────────────────────────


def _extract(label, ac, obj_attr):
    """Return a results dict from a solved ACOPF instance."""
    m = ac.model
    net = ac.net
    # clean label for table display (plot uses "short" label separately)
    display_label = label.replace("\n", " ").strip()

    obj_val = pyo.value(getattr(m, obj_attr))

    eg_p = net.res_ext_grid.p_mw.sum()
    eg_q = net.res_ext_grid.q_mvar.sum()

    slack_bus = net.ext_grid.bus.iloc[0]
    vm = net.res_bus.vm_pu.drop(index=slack_bus)

    pv_actual = net.res_sgen.p_mw.sum()
    pv_rated = net.sgen["p_mw"].sum()
    curtailment_pct = (
        100.0 * (1.0 - pv_actual / pv_rated) if pv_rated > 0 else 0.0
    )

    return {
        "objective": display_label,
        "obj_value": round(obj_val, 6),
        "eg_p_mw": round(eg_p, 4),
        "eg_q_mvar": round(eg_q, 4),
        "vm_min_pu": round(vm.min(), 4),
        "vm_max_pu": round(vm.max(), 4),
        "max_loading_%": round(net.res_line.loading_percent.max(), 2),
        "pv_mw": round(pv_actual, 4),
        "pv_curtail_%": round(curtailment_pct, 1),
    }


# ── objective catalogue ───────────────────────────────────────────────────────

_SPECS = [
    {
        "label": "1 · min voltage\n  deviation",
        "short": "V-dev",
        "setup": lambda ac: ac.add_voltage_deviation_objective(),
        "attr": "obj_v_deviation",
    },
    {
        "label": "2 · min reactive\n  generation",
        "short": "Q-gen",
        "setup": lambda ac: ac.add_reactive_power_flow_objective(),
        "attr": "obj_reactive",
    },
    {
        "label": "3 · min active\n  import",
        "short": "P-import",
        "setup": lambda ac: _add_import_objective(ac.model),
        "attr": "obj_import",
    },
    {
        "label": "4 · min network\n  losses",
        "short": "Losses",
        "setup": lambda ac: _add_losses_objective(ac.model),
        "attr": "obj_losses",
    },
]

# ── plot ──────────────────────────────────────────────────────────────────────


def _plot(df, save_path):
    """Four-panel bar chart comparing key operating-point metrics."""
    x = np.arange(len(df))
    colors = ["#00549F", "#E30066", "#57AB27", "#F6A800"]

    fig, axes = plt.subplots(2, 2, figsize=(10, 7))
    fig.suptitle(
        f"AC-OPF objective tradeoffs  ·  {NET_NAME}  (PV × {PV_SCALE}, Vmax = {VM_MAX})",
        fontsize=11,
        fontweight="bold",
    )

    panels = [
        (axes[0, 0], "eg_p_mw", "Ext-grid P [MW]", "Active import"),
        (axes[0, 1], "eg_q_mvar", "Ext-grid Q [MVAr]", "Reactive import"),
        (axes[1, 0], "pv_curtail_%", "PV curtailment [%]", "PV curtailment"),
        (
            axes[1, 1],
            "max_loading_%",
            "Max line loading [%]",
            "Max line loading",
        ),
    ]

    for ax, col, ylabel, title in panels:
        bars = ax.bar(x, df[col], color=colors, width=0.55, edgecolor="white")
        ax.set_xticks(x)
        ax.set_xticklabels(df["short"], fontsize=9)
        ax.set_ylabel(ylabel, fontsize=9)
        ax.set_title(title, fontsize=10, fontweight="bold")
        ax.grid(axis="y", linestyle="--", alpha=0.5)
        for bar in bars:
            h = bar.get_height()
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                h + 0.001 * abs(h + 1e-9),
                f"{h:.3f}",
                ha="center",
                va="bottom",
                fontsize=8,
            )

    # add vm band overlay on curtailment panel (secondary)
    ax2 = axes[1, 0].twinx()
    ax2.plot(x, df["vm_max_pu"], "k^--", ms=6, label="vm_max")
    ax2.plot(x, df["vm_min_pu"], "kv--", ms=6, label="vm_min")
    ax2.axhline(VM_MAX, color="red", lw=0.8, ls=":")
    ax2.axhline(VM_MIN, color="red", lw=0.8, ls=":")
    ax2.set_ylabel("Voltage [p.u.]", fontsize=9)
    ax2.set_ylim(0.93, 1.09)
    ax2.legend(fontsize=8, loc="upper right")

    fig.tight_layout()
    os.makedirs(RESULTS_DIR, exist_ok=True)
    fig.savefig(save_path, dpi=150, bbox_inches="tight")
    return fig


# ── main ─────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    # ── 1. scenario ──────────────────────────────────────────────────────────
    net = _build_scenario()

    print(f"Network  : {NET_NAME}")
    print(
        f"Scenario : lW loadcase  |  PV × {PV_SCALE}  |  Vmax = {VM_MAX} p.u."
    )
    print(
        f"PV total : {net.sgen['p_mw'].sum():.4f} MW  "
        f"(rated {net.sgen['p_mw'].sum() / PV_SCALE:.4f} MW × {PV_SCALE})"
    )
    print(f"Load total: {net.load.p_mw.sum():.4f} MW")

    # ── 2. pandapower baseline (unconstrained) ────────────────────────────────
    net_pp = copy.deepcopy(net)
    pp.runpp(net_pp, voltage_depend_loads=False)
    slack_bus = net_pp.ext_grid.bus.iloc[0]

    print(
        f"\nPandapower baseline (no OPF):  "
        f"vm_max = {net_pp.res_bus.vm_pu.drop(index=slack_bus).max():.4f} pu  "
        f"(Vmax = {VM_MAX})"
    )
    if net_pp.res_bus.vm_pu.drop(index=slack_bus).max() > VM_MAX:
        print("  → voltage violation present; OPF curtailment is necessary")
    print()

    # ── 3. solve four objectives ──────────────────────────────────────────────
    rows = []
    shorts = []

    for spec in _SPECS:
        label = spec["label"]
        print(
            f"Solving: {label.replace(chr(10), ' ')} ...", end="  ", flush=True
        )

        ac = _build_acopf(net)
        spec["setup"](ac)
        result = ac.solve(solver=SOLVER, print_solver_output=False)

        converged = (
            result is not None
            and result.solver.termination_condition.value == "optimal"
        )
        status = (
            "OK"
            if converged
            else f"FAILED ({result.solver.termination_condition.value})"
        )
        print(status)

        if converged:
            row = _extract(label, ac, spec["attr"])
            row["short"] = spec["short"]
            rows.append(row)
            shorts.append(spec["short"])

    # ── 4. results table ──────────────────────────────────────────────────────
    df = pd.DataFrame(rows)
    display_cols = [
        "objective",
        "eg_p_mw",
        "eg_q_mvar",
        "vm_min_pu",
        "vm_max_pu",
        "max_loading_%",
        "pv_mw",
        "pv_curtail_%",
    ]

    pd.set_option("display.float_format", "{:.4f}".format)
    pd.set_option("display.width", 120)
    pd.set_option("display.max_colwidth", 30)

    print("\n" + "─" * 90)
    print("RESULTS COMPARISON")
    print("─" * 90)
    print(df[display_cols].to_string(index=False))
    print("─" * 90)

    # per-objective objective value (different units — shown separately)
    print("\nObjective values (native units per formulation):")
    for _, row in df.iterrows():
        print(f"  {row['short']:10s}  {row['obj_value']:.6f}")

    # ── 5. delta table (max − min across objectives) ──────────────────────────
    numeric = df[
        [
            "eg_p_mw",
            "eg_q_mvar",
            "vm_min_pu",
            "vm_max_pu",
            "max_loading_%",
            "pv_curtail_%",
        ]
    ]
    spread = (numeric.max() - numeric.min()).rename("range (max − min)")
    print(
        "\nSpread across objectives (shows sensitivity to objective choice):"
    )
    print(spread.to_string())

    # ── 6. plot ───────────────────────────────────────────────────────────────
    save_path = os.path.join(RESULTS_DIR, "objective_tradeoff_lv_rural1.png")
    _plot(df, save_path)
    plt.close("all")
    print(f"\nPlot saved → {save_path}")

    # ── 7. learning summary ───────────────────────────────────────────────────
    print("""
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
WHAT THIS DEMO SHOWS
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

All four solves use the same network, the same constraints (voltage
bounds, line ratings, PV P/Q limits), and the same IPOPT solver.
Only the objective function changes.  Yet the operating points differ:

  • min voltage deviation  — curtails PV symmetrically to keep all
    bus voltages close to 1.0 pu.  Voltages cluster near 1.0,
    but this may sacrifice extra reactive compensation.

  • min reactive generation — minimises reactive power from PV
    inverters.  The OPF handles voltage deviations with active-power
    curtailment instead of reactive compensation, so vm_max may drift
    further from 1.0 than in the voltage-deviation case.

  • min active import — maximises use of local PV generation (allows
    the grid to export surplus).  Curtailment is minimal; bus voltages
    are pushed to (but not past) Vmax.  Highest vm_max of all cases.

  • min network losses — reveals a structural degeneracy: minimum
    I²R losses = minimum current = curtail nearly all PV.  With
    ~99 % PV curtailment the feeder is almost dead (max line loading
    < 2 %), so losses approach zero.  In practice this objective
    must be combined with a minimum-generation or minimum-curtailment
    penalty to avoid the trivially null solution.

Key takeaway: in distribution-grid OPF the choice of objective is a
policy decision, not just a numerical preference.  "Minimise losses"
and "minimise curtailment" are structurally opposed in high-PV grids.
Practitioners must choose—or combine—objectives deliberately.
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
""")
