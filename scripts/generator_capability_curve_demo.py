"""Generator reactive capability curve demo for potpourri.

Shows how reactive power limits affect AC-OPF dispatch and voltages on a
small radial LV feeder with three generator types:

  PV inverter  — rectangular Q bounds derived from a power-factor rule
                 (cos φ ≥ pf_pv, typical grid-code requirement)
  Wind inverter — tighter rectangular bounds  (cos φ ≥ pf_wind)
  Battery       — circular capability: P² + Q² ≤ S²  (built-in
                  ``stor_inverter_cap`` constraint, always active)

Two scenarios are compared:

  Scenario A – wide Q range  (pf ≥ 0.70 for PV and wind)
  Scenario B – grid-code     (PV pf ≥ 0.90, wind pf ≥ 0.95)

Both scenarios use the same voltage-deviation objective  Σ(v−1)²  so the
solver exploits all available reactive freedom to keep voltages near 1 p.u.
The tighter limits in Scenario B leave residual voltage deviations and
shift reactive burden toward the external grid.

Network topology (0.4 kV, 1 MVA base):

    ext_grid ── bus0 ── bus1 ── bus2 ─── bus3 ─── bus4
                        load   load       load      load
                               PV 0.10   wind 0.08  bat 50 kVA
                               MW        MW

Results reported:
  • Q dispatch per generator
  • bus voltage profile
  • ext-grid reactive import
  • battery (P, Q) operating point vs. its S² capability circle

Requires: IPOPT (via conda or environment.yaml).

Institut für Elektrische Anlagen und Netze, Digitalisierung und
Energiewirtschaft (IAEW)
(c) 2023, Steffen Kortmann
"""

import copy
import math
import warnings

import pandas as pd
import pandapower as pp
import pyomo.environ as pyo

from potpourri.models.ACOPF_base import ACOPF

warnings.filterwarnings("ignore")

SOLVER = "ipopt"

# ── power-factor scenario parameters ─────────────────────────────────────────
PF_WIDE = 0.70  # Scenario A: loose limit  (cos φ ≥ 0.70)
PF_PV = 0.90  # Scenario B: VDE-AR-N 4105 / EN 50549 for PV
PF_WIND = 0.95  # Scenario B: stricter for wind generators


# ── helpers ───────────────────────────────────────────────────────────────────


def q_lim_from_pf(p_mw: float, pf: float) -> float:
    """Return the Q limit [Mvar] for rated P [MW] at minimum power factor pf."""
    return abs(p_mw) * math.tan(math.acos(pf))


# ── network builder ───────────────────────────────────────────────────────────


def build_feeder() -> pp.pandapowerNet:
    """Return a 5-bus radial 0.4 kV feeder with PV, wind, and a battery.

    Line impedances are chosen so that the ±5 % voltage band is meaningfully
    stressed by the 10–15 % net renewable injection surplus, making reactive
    dispatch clearly visible in the voltage profile.
    """
    net = pp.create_empty_network(sn_mva=1.0)

    # buses (all 0.4 kV LV)
    kw = dict(max_vm_pu=1.05, min_vm_pu=0.95)
    b0 = pp.create_bus(net, vn_kv=0.4, name="slack", **kw)
    b1 = pp.create_bus(net, vn_kv=0.4, name="bus1", **kw)
    b2 = pp.create_bus(net, vn_kv=0.4, name="bus2-PV", **kw)
    b3 = pp.create_bus(net, vn_kv=0.4, name="bus3-wind", **kw)
    b4 = pp.create_bus(net, vn_kv=0.4, name="bus4-bat", **kw)

    pp.create_ext_grid(
        net, bus=b0, vm_pu=1.02, max_q_mvar=1.0, min_q_mvar=-1.0
    )

    # radial lines: r/x typical for LV NAYY cable, short spans so Z_pu ≈ 0.3–0.5
    lp = dict(
        r_ohm_per_km=0.5, x_ohm_per_km=0.15, c_nf_per_km=200, max_i_ka=0.2
    )
    pp.create_line_from_parameters(
        net, b0, b1, length_km=0.12, **lp, name="l01"
    )
    pp.create_line_from_parameters(
        net, b1, b2, length_km=0.15, **lp, name="l12"
    )
    pp.create_line_from_parameters(
        net, b2, b3, length_km=0.15, **lp, name="l23"
    )
    pp.create_line_from_parameters(
        net, b3, b4, length_km=0.12, **lp, name="l34"
    )

    # loads (pf ≈ 0.90 inductive) — total 0.17 MW
    pp.create_load(net, b1, p_mw=0.040, q_mvar=0.019, name="load1")
    pp.create_load(net, b2, p_mw=0.045, q_mvar=0.022, name="load2")
    pp.create_load(net, b3, p_mw=0.045, q_mvar=0.022, name="load3")
    pp.create_load(net, b4, p_mw=0.040, q_mvar=0.019, name="load4")

    # PV at bus2 — active power fixed at 0.10 MW; Q free within limits
    pp.create_sgen(
        net,
        b2,
        p_mw=0.10,
        q_mvar=0.0,
        name="PV",
        type="PV",
        controllable=True,
        max_p_mw=0.10,
        min_p_mw=0.10,
    )

    # wind at bus3 — active power fixed at 0.08 MW
    pp.create_sgen(
        net,
        b3,
        p_mw=0.08,
        q_mvar=0.0,
        name="Wind",
        type="WP",
        controllable=True,
        max_p_mw=0.08,
        min_p_mw=0.08,
    )

    # battery at bus4 — inverter rated at 50 kVA; stor_inverter_cap enforces
    # the P²+Q² ≤ S² circle automatically when net.storage is non-empty
    pp.create_storage(
        net,
        b4,
        p_mw=0.0,
        max_e_mwh=0.10,
        sn_mva=0.05,
        scaling=1.0,
        efficiency_percent=95.0,
        soc_percent=50.0,
        name="Battery",
    )

    net.line["max_loading_percent"] = 100.0
    return net


def set_q_limits(net: pp.pandapowerNet, pf_pv: float, pf_wind: float) -> None:
    """Set max_q_mvar / min_q_mvar on every sgen from a power-factor rule."""
    for idx, row in net.sgen.iterrows():
        pf = pf_pv if row.get("type", "") == "PV" else pf_wind
        ql = q_lim_from_pf(row.p_mw, pf)
        net.sgen.at[idx, "max_q_mvar"] = ql
        net.sgen.at[idx, "min_q_mvar"] = -ql


# ── ACOPF solver ──────────────────────────────────────────────────────────────


def solve_scenario(net: pp.pandapowerNet, label: str) -> dict:
    """Build and solve an ACOPF with voltage-deviation objective.

    Calls add_OPF() (which enforces the Q limits and voltage bounds) followed
    by add_voltage_deviation_objective().  The battery's circular capability
    constraint (stor_inverter_cap) is added automatically by Basemodel when
    net.storage is non-empty.
    """
    opf = ACOPF(net)
    opf.add_OPF()
    opf.add_voltage_deviation_objective()

    result = opf.solve(solver=SOLVER, print_solver_output=False)
    if result is None or not pyo.check_optimal_termination(result):
        raise RuntimeError(
            f"[{label}] solver did not find an optimal solution"
        )

    m = opf.model
    base = pyo.value(m.baseMVA)

    # sgen indices from net (sgen[0]=PV, sgen[1]=Wind, storage[0]=Battery)
    pv_idx = net.sgen.index[net.sgen.name == "PV"].item()
    wind_idx = net.sgen.index[net.sgen.name == "Wind"].item()
    bat_idx = net.storage.index[net.storage.name == "Battery"].item()

    q_pv = pyo.value(m.qsG[pv_idx]) * base
    q_wind = pyo.value(m.qsG[wind_idx]) * base
    p_bat = pyo.value(m.pSTOR[bat_idx]) * base  # Pchg − Pdis
    q_bat = pyo.value(m.qSTOR[bat_idx]) * base
    q_ext = sum(pyo.value(m.qG[g]) for g in m.eG) * base
    p_ext = sum(pyo.value(m.pG[g]) for g in m.eG) * base

    vm = {b: pyo.value(m.v[b]) for b in m.B}

    return {
        "label": label,
        "q_pv": q_pv,
        "q_wind": q_wind,
        "p_bat": p_bat,
        "q_bat": q_bat,
        "p_ext": p_ext,
        "q_ext": q_ext,
        "vm": vm,
        "vm_min": min(vm.values()),
        "vm_max": max(vm.values()),
        "obj": pyo.value(m.obj_v_deviation),
    }


# ── main ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    base_net = build_feeder()

    # ── Scenario A: wide Q limits (pf ≥ 0.70) ────────────────────────────
    net_a = copy.deepcopy(base_net)
    set_q_limits(net_a, pf_pv=PF_WIDE, pf_wind=PF_WIDE)

    q_pv_lim_a = q_lim_from_pf(net_a.sgen.at[0, "p_mw"], PF_WIDE)
    q_wind_lim_a = q_lim_from_pf(net_a.sgen.at[1, "p_mw"], PF_WIDE)

    print("Solving Scenario A (wide limits, pf ≥ 0.70) …", flush=True)
    res_a = solve_scenario(net_a, "A – wide (pf ≥ 0.70)")

    # ── Scenario B: grid-code limits ──────────────────────────────────────
    net_b = copy.deepcopy(base_net)
    set_q_limits(net_b, pf_pv=PF_PV, pf_wind=PF_WIND)

    q_pv_lim_b = q_lim_from_pf(net_b.sgen.at[0, "p_mw"], PF_PV)
    q_wind_lim_b = q_lim_from_pf(net_b.sgen.at[1, "p_mw"], PF_WIND)

    print(
        "Solving Scenario B (grid code, PV pf ≥ 0.90 / wind pf ≥ 0.95) …",
        flush=True,
    )
    res_b = solve_scenario(net_b, "B – grid code")

    # ── Q-limit summary ───────────────────────────────────────────────────
    bat_s_rated = base_net.storage.sn_mva.iloc[0]  # MVA

    print("\n── Q limits [Mvar] ─────────────────────────────────────────────")
    print(f"  {'':30s}  {'PV':>8}  {'Wind':>8}  {'Battery':>12}")
    print(
        f"  {'Scenario A (pf ≥ 0.70)':30s}  "
        f"±{q_pv_lim_a:.4f}  ±{q_wind_lim_a:.4f}  "
        f"circle S={bat_s_rated:.3f} MVA"
    )
    print(
        f"  {'Scenario B (grid code)':30s}  "
        f"±{q_pv_lim_b:.4f}  ±{q_wind_lim_b:.4f}  "
        f"circle S={bat_s_rated:.3f} MVA"
    )

    # ── Dispatch and voltage comparison ───────────────────────────────────
    rows = []
    for res in (res_a, res_b):
        rows.append(
            {
                "Scenario": res["label"],
                "Q PV [Mvar]": f"{res['q_pv']:+.4f}",
                "Q Wind [Mvar]": f"{res['q_wind']:+.4f}",
                "P Battery [MW]": f"{res['p_bat']:+.4f}",
                "Q Battery [Mvar]": f"{res['q_bat']:+.4f}",
                "P ext-grid [MW]": f"{res['p_ext']:+.4f}",
                "Q ext-grid [Mvar]": f"{res['q_ext']:+.4f}",
                "vm_min [p.u.]": f"{res['vm_min']:.5f}",
                "vm_max [p.u.]": f"{res['vm_max']:.5f}",
                "Obj Σ(v−1)²": f"{res['obj']:.6f}",
            }
        )

    df = pd.DataFrame(rows).set_index("Scenario")
    print("\n── Dispatch and voltage summary ────────────────────────────────")
    print(df.T.to_string())

    # ── Per-bus voltage profile ───────────────────────────────────────────
    print("\n── Bus voltages [p.u.] ─────────────────────────────────────────")
    bus_names = {b: base_net.bus.at[b, "name"] for b in sorted(res_a["vm"])}
    vm_df = pd.DataFrame(
        {
            "bus name": [bus_names[b] for b in sorted(res_a["vm"])],
            "Scen A": [res_a["vm"][b] for b in sorted(res_a["vm"])],
            "Scen B": [res_b["vm"][b] for b in sorted(res_b["vm"])],
            "ΔV (A−B)": [
                res_a["vm"][b] - res_b["vm"][b] for b in sorted(res_a["vm"])
            ],
        },
        index=sorted(res_a["vm"]),
    )
    vm_df.index.name = "bus"
    print(vm_df.to_string(float_format="{:.5f}".format))

    # ── Battery operating point ───────────────────────────────────────────
    print("\n── Battery operating point vs. capability circle ───────────────")
    print(f"  Battery inverter rating: S = {bat_s_rated * 1000:.0f} kVA")
    for res in (res_a, res_b):
        s_util = math.hypot(res["p_bat"], res["q_bat"]) / bat_s_rated
        print(
            f"  {res['label']:<35}  "
            f"P = {res['p_bat'] * 1e3:+6.1f} kW  "
            f"Q = {res['q_bat'] * 1e3:+6.1f} kvar  "
            f"|S|/S_rated = {s_util:.3f}"
        )

    print(
        "\nKey takeaway: Tighter Q limits (Scenario B) reduce the optimizer's"
        "\nability to absorb reactive power at PV/wind buses, leaving voltages"
        "\nfurther from 1.0 p.u. and pushing more reactive burden onto the"
        "\next-grid.  The battery's circular P²+Q²≤S² constraint is always"
        "\nactive; its Q headroom depends on its real-power operating point."
    )

    # ── Optional plot ─────────────────────────────────────────────────────
    try:
        import numpy as np
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
        fig.suptitle("Generator Reactive Capability Curve Demo", fontsize=13)

        bus_ids = sorted(res_a["vm"])
        bus_labels = [bus_names[b] for b in bus_ids]

        # 1 — voltage profile
        ax = axes[0]
        ax.plot(
            bus_labels,
            [res_a["vm"][b] for b in bus_ids],
            "o-",
            label="A – wide (pf≥0.70)",
        )
        ax.plot(
            bus_labels,
            [res_b["vm"][b] for b in bus_ids],
            "s--",
            label="B – grid code",
        )
        ax.axhline(1.05, color="r", lw=0.8, ls=":", label="V limits")
        ax.axhline(0.95, color="r", lw=0.8, ls=":")
        ax.set_ylim(0.93, 1.07)
        ax.set_ylabel("Voltage [p.u.]")
        ax.set_title("Voltage profile")
        ax.legend(fontsize=8)
        ax.tick_params(axis="x", labelrotation=20)

        # 2 — Q dispatch bar chart
        ax = axes[1]
        labels = ["PV", "Wind", "Battery", "Ext-grid"]
        vals_a = [
            res_a["q_pv"],
            res_a["q_wind"],
            res_a["q_bat"],
            res_a["q_ext"],
        ]
        vals_b = [
            res_b["q_pv"],
            res_b["q_wind"],
            res_b["q_bat"],
            res_b["q_ext"],
        ]
        x = np.arange(len(labels))
        w = 0.35
        ax.bar(x - w / 2, [v * 1e3 for v in vals_a], w, label="A – wide")
        ax.bar(x + w / 2, [v * 1e3 for v in vals_b], w, label="B – grid code")
        ax.axhline(0, color="k", lw=0.6)
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.set_ylabel("Q [kvar]")
        ax.set_title("Reactive dispatch")
        ax.legend(fontsize=8)

        # 3 — battery capability circle
        ax = axes[2]
        theta = np.linspace(0, 2 * np.pi, 400)
        r_kva = bat_s_rated * 1e3
        ax.plot(
            r_kva * np.cos(theta),
            r_kva * np.sin(theta),
            "k--",
            lw=1.2,
            label=f"S = {r_kva:.0f} kVA",
        )
        ax.plot(
            res_a["p_bat"] * 1e3,
            res_a["q_bat"] * 1e3,
            "o",
            ms=10,
            label="A – wide",
        )
        ax.plot(
            res_b["p_bat"] * 1e3,
            res_b["q_bat"] * 1e3,
            "s",
            ms=10,
            label="B – grid code",
        )
        ax.axhline(0, color="k", lw=0.5)
        ax.axvline(0, color="k", lw=0.5)
        ax.set_xlabel("P [kW]")
        ax.set_ylabel("Q [kvar]")
        ax.set_title("Battery capability circle")
        ax.legend(fontsize=8)
        ax.set_aspect("equal")

        plt.tight_layout()
        out = "generator_capability_curve_demo.png"
        plt.savefig(out, dpi=120)
        print(f"\nPlot saved to {out}")
        plt.show()

    except ImportError:
        print("\n(matplotlib not available — skipping plot)")
