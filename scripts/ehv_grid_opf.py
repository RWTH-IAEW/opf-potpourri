# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""OPF feasibility test: 1-EHVHV-mixed-all-0-no_sw (large EHV/HV grid).

Tests whether potpourri can compute both DC-OPF and AC-OPF on the largest
standard simbench grid — the mixed EHV/HV network spanning 380 kV and 110 kV
voltage levels with several hundred buses.

Objective (both formulations): minimise active power import from the external
grid(s).  For DC the OPF reduces to a pure LP; for AC it is a non-convex NLP.

Solvers: GLPK for DC-OPF, IPOPT for AC-OPF.

Author: Steffen Kortmann (2026)
"""

from __future__ import annotations

import copy
import time
import warnings

import pyomo.environ as pyo
import simbench as sb

from potpourri.models.ACOPF_base import ACOPF
from potpourri.models.DCOPF import DCOPF

warnings.filterwarnings("ignore")

# ── Configuration ─────────────────────────────────────────────────────────────
NET_NAME = "1-EHVHV-mixed-all-0-no_sw"
PROFILE_IDX = 672
DC_SOLVER = "glpk"
AC_SOLVER = "ipopt"
# ──────────────────────────────────────────────────────────────────────────────


def configure_for_opf(net) -> None:
    """Set OPF bounds on the network in-place."""
    # Voltage limits: ±10 % on all buses
    net.bus["max_vm_pu"] = 1.1
    net.bus["min_vm_pu"] = 0.9

    # Thermal limits: use 100 % of rated capacity
    net.line["max_loading_percent"] = 100.0
    net.trafo["max_loading_percent"] = 100.0
    if "trafo3w" in net and len(net.trafo3w):
        net.trafo3w["max_loading_percent"] = 100.0

    # External grid: wide P/Q bounds (slack)
    net.ext_grid["max_p_mw"] = 1.0e6
    net.ext_grid["min_p_mw"] = -1.0e6
    net.ext_grid["max_q_mvar"] = 1.0e6
    net.ext_grid["min_q_mvar"] = -1.0e6

    # Controllable synchronous generators — allow full redispatch within
    # their rated P range; fix Q bounds if not already set.
    if len(net.gen):
        net.gen["controllable"] = True
        if "max_p_mw" not in net.gen.columns:
            net.gen["max_p_mw"] = net.gen["p_mw"].clip(lower=0.0)
        if "min_p_mw" not in net.gen.columns:
            net.gen["min_p_mw"] = 0.0
        if "max_q_mvar" not in net.gen.columns:
            net.gen["max_q_mvar"] = net.gen["p_mw"].abs() * 0.5
        if "min_q_mvar" not in net.gen.columns:
            net.gen["min_q_mvar"] = -net.gen["p_mw"].abs() * 0.5

    # Static generators (wind / PV) — mark as controllable with curtailment
    if len(net.sgen):
        net.sgen["controllable"] = True
        if "max_p_mw" not in net.sgen.columns:
            net.sgen["max_p_mw"] = net.sgen["p_mw"]
        net.sgen["max_p_mw"] = net.sgen["p_mw"]
        net.sgen["min_p_mw"] = 0.0
        sn = net.sgen["sn_mva"].fillna(net.sgen["p_mw"].abs() * 1.1)
        net.sgen["sn_mva"] = sn
        net.sgen["max_q_mvar"] = 0.4 * sn
        net.sgen["min_q_mvar"] = -0.4 * sn


def print_network_stats(net) -> None:
    print(f"  Buses       : {len(net.bus)}")
    print(f"  Lines       : {len(net.line)}")
    print(f"  Transformers: {len(net.trafo)}")
    print(f"  Ext grids   : {len(net.ext_grid)}")
    print(f"  Generators  : {len(net.gen)}")
    print(f"  Static gens : {len(net.sgen)}")
    print(f"  Loads       : {len(net.load)}")
    kv_levels = sorted(net.bus.vn_kv.unique())
    print(f"  Voltage levels (kV): {kv_levels}")


def run_dcopf(net_template) -> dict:
    """Solve DC-OPF and return a result summary dict."""
    net = copy.deepcopy(net_template)

    t0 = time.perf_counter()
    dcopf = DCOPF(net)
    dcopf.add_OPF()

    # Minimise active power drawn from external grid(s)
    dcopf.model.obj = pyo.Objective(
        expr=sum(dcopf.model.pG[g] for g in dcopf.model.G),
        sense=pyo.minimize,
    )

    res = dcopf.solve(solver=DC_SOLVER, print_solver_output=False)
    elapsed = time.perf_counter() - t0

    ok = res is not None and pyo.check_optimal_termination(res)
    base = dcopf.model.baseMVA

    ext_p_mw = (
        sum(pyo.value(dcopf.model.pG[g]) for g in dcopf.model.G) * base
        if ok
        else float("nan")
    )

    # Five most-loaded lines
    loadings = {}
    if ok:
        for line in dcopf.model.L:
            p_from = abs(pyo.value(dcopf.model.pLfrom[line]))
            s_max = pyo.value(dcopf.model.SLmax[line]) or 1e-9
            loadings[line] = p_from / s_max * 100

    return {
        "ok": ok,
        "time_s": elapsed,
        "termination": str(res.solver.termination_condition)
        if res
        else "error",
        "ext_p_mw": ext_p_mw,
        "top5_loading": sorted(loadings.items(), key=lambda x: -x[1])[:5],
    }


def run_acopf(net_template) -> dict:
    """Solve AC-OPF and return a result summary dict."""
    net = copy.deepcopy(net_template)

    t0 = time.perf_counter()
    ac = ACOPF(net)
    ac.add_OPF(thermal_limit="current", free_slack_vm=True)
    ac.add_voltage_deviation_objective()

    res = ac.solve(solver=AC_SOLVER, print_solver_output=False)
    elapsed = time.perf_counter() - t0

    ok = res is not None and pyo.check_optimal_termination(res)

    vm_min = ac.net.res_bus.vm_pu.min() if ok else float("nan")
    vm_max = ac.net.res_bus.vm_pu.max() if ok else float("nan")
    ext_p_mw = ac.net.res_ext_grid.p_mw.sum() if ok else float("nan")
    losses_mw = ac.net.res_line.pl_mw.sum() if ok else float("nan")

    return {
        "ok": ok,
        "time_s": elapsed,
        "termination": str(res.solver.termination_condition)
        if res
        else "error",
        "vm_min": vm_min,
        "vm_max": vm_max,
        "ext_p_mw": ext_p_mw,
        "losses_mw": losses_mw,
    }


def main() -> None:
    print(f"Network: {NET_NAME}")
    print(f"Snapshot: profile index {PROFILE_IDX}  (summer midday)\n")

    # ── Load and configure ────────────────────────────────────────────────
    net_template = sb.get_simbench_net(NET_NAME)
    profiles = sb.get_absolute_values(
        net_template, profiles_instead_of_study_cases=True
    )

    if ("sgen", "p_mw") in profiles:
        net_template.sgen["p_mw"] = profiles[("sgen", "p_mw")].iloc[
            PROFILE_IDX
        ]
    if ("load", "p_mw") in profiles:
        net_template.load["p_mw"] = profiles[("load", "p_mw")].iloc[
            PROFILE_IDX
        ]
    if ("load", "q_mvar") in profiles:
        net_template.load["q_mvar"] = profiles[("load", "q_mvar")].iloc[
            PROFILE_IDX
        ]

    configure_for_opf(net_template)
    print("Network statistics:")
    print_network_stats(net_template)

    total_load = net_template.load["p_mw"].sum()
    total_sgen = (
        net_template.sgen["p_mw"].sum() if len(net_template.sgen) else 0.0
    )
    total_gen = (
        net_template.gen["p_mw"].sum() if len(net_template.gen) else 0.0
    )
    print(f"\n  Load total  : {total_load:.1f} MW")
    print(f"  Gen total   : {total_gen:.1f} MW  (synchronous)")
    print(f"  sGen total  : {total_sgen:.1f} MW  (wind / PV)")

    # ── DC OPF ───────────────────────────────────────────────────────────
    print(f"\n{'─' * 60}")
    print(f"DC-OPF  (solver: {DC_SOLVER})")
    print(f"{'─' * 60}")
    dc = run_dcopf(net_template)
    tag = "✓  OPTIMAL" if dc["ok"] else f"✗  {dc['termination']}"
    print(f"  Status      : {tag}")
    print(f"  Solve time  : {dc['time_s']:.2f} s")
    print(f"  Ext-grid P  : {dc['ext_p_mw']:+.1f} MW")
    if dc["top5_loading"]:
        print("  Top-5 line loadings:")
        for line, pct in dc["top5_loading"]:
            print(f"    Line {line:>4}: {pct:.1f} %")

    # ── AC OPF ───────────────────────────────────────────────────────────
    print(f"\n{'─' * 60}")
    print(f"AC-OPF  (solver: {AC_SOLVER})")
    print(f"{'─' * 60}")
    ac = run_acopf(net_template)
    tag = "✓  OPTIMAL" if ac["ok"] else f"✗  {ac['termination']}"
    print(f"  Status      : {tag}")
    print(f"  Solve time  : {ac['time_s']:.1f} s")
    print(f"  Voltage     : {ac['vm_min']:.4f} – {ac['vm_max']:.4f} p.u.")
    print(f"  Ext-grid P  : {ac['ext_p_mw']:+.1f} MW")
    print(f"  Line losses : {ac['losses_mw']:.2f} MW")

    # ── Summary ──────────────────────────────────────────────────────────
    print(f"\n{'=' * 60}")
    print("Summary")
    print(f"{'=' * 60}")
    both_ok = dc["ok"] and ac["ok"]
    if both_ok:
        print("Both DC-OPF and AC-OPF solved to optimality.")
        print(
            f"DC took {dc['time_s']:.1f} s,  AC took {ac['time_s']:.1f} s "
            f"({ac['time_s'] / dc['time_s']:.0f}× slower)."
        )
    else:
        if not dc["ok"]:
            print(f"DC-OPF failed: {dc['termination']}")
        if not ac["ok"]:
            print(f"AC-OPF failed: {ac['termination']}")


if __name__ == "__main__":
    main()
