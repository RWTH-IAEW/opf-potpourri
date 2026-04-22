"""Multi-period AC OPF example for potpourri.

Demonstrates how to solve a 24-hour AC OPF over a SimBench low-voltage
network.  Each of the 96 time steps (15-minute intervals) is coupled through
the shared voltage and line-loading constraints; the solver finds the optimal
reactive-power dispatch for the entire day simultaneously.

Workflow:
  1. Load a SimBench LV network — profiles are attached automatically.
  2. Set operational limits (voltage bounds, line loading, sgen curtailment).
  3. Build ACOPF_multi_period(net, toT=96).
  4. Call add_OPF() and add_voltage_deviation_objective().
  5. Solve with IPOPT and inspect time-series results.

Key difference from snapshot OPF: all time steps are optimised jointly.
There is no temporal coupling of energy (no storage here), but the shared
network constraints make the problem a single large NLP.

Network: 1-LV-urban6--0-sw (urban low-voltage, diverse load/PV mix).

Institut für Elektrische Anlagen und Netze, Digitalisierung und
Energiewirtschaft (IAEW)
(c) 2023, Steffen Kortmann
"""

import warnings

import pyomo.environ as pe
import simbench as sb

from potpourri.models_multi_period.ACOPF_multi_period import ACOPF_multi_period

warnings.filterwarnings("ignore")

SOLVER = "ipopt"
NET_NAME = "1-LV-urban6--0-sw"
FROM_T = 0
TO_T = 96  # 96 × 15 min = 24 hours


if __name__ == "__main__":
    # ── 1. Load network ───────────────────────────────────────────────────
    net = sb.get_simbench_net(NET_NAME)

    # Operational limits
    net.bus["max_vm_pu"] = 1.05
    net.bus["min_vm_pu"] = 0.95
    net.line["max_loading_percent"] = 80.0
    net.sgen["controllable"] = True
    net.sgen["max_p_mw"] = net.sgen["p_mw"]
    net.sgen["min_p_mw"] = 0.0
    net.ext_grid["max_q_mvar"] = 500.0
    net.ext_grid["min_q_mvar"] = -500.0

    print(f"Network  : {NET_NAME}")
    print(f"Horizon  : t={FROM_T} … {TO_T - 1}  ({TO_T} × 15 min = 24 h)")
    print(f"Solver   : {SOLVER}\n")

    # ── 2. Build and solve the multi-period model ─────────────────────────
    opf = ACOPF_multi_period(net, toT=TO_T, fromT=FROM_T)
    opf.add_OPF()
    opf.add_voltage_deviation_objective()

    print("Solving ...", flush=True)
    result = opf.solve(solver=SOLVER, print_solver_output=False)

    term = (
        result.solver.termination_condition.value
        if result is not None
        else "no_result"
    )
    print(f"Termination: {term}\n")

    if term != "optimal":
        print("Solver did not converge – aborting.")
        raise SystemExit(1)

    # ── 3. Inspect results ────────────────────────────────────────────────
    base = opf.model.baseMVA

    # Hourly ext-grid import (sample every 4 time steps = 1 hour)
    print("Ext-grid active power — hourly summary (MW):")
    print(f"  {'Hour':>6}  {'P ext-grid (MW)':>18}")
    for t in list(opf.model.T)[::4]:
        hour = (t - FROM_T) / 4
        p_mw = sum(pe.value(opf.model.pG[g, t]) * base for g in opf.model.eG)
        print(f"  {hour:>5.0f}h  {p_mw:>18.4f}")

    # Voltage band over the day
    v_max_by_t = [
        max(pe.value(opf.model.v[b, t]) for b in opf.model.B)
        for t in opf.model.T
    ]
    v_min_by_t = [
        min(pe.value(opf.model.v[b, t]) for b in opf.model.B)
        for t in opf.model.T
    ]
    print(
        f"\nVoltage band over full horizon:\n"
        f"  vm_max = {max(v_max_by_t):.4f} p.u.  "
        f"  vm_min = {min(v_min_by_t):.4f} p.u."
    )

    # Objective value
    obj_val = pe.value(opf.model.obj_v_deviation)
    print(f"\nObjective Σ(v−1)² over all buses and steps: {obj_val:.6f}")

    print(
        "\nKey takeaway: ACOPF_multi_period solves all 96 time steps as one "
        "large NLP.  The reactive power dispatch at each step is optimised "
        "jointly — unlike snapshot OPF where each step is independent."
    )
