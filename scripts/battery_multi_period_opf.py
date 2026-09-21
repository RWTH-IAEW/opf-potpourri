# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Multi-period AC OPF with battery storage for potpourri.

Battery storage shifts energy in time: charge when generation is high or
demand is low, discharge when the reverse is true.  The converter also
supplies reactive power, bounded by its apparent-power circle.  This example
adds ``Battery_multi_period`` to a 24-hour AC OPF and shows how batteries
affect the voltage profile over the day.

Workflow:
  1. Solve a 24-hour AC OPF without batteries (baseline).
  2. Attach Battery_multi_period, re-solve, and compare Σ(v−1)².
  3. Inspect per-battery SOC and power trajectories.
  4. Reactive support: widen the converter and add a power-factor floor.

Battery constructor parameters:
  - ``scenario``:          preset penetration level (0–3)
  - ``penetration``:       % of non-slack buses to equip (overrides scenario)
  - ``power_pu``:          symmetric charge/discharge limit [p.u.]
  - ``capacity_pu_h``:     energy capacity [p.u.·h]
  - ``efficiency``:        one-way efficiency (0–1); round trip is η²
  - ``soc_min``:           minimum SOC fraction
  - ``initial_soc_fraction``: SOC at t=0 as fraction of soc_max
  - ``terminal_soc``:      SOC at the last step ("cyclic", a float, or None)
  - ``s_inv_pu``:          converter apparent-power rating [p.u.]; the S²
                           circle bounds P and Q jointly.  Defaults to
                           ``power_pu``, which leaves no reactive headroom at
                           full active power
  - ``cos_phi_min``:       power-factor floor on the reactive dispatch
  - ``q_control``:         grid-code capability area ("qp"/"qu"/"both")

Network: 1-LV-rural1--0-sw  (96 time steps = 1 day at 15-min resolution).

.. note::

   Battery placement is seeded (``DEFAULT_PLACEMENT_SEED``), so re-running this
   script equips the same buses and reproduces the same numbers.  Pass
   ``seed=...`` to a device constructor to sample a different placement.

Institut für Elektrische Anlagen und Netze, Digitalisierung und
Energiewirtschaft (IAEW)
Author: Steffen Kortmann (2023)
"""

import logging
import warnings

import pyomo.environ as pyo
import simbench as sb

from potpourri.models_multi_period.ACOPF_multi_period import ACOPF_multi_period
from potpourri.technologies.battery import Battery_multi_period

warnings.filterwarnings("ignore")
# Step 4 deliberately solves two cases that have no feasible dispatch, and
# Pyomo logs a multi-line warning whenever it loads such a result. The status
# is reported explicitly below, so silence the raw log to keep the output
# readable.
logging.getLogger("pyomo.core").setLevel(logging.ERROR)
logging.getLogger("pyomo.solvers").setLevel(logging.ERROR)

# ── Configuration ─────────────────────────────────────────────────────────────
SOLVER = "ipopt"
NET_NAME = "1-LV-rural1--0-sw"
FROM_T = 0
TO_T = 96  # 96 × 15 min = 1 day
# ──────────────────────────────────────────────────────────────────────────────


if __name__ == "__main__":
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

    # ── 1. Baseline — no battery storage ─────────────────────────────────
    print("Solving baseline (no batteries) ...", flush=True)
    opf_base = ACOPF_multi_period(net, toT=TO_T, fromT=FROM_T)
    opf_base.add_OPF()
    opf_base.add_voltage_deviation_objective()
    opf_base.solve(solver=SOLVER, print_solver_output=False)

    # Read the model's own objective rather than re-deriving a sum of squared
    # deviations. add_voltage_deviation_objective() excludes the slack from the
    # (v − 1)² term and instead penalises its distance from the base-case
    # magnitude v_b0, so a hand-rolled Σ(v−1)² over *all* buses measures a
    # different quantity — and can move in the opposite direction to the value
    # actually being minimised.
    v_dev_base = pyo.value(opf_base.model.obj_v_deviation)
    print(f"  Baseline  objective = {v_dev_base:.6f}")

    # ── 2. With battery storage — scenario 1 (≈ 8 % penetration) ─────────
    print("\nSolving with battery storage (scenario=1) ...", flush=True)
    opf_bat = ACOPF_multi_period(net, toT=TO_T, fromT=FROM_T)

    battery = Battery_multi_period(opf_bat.net, T=TO_T - FROM_T, scenario=1)
    battery.get_all(opf_bat.model)

    opf_bat.add_OPF()
    opf_bat.add_voltage_deviation_objective()
    opf_bat.solve(solver=SOLVER, print_solver_output=False)

    v_dev_bat = pyo.value(opf_bat.model.obj_v_deviation)
    change = (v_dev_base - v_dev_bat) / v_dev_base * 100
    print(f"  Batteries objective = {v_dev_bat:.6f}")
    print(f"  Objective reduction: {change:+.2f} %")

    # ── 3. SOC profile for each battery ──────────────────────────────────
    print("\nBattery SOC (fraction) — every hour:")
    bats = list(opf_bat.model.BAT)
    header = f"  {'hour':>6}" + "".join(f"  bat{b:>2}" for b in bats)
    print(header)
    for t in list(opf_bat.model.T)[::4]:
        hour = (t - FROM_T) / 4
        row = f"  {hour:>5.0f}h"
        for b in bats:
            row += f"  {pyo.value(opf_bat.model.BAT_SOC[b, t]):>5.2f}"
        print(row)

    # ── 4. Reactive support from the converter ────────────────────────────
    # Midday PV pushes this feeder against its 1.05 p.u. ceiling, and holding
    # it needs reactive absorption. The three cases below differ only in how
    # much reactive capability the converter is allowed:
    #
    #   * s_inv_pu defaults to power_pu, so the S² circle leaves no room for Q
    #     once the battery is charging or discharging at its power limit.
    #   * Oversizing the converter buys that headroom back.
    #   * A cos φ floor ties Q to the active power, so it takes the headroom
    #     away again whenever P is small — which is exactly when voltage
    #     support is wanted.
    #
    # Some of these cases may have no feasible dispatch at all: with little
    # reactive capability the band cannot always be held, and whether a given
    # case converges also depends on which buses the batteries landed on. The
    # placement is seeded, so a given run is repeatable, but the status is
    # reported rather than an objective value read off an infeasible iterate.
    base = net.sn_mva
    print("\nReactive support (converter sizing vs power-factor floor):")
    for label, kwargs in [
        ("S = P (no Q headroom)", {}),
        ("S = 1.5 P", {"s_inv_pu": 0.018 / base}),
        (
            "S = 1.5 P, cos φ ≥ 0.95",
            {"s_inv_pu": 0.018 / base, "cos_phi_min": 0.95},
        ),
    ]:
        opf_q = ACOPF_multi_period(net, toT=TO_T, fromT=FROM_T)
        bat_q = Battery_multi_period(
            opf_q.net,
            T=TO_T - FROM_T,
            penetration=20.0,  # 20 % of non-slack buses
            power_pu=0.012 / base,  # 12 kW charge/discharge limit
            capacity_pu_h=0.030 / base,  # 30 kWh capacity
            efficiency=0.95,
            soc_min=0.1,
            initial_soc_fraction=0.4,
            **kwargs,
        )
        bat_q.get_all(opf_q.model)
        opf_q.add_OPF()
        opf_q.add_voltage_deviation_objective()
        results = opf_q.solve(solver=SOLVER, print_solver_output=False)

        if not pyo.check_optimal_termination(results):
            print(
                f"  {label:<24s} no feasible dispatch — the voltage band "
                f"cannot be held with this reactive capability"
            )
            continue

        dev = pyo.value(opf_q.model.obj_v_deviation)
        q_peak = max(
            abs(pyo.value(opf_q.model.BAT_Q[b, t]))
            for b in opf_q.model.BAT
            for t in opf_q.model.T
        )
        print(
            f"  {label:<24s} objective = {dev:.6f}   "
            f"peak |Q| = {q_peak * base * 1000:.2f} kvar   "
            f"(S_inv = {bat_q.bat_s_inv * base * 1000:.0f} kVA)"
        )

    print(
        "\nKey takeaway: the battery shifts energy in time and supports "
        "voltage with reactive power, both bounded by the converter's S² "
        "circle.  Sizing the converter above the active-power limit is what "
        "buys reactive headroom at full charge or discharge — and a cos φ "
        "floor gives it back up, because it ties Q to P just when P is low."
        "\n\nAll figures above are the model's own objective.  Comparing a "
        "hand-rolled\nΣ(v−1)² instead would compare a different quantity: the "
        "objective leaves the\nslack out of that term and penalises its "
        "distance from the base-case magnitude."
    )
