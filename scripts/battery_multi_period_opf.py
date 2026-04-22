"""Multi-period AC OPF with battery storage for potpourri.

Battery storage shifts energy in time: charge when generation is high or
demand is low, discharge when the reverse is true.  This example adds
``Battery_multi_period`` to a 24-hour AC OPF and shows how batteries affect
the voltage profile and line loading over the day.

Workflow:
  1. Solve a 24-hour AC OPF without batteries (baseline).
  2. Attach Battery_multi_period and re-solve.
  3. Compare total voltage deviation Σ(v−1)² between both cases.
  4. Inspect per-battery SOC and power trajectories.

Battery constructor parameters:
  - ``scenario``:          preset penetration level (0–3)
  - ``penetration``:       % of non-slack buses to equip (overrides scenario)
  - ``power_pu``:          symmetric charge/discharge limit [p.u.]
  - ``capacity_pu_h``:     energy capacity [p.u.·h]
  - ``efficiency``:        one-way efficiency (0–1)
  - ``soc_min``:           minimum SOC fraction
  - ``initial_soc_fraction``: SOC at t=0 as fraction of soc_max

Network: 1-LV-urban6--0-sw  (96 time steps = 1 day at 15-min resolution).

Institut für Elektrische Anlagen und Netze, Digitalisierung und
Energiewirtschaft (IAEW)
(c) 2023, Steffen Kortmann
"""

import warnings

import pyomo.environ as pe
import simbench as sb

from potpourri.models_multi_period.ACOPF_multi_period import ACOPF_multi_period
from potpourri.technologies.battery import Battery_multi_period

warnings.filterwarnings("ignore")

SOLVER = "ipopt"
NET_NAME = "1-LV-urban6--0-sw"
FROM_T = 0
TO_T = 96  # 96 × 15 min = 1 day


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

    v_dev_base = sum(
        (pe.value(opf_base.model.v[b, t]) - 1.0) ** 2
        for b in opf_base.model.B
        for t in opf_base.model.T
    )
    print(f"  Baseline  Σ(v−1)² = {v_dev_base:.4f}")

    # ── 2. With battery storage — scenario 1 (≈ 8 % penetration) ─────────
    print("\nSolving with battery storage (scenario=1) ...", flush=True)
    opf_bat = ACOPF_multi_period(net, toT=TO_T, fromT=FROM_T)

    battery = Battery_multi_period(opf_bat.net, T=TO_T - FROM_T, scenario=1)
    battery.get_all(opf_bat.model)

    opf_bat.add_OPF()
    opf_bat.add_voltage_deviation_objective()
    opf_bat.solve(solver=SOLVER, print_solver_output=False)

    v_dev_bat = sum(
        (pe.value(opf_bat.model.v[b, t]) - 1.0) ** 2
        for b in opf_bat.model.B
        for t in opf_bat.model.T
    )
    improvement = (v_dev_base - v_dev_bat) / v_dev_base * 100
    print(f"  Batteries Σ(v−1)² = {v_dev_bat:.4f}")
    print(f"  Voltage-deviation reduction: {improvement:.1f} %")

    # ── 3. SOC profile for each battery ──────────────────────────────────
    print("\nBattery SOC (fraction) — every hour:")
    bats = list(opf_bat.model.BAT)
    header = f"  {'hour':>6}" + "".join(f"  bat{b:>2}" for b in bats)
    print(header)
    for t in list(opf_bat.model.T)[::4]:
        hour = (t - FROM_T) / 4
        row = f"  {hour:>5.0f}h"
        for b in bats:
            row += f"  {pe.value(opf_bat.model.BAT_SOC[b, t]):>5.2f}"
        print(row)

    # ── 4. Custom battery sizing ──────────────────────────────────────────
    base = net.sn_mva
    opf_custom = ACOPF_multi_period(net, toT=TO_T, fromT=FROM_T)
    battery_custom = Battery_multi_period(
        opf_custom.net,
        T=TO_T - FROM_T,
        penetration=20.0,  # 20 % of non-slack buses
        power_pu=0.012 / base,  # 12 kW charge/discharge limit
        capacity_pu_h=0.030 / base,  # 30 kWh capacity
        efficiency=0.95,
        soc_min=0.1,
        initial_soc_fraction=0.4,
    )
    battery_custom.get_all(opf_custom.model)
    print(
        f"\nCustom battery: {len(battery_custom.random_indexes)} units, "
        f"power={battery_custom.bat_power:.4f} p.u., "
        f"cap={battery_custom.bat_cap:.4f} p.u.·h, "
        f"η={battery_custom.bat_efficiency}"
    )

    print(
        "\nKey takeaway: Battery_multi_period attaches SOC dynamics as Pyomo "
        "constraints.  The optimizer shifts energy temporally to smooth "
        "voltages — the effect grows with battery penetration and capacity."
    )
