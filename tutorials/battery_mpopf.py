# %% md
# # Multi-Period OPF with Battery Storage
#
# Battery storage shifts energy in time: charge when generation is high or
# demand is low, discharge when the reverse is true.  This tutorial adds
# `Battery_multi_period` to a 24-hour AC OPF and shows how the batteries
# affect voltage profiles and line loading over the day.
#
# We will:
# 1. Solve a 24-hour AC OPF without batteries as a baseline
# 2. Add battery storage and re-solve
# 3. Compare voltage deviation between both cases
# 4. Demonstrate the new flexible constructor parameters

# %% md
# ## 1. Imports

# %%
import pyomo.environ as pe
import simbench as sb
import warnings

warnings.filterwarnings("ignore", category=FutureWarning)

from potpourri.models_multi_period.ACOPF_multi_period import ACOPF_multi_period
from potpourri.technologies.battery import Battery_multi_period

# %% md
# ## 2. Network and horizon
#
# We optimise over 96 time steps (one full day at 15-minute resolution).
# SimBench provides the load and generation profiles automatically.

# %%
net = sb.get_simbench_net("1-LV-urban6--0-sw")

fromT = 0
toT = 96  # 96 × 15 min = 24 hours

# Operational limits
net.bus["max_vm_pu"] = 1.05
net.bus["min_vm_pu"] = 0.95
net.line["max_loading_percent"] = 80.0

net.sgen["controllable"] = True
net.sgen["max_p_mw"] = net.sgen["p_mw"]
net.sgen["min_p_mw"] = 0.0

net.ext_grid["max_q_mvar"] = 500.0
net.ext_grid["min_q_mvar"] = -500.0

# %% md
# ## 3. Baseline — no battery storage

# %%
opf_base = ACOPF_multi_period(net, toT=toT, fromT=fromT)
opf_base.add_OPF()
opf_base.add_voltage_deviation_objective()
opf_base.solve(solver="ipopt", print_solver_output=False)

v_dev_base = [
    sum(
        (pe.value(opf_base.model.v[b, t]) - 1.0) ** 2 for b in opf_base.model.B
    )
    for t in opf_base.model.T
]
print(f"Baseline  |  total Σ(v-1)²: {sum(v_dev_base):.4f}")

# %% md
# ## 4. With battery storage — scenario-based placement
#
# `Battery_multi_period` places batteries randomly at a fraction of buses
# controlled by the `scenario` argument:
#
# | scenario | penetration |
# |----------|-------------|
# | 0        | 1 %         |
# | 1        | 7.9 %       |
# | 2        | 9.9 %       |
# | 3        | 10.6 %      |
#
# We use `scenario=1` (≈ 8 % of buses).

# %%
opf_bat = ACOPF_multi_period(net, toT=toT, fromT=fromT)

battery = Battery_multi_period(opf_bat.net, T=toT - fromT, scenario=1)
battery.get_all(opf_bat.model)

opf_bat.add_OPF()
opf_bat.add_voltage_deviation_objective()
opf_bat.solve(solver="ipopt", print_solver_output=False)

v_dev_bat = [
    sum((pe.value(opf_bat.model.v[b, t]) - 1.0) ** 2 for b in opf_bat.model.B)
    for t in opf_bat.model.T
]
print(f"With batteries  |  total Σ(v-1)²: {sum(v_dev_bat):.4f}")
improvement = (sum(v_dev_base) - sum(v_dev_bat)) / sum(v_dev_base) * 100
print(f"Voltage deviation reduction: {improvement:.1f} %")

# %% md
# ## 5. Explicit constructor parameters
#
# Instead of a scenario, every parameter can be set directly.  This is
# useful when sizing batteries for a specific grid or study:
#
# | parameter          | meaning                                         |
# |--------------------|-------------------------------------------------|
# | `penetration`      | % of non-slack buses to equip                   |
# | `power_pu`         | symmetric charge/discharge limit (p.u.)         |
# | `capacity_pu_h`    | energy capacity (p.u. power × hours)            |
# | `efficiency`       | one-way charge/discharge efficiency (0–1)       |
# | `soc_min`          | minimum state of charge                         |
# | `initial_soc_fraction` | SOC at t=0 as fraction of `soc_max`         |
#
# Convert from physical units: `power_pu = power_MW / net.sn_mva`

# %%
baseMVA = net.sn_mva

opf_custom = ACOPF_multi_period(net, toT=toT, fromT=fromT)
battery_custom = Battery_multi_period(
    opf_custom.net,
    T=toT - fromT,
    penetration=20.0,  # 20 % of non-slack buses
    power_pu=0.012 / baseMVA,  # 12 kW charge/discharge limit
    capacity_pu_h=0.030 / baseMVA,  # 30 kWh capacity
    efficiency=0.95,
    soc_min=0.1,
    initial_soc_fraction=0.4,
)
battery_custom.get_all(opf_custom.model)
print(
    f"Custom battery: {len(battery_custom.random_indexes)} batteries placed, "
    f"power={battery_custom.bat_power:.4f} p.u., "
    f"cap={battery_custom.bat_cap:.4f} p.u.·h, "
    f"η={battery_custom.bat_efficiency}"
)

# %% md
# ## 6. Battery SOC profiles
#
# Inspect the state-of-charge trajectory for the scenario=1 case.

# %%
baseMVA = opf_bat.model.baseMVA

print("\nBattery SOC profiles (every 4 time steps = 1 hour):")
print(f"{'time':>6}", end="")
for b in opf_bat.model.BAT:
    print(f"  bat{b:>2}", end="")
print()

for t in list(opf_bat.model.T)[::4]:  # every hour
    hour = (t - fromT) / 4
    print(f"{hour:>5.1f}h", end="")
    for b in opf_bat.model.BAT:
        soc = pe.value(opf_bat.model.BAT_SOC[b, t])
        print(f"  {soc:>5.2f}", end="")
    print()

# %% md
# ## 7. Battery power dispatch
#
# Positive `BAT_P` = charging (consuming grid power), negative = discharging.

# %%
print("\nBattery power (MW) — first battery (bat 0), every hour:")
print(f"{'hour':>6}  {'BAT_P (MW)':>12}  {'action':>12}")
bat0 = list(opf_bat.model.BAT)[0]
for t in list(opf_bat.model.T)[::4]:
    hour = (t - fromT) / 4
    p_mw = pe.value(opf_bat.model.BAT_P[bat0, t]) * baseMVA
    action = (
        "charging"
        if p_mw > 1e-4
        else ("discharging" if p_mw < -1e-4 else "idle")
    )
    print(f"{hour:>6.1f}  {p_mw:>12.4f}  {action:>12}")

# %% md
# ## Conclusion
#
# The refactored `Battery_multi_period` accepts explicit parameters for every
# physical property of the battery (power limit, capacity, efficiency, SOC
# bounds, initial SOC), in addition to the scenario-based shorthand.  The
# `penetration=` keyword bypasses the fixed scenario lookup entirely, enabling
# sensitivity studies over battery sizing and deployment density.
