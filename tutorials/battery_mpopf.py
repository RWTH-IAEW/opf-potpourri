# %% md
# # Multi-Period OPF with Battery Storage
#
# Battery storage can shift energy in time: charge when generation is high or
# demand is low, discharge when the reverse is true. This tutorial adds
# `Battery_multi_period` to a 24-hour AC OPF and shows how the batteries affect
# voltage profiles and line loading over the day.
#
# We will:
# 1. Solve a 24-hour AC OPF without batteries as a baseline
# 2. Add battery storage and re-solve
# 3. Compare voltage deviation and peak line loading between both cases

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
opf_base.solve(solver="gurobi_direct", print_solver_output=False)

# Collect voltage deviation and max line loading across all time steps
v_dev_base = []
loading_base = []
baseMVA = opf_base.model.baseMVA

for t in opf_base.model.T:
    v_dev_base.append(
        sum(
            (pe.value(opf_base.model.v[b, t]) - 1.0) ** 2
            for b in opf_base.model.B
        )
    )

print(f"Baseline  |  total Σ(v-1)²: {sum(v_dev_base):.4f}")

# %% md
# ## 4. With battery storage
#
# `Battery_multi_period` places batteries randomly at a percentage of buses
# controlled by the `scenario` argument:
#
# | scenario | penetration |
# |----------|-------------|
# | 0        | 1 %         |
# | 1        | 7.9 %       |
# | 2        | 9.9 %       |
# | 3        | 10.6 %      |
#
# We use scenario=1 (≈ 8% of buses) and pass `T` as the number of time steps
# so the battery object can build its internal time range.

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
print(f"Voltage deviation reduction: {improvement:.1f}%")

# %% md
# ## 5. Battery SOC profiles
#
# Inspect the state-of-charge trajectory for each battery over the day.

# %%
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
# ## 6. Battery power dispatch
#
# Positive BAT_P = charging (consuming grid power), negative = discharging.

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
# Adding battery storage allowed the multi-period OPF to shift charging to
# periods of low demand or high PV generation, reducing the total voltage
# deviation across the day. The SOC profile confirms the expected charge-
# during-day / discharge-in-evening behaviour typical of residential PV+battery
# systems.
