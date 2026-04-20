# %% md
# # EV Charging Optimisation with Day-Ahead Prices
#
# Electric vehicles connected to the grid can be smart-charged: the optimiser
# decides *when* to charge each vehicle to minimise the total electricity cost,
# given day-ahead market prices, while ensuring every vehicle has enough charge
# for its planned trips.
#
# The `EV_multi_period` module supports:
# - **V1G** — unidirectional charging only
# - **V2G** — bidirectional charging/discharging (vehicle-to-grid)
# - **Non-controllable** — fixed charging profile (reference case)
#
# The `add_arbitrage_objective()` minimises the sum of charging cost over the
# horizon: `Σ_v Σ_t  p_opf[v,t] · price[t] · Δt`
#
# We will:
# 1. Set up a 24-hour OPF with 10 EVs and the arbitrage objective
# 2. Compare EV-managed vs. uncontrolled charging profiles
# 3. Examine per-vehicle SOC trajectories

# %% md
# ## 1. Imports

# %%
import pyomo.environ as pe
import simbench as sb
import warnings

warnings.filterwarnings("ignore", category=FutureWarning)

from potpourri.models_multi_period.ACOPF_multi_period import ACOPF_multi_period

# %% md
# ## 2. Network and time horizon

# %%
net = sb.get_simbench_net("1-LV-urban6--0-sw")

fromT = 0
toT = 96  # one full day, 15-minute intervals

net.bus["max_vm_pu"] = 1.05
net.bus["min_vm_pu"] = 0.95
net.line["max_loading_percent"] = 100.0
net.ext_grid["max_q_mvar"] = 500.0
net.ext_grid["min_q_mvar"] = -500.0

net.sgen["controllable"] = True
net.sgen["max_p_mw"] = net.sgen["p_mw"]
net.sgen["min_p_mw"] = 0.0

# %% md
# ## 3. Build the multi-period model with EVs
#
# Passing `num_vehicles=10` loads 10 EV profiles from `data/emob_profiles.pkl`
# and randomly assigns the corresponding charging points to network buses.
# Day-ahead prices are loaded automatically from `data/da_prices_hourly_2022.xlsx`.

# %%
opf = ACOPF_multi_period(net, toT=toT, fromT=fromT, num_vehicles=10)
opf.add_OPF()
opf.add_arbitrage_objective()  # minimise Σ p_opf · price · Δt

print(f"Number of vehicles: {len(list(opf.model.veh))}")
print(f"V2G vehicles:       {len(list(opf.model.veh_v2g))}")
print(f"V1G vehicles:       {len(list(opf.model.veh_v1g))}")
print(f"Non-controllable:   {len(list(opf.model.veh_nc))}")

# %% md
# ## 4. Solve

# %%
opf.solve(solver="ipopt", print_solver_output=False, time_limit=600)

baseMVA = opf.model.baseMVA

# %% md
# ## 5. Day-ahead price profile
#
# The optimiser has access to the full 24-hour price signal. Prices are in €/MWh
# at hourly resolution, repeated 4× for 15-minute time steps.

# %%
# Retrieve price profile from EV object
ev_obj = next(
    f for f in opf.flexibilities if f.__class__.__name__ == "EV_multi_period"
)
prices = [pe.value(opf.model.p_da[t]) for t in opf.model.T]

print("Day-ahead prices (€/MWh), first 24 time steps:")
print(f"{'t':>4}  {'hour':>6}  {'price':>8}")
for t in list(opf.model.T)[:24]:
    hour = (t - fromT) / 4
    print(f"{t:>4}  {hour:>6.2f}  {prices[t - fromT]:>8.2f}")

# %% md
# ## 6. Aggregate EV charging profile

# %%
print(
    "\nAggregate EV power (MW) over the day (controllable vs. uncontrolled):"
)
print(f"{'hour':>6}  {'optimised (MW)':>16}  {'reference (MW)':>16}")

for t in list(opf.model.T)[::4]:  # every hour
    hour = (t - fromT) / 4

    p_opt = sum(
        pe.value(opf.model.p_opf[v, t]) * baseMVA for v in opf.model.veh
    )

    # Reference: use the fixed uncontrolled profile from p_req_ev if available
    p_ref = sum(
        pe.value(opf.model.p_req_ev[v, t]) * baseMVA
        for v in opf.model.veh
        if (v, t) in opf.model.p_req_ev
    )

    print(f"{hour:>6.1f}  {p_opt:>16.3f}  {p_ref:>16.3f}")

# %% md
# ## 7. Per-vehicle SOC trajectories
#
# Each vehicle must meet its departure SOC target. The `buffer_soc` variable
# provides soft relaxation so the optimiser can slightly miss the target rather
# than becoming infeasible.

# %%
print("\nSOC trajectories for first 3 vehicles (every 4 time steps = 1 hour):")
vehs_to_show = list(opf.model.veh)[:3]

header = f"{'hour':>6}" + "".join(
    f"  {'veh' + str(v):>7}" for v in vehs_to_show
)
print(header)

for t in list(opf.model.T)[::4]:
    hour = (t - fromT) / 4
    row = f"{hour:>6.1f}"
    for v in vehs_to_show:
        soc = pe.value(opf.model.soc[v, t])
        row += f"  {soc:>7.3f}"
    print(row)

# %% md
# ## 8. Total charging cost comparison

# %%
# Optimised cost
cost_opt = sum(
    pe.value(opf.model.p_opf[v, t])
    * baseMVA
    * pe.value(opf.model.p_da[t])
    * (opf.model.deltaT / 60.0 if hasattr(opf.model, "deltaT") else 0.25)
    for v in opf.model.veh
    for t in opf.model.T
)

# Uncontrolled reference cost
cost_ref = sum(
    pe.value(opf.model.p_req_ev[v, t])
    * baseMVA
    * pe.value(opf.model.p_da[t])
    * 0.25
    for v in opf.model.veh
    for t in opf.model.T
    if (v, t) in opf.model.p_req_ev
)

print(f"Optimised charging cost:    €{cost_opt:.2f}")
print(f"Uncontrolled charging cost: €{cost_ref:.2f}")
if cost_ref > 0:
    print(
        f"Cost saving:                €{cost_ref - cost_opt:.2f}  "
        f"({(cost_ref - cost_opt) / cost_ref * 100:.1f}%)"
    )

# %% md
# ## Conclusion
#
# The arbitrage-optimal EV charging schedule shifts charging to low-price hours
# (typically overnight and early morning) and away from the evening peak.
# V2G-capable vehicles discharge during high-price periods, further reducing
# net cost. The AC power flow constraints in the OPF ensure that this flexible
# charging does not violate grid limits at any time step.
