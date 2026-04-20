# %% md
# # Hosting Capacity Analysis with HC_ACOPF
#
# The hosting capacity (HC) of a distribution grid is the maximum amount of
# distributed generation that can be connected without violating operational
# limits. `HC_ACOPF` extends the AC OPF with:
#
# - Binary variables `y[w] ∈ {0, 1}` — whether wind generator `w` is active
# - An apparent power envelope per generator: `SWmin·y ≤ S² ≤ SWmax·y`
# - Q-P and Q-U droop curves enforced by grid code constraints
# - A default objective that maximises total wind generation minus network losses
#
# If the network's `sgen` table has no `wind_hc` column, `HC_ACOPF` automatically
# places one candidate wind generator at every bus that is not an external grid bus.
#
# We will:
# 1. Run a baseline AC power flow
# 2. Solve the HC optimisation and identify which buses can host wind power
# 3. Sweep the `eps` parameter to trade off wind generation against losses

# %% md
# ## 1. Imports

# %%
import os

os.environ["NEOS_EMAIL"] = "test-email@test.com"
import copy
import pyomo.environ as pe
import pandapower as pp
import simbench as sb
import warnings

warnings.filterwarnings("ignore", category=FutureWarning)

from potpourri.models.HC_ACOPF import HC_ACOPF

# %% md
# ## 2. Load network and set limits
#
# We use a low-voltage rural network under a high-load, low-wind scenario so
# that the hosting capacity is constrained by line loading and voltage limits.

# %%
net = sb.get_simbench_net("1-LV-rural1--0-sw")

profiles = sb.get_absolute_values(net, profiles_instead_of_study_cases=True)
select_profile_idx = 1190  # high-load winter evening

net.sgen["p_mw"] = profiles[("sgen", "p_mw")].iloc[select_profile_idx]
net.load["p_mw"] = profiles[("load", "p_mw")].iloc[select_profile_idx]
net.load["q_mvar"] = profiles[("load", "q_mvar")].iloc[select_profile_idx]

# Operational limits
net.bus["max_vm_pu"] = 1.05
net.bus["min_vm_pu"] = 0.95
net.line["max_loading_percent"] = 80.0
net.ext_grid["max_q_mvar"] = 500.0
net.ext_grid["min_q_mvar"] = -500.0

pp.runpp(net, voltage_depend_loads=False)
print(f"Baseline max line loading:  {net.res_line.loading_percent.max():.1f}%")
print(
    f"Baseline voltage range:     {net.res_bus.vm_pu.min():.3f} – {net.res_bus.vm_pu.max():.3f} p.u."
)

# %% md
# ## 3. Solve HC-OPF
#
# `HC_ACOPF` places candidate wind generators at every non-slack bus. The solver
# finds the binary assignment `y[w]` and the operating point `(psG[w], qsG[w])`
# that maximises total injected wind power while keeping all constraints feasible.

# %%
hc = HC_ACOPF(net)
hc.add_OPF()
hc.solve(solver="neos", print_solver_output=False)

print(
    f"\nObjective (wind generation − losses): {pe.value(hc.model.obj):.4f} p.u."
)

# %% md
# ## 4. Examine results
#
# The binary variable `y[w]` tells us which buses can host a wind turbine.

# %%
baseMVA = hc.model.baseMVA

active_wind = []
inactive_wind = []

for w in hc.model.WIND_HC:
    y = pe.value(hc.model.y[w])
    p_mw = pe.value(hc.model.psG[w]) * baseMVA
    q_mvar = pe.value(hc.model.qsG[w]) * baseMVA
    bus = hc.net.sgen.bus.iloc[w]
    if y and y > 0.5:
        active_wind.append((w, bus, p_mw, q_mvar))
    else:
        inactive_wind.append((w, bus))

print(f"\nActive wind generators ({len(active_wind)}):")
print(f"  {'sgen':>5}  {'bus':>5}  {'P (MW)':>8}  {'Q (Mvar)':>10}")
for w, bus, p, q in active_wind:
    print(f"  {w:>5}  {bus:>5}  {p:>8.3f}  {q:>10.3f}")

print(
    f"\nInactive buses ({len(inactive_wind)}) — hosting would violate limits:"
)
print(
    "  " + ", ".join(f"bus {b}" for _, b in inactive_wind[:10]),
    "..." if len(inactive_wind) > 10 else "",
)

total_wind_mw = sum(p for _, _, p, _ in active_wind)
print(f"\nTotal hosted wind capacity: {total_wind_mw:.2f} MW")

# %% md
# ## 5. Apparent power envelope sweep
#
# By tightening the minimum apparent power `SWmin`, we can force the solver to
# commit to fewer but larger wind installations. This corresponds to a practical
# constraint on minimum generator size.

# %%
swmin_values = [0.0, 0.01, 0.02, 0.05]  # in p.u. (multiply by baseMVA for MW)

print(f"{'SWmin (MW)':>12}  {'Active sites':>14}  {'Total wind (MW)':>17}")
for swmin_pu in swmin_values:
    hc_sweep = HC_ACOPF(copy.deepcopy(net))
    hc_sweep._calc_opf_parameters(SWmin=swmin_pu * baseMVA)
    hc_sweep.add_OPF()
    hc_sweep.solve(solver="neos", print_solver_output=False)

    n_active = sum(
        1
        for w in hc_sweep.model.WIND_HC
        if (v := pe.value(hc_sweep.model.y[w])) and v > 0.5
    )
    total_mw = sum(
        pe.value(hc_sweep.model.psG[w]) * baseMVA
        for w in hc_sweep.model.WIND_HC
        if (v := pe.value(hc_sweep.model.y[w])) and v > 0.5
    )
    print(f"{swmin_pu * baseMVA:>12.2f}  {n_active:>14}  {total_mw:>17.2f}")

# %% md
# ## 6. Weighted wind-vs-loss objective
#
# `add_loss_obj()` replaces the default objective with a weighted combination:
#
#     max  ε · Σ psG[w]  +  (1 − ε) · (− Σ losses)
#
# `eps = 1.0` is equivalent to the default (pure wind maximisation). As `eps`
# decreases, the solver increasingly prioritises loss reduction over wind hosting.

# %%
eps_values = [1.0, 0.8, 0.5]

print(
    f"{'eps':>6}  {'Active sites':>14}  {'Total wind (MW)':>17}  {'Losses (MW)':>13}"
)
for eps in eps_values:
    hc_w = HC_ACOPF(copy.deepcopy(net))
    hc_w.add_OPF()
    hc_w.add_loss_obj()
    hc_w.model.eps.set_value(eps)
    hc_w.solve(solver="neos", print_solver_output=False)

    n_active = sum(
        1
        for w in hc_w.model.WIND_HC
        if (v := pe.value(hc_w.model.y[w])) and v > 0.5
    )
    total_mw = sum(
        pe.value(hc_w.model.psG[w]) * baseMVA
        for w in hc_w.model.WIND_HC
        if (v := pe.value(hc_w.model.y[w])) and v > 0.5
    )
    losses = sum(
        (pe.value(hc_w.model.pLfrom[l]) + pe.value(hc_w.model.pLto[l]))
        * baseMVA
        for l in hc_w.model.L
    )
    print(f"{eps:>6.1f}  {n_active:>14}  {total_mw:>17.2f}  {losses:>13.4f}")

# %% md
# ## Conclusion
#
# `HC_ACOPF` automates the hosting capacity problem: it finds which buses can
# actively host wind power and how much, while enforcing the full AC power flow
# physics and grid-code reactive power requirements. The `eps` parameter and the
# `SWmin` apparent power floor allow sensitivity analysis without rebuilding the
# model from scratch.
