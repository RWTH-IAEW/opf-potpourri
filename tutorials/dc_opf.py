#%% md
# # DC Optimal Power Flow
#
# This tutorial introduces the linearised DC power flow formulation. The DC
# approximation ignores voltage magnitudes and reactive power, reducing the
# power flow equations to a linear system in voltage angles. It runs faster
# than AC OPF and is useful for initial screening of operational limits.
#
# We will:
# 1. Run a pure DC power flow and verify it against pandapower's built-in solver
# 2. Set up and solve a DC OPF that maximises use of local renewable generation
# 3. Inspect line flows and generator dispatch

#%% md
# ## 1. Imports

#%%
import os
os.environ['NEOS_EMAIL'] = "test-email@test.com"
import numpy as np
import pyomo.environ as pe
import pandapower as pp
import simbench as sb
import warnings
warnings.filterwarnings('ignore', category=FutureWarning)

from potpourri.models.DC import DC
from potpourri.models.DCOPF import DCOPF

#%% md
# ## 2. Load network
#
# We use a SimBench low-voltage rural network. After loading the absolute
# generation and load profiles we pick a single time step for the static
# optimisation.

#%%
net = sb.get_simbench_net("1-LV-rural1--0-sw")

profiles = sb.get_absolute_values(net, profiles_instead_of_study_cases=True)
select_profile_idx = 672  # midday in summer

net.sgen["p_mw"]    = profiles[("sgen", "p_mw")].iloc[select_profile_idx]
net.load["p_mw"]    = profiles[("load", "p_mw")].iloc[select_profile_idx]
net.load["q_mvar"]  = profiles[("load", "q_mvar")].iloc[select_profile_idx]

#%% md
# ## 3. Pure DC power flow
#
# Instantiating `DC` runs `pp.runpp()` internally and builds the Pyomo model.
# Calling `solve()` without `add_OPF()` solves a pure power flow (no
# optimisation — all variables are fixed to their setpoints).

#%%
pp.rundcpp(net)                   # pandapower reference

pf = DC(net)
pf.solve(solver="neos", neos_opt="cplex", print_solver_output=False)

# Compare angle results at each bus
print("Bus voltage angles — pandapower vs Pyomo DC:")
print(f"{'Bus':>5}  {'pp va (°)':>12}  {'pyo va (°)':>12}  {'diff':>10}")
for b in pf.model.B:
    pp_va  = net.res_bus.va_degree.iloc[b]
    pyo_va = pe.value(pf.model.delta[b]) * 180 / np.pi
    print(f"{b:>5}  {pp_va:>12.4f}  {pyo_va:>12.4f}  {pp_va - pyo_va:>10.6f}")

#%% md
# The angles should agree to within numerical tolerance. Small differences arise
# from the per-unit normalisation convention.

#%% md
# ## 4. DC OPF — maximise local generation
#
# We now set renewable generators as controllable and find the dispatch that
# minimises the power drawn from the external grid (slack bus), effectively
# maximising local generation. This is a linear programme solved by GLPK.

#%%
import copy

net_opf = copy.deepcopy(net)

# Mark PV/wind generators as controllable with realistic limits
net_opf.sgen["controllable"] = True
net_opf.sgen["max_p_mw"]     = net_opf.sgen["p_mw"]
net_opf.sgen["min_p_mw"]     = net_opf.sgen["p_mw"] * 0.0   # can curtail down to zero

# External grid: allow both import and export
net_opf.ext_grid["max_p_mw"] = 10000.
net_opf.ext_grid["min_p_mw"] = -10000.

# Line loading limit
net_opf.line["max_loading_percent"] = 80.

#%%
dcopf = DCOPF(net_opf)
dcopf.add_OPF()

# Add a custom objective: minimise power imported from the external grid.
# Positive pG[g] = import; negative = export.
dcopf.model.obj = pe.Objective(
    expr=sum(dcopf.model.pG[g] for g in dcopf.model.G),
    sense=pe.minimize
)

dcopf.solve(solver='gurobi_direct', print_solver_output=False)

#%% md
# ## 5. Results

#%%
baseMVA = dcopf.model.baseMVA

print("== External grid dispatch ==")
for g in dcopf.model.G:
    p_mw = pe.value(dcopf.model.pG[g]) * baseMVA
    print(f"  Generator {g}: {p_mw:+.3f} MW  (+ = import, - = export)")

print("\n== Static generator dispatch ==")
for g in dcopf.model.sG:
    p_mw = pe.value(dcopf.model.psG[g]) * baseMVA
    print(f"  sgen {g}: {p_mw:.3f} MW")

print("\n== Line loading (top 5 most loaded) ==")
line_loading = {}
for l in dcopf.model.L:
    p_from = abs(pe.value(dcopf.model.pLfrom[l]) * baseMVA)
    p_to   = abs(pe.value(dcopf.model.pLto[l])   * baseMVA)
    smax   = pe.value(dcopf.model.SLmax[l])        * baseMVA
    # DC: S ≈ P (no reactive power)
    loading = max(p_from, p_to) / (smax if smax > 0 else 1) * 100
    line_loading[l] = loading

for l, loading in sorted(line_loading.items(), key=lambda x: -x[1])[:5]:
    print(f"  Line {l:>3}: {loading:.1f}% loading")

#%% md
# ## Conclusion
#
# The DC OPF dispatched the controllable generators to minimise external grid
# import while respecting the 80% line loading constraint. Because the DC
# formulation is linear, GLPK solves it in milliseconds — useful for fast
# pre-screening before running the more accurate but slower AC OPF.
