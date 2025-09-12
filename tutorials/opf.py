#%% md
# # Power System Analysis with potpourri and ACOPF
# 
# This notebook demonstrates how to use the potpourri package for power system analysis, specifically focusing on the solution of an Optimal Power Flow (ACOPF) problem. We will walk through the process of setting up the network, configuring the data, and solving the ACOPF problem using the potpourri.models.ACOPF_base class.
#%% md
# 1. Import Libraries
# 
# We begin by importing the required libraries for the analysis, including simbench, pyomo, and pandapower.
#%%
import copy
import numpy as np
import simbench as sb
import pyomo.environ as pe
import pandapower as pp
import warnings
warnings.filterwarnings('ignore', category=FutureWarning)
#%% md
# The ACOPF model from potpourri is also imported to set up and solve the optimal power flow problem.
#%%
from potpourri.models.ACOPF_base import ACOPF
#%% md
# 2. Load Network Data
# 
# In this section, we load the network model using simbench. We can use different network scenarios provided by Simbench; in this case, we'll use a low-voltage rural network.
# We define case keys for different load profiles (lW and hL). These profiles will be used to simulate different operational conditions in the network.
# The profiles variable is populated with absolute values from the Simbench data, which provides the power consumption and generation profiles.
# We select a specific profile index (e.g., select_profile_idx = 1190) to assign power generation and load data to the network.
#%%
# net = sb.get_simbench_net("1-HV-mixed--0-no_sw")
net = sb.get_simbench_net("1-LV-rural1--0-sw")

case_keys = ['lW', 'hL']

hcs = []
obj = []

# Load net data
profiles = sb.get_absolute_values(net, profiles_instead_of_study_cases=True)

select_profile_idx = 1190

net.sgen["p_mw"] = profiles[('sgen', 'p_mw')].iloc[select_profile_idx]
net.load["p_mw"] = profiles[('load', 'p_mw')].iloc[select_profile_idx]
net.load["q_mvar"] = profiles[('load', 'q_mvar')].iloc[select_profile_idx]
#%% md
# 4. Set Maximum and Minimum Reactive Power for Generators
# 
# Next, we calculate the maximum and minimum reactive power for each generator based on a fixed power factor (0.95 in this case). This step is crucial for setting up reactive power limits for the optimization.
#%%
power_factor = 0.95
q_max = np.sqrt((net.sgen["p_mw"] / power_factor) ** 2 - net.sgen["p_mw"] ** 2)
net.sgen["max_q_mvar"] = q_max
net.sgen["min_q_mvar"] = -q_max
#%% md
# 5. Set Generator Power Limits and Run Power Flow
# 
# We then set the power limits for each generator (max_p_mw and min_p_mw) and mark them as controllable. After this setup, we perform a power flow calculation using pandapower to get the initial state of the grid.
#%%
net.sgen["max_p_mw"] = net.sgen["p_mw"]
net.sgen["min_p_mw"] = net.sgen["p_mw"]
net.sgen['controllable'] = True
pp.runpp(net)
#%% md
# 6. Create a Copy of the Network for OPF
# 
# We create a deep copy of the network data to set up the optimal power flow (OPF) scenario. The external grid (ext_grid) power limits are updated based on the power flow results.
#%%
net_case_opf = copy.deepcopy(net)

net_case_opf.ext_grid["max_p_mw"] = net_case_opf.res_ext_grid["p_mw"]
net_case_opf.ext_grid["min_p_mw"] = net_case_opf.res_ext_grid["p_mw"]
net_case_opf.ext_grid["max_q_mvar"] = 0.05
net_case_opf.ext_grid["min_q_mvar"] = -0.05
#%% md
# 7. Set Up ACOPF Model
# 
# Now, we instantiate the ACOPF model from potpourri and add the optimization problem. In this case, we are minimizing the reactive power flow in addition to other power flow objectives. Before solving the OPF problem, we print the grid state (power generation and reactive power) for the generators:
#%%
hc = ACOPF(net_case_opf)
hc.add_OPF()
hc.add_reactive_power_flow_objective()
print("Grid-state before OPF:")
print(net.sgen[["p_mw", "q_mvar"]])
print("Grid-state before OPF:")
print("SGEN P:")
for g in hc.model.sG:
    print(pe.value(hc.model.PsG[g]))
print("SGEN Q:")
for g in hc.model.sG:
    print(pe.value(hc.model.QsG[g]))
#%% md
# 8. Solve the ACOPF Problem
# 
# We solve the ACOPF optimization problem using the ipopt solver. The solver adjusts the power generation and reactive power to minimize the objective function while satisfying the grid constraints. After solving the OPF problem, we print the updated grid state and generator outputs:
#%%
hc.solve(solver='neos', print_solver_output=False)
# hc.solve(solver='ipopt', print_solver_output=False)

# Print summary of changes
print("Grid-state after OPF:")
print(net.sgen[["p_mw", "q_mvar"]])
print("Grid-state after OPF:")
print("SGEN P:")
for g in hc.model.sG:
    print(pe.value(hc.model.psG[g]))
print("SGEN Q:")
for g in hc.model.sG:
    print(pe.value(hc.model.qsG[g]))
#%% md
# 9. Run OPF for Different Load Cases
# 
# In this section, we run the OPF for different load cases (lW and hL). We update the load values (p_mw and q_mvar) for each case and perform the optimization.
#%%
    for case in case_keys:
        net_case = copy.deepcopy(net)

        factors = net_case.loadcases.loc[case]

        net_case.load.p_mw *= factors['pload']
        net_case.load.q_mvar *= factors['qload']

        net_case.sgen.loc[net_case.sgen.type == 'Wind', 'scaling'] = factors['Wind_p']
        net_case.sgen.loc[net_case.sgen.type == 'PV', 'scaling'] = factors['PV_p']
        net_case.sgen.loc[(net_case.sgen.type != 'Wind') & (net_case.sgen.type != 'PV'), 'scaling'] = factors['RES_p']

        net_case.ext_grid.vm_pu = factors['Slack_vm']

        hc = ACOPF(net_case)
        hc.add_voltage_deviation_objective()
        hc.solve(solver='neos')

        hcs.append(copy.deepcopy(hc))
        obj.append(pe.value(hc.model.obj_v_deviation))

        for g in hc.model.sG:
           print(f"Static Gen P[MW] {g}: {pe.value(hc.model.psG[g])}")
           print(f"Static Gen Q[MVar] {g}: {pe.value(hc.model.qsG[g])}")

        line_flow = pe.value(hc.model.pLfrom[0])
        print(f"Line flow over first line P[MW]: {line_flow}")
#%% md
# Conclusion
# 
# In this notebook, we have demonstrated how to use the potpourri package to solve an Optimal Power Flow (ACOPF) problem. We started by loading network data from Simbench, set up the ACOPF optimization problem, solved it using Pyomo, and analyzed the results. This workflow can be extended to other power system models and optimization objectives.