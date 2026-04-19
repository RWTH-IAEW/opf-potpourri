#%% md
# # Multi-period Optimal Power Flow (MPOPF)
# ---
# 
# 
# 
# The multiperiod ACOPF model from potpourri is imported to set up and solve the multiperiod OPF problem.
# 
# 
#%%
import os
os.environ['NEOS_EMAIL'] = "test-email@test.com"
import simbench as sb
import pyomo.environ as pe
import warnings
warnings.filterwarnings('ignore', category=FutureWarning)
#%%
from potpourri.models_multi_period.ACOPF_multi_period import ACOPF_multi_period
#%% md
# 1. Load network data
# 
# A low-voltage rural SimBench network model is used in this example.
# 
#%%
net = sb.get_simbench_net("1-LV-urban6--0-sw")
#%% md
# 
# 2. Set up optimization horizon
# 
# The optimization horizon spans from a starting time (fromT) to an ending time (toT). In this case, we set the horizon to 96 time steps (1 day) with 15-minute intervals. The standard timeseries data used for load and generation elements in the grid are the Simbench profiles for each grid element over the defined time horizon.
#%%
fromT = 0
toT = fromT + 96 # 1 day for 15 min intervals
#%% md
# 3. Customize Generator Power Limits
# 
# We then set the power limits for each generator (max_p_mw and min_p_mw) and mark them as controllable. If these limits are not explicitly set, the default max. values provided by SimBench are used.
#%%
net.sgen["max_p_mw"] = net.sgen["p_mw"]
net.sgen["min_p_mw"] = net.sgen["p_mw"]
net.sgen['controllable'] = True
#%% md
# 
# 4. Set up ACOPF Model
# 
# Now, we instantiate the multiperiod ACOPF model from potpourri and add the optimization problem. An optional input for the instance of the OPF is the number of electric vehicles 'num_vehicles' that can offer additional flexibility by optimizing their charging behaviour while taking their mobility needs into account. Currently, the optimization of up to 100 EVs is allowed. To take into account is that the addition of EVs increases the computation time of the optimization problem significantly.
# The objective function can be set seperately. In this case, we are minimizing the voltage deviation.
#%%
opf = ACOPF_multi_period(net, toT, fromT, num_vehicles=5)
opf.add_OPF()
opf.add_voltage_deviation_objective()
#%% md
# 4. Solve Multi-period OPF
# 
# We solve the ACOPF optimization problem using the ipopt solver. The solver adjusts the power generation and reactive power over the entire time horizon to minimize the objective function while satisfying the grid constraints. After solving the OPF problem, we access some results, such as the power generation of an exemplary pv generation over the entire optimization horizon.
# 
#%%
opf.solve(solver="neos", neos_opt='knitro', print_solver_output=False)

# example of how to access the results -> sg output throughout the time horizon
for t in opf.model.T:
    print(pe.value(opf.model.psG[0, t]))

#%% md
# Conclusion
# 
# In this notebook, we have demonstrated how to use the potpourri package to solve a Multi-Period Optimal Power Flow (MP-ACOPF) problem. We started by loading network data from Simbench, set up the multiperiod ACOPF optimization horizon and optimization problem, solved it using Pyomo, and showed an example of accessing the results. This workflow can be extended to other power system models and optimization objectives.
# 
# 