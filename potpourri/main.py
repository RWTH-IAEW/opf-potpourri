import os
import config
import potpourri.scripts.run as potpourri
from potpourri.scripts.utils import showParameter

showParameter(testcase=config.test_case)
    
# Check if the 'WSL_DISTRO_NAME' environment variable is present
if 'WSL_DISTRO_NAME' in os.environ:
    from potpourri.scripts.set_gurobi_key import set_gurobi_key, retrieve_wls_gurobi_license
    gurobi_license = retrieve_wls_gurobi_license()
    set_gurobi_key(gurobi_license)
else:
    print("System is not running on Windows Subsystem for Linux (WSL).")

potpourri.dcopf(tc=config.test_case, solver="glpk", print_output=True, print_solver_output=True, print_model=False)
potpourri.dcopf(tc=config.test_case, solver="gurobi", print_output=True, print_solver_output=True, print_model=False,
                objective="quadratic")
potpourri.acopf(tc=config.test_case, solver="ipopt", print_output=True, print_solver_output=True, print_model=False)
potpourri.acopf(tc=config.test_case, solver="ipopt", print_output=True, print_solver_output=True, print_model=False,
                objective="min_gen")
potpourri.acopf(tc=config.test_case, solver="ipopt", print_output=True, print_solver_output=True, print_model=False,
                objective="max_cost")
potpourri.scopf(tc=config.test_case, solver="gurobi", print_output=True, print_solver_output=True, print_model=False)
