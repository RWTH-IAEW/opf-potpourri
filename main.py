import os
import potpourri.scripts.run as potpourri
from potpourri.scripts.utils import showParameter

test_case_folder = "testcases/matpower/"
test_case = "case14_battery2.xlsx" # "case14_battery.xlsx" #"case24_ieee_rts.xlsx"
test_case = os.path.join(test_case_folder, test_case)
showParameter(testcase=test_case)
    
# Check if the 'WSL_DISTRO_NAME' environment variable is present
if 'WSL_DISTRO_NAME' in os.environ:
    from potpourri.scripts.set_gurobi_key import set_gurobi_key, retrieve_wls_gurobi_license
    gurobi_license = retrieve_wls_gurobi_license()
    set_gurobi_key(gurobi_license)
else:
    print("System is not running on Windows Subsystem for Linux (WSL).")

#potpourri.dcopf(tc=test_case, solver="glpk", print_output=True, print_solver_output=True, print_model=False)
# potpourri.dcopf(tc=test_case, solver="gurobi", print_output=True, print_solver_output=True, print_model=False, objective="quadratic")
#potpourri.acopf(tc=test_case, solver="ipopt", print_output=True, print_solver_output=True, print_model=False)
# potpourri.acopf(tc=test_case, solver="ipopt", print_output=True, print_solver_output=True, print_model=False, objective="min_gen")
#potpourri.acopf(tc=test_case, solver="ipopt", print_output=True, print_solver_output=True, print_model=False, objective="multiperiod")
#potpourri.uc(tc=test_case, solver="glpk", print_output=True, print_solver_output=True, print_model=False)
# potpourri.scopf(tc=test_case, solver="gurobi", print_output=True, print_solver_output=True, print_model=False)
potpourri.acopf(tc=test_case, solver="ipopt", print_output=True, print_solver_output=True, print_model=False, objective="multiperiod_battery")
#potpourri.acopf(tc=test_case, solver="mindtpy", print_output=True, print_solver_output=True, print_model=False, objective="multiperiod_battery")
#potpourri.acopf(tc=test_case, solver="apopt.py", print_output=True, print_solver_output=True, print_model=False, objective="multiperiod_battery")
#