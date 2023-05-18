import os
import potpourri.scripts.run as potpourri
from potpourri.scripts.utils import showParameter

test_case_folder = "testcases/matpower/"
test_case = "case24_ieee_rts.xlsx"
test_case = os.path.join(test_case_folder, test_case)
showParameter(testcase=test_case)

potpourri.dcopf(tc=test_case, solver="gurobi", print_output=True, print_solver_output=True, print_model=False)
potpourri.dcopf(tc=test_case, solver="gurobi", print_output=True, print_solver_output=True, print_model=False, objective="quadratic")
potpourri.acopf(tc=test_case, solver="ipopt", print_output=True, print_solver_output=True, print_model=False)
potpourri.acopf(tc=test_case, solver="ipopt", print_output=True, print_solver_output=True, print_model=False, objective="min_gen")
potpourri.acopf(tc=test_case, solver="ipopt", print_output=True, print_solver_output=True, print_model=False, objective="max_cost")
potpourri.scopf(tc=test_case, solver="gurobi", print_output=True, print_solver_output=True, print_model=False)