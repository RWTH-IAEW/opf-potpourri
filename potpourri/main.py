import os
os.system('clear')

import scripts.run as potpourri

test_case_folder = "potpourri/testcases/"

for i in os.listdir(test_case_folder):
    if i.endswith(".xlsx") and i != "tempelate.xlsx":
        test_case = os.path.join(test_case_folder, i)
        print(test_case)
        i = []
        # potpourri.dcopf(tc=test_case, solver="glpk", print_output=False, print_solver_output=True, print_model=False)
        try:
            potpourri.dcopf(tc=test_case, solver="gurobi", print_output=False, print_solver_output=False, print_model=False)
        except:
            print(f"Model {test_case} has not succesfully been solved")
            i.append(test_case)
        print(i)
        # potpourri.dcopf(tc=test_case, solver="gurobi", print_output=False, print_solver_output=True, print_model=False, objective="quadratic")
        # potpourri.acopf(tc=test_case, solver="ipopt", print_output=False, print_solver_output=True, print_model=False)
        # potpourri.acopf(tc=test_case, solver="ipopt", print_output=False, print_solver_output=True, print_model=False, objective="min_gen")
        # potpourri.acopf(tc=test_case, solver="ipopt", print_output=False, print_solver_output=True, print_model=False, objective="max_cost")
        # potpourri.scopf(tc=test_case, solver="gurobi", print_output=False, print_solver_output=True, print_model=False)