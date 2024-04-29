import os
os.system('clear')

import potpourri.scripts.run as potpourri

test_case_folder = "testcases/matpower"

infeas_list = []
success_list = []

for i in os.listdir(test_case_folder):
    if i.endswith(".xlsx") and i != "tempelate.xlsx":
        test_case = os.path.join(test_case_folder, i)
        print(test_case)
        
        try:
            potpourri.dcopf(tc=test_case, solver="gurobi", print_output=False, print_solver_output=False, print_model=False)
            success_list.append(test_case)
        except:
            print(f"Model {test_case} has not succesfully been solved")
            infeas_list.append(test_case)
            
        # potpourri.dcopf(tc=test_case, solver="glpk", print_output=False, print_solver_output=True, print_model=False)
        # potpourri.dcopf(tc=test_case, solver="gurobi", print_output=False, print_solver_output=True, print_model=False, objective="quadratic")
        # potpourri.acopf(tc=test_case, solver="ipopt", print_output=False, print_solver_output=True, print_model=False)
        # potpourri.acopf(tc=test_case, solver="ipopt", print_output=False, print_solver_output=True, print_model=False, objective="min_gen")
        # potpourri.acopf(tc=test_case, solver="ipopt", print_output=False, print_solver_output=True, print_model=False, objective="max_cost")
        # potpourri.scopf(tc=test_case, solver="gurobi", print_output=False, print_solver_output=True, print_model=False)
        
# Open the log file in write mode
with open('test_networks.txt', 'w') as log_file:
    log_file.write("Infeasible networks: \n")
    for item in infeas_list:
        log_file.write("- " + str(item) + '\n')
    log_file.write("Feasible networks: \n")
    for item in success_list:
        log_file.write("- " + str(item) + '\n')