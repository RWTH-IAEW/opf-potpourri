<<<<<<< HEAD
#==================================================================
# runfile.py
# This is a top level OATS script. Simulations can be ran using this script
# ---Author---
# W. Bukhsh,
# wbukhsh@gmail.com
# OATS
# Copyright (c) 2017 by W. Bukhsh, Glasgow, Scotland
# OATS is distributed under the GNU GENERAL PUBLIC LICENSE v3. (see LICENSE file for details).
#==================================================================
import logging
import datetime
import os
from .runcase import runcase

# Get the current date and time
current_datetime = datetime.datetime.now()
# Format the datetime as desired
datetime_str = current_datetime.strftime("%Y-%m-%d_%H-%M-%S")

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)-8s %(message)s',
                    datefmt='%a, %d %b %Y %H:%M:%S',
                    filename=f'logs/log_{datetime_str}.log',
                    filemode='w')
logging.info("OATS log file")
logging.info("Program started")

oats_dir = os.path.dirname(os.path.realpath(__file__))
default_testcase = "/home/skortmann/potpourri/testcases/matpower/case24_ieee_rts.xlsx"

#----------------------------------------------------------------------
# DC Load flow
def dclf(tc='default',solver='ipopt',neos=True,out=0):
    """
    Solves DC load flow problem

    ARGUMENTS:
        **tc** (*.xlsx file)  - OATS test case. See OATS data format for details

        **solver** (str)  - name of a solver. Defualt is 'ipopt'

        **neos** (bool) - If True, the problem is solved using NEOS otherwise using a localy install solver.

        **out** (bool) - If True, the output is displayed on screen.
    """

    if tc == 'default':
        tc = default_testcase
    #options
    opt=({'neos':neos,\
    'solver':solver,'out':out})
    testcase = tc
    model ='DCLF'
    # ==log==
    logging.info("Solver selected: "+opt['solver'])
    logging.info("Testcase selected: "+testcase)
    logging.info("Model selected: "+model)
    runcase(testcase,model,opt)
    logging.info("Done!")
    
# AC Load flow
def aclf(tc='default',solver='ipopt',neos=True,out=0):
    """
    Solves AC load flow problem

    ARGUMENTS:
        **tc** (*.xlsx file)  - OATS test case. See OATS data format for details

        **solver** (str)  - name of a solver. Defualt is 'ipopt'

        **out** (bool) - If True, the output is displayed on screen.
    """
    if tc == 'default':
        tc = default_testcase
    opt=({'solver':solver,'out':out})
    testcase = tc
    model ='ACLF'
    # ==log==
    logging.info("Solver selected: "+opt['solver'])
    logging.info("Testcase selected: "+testcase)
    logging.info("Model selected: "+model)
    runcase(testcase,model,opt)
    logging.info("Done!")
    
# DC optimal power flow problem
def dcopf(tc='default',
          solver='ipopt', 
          print_output:bool = False,
          print_solver_output:bool = False,
          print_model:bool = False,
          objective:str = None):
    """
    Solves DC optimal power flow problem

    ARGUMENTS:
        **tc** (*.xlsx file)  - OATS test case. See OATS data format for details

        **solver** (str)  - name of a solver. Defualt is 'ipopt'
        
        **out** (bool) - If True, the output is displayed on screen.
    """

    if tc == 'default':
        tc = default_testcase
    #options
    opt=({'solver':solver,'print_output':print_output,'print_solver_output':print_solver_output,'print_model':print_model})
    testcase = tc
    if objective == None:
        model ='DCOPF'
    elif objective == "quadratic":
        model = 'DCOPF_quadratic'
    # ==log==
    logging.info("Solver selected: "+opt['solver'])
    logging.info("Testcase selected: "+testcase)
    logging.info("Model selected: "+model)
    runcase(testcase,model,opt)
    logging.info("Done!")

# AC optimal power flow problem
def acopf(tc:str = 'default', 
          solver:str = 'ipopt', 
          print_output:bool = False, 
          print_solver_output:bool = False,
          print_model:bool = False, 
          objective:str = None):
    """
    Solves AC optimal power flow problem

    ARGUMENTS:
        **tc** (*.xlsx file)  - OATS test case. See OATS data format for details

        **solver** (str)  - name of a solver. Defualt is 'ipopt'

        **out** (bool) - If True, the output is displayed on screen.
    """

    if tc == 'default':
        tc = default_testcase
    #options
    opt=({'solver':solver,'print_output':print_output,'print_solver_output':print_solver_output,'print_model':print_model})
    testcase = tc
    model = "ACOPF"
    if objective == "min_gen":
        model = "ACOPF_min_gen"
    elif objective == "max_cost":
        model = "ACOPF_max_cost"
    # ==log==
    logging.info("Solver selected: "+opt['solver'])
    logging.info("Testcase selected: "+testcase)
    logging.info("Model selected: "+model)
    runcase(testcase,model,opt)
    logging.info("Done!")

# security constrained optimal power flow problem
def scopf(tc:str = 'default', 
          solver:str = 'ipopt', 
          print_output:bool = False, 
          print_solver_output:bool = False,
          print_model:bool = False, 
          objective:str = None):
    """
    Solves security constrained optimal power flow problem

    ARGUMENTS:
        **tc** (*.xlsx file)  - OATS test case. See OATS data format for details

        **solver** (str)  - name of a solver. Defualt is 'ipopt'

        **out** (bool) - If True, the output is displayed on screen.
    """

    if tc == 'default':
        tc = default_testcase
    #options
    opt=({'solver':solver,'print_output':print_output,'print_solver_output':print_solver_output,'print_model':print_model})
    testcase = tc
    model ='SCOPF'
    # ==log==
    logging.info("Solver selected: "+opt['solver'])
    logging.info("Testcase selected: "+testcase)
    logging.info("Model selected: "+model)
    runcase(testcase,model,opt)
    logging.info("Done!")

# unit commitment problem
def uc(tc:str = 'default', 
          solver:str = 'ipopt', 
          print_output:bool = False, 
          print_solver_output:bool = False,
          print_model:bool = False, 
          objective:str = None):
    """
    Solves unit commitment problem

    ARGUMENTS:
        **tc** (*.xlsx file)  - OATS test case. See OATS data format for details

        **solver** (str)  - name of a solver. Defualt is 'ipopt'

        **out** (bool) - If True, the output is displayed on screen.
    """

    if tc == 'default':
        tc = default_testcase
    #options
    opt=({'solver':solver,'print_output':print_output,'print_solver_output':print_solver_output,'print_model':print_model})
    testcase = tc
    model ='UC'
    # ==log==
    logging.info("Solver selected: "+opt['solver'])
    logging.info("Testcase selected: "+testcase)
    logging.info("Model selected: "+model)
    runcase(testcase,model,opt)
    logging.info("Done!")

# user defined model
def model(model:str = 'ACOPF', 
          tc:str = 'default', 
          solver:str = 'ipopt', 
          print_output:bool = False, 
          print_solver_output:bool = False,
          print_model:bool = False, 
          objective:str = None):
    """
    Solves a user defined model

    ARGUMENTS:
        **model** (*.py file)  - A model file that describes the optimisation problem using PYOMO modelling language.

        **tc** (*.xlsx file)  - OATS test case. See OATS data format for details

        **solver** (str)  - name of a solver. Defualt is 'ipopt'

        **out** (bool) - If True, the output is displayed on screen.
    """

    if tc == 'default':
        tc = default_testcase
    #options
    opt=({'solver':solver,'print_output':print_output,'print_solver_output':print_solver_output,'print_model':print_model,'user_def_model':True})
    testcase = tc
    model =model
    # ==log==
    logging.info("Solver selected: "+opt['solver'])
    logging.info("Testcase selected: "+testcase)
    logging.info("Model selected: "+model)
    runcase(testcase,model,opt)
    logging.info("Done!")
=======
#==================================================================
# runfile.py
# This is a top level OATS script. Simulations can be ran using this script
# ---Author---
# W. Bukhsh,
# wbukhsh@gmail.com
# OATS
# Copyright (c) 2017 by W. Bukhsh, Glasgow, Scotland
# OATS is distributed under the GNU GENERAL PUBLIC LICENSE v3. (see LICENSE file for details).
#==================================================================
import logging
import datetime
import os
from .runcase import runcase

# Get the current date and time
current_datetime = datetime.datetime.now()
# Format the datetime as desired
datetime_str = current_datetime.strftime("%Y-%m-%d_%H-%M-%S")

import logging

# Configure the logger
logger = logging.getLogger()
logger.setLevel(logging.INFO)

# Create a file handler to write logs to a file
log_file = f'potpourri/logs/log_{datetime_str}.log'
file_handler = logging.FileHandler(log_file, mode='w')
file_handler.setLevel(logging.INFO)

# Create a console handler to display logs in the terminal/console
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)

# Define the log format
formatter = logging.Formatter('%(asctime)s %(levelname)-8s %(message)s', datefmt='%a, %d %b %Y %H:%M:%S')
file_handler.setFormatter(formatter)
console_handler.setFormatter(formatter)

# Add the handlers to the logger
logger.addHandler(file_handler)
logger.addHandler(console_handler)

# Log messages
logger.info("OATS log file")
logger.info("Program started")

oats_dir = os.path.dirname(os.path.realpath(__file__))
current_path = os.getcwd()
default_testcase = current_path + "/testcases/matpower/case24_ieee_rts.xlsx"

#----------------------------------------------------------------------
# DC Load flow
def dclf(tc='default',solver='ipopt',neos=True,out=0):
    """
    Solves DC load flow problem

    ARGUMENTS:
        **tc** (*.xlsx file)  - OATS test case. See OATS data format for details

        **solver** (str)  - name of a solver. Defualt is 'ipopt'

        **neos** (bool) - If True, the problem is solved using NEOS otherwise using a localy install solver.

        **out** (bool) - If True, the output is displayed on screen.
    """

    if tc == 'default':
        tc = default_testcase
    #options
    opt=({'neos':neos,\
    'solver':solver,'out':out})
    testcase = tc
    model ='DCLF'
    # ==log==
    logging.info("Solver selected: "+opt['solver'])
    logging.info("Testcase selected: "+testcase)
    logging.info("Model selected: "+model)
    runcase(testcase,model,opt)
    logging.info("Done!")
    
# AC Load flow
def aclf(tc='default',solver='ipopt',neos=True,out=0):
    """
    Solves AC load flow problem

    ARGUMENTS:
        **tc** (*.xlsx file)  - OATS test case. See OATS data format for details

        **solver** (str)  - name of a solver. Defualt is 'ipopt'

        **out** (bool) - If True, the output is displayed on screen.
    """
    if tc == 'default':
        tc = default_testcase
    opt=({'solver':solver,'out':out})
    testcase = tc
    model ='ACLF'
    # ==log==
    logging.info("Solver selected: "+opt['solver'])
    logging.info("Testcase selected: "+testcase)
    logging.info("Model selected: "+model)
    runcase(testcase,model,opt)
    logging.info("Done!")
    
# DC optimal power flow problem
def dcopf(tc='default',
          solver='ipopt', 
          print_output:bool = False,
          print_solver_output:bool = False,
          print_model:bool = False,
          objective:str = None):
    """
    Solves DC optimal power flow problem

    ARGUMENTS:
        **tc** (*.xlsx file)  - OATS test case. See OATS data format for details

        **solver** (str)  - name of a solver. Defualt is 'ipopt'
        
        **out** (bool) - If True, the output is displayed on screen.
    """

    if tc == 'default':
        tc = default_testcase
    #options
    opt=({'solver':solver,'print_output':print_output,'print_solver_output':print_solver_output,'print_model':print_model})
    testcase = tc
    if objective == None:
        model ='DCOPF'
    elif objective == "quadratic":
        model = 'DCOPF_quadratic'
    # ==log==
    logging.info("Solver selected: "+opt['solver'])
    logging.info("Testcase selected: "+testcase)
    logging.info("Model selected: "+model)
    runcase(testcase,model,opt)
    logging.info("Done!")

# AC optimal power flow problem
def acopf(tc:str = 'default', 
          solver:str = 'ipopt', 
          print_output:bool = False, 
          print_solver_output:bool = False,
          print_model:bool = False, 
          objective:str = None):
    """
    Solves AC optimal power flow problem

    ARGUMENTS:
        **tc** (*.xlsx file)  - OATS test case. See OATS data format for details

        **solver** (str)  - name of a solver. Defualt is 'ipopt'

        **out** (bool) - If True, the output is displayed on screen.
    """

    if tc == 'default':
        tc = default_testcase
    #options
    opt=({'solver':solver,'print_output':print_output,'print_solver_output':print_solver_output,'print_model':print_model})
    testcase = tc
    model = "ACOPF"
    if objective == "min_gen":
        model = "ACOPF_min_gen"
    elif objective == "max_cost":
        model = "ACOPF_max_cost"
    elif objective== "multiperiod":
        model = "ACOPF_multiperiod"    
    elif objective == "multiperiod_battery":
        model = "ACOPF_multiperiod_battery" ###burasi ne yapiyor
    elif objective == "multiperiod_battery_exact":
        model = "ACOPF_multiperiod_battery_exact"
    elif objective == "multiperiod_battery_simplified":
        model = "ACOPF_multiperiod_battery_simplified"
    elif objective == "multiperiod_battery_relaxed":
        model = "ACOPF_multiperiod_battery_relaxed"
    # ==log==
    logging.info("Solver selected: "+opt['solver'])
    logging.info("Testcase selected: "+testcase)
    logging.info("Model selected: "+model)
    runcase(testcase,model,opt)
    logging.info("Done!")

# security constrained optimal power flow problem
def scopf(tc:str = 'default', 
          solver:str = 'ipopt', 
          print_output:bool = False, 
          print_solver_output:bool = False,
          print_model:bool = False, 
          objective:str = None):
    """
    Solves security constrained optimal power flow problem

    ARGUMENTS:
        **tc** (*.xlsx file)  - OATS test case. See OATS data format for details

        **solver** (str)  - name of a solver. Defualt is 'ipopt'

        **out** (bool) - If True, the output is displayed on screen.
    """

    if tc == 'default':
        tc = default_testcase
    #options
    opt=({'solver':solver,'print_output':print_output,'print_solver_output':print_solver_output,'print_model':print_model})
    testcase = tc
    model ='SCOPF'
    # ==log==
    logging.info("Solver selected: "+opt['solver'])
    logging.info("Testcase selected: "+testcase)
    logging.info("Model selected: "+model)
    runcase(testcase,model,opt)
    logging.info("Done!")

# unit commitment problem
def uc(tc:str = 'default', 
          solver:str = 'ipopt', 
          print_output:bool = False, 
          print_solver_output:bool = False,
          print_model:bool = False, 
          objective:str = None):
    """
    Solves unit commitment problem

    ARGUMENTS:
        **tc** (*.xlsx file)  - OATS test case. See OATS data format for details

        **solver** (str)  - name of a solver. Defualt is 'ipopt'

        **out** (bool) - If True, the output is displayed on screen.
    """

    if tc == 'default':
        tc = default_testcase
    #options
    opt=({'solver':solver,'print_output':print_output,'print_solver_output':print_solver_output,'print_model':print_model})
    testcase = tc
    model ='UC'
    # ==log==
    logging.info("Solver selected: "+opt['solver'])
    logging.info("Testcase selected: "+testcase)
    logging.info("Model selected: "+model)
    runcase(testcase,model,opt)
    logging.info("Done!")

# user defined model
def model(model:str = 'ACOPF', 
          tc:str = 'default', 
          solver:str = 'ipopt', 
          print_output:bool = False, 
          print_solver_output:bool = False,
          print_model:bool = False, 
          objective:str = None):
    """
    Solves a user defined model

    ARGUMENTS:
        **model** (*.py file)  - A model file that describes the optimisation problem using PYOMO modelling language.

        **tc** (*.xlsx file)  - OATS test case. See OATS data format for details

        **solver** (str)  - name of a solver. Defualt is 'ipopt'

        **out** (bool) - If True, the output is displayed on screen.
    """

    if tc == 'default':
        tc = default_testcase
    #options
    opt=({'solver':solver,'print_output':print_output,'print_solver_output':print_solver_output,'print_model':print_model,'user_def_model':True})
    testcase = tc
    model =model
    # ==log==
    logging.info("Solver selected: "+opt['solver'])
    logging.info("Testcase selected: "+testcase)
    logging.info("Model selected: "+model)
    runcase(testcase,model,opt)
    logging.info("Done!")
>>>>>>> ma_lohse
