#==================================================================
# runcase.py
# This is a second level OATS script, which is called by the runfile.
# This script receives model,testcase, options as an input to run simulation
# ---Author---
# W. Bukhsh,
# wbukhsh@gmail.com
# April 2018
# OATS
# Copyright (c) 2017 by W. Bukhsh, Glasgow, Scotland
# OATS is distributed under the GNU GENERAL PUBLIC LICENSE v3. (see LICENSE file for details).
#==================================================================

#===============Import===============
from __future__ import division
import os
from pyomo.environ import *
from pyomo.opt import SolverFactory
from pyomo.opt import SolverStatus, TerminationCondition
import logging
from .selecttestcase import selecttestcase
from .printdata import printdata
from .printoutput import printoutput
import imp
#====================================

def runcase(testcase,mod,opt=None):
    
    oats_dir = os.path.dirname(os.path.realpath(__file__))
    if 'user_def_model' in opt:
        modelf = imp.load_source(mod, mod+'.py')
        model = modelf.model
    else:
        try:
            modelf = imp.load_source(mod, 'potpourri/models/' +mod+ '.py')
            model = modelf.model
            logging.info("Given model file found and selected from the models library")
        except Exception:
            logging.error("Given model file not found in the 'models' library", exc_info=False)
            raise
    try:
        ptc = selecttestcase(testcase) #read test case
        logging.info("Given testcase file found and selected from the testcase library")
    except Exception:
        logging.error("Given testcase  not found in the 'testcases' library", exc_info=False)
        raise

    datfile = 'datafile.dat'
    r = printdata(datfile,ptc,mod,opt)
    r.reducedata()
    r.printheader()
    if 'UC' in mod:
        r.printUCdat()
    else:
        r.printkeysets()
        r.printnetwork()
        #'OPF' or 'LF'
        if 'OPF' in mod:
            r.printOPF()
        elif 'LF' in mod:
            r.printLF()
        #'AC' or 'DC'
        if 'AC' in mod:
            r.printAC()
        elif 'DC' or 'SC' in mod:
            r.printDC()

        if mod=='ACLF':
            r.printACLF()
        elif mod=='DCOPF':
            r.printDCOPF()
        elif 'ACOPF' in mod:
            r.printACOPF()
            if 'multiperiod' in mod:
                r.printACOPF_multiperiod()
                if 'battery' in mod:
                    r.printACOPF_multiperiod_battery()
        elif mod=='SCOPF' or mod=='SCOPF_BM':
            r.printSCdat()
        if mod=='ACOPF_BM':
            r.printBM()
        if mod=='DCOPF_BM' or mod=='SCOPF_BM':
            r.printDCBM()



    ###############Solver settings####################
    optimise = SolverFactory(opt['solver'])
    #opt.options['mipgap'] = 0.1
    #################################################

    ############Solve###################
    
    instance = model.create_instance(datfile)
    instance.dual = Suffix(direction=Suffix.IMPORT)
    print("_______________________________________________")
    #instance.pprint(filename='model.lp')
    #filename = os.path.join(os.path.dirname(__file__), 'model.lp')
    #instance.write(filename, io_options={'symbolic_solver_labels': True})
    with open("/home/bengisu/potpourri/testcases/matpower/doc", 'w') as output_file:    
        instance.pprint(output_file)

    
    if not (opt['print_solver_output']):
        results = optimise.solve(instance,tee=False)
    else:
        results = optimise.solve(instance,tee=True)
    instance.solutions.load_from(results)
    
    # ##################################
    #
    # ############Output###################
    #print('_______________________________________Test_________________________________________________________')
    #print(instance.pChar[g[1],t])
    #for t in instance.T:
    #   for g in instance.Bbs:
    #       print(instance.e[g[1],t].value*instance.baseMVA)
    #(print(g) for g in instance.Bbs)

    print('________________________________________________________________________________________')
    
    o = printoutput(results, instance,mod)
    
    if (opt['print_output']):
        o.solutionstatus()
        #if (opt['print_model']):
        #instance.display()
        o.printsummary()
        
        if 'multiperiod' in mod:
            print("_________________________________Multiperiod_________________________________")
            o.printACOPF_multiperiod(testcase)
            #if 'battery' in mod:
            #    print("---------------battery---------------")
            #    o.print_Battery(testcase)    
        else:
            o.printoutputxls(testcase)

    else:
        o.greet()
        o.solutionstatus()
        
    if 'UC' in mod:
        o.printUC()
    

    # elif 'ACOPF_multiperiod' in mod:
    #     o.printACOPF_multiperiod()