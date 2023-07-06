#==================================================================
# A Python script that reads Simbench/Pandapower test case and converts into OATS test case
# ---Author---
# S. Kortmann,
# s.kortmann@iaew.rwth-aachen.de
# POTPOURRI
# Copyright (c) 2023 by S.Kortmann, Aachen, Germany
# POTPOURRI is distributed under the GNU GENERAL PUBLIC LICENSE v3. (see LICENSE file for details).
#
# TODO:
# https://oats.readthedocs.io/en/latest/dataformat.html
# Use information from OATS and fill it with Pandapower network structure
# Current problems:
# Not all parameters can be retrieved from Pandapower
# Potential solution:
# Take ppc into account
#==================================================================

import pandas as pd
import pandapower as pp
import simbench as sb
import random
import math

def convert_pp_net(net: pp.pandapowerNet, sb_code1:str, value_of_lost_load: int = 100000, time_period: int = 96) -> None:
    
    # Create ppc from net
    ppc = pp.converter.to_ppc(net)
    
    #value of lost load
    VOLL = str(value_of_lost_load)

    zneNTCcol = ['interconnection_ID', 'from_zone', 'to_zone', 'TransferCapacityTo(MW)', 'TransferCapacityFr(MW)']
    znecol    = ['zone', 'reserve(MW)']
    
    dfwnd    = pd.DataFrame(columns=['busname','name','stat','PG','QG','PGLB','PGUB','QGLB','QGUB','VS','contingency','probability'])
    dfzne    = pd.DataFrame(columns=znecol)
    dfzneNTC = pd.DataFrame(columns=znecol)
    baseMVA = pd.DataFrame({'baseMVA': [100]})
    
    # Create bus
    dfbus    = pd.DataFrame(columns=['name','baseKV','type','zone','VM','VA','VNLB','VNUB','VELB','VEUB'])
    dfbus['name'] = net.bus.index
    dfbus['baseKV'] = net.bus.vn_kv
    # TODO
    # bugfix the type in ppc
    dfbus['type'] = ppc["bus"][:, 1][:len(dfbus)]
    dfbus['zone'] = net.bus.zone
    dfbus['VM'] = 1.0 #net.res_bus.vm_pu
    dfbus['VA'] = 0 #net.res_bus.va_degree
    dfbus['VNLB'] = net.bus.min_vm_pu
    dfbus['VNUB'] = net.bus.max_vm_pu
    dfbus['VELB'] = net.bus.min_vm_pu
    dfbus['VEUB'] = net.bus.max_vm_pu
    
    dfdem    = pd.DataFrame(columns=['name','busname','real','reactive','stat','VOLL'])
    dfdem['name'] = net.load.name.str.replace(" ", "_")
    dfdem['busname'] = net.load.bus
    dfdem['real'] = net.load.p_mw * 1000
    dfdem['reactive'] = net.load.q_mvar * 1000
    dfdem['stat'] = str(1)
    dfdem['VOLL'] = VOLL
    
    dfsht    = pd.DataFrame(columns=['busname','name','GL','BL','stat'])
    dfsht['busname'] = None
    dfsht['name'] = None    
    dfsht['GL'] = None
    dfsht['BL'] = None
    dfsht['stat'] = None
    
    dfbrn    = pd.DataFrame(columns=['name','from_busname','to_busname','stat','r','x','b','ShortTermRating','ContinousRating','angLB','angUB','contingency','probability'])
    dfbrn['name'] = net.line.name.str.replace(" ", "_")
    dfbrn['from_busname'] =  net.line.from_bus
    dfbrn['to_busname'] =  net.line.to_bus
    dfbrn['stat'] =  str(1)
    dfbrn['r'] = net.line.r_ohm_per_km*net.line.length_km
    dfbrn['x'] = net.line.x_ohm_per_km*net.line.length_km
        # TODO
        # update 30.06.2023 : no idea about b :-D 
        # Check if b == c???
    dfbrn['b'] = 0.1 #-(net.line.x_ohm_per_km*net.line.length_km)/((net.line.x_ohm_per_km*net.line.length_km)**2 + (net.line.r_ohm_per_km*net.line.length_km)**2)
    dfbrn['ShortTermRating'] = str(9999)
    dfbrn['ContinousRating'] = str(9999)
    dfbrn['angLB'] = -360
    dfbrn['angUB'] = 360
    dfbrn['contingency'] = str(1)
    dfbrn['probability'] = str(0.0001)
    
    dftrn    = pd.DataFrame(columns=['name','from_busname','to_busname','stat','type','r','x','ShortTermRating','ContinousRating','angLB','angUB','PhaseShift','TapRatio','TapLB','TapUB','contingency','probability'])
    dftrn['name'] = net.trafo.name.str.replace(" ", "_")
    dftrn['from_busname'] = net.trafo.hv_bus
    dftrn['to_busname'] = net.trafo.lv_bus
    # TODO
    # Check meaning of stat
    dftrn['stat'] = str(1)
    dftrn['type'] = str(1)
    # TODO
    # Calculate from given formulas:
    # https://pandapower.readthedocs.io/en/v2.13.1/elements/trafo.html
    r_k = net.trafo.vkr_percent/100*(net.sn_mva/net.trafo.sn_mva)
    dftrn['r'] = r_k # 0.1
    z_k = net.trafo.vk_percent/100*(net.sn_mva/net.trafo.sn_mva)
    dftrn['x'] = 0.1 #math.sqrt(z_k**2-r_k**2)
    dftrn['ShortTermRating'] = str(9999)
    dftrn['ContinousRating'] = str(9999)        
    dftrn['angLB'] = -360
    dftrn['angUB'] = 360
    dftrn['PhaseShift'] = 0
    dftrn['TapRatio'] = 0.98
    dftrn['TapLB'] = str(0.98*(1-0.05))
    dftrn['TapUB'] = str(0.98*(1+0.05))
    dftrn['contingency'] = str(1)
    dftrn['probability'] = str(0.0001)
                    
    column_names = ['gen', 'start', 'shut', 'c2', 'c1', 'c0']
    costdat = pd.DataFrame(columns=column_names)
    costdat['gen'] = net.sgen.name
    costdat['start'] = str(1+ random.randint(0,100))
    costdat['shut'] = str(1 + random.randint(0,100))
    costdat['c2'] = str(1 + random.randint(0,100))
    costdat['c1'] = str(1+ random.randint(0,100))
    costdat['c0'] = str(1+ random.randint(0,100))
    
    dfgen    = pd.DataFrame(columns=['busname','name','stat','type','PG','QG','PGLB','PGUB','QGLB','QGUB','VS','RampDown(MW/hr)','RampUp(MW/hr)','MinDownTime(hr)','MinUpTime(hr)','FuelType','contingency','probability','startup','shutdown','costc2','costc1','costc0'])
    dfgen['busname'] = net.sgen.bus
    dfgen['name'] = net.sgen.name.str.replace(" ", "_")
        # TODO
         # Check meaning of stat
    dfgen['stat'] = str(1)
    dfgen['type'] = str(1)
    dfgen['PG'] = net.sgen.p_mw * 1000
    dfgen['QG'] = net.sgen.q_mvar * 1000
    dfgen['PGLB'] = 0
    dfgen['PGUB'] = net.sgen.p_mw * 1000
    dfgen['QGLB'] = -net.sgen.q_mvar * 1000
    dfgen['QGUB'] = net.sgen.q_mvar * 1000
    dfgen['VS'] = str(1)
        # TODO
        # Look up ramping times
    dfgen['RampDown(MW/hr)'] = 100
    dfgen['RampUp(MW/hr)'] = 100
    dfgen['MinDownTime(hr)'] = str(1)
    dfgen['MinUpTime(hr)'] = str(1)
    dfgen['FuelType'] = 'NA'
    dfgen['contingency'] = str(0)
    dfgen['probability'] = str(0.0001)
    dfgen['startup'] = costdat['start']
    dfgen['shutdown'] = costdat['shut']
    dfgen['costc2'] = costdat['c2']
    dfgen['costc1'] = costdat['c1']
    dfgen['costc0'] = costdat['c0']
    
    dfts     = pd.DataFrame(columns=["timeperiod"])
    
    # calculate absolute profiles
    profiles = sb.get_absolute_values(net, profiles_instead_of_study_cases=True)
    
    for (i, name) in enumerate(net.load.name):
        new_name = name.replace(" ", "_")
        dfts[f"{new_name}"] = profiles[("load", "p_mw")][i] * 1000
        
    dfts["timeperiod"] = list(range(len(profiles[("load", "p_mw")])))
    
    # Create the MultiIndex
    new_columns = [("Demand", column) if column != "timeperiod" else (" ","timeperiod") for column in dfts.columns]
    dfts.columns = pd.MultiIndex.from_tuples(new_columns)
    
    dfts[(" ","timeperiod")] = list(range(len(profiles[("load", "p_mw")])))
    
    #reduce for first 24 time steps
    #dfts = dfts.resample("1H").sum()
    dfts = dfts.head(time_period)
    
    print(dfts.columns)
    
    dfbat     = pd.DataFrame(columns=["name", "busname", "PCMAX", "PCMIN", "PDMAX", "PDMIN", "nchar", "ndis", "EMAX", "EMIN", "E0", "c3"])
    dfbat["name"] = net.load.name.str.replace(" ", "_") + "_bat"
    dfbat["busname"] = net.load.index
    dfbat["PCMAX"] = 20
    dfbat["PCMIN"] = 0
    dfbat["PDMAX"] = 20
    dfbat["PDMIN"] = 0
    dfbat["nchar"] = 0.8
    dfbat["ndis"] = 0.9
    dfbat["EMAX"] = 200
    dfbat["EMIN"] = 0
    dfbat["E0"] = 0
    dfbat["c3"] = 1

    # Define the MultiIndex
    multi_index = [(" ", "timeperiod"), ("busname", "1"), ("busname", "2"), ("PS", "1"), ("PS", "2"),
                ("QS", "1"), ("QS", "2"), ("c4", "1"), ("c4", "2")]

    # Initialize the DataFrame with the MultiIndex
    dfsolar = pd.DataFrame(columns=multi_index)
    
    dfsolar.columns = pd.MultiIndex.from_tuples(multi_index)

    dfsolar[(" ","timeperiod")] = list(range(time_period))
    dfsolar[("busname","1")] = 1
    dfsolar[("busname","2")] = 2
    #T_----odo
    dfsolar[("PS","1")] = 0
    dfsolar[("PS","2")] = 0
    dfsolar[("QS","1")] = 0
    dfsolar[("QS","2")] = 0
    dfsolar[("c4","1")] = 11
    dfsolar[("c4","2")] = 12

    name = sb_code1 + '.xlsx'
     
    with pd.ExcelWriter(name) as writer:  
        dfbus.to_excel(writer, sheet_name = 'bus',index=False , header=True)
        dfdem.to_excel(writer, sheet_name = 'demand',index=False, header=True)
        dfbrn.to_excel(writer, sheet_name = 'branch',index=False, header=True)
        dftrn.to_excel(writer, sheet_name = 'transformer',index=False, header=True)
        dfwnd.to_excel(writer, sheet_name = 'wind',index=False, header=True)
        dfsht.to_excel(writer, sheet_name = 'shunt',index=False, header=True)
        dfzne.to_excel(writer, sheet_name = 'zone',index=False, header=True)
        dfzneNTC.to_excel(writer, sheet_name = 'zonalNTC',index=False, header=True)
        dfgen.to_excel(writer, sheet_name = 'generator',index=False, header=True)
        dfts.to_excel(writer, sheet_name = 'timeseries',index=True, header=True)
        dfbat.to_excel(writer, sheet_name = 'Battery',index=False, header = True)
        baseMVA.to_excel(writer, sheet_name = 'baseMVA',index=False, header=True)
        dfsolar.to_excel(writer, sheet_name = 'solar',index=True, header = True)
    
    return

if __name__ == "__main__":
    
    sb_code1 = "1-LV-rural1--1-no_sw" #"1-LV-semiurb5--0-no_sw"  # rural MV grid of scenario 0 with full switchs
    net = sb.get_simbench_net(sb_code1)
    pp.runpp(net)
    
    convert_pp_net(net, sb_code1)
