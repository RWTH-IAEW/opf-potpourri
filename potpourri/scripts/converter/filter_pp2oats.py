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

def convert_pp_net(net: pp.pandapowerNet, value_of_lost_load: int = 100000) -> None:
    
    # Create ppc from net
    ppc = pp.converter.to_ppc(net)
    
    #value of lost load
    VOLL = str(value_of_lost_load)

    zneNTCcol = ['interconnection_ID', 'from_zone', 'to_zone', 'TransferCapacityTo(MW)', 'TransferCapacityFr(MW)']
    znecol    = ['zone', 'reserve(MW)']
    
    dfwnd    = pd.DataFrame(columns=['busname','name','stat','PG','QG','PGLB','PGUB','QGLB','QGUB','VS','contingency','probability'])
    dfzne    = pd.DataFrame(columns=znecol)
    dfzneNTC = pd.DataFrame(columns=znecol)
    dfts     = pd.DataFrame()
    baseMVA = pd.DataFrame({'baseMVA': [100]})
    
    # Create bus
    dfbus    = pd.DataFrame(columns=['name','baseKV','type','zone','VM','VA','VNLB','VNUB','VELB','VEUB'])
    dfbus['name'] = net.bus.name
    dfbus['baseKV'] = net.bus.vn_kv
    # TODO
    # bugfix the type in ppc
    dfbus['type'] = ppc["bus"][:, 1][:len(dfbus)]
    dfbus['zone'] = net.bus.zone
    dfbus['VM'] = net.bus.vn_kv
    dfbus['VA'] = 0
    dfbus['VNLB'] = net.bus.min_vm_pu
    dfbus['VNUB'] = net.bus.max_vm_pu
    dfbus['VELB'] = net.bus.min_vm_pu
    dfbus['VEUB'] = net.bus.max_vm_pu
    
    dfdem    = pd.DataFrame(columns=['name','busname','real','reactive','stat','VOLL'])
    dfdem['name'] = 'D' + net.load.name
    dfdem['busname'] = net.load.bus
    dfdem['real'] = net.load.p_mw
    dfdem['reactive'] = net.load.q_mvar
    dfdem['stat'] = str(1)
    dfdem['VOLL'] = VOLL
    
    dfsht    = pd.DataFrame(columns=['busname','name','GL','BL','stat'])
    dfsht['busname'] = None
    dfsht['name'] = None    
    dfsht['GL'] = None
    dfsht['BL'] = None
    dfsht['stat'] = None
    
    dfbrn    = pd.DataFrame(columns=['name','from_busname','to_busname','stat','r','x','b','ShortTermRating','ContinousRating','angLB','angUB','contingency','probability'])
    dfbrn['name'] = 'L'+ net.line.name
    dfbrn['from_busname'] =  net.line.from_bus
    dfbrn['to_busname'] =  net.line.to_bus
    dfbrn['stat'] =  str(1)
    dfbrn['r'] = net.line.r_ohm_per_km*net.line.length_km
    dfbrn['x'] = net.line.x_ohm_per_km*net.line.length_km
    # TODO
    # Check if b == c???
    dfbrn['b'] = net.line.c_nf_per_km*net.line.length_km
    dfbrn['ShortTermRating'] = str(9999)
    dfbrn['ContinousRating'] = str(9999)
    dfbrn['angLB'] = -360
    dfbrn['angUB'] = 360
    dfbrn['contingency'] = str(1)
    dfbrn['probability'] = str(0.0001)
    
    dftrn    = pd.DataFrame(columns=['name','from_busname','to_busname','stat','type','r','x','ShortTermRating','ContinousRating','angLB','angUB','PhaseShift','TapRatio','TapLB','TapUB','contingency','probability'])
    dftrn['name'] = 'T'+net.trafo.name
    dftrn['from_busname'] = net.trafo.hv_bus
    dftrn['to_busname'] = net.trafo.lv_bus
    # TODO
    # Check meaning of stat
    dftrn['stat'] = str(1)
    dftrn['type'] = str(1)
    # TODO
    # Calculate from given formulas:
    # https://pandapower.readthedocs.io/en/v2.13.1/elements/trafo.html
    dummy_net_sn_mva = 1
    dummy_sn_mva = 10
    dftrn['r'] = (net.trafo.vkr_percent/100)*(dummy_net_sn_mva/dummy_sn_mva)
    dftrn['x'] = None
    dftrn['ShortTermRating'] = str(9999)
    dftrn['ContinousRating'] = str(9999)       
    # TODO
    # Check values 
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
    dfgen['name'] = 'G' + net.sgen.name
    # TODO
    # Check meaning of stat
    dfgen['stat'] = str(1)
    dfgen['type'] = str(1)
    dfgen['PG'] = net.sgen.p_mw
    dfgen['QG'] = net.sgen.q_mvar
    dfgen['PGLB'] = 0
    dfgen['PGUB'] = net.sgen.p_mw
    dfgen['QGLB'] = -net.sgen.q_mvar
    dfgen['QGUB'] = net.sgen.q_mvar
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
    
    with pd.ExcelWriter('test_simbench.xlsx') as writer:  
        dfbus.to_excel(writer, sheet_name = 'bus',index=False , header=True)
        dfdem.to_excel(writer, sheet_name = 'demand',index=False, header=True)
        dfbrn.to_excel(writer, sheet_name = 'branch',index=False, header=True)
        dftrn.to_excel(writer, sheet_name = 'transformer',index=False, header=True)
        dfwnd.to_excel(writer, sheet_name = 'wind',index=False, header=True)
        dfsht.to_excel(writer, sheet_name = 'shunt',index=False, header=True)
        dfzne.to_excel(writer, sheet_name = 'zone',index=False, header=True)
        dfzneNTC.to_excel(writer, sheet_name = 'zonalNTC',index=False, header=True)
        dfgen.to_excel(writer, sheet_name = 'generator',index=False, header=True)
        dfts.to_excel(writer, sheet_name = 'timeseries',index=False, header=True)
        baseMVA.to_excel(writer, sheet_name = 'baseMVA',index=False, header=True)
    
    return

if __name__ == "__main__":
    
    sb_code1 = "1-MV-rural--0-sw"  # rural MV grid of scenario 0 with full switchs
    net = sb.get_simbench_net(sb_code1)
    pp.runpp(net)
    
    convert_pp_net(net)
