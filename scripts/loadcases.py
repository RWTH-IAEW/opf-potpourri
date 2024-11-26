import copy

import simbench as sb
import pyomo.environ as pe

from potpourri.models.ACOPF_base import ACOPF

if __name__ == '__main__':
    net = sb.get_simbench_net("1-HV-mixed--0-no_sw")

    case_keys = ['lW', 'hL']

    hcs = []
    obj = []

    for case in case_keys:
        net_case = copy.deepcopy(net)
        factors = net_case.loadcases.loc[case]

        net_case.load.p_mw *= factors['pload']
        net_case.load.q_mvar *= factors['qload']

        net_case.sgen.scaling[net_case.sgen.type == 'Wind'] = factors['Wind_p']
        net_case.sgen.scaling[net_case.sgen.type == 'PV'] = factors['PV_p']
        net_case.sgen.scaling[(net_case.sgen.type != 'Wind') & (net_case.sgen.type != 'PV')] = factors['RES_p']

        net_case.ext_grid.vm_pu = factors['Slack_vm']

        hc = ACOPF(net_case)
        hc.add_voltage_deviation_objective()
        hc.solve(solver='neos')

        hcs.append(copy.deepcopy(hc))
        obj.append(pe.value(hc.model.obj))

        for g in hc.model.sG:
            print(f"Gen {g}: {pe.value(hc.model.pG[g])}")
            print(f"Gen {g}: {pe.value(hc.model.qG[g])}")


        line_flow = pe.value(hc.model.pLfrom[0])