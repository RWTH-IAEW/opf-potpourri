

def apply_loadcase_to_sb_net(net, case):
    # apply loadcase 'lW' to net
    factors = net.loadcases.loc[case]
    net.load.p_mw *= factors['pload']
    net.load.q_mvar *= factors['qload']
    net.sgen.scaling[net.sgen.type == 'Wind'] = factors['Wind_p']
    net.sgen.scaling[net.sgen.type == 'PV'] = factors['PV_p']
    net.sgen.scaling[(net.sgen.type != 'Wind') & (net.sgen.type != 'Solar')] = factors['RES_p']
    net.ext_grid.vm_pu = factors['Slack_vm']

    return net


def add_regulatory_q_control_to_wind(net, variant):
    # add wind control variant to existing wind generators
    net.sgen['controllable'] = False
    net.sgen['controllable'][net.sgen.type == 'Wind'] = True
    net.sgen['p_inst_mw'] = net.sgen['p_mw']
    net.sgen['var_q'] = None
    net.sgen['var_q'][net.sgen.type == 'Wind'] = variant

    return net