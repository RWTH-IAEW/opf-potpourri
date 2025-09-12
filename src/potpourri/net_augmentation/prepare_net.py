import pandas as pd
import pandapower as pp

def apply_loadcase_to_sb_net(net, case):
    # apply loadcase 'lW' to net
    factors = net.loadcases.loc[case]
    net.load.p_mw *= factors['pload']
    net.load.q_mvar *= factors['qload']
    net.sgen.loc[net.sgen.type == 'Wind', 'scaling'] = factors['Wind_p']
    net.sgen.loc[net.sgen.type == 'PV', 'scaling'] = factors['PV_p']
    net.sgen.loc[(net.sgen.type != 'Wind') & (net.sgen.type != 'PV'), 'scaling'] = factors['RES_p']
    net.ext_grid.loc[:, 'vm_pu'] = factors['Slack_vm']

    return net


def add_regulatory_q_control_to_wind(net, variant):
    # add wind control variant to existing wind generators
    net.sgen['controllable'] = False
    net.sgen['controllable'][net.sgen.type == 'Wind'] = True
    net.sgen['p_inst_mw'] = net.sgen['p_mw']
    net.sgen['var_q'] = None
    net.sgen['var_q'][net.sgen.type == 'Wind'] = variant

    return net

def upgrade_pandapower_net(old_net):
    # Create a new empty pandapowerNet with the current version of pandapower
    new_net = pp.create_empty_network()

    # Iterate through each key in the old net and copy it to the new net
    for key, value in old_net.items():
            # If the attribute exists in the new net, copy the data
            if isinstance(value, pd.DataFrame):
                new_net[key] = value.copy()  # Copy DataFrame
            else:
                new_net[key] = value  # Copy other data

    return new_net
