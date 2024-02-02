from potpourri.models.class_based.HC_ACOPF import HC_ACOPF

import pandapower as pp
import pickle

def add_scenario_load_sgen_to_net(net, busses):
    net.load.in_service = False
    net.sgen.in_service = False

    pp.create_loads(net, busses, p_mw=0, name='load_scenario')
    pp.create_sgens(net, busses, p_mw=0, name='sgen_scenario')


if __name__ == '__main__':
    with open("../../data/scenarios/1-HV-mixed--0-no_sw_scenario_types.pkl", "rb") as f:
        net = pickle.load(f)

    with open("../../data/scenarios/1-HV-mixed_loadcases_assumption_scenarios_with_q.pkl", "rb") as f:
        scenarios = pickle.load(f)

    add_scenario_load_sgen_to_net(net, scenarios['load_sgen_busses'])

    net.load.p_mw[net.load.name == 'load_scenario'] = scenarios['load'][0][:, 0]
    net.load.q_mvar[net.load.name == 'load_scenario'] = scenarios['load_q'][0][:, 0]

    net.sgen.p_mw[net.sgen.name == 'sgen_scenario'] = scenarios['sgen'][0][:, 0]
    net.sgen.q_mvar[net.sgen.name == 'sgen_scenario'] = scenarios['sgen_q'][0][:, 0]

    hc = HC_ACOPF(net)
    hc.solve()
    hc.add_OPF()


    # for i in range(len(scenarios['load'][0][0])):
    #     net.load.p_mw[net.load.name == 'load_scenario'] = scenarios['load'][0][:, i]
    #     net.load.q_mvar[net.load.name == 'load_scenario'] = scenarios['load_q'][0][:, i]
    #
    #     net.sgen.p_mw[net.sgen.name == 'sgen_scenario'] = scenarios['sgen'][0][:, i]
    #     net.sgen.q_mvar[net.sgen.name == 'sgen_scenario'] = scenarios['sgen_q'][0][:, i]

