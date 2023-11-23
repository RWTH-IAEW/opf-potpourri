import pandapower as pp
from potpourri.models.class_based.ACOPF_base import ACOPF
from potpourri.models.class_based.DCOPF import DCOPF
from potpourri.models.class_based.DCLF import DCLF
from potpourri.models.class_based.ACLF import ACLF
from potpourri.models.class_based.basemodel import Basemodel
from potpourri.scripts.printoutput import printoutput
from potpourri.models.class_based.HC_OPF import HCOPF

def create_testnet():
    net = pp.create_empty_network()

    pp.create_buses(net, 3, name=["b1", "b2", "b3"], vn_kv=110, max_vm_pu=1.118, min_vm_pu=0.9)
    pp.create_lines(net, [0, 0, 1], [1, 2, 2], name=["l1", "l2", "l3"], length_km=[2, 5, 7],
                    std_type='305-AL1/39-ST1A 110.0')
    pp.create_sgens(net, [2, 1, 0], name=["gen2", "gen1", "gen0"], p_mw=[3, 10, 15])
    pp.create_ext_grid(net, 0)
    pp.create_loads(net, [0, 1, 2], name=["load1", "load2", "load3"], p_mw=[2, 12, 16])
    # pp.create_sgen(net, 0, p_mw=1, controllable=True)
    # pp.create_sgen(net, 0, p_mw=1, controllable=True, max_p_mw=10)
    # pp.create_load(net, 1, p_mw=1, controllable=True)
    # pp.create_load(net, 2, p_mw=2, controllable=True, max_p_mw=5)

    pp.create_sgens(net, net.bus.index, p_mw=0, wind_hc=True)

    pp.create_shunt(net, 2, p_mw=3, q_mvar=5)
    return net

if __name__ == '__main__':
    net = create_testnet()

    dcopf = DCOPF(net)
    dclf = DCLF(net)
    aclf = ACLF(net)
#    acopf = ACOPF(net)
    base = Basemodel(net)
    hc = HCOPF(net)

    models = [dclf, aclf, dcopf, hc]
    o = []

    for sp in models:
        sp.solve()
        o.append(printoutput(sp.results, sp.model, sp.model.name))
        sp.model.pG.pprint()

    # wind = []
    # for w in sp_hc.model.W:
    #     wind.append({'bus': w, 'p_pu': pe.value(sp_hc.model.pW[w]), 'q_pu': pe.value(sp_hc.model.qW[w])})
    # wind = pd.DataFrame(wind)
    # pp.create_sgens(net, wind.bus, p_mw=wind.p_pu * sp_hc.model.baseMVA, q_mvar=wind.q_pu * sp_hc.model.baseMVA)
    # pp.runpp(net)
