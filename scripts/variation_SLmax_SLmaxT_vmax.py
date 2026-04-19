import pandapower as pp
import simbench as sb
import matplotlib.pyplot as plt

from src.potpourri.models.HC_ACOPF import HC_ACOPF
from src.potpourri.plotting.plot_functions import plot_wind_hc_results, plot_qu_res

if __name__ == '__main__':
    # net = pp.networks.simple_four_bus_system()

    net = sb.get_simbench_net("1-HV-mixed--0-no_sw")

    hc = HC_ACOPF(net, SWmax=10000, SWmin=10, peGmax=1000000)
    hc.solve()
    hc.add_OPF()
    hc.solve(solver="mindtpy", mip_solver='gurobi')

    hc_t = HC_ACOPF(net, SWmax=10000, SWmin=10, peGmax=1000000)
    hc_t.solve()
    # hc_t.SLmaxT_data *= 2
    hc_t.add_OPF()
    hc_t.model.transf_lim1.deactivate()
    hc_t.model.transf_lim2.deactivate()
    hc_t.solve(solver="mindtpy", mip_solver='gurobi')

    hc_l = HC_ACOPF(net, SWmax=10000, SWmin=10, peGmax=1000000)
    hc_l.solve()
    # hc_l.SLmax_data *= 2
    hc_l.add_OPF()
    hc_l.model.line_lim_from.deactivate()
    hc_l.model.line_lim_to.deactivate()
    hc_l.solve(solver="mindtpy", mip_solver='gurobi')

    # hc_lt = HC_ACOPF(net, SWmax=10000, SWmin=10, peGmax=1000000)
    # hc_lt.solve()
    # hc_lt.SLmax_data *= 2
    # hc_lt.SLmaxT_data *= 2
    # hc_lt.add_OPF()
    # hc_lt.solve(solver="mindtpy", mip_solver='gurobi')

    hc_u = HC_ACOPF(net, SWmax=10000, SWmin=10, peGmax=1000000)
    hc_u.solve()
    hc_u.add_OPF()
    hc_u.model.v_constraint.deactivate()
    hc_u.solve(solver="mindtpy", mip_solver='gurobi')

    hc_q = HC_ACOPF(net, SWmax=10000, SWmin=10, peGmax=1000000)
    hc_q.solve()
    hc_q.add_OPF()
    hc_q.model.QW_max_constraint.deactivate()
    hc_q.model.QU_max_hc_constraint.deactivate()
    hc_q.solve(solver="mindtpy", mip_solver='gurobi')

    nets = [hc.net, hc_t.net, hc_l.net, hc_u.net, hc_q.net]
    for net in nets:
        plot_wind_hc_results([net])
        plot_qu_res(net)

    # array with objective values rounded to 2 decimals
    objectives = [round(hc.model.obj_hc(), 2), round(hc_t.model.obj_hc(), 2), round(hc_l.model.obj_hc(), 2),
                  round(hc_u.model.obj_hc(), 2), round(hc_q.model.obj_hc(), 2)]

    plt.bar(range(len(objectives)), objectives, tick_label=objectives)
