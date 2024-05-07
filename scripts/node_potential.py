import pickle
import pyomo.environ as pe
import pandapower as pp
from tqdm import tqdm

from potpourri.models.HC_ACOPF import HC_ACOPF
from scripts.prepare_net import apply_loadcase_to_sb_net
from rwth_colors import colors

def write_wind_node_potential_to_buses(opf, node_potentials):
    # Step 1: Initialize a new column in hc.net.bus named 'node_potential' with default values as NaN.
    opf.net.bus['node_potential_p_mw'] = float('NaN')

    # Step 2: Iterate over the node_potentials list
    for sgen_index, node_potential in node_potentials:
        # Step 3: Find the corresponding bus index in hc.net.sgen.bus
        bus_index = opf.static_generation_data.bus[sgen_index]
        # Step 4: Assign the node potential value to the 'node_potential' column of hc.net.bus
        opf.net.bus['node_potential_p_mw'].iloc[bus_index] = node_potential


def plot_node_potential_and_hc(net, wind_pot=False):
    wind_index = net.sgen[net.sgen.wind_hc].index
    sgen_wind_bus = net.sgen.bus[wind_index]
    sgen_wind_p = net.sgen.p_mw[wind_index]
    net.res_bus.loc[sgen_wind_bus, 'p_wind_mw'] = sgen_wind_p.where(sgen_wind_p >= 1e-3, 0).values.astype(float)

    # create marker trace with marker size scaled according to wind generation
    hosting_capacity_trace = pp.plotting.create_weighted_marker_trace(net, "res_bus",
                                                                      elm_ids=net.res_bus.index[
                                                                          net.res_bus.p_wind_mw > 0],
                                                                      column_to_plot="p_wind_mw",
                                                                      trace_name='Maximum hosting capacity',
                                                                      scale_legend_unit='MW',
                                                                      color=colors['red'])

    node_potential_trace = pp.plotting.create_weighted_marker_trace(net, "bus", elm_ids=net.bus.index[
        net.bus.node_potential_p_mw.notna()],
                                                                    column_to_plot="node_potential_p_mw",
                                                                    trace_name='Node potential', color=colors['blue'],
                                                                    scale_legend_unit='MW')

    if wind_pot:
        wind_potential_trace = pp.plotting.create_weighted_marker_trace(net, "bus", elm_ids=net.bus.index[
            net.bus.windpot_p_mw.notna()],
                                                                        column_to_plot="windpot_p_mw",
                                                                        trace_name='Wind potential', color=colors['green'],
                                                                        scale_legend_unit='MW')
        fig = pp.plotting.simple_plotly(net, bus_size=3, bus_color='black',
                                  additional_traces=[wind_potential_trace, node_potential_trace, hosting_capacity_trace], showlegend=True, auto_open=False)
    else:
        fig = pp.plotting.simple_plotly(net, bus_size=3, bus_color='black',
                              additional_traces=[node_potential_trace, hosting_capacity_trace], showlegend=True, auto_open=False, ext_grid_color=colors['green'])
    fig.update_layout(
            {
                "paper_bgcolor": "rgba(0, 0, 0, 0)",
                "plot_bgcolor": "rgba(0, 0, 0, 0)",
            }
        )

    fig.data[-2].marker.color = colors['blue']
    fig.data[-1].marker.color = colors['red']
    fig.data[-1].marker.opacity = 1.0
    fig.data[-2].marker.opacity = 1.0

    fig.show(renderer="browser")
    return fig

if __name__ == "__main__":
    net = pickle.load(
        open('potpourri/data/windpot/sb_hv_grid_with_potential_3MW_230m.pkl', 'rb'))

    case = 'lW'
    net = apply_loadcase_to_sb_net(net, case)
    # add_regulatory_q_control_to_wind(net, 1)

    hc_node = HC_ACOPF(net)
    hc_node.solve()
    hc_node.add_OPF()
    hc_node.add_tap_changer_linear()

    # for w in hc_node.model.WINDc:
    #     hc_node.model.psG[w].fix()

    hc_node.fix_vars('y', 0)
    for w in hc_node.model.WIND_HC:
        hc_node.model.psG[w].fix(0.)
        hc_node.model.qsG[w].fix(0.)

    hc_node.model.obj.deactivate()

    def obj_node(model):
        return sum(model.psG[w] for w in model.WIND_HC)

    hc_node.model.obj_node = pe.Objective(rule=obj_node, sense=pe.maximize)

    node_potentials = []

    # calculate node potential for each node individually (each wind hc static generator)
    for w in tqdm(hc_node.model.WIND_HC):
        # if wind potential is 0, node potential is 0, skip calculation
        if hc_node.model.pWmax[w]() == 0:
            node_potentials.append((w, 0.))
            continue

        hc_node.model.y[w] = 1
        hc_node.model.psG[w].unfix()
        hc_node.model.qsG[w].unfix()

        hc_node.solve()

        node_potentials.append((w, pe.value(hc_node.model.obj_node)))

        hc_node.model.y[w] = 0
        hc_node.model.psG[w].fix(0.)
        hc_node.model.qsG[w].fix(0.)

    hc_node.model.obj_node.deactivate()
    hc_node.unfix_vars('y', 0)
    for w in hc_node.model.WIND_HC:
        hc_node.model.psG[w].unfix()
        hc_node.model.qsG[w].unfix()

    # calculate maximum hosting capacity
    hc_node.change_vals('SWmin', 10)
    hc_node.model.obj.activate()
    hc_node.solve(solver='mindtpy')

    write_wind_node_potential_to_buses(hc_node, node_potentials)
    plot_node_potential_and_hc(hc_node.net)

    results_dir= '../potpourri/results/node_potential'