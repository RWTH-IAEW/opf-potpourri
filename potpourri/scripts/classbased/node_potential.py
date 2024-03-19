import pickle
import pyomo.environ as pe
import pandapower as pp
from potpourri.models.class_based.HC_ACOPF import HC_ACOPF
from potpourri.scripts.classbased.plot_functions import plot_wind_hc_results


def write_wind_node_potential_to_buses(opf, node_potentials):
    # Step 1: Initialize a new column in hc.net.bus named 'node_potential' with default values as NaN.
    opf.net.bus['node_potential_mw'] = float('NaN')

    # Step 2: Iterate over the node_potentials list
    for sgen_index, node_potential in node_potentials:
        # Step 3: Find the corresponding bus index in hc.net.sgen.bus
        bus_index = opf.static_generation_data.bus[sgen_index]
        # Step 4: Assign the node potential value to the 'node_potential' column of hc.net.bus
        opf.net.bus['node_potential_mw'].iloc[bus_index] = node_potential


def plot_node_potential_and_hc(net):
    wind_index = net.sgen[net.sgen.wind_hc].index
    sgen_wind_bus = net.sgen.bus[wind_index]
    sgen_wind_p = net.res_sgen.p_mw[wind_index]
    net.res_bus.loc[sgen_wind_bus, 'p_wind_mw'] = sgen_wind_p.values.astype(float)

    # create marker trace with marker size scaled according to wind generation
    hosting_capacity_trace = pp.plotting.create_weighted_marker_trace(net, "res_bus",
                                                                      elm_ids=net.res_bus.index[
                                                                          net.res_bus.p_wind_mw > 0],
                                                                      column_to_plot="p_wind_mw", marker_scaling=2,
                                                                      trace_name='Maximum hosting capacity',
                                                                      scale_legend_unit='MW')

    node_potential_trace = pp.plotting.create_weighted_marker_trace(net, "bus", elm_ids=net.bus.index[
        net.bus.node_potential.notna()],
                                                                    column_to_plot="node_potential_mw",
                                                                    marker_scaling=2,
                                                                    trace_name='Node potential', color='blue',
                                                                    scale_legend_unit='MW')

    pp.plotting.simple_plotly(net, bus_size=3, bus_color='black',
                              additional_traces=[node_potential_trace, hosting_capacity_trace], showlegend=True)



if __name__ == "__main__":
    net = pickle.load(
        open('C:\\Users\\f.lohse\PycharmProjects\potpourri\potpourri\data\simbench_hv_grid_with_potential.pkl', 'rb'))

    hc_node = HC_ACOPF(net)
    hc_node.solve()
    hc_node.add_OPF()

    hc_node.fix_vars('y', 0)
    hc_node.model.obj.deactivate()

    node_potentials = []

    for w in hc_node.model.WIND_HC:
        # calculate node potential for each node individually (each wind hc static generator)
        hc_node.model.y[w].unfix()
        hc_node.model.y[w] = 1

        def obj_node(model):
            return model.psG[w]

        hc_node.model.obj_node = pe.Objective(rule=obj_node, sense=pe.maximize)

        hc_node.solve()

        node_potentials.append((w, pe.value(hc_node.model.obj_node)))

        print(pe.value(hc_node.model.obj_node))

        hc_node.model.y[w].fix(0)
        hc_node.model.psG[w] = 0.
        hc_node.model.qsG[w] = 0.

    hc_node.model.obj_node.deactivate()
    hc_node.unfix_vars('y', 0)
    hc_node.model.obj.activate()
    hc_node.solve(solver='mindtpy')

    write_wind_node_potential_to_buses(hc_node, node_potentials)
    plot_node_potential_and_hc(hc_node.net)
