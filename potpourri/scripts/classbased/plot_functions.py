import pandapower as pp
import matplotlib.pyplot as plt
import numpy as np


def plot_wind_hc_results(nets):
    # marker_trace = pp.plotting.create_weighted_marker_trace(net, "sgen", elm_ids=net.sgen.index[net.sgen.wind_hc],
    #                                                         column_to_plot="p_mw", marker_scaling=2, color="green")
    # pp.plotting.simple_plotly(net, bus_size=0.25, additional_traces=[marker_trace])

    # get maximum wind generation from all wind generator from all nets
    max_wind_gen = max([net.sgen[net.sgen.wind_hc].p_mw.max() for net in nets])

    for net in nets:
        sgen_wind = net.sgen[net.sgen.wind_hc]
        net.res_bus.loc[sgen_wind.bus, 'p_wind_mw'] = sgen_wind['p_mw'].values

        # create marker trace with marker size scaled according to wind generation
        marker_trace = pp.plotting.create_weighted_marker_trace(net, "res_bus",
                                                                elm_ids=net.res_bus.index[net.res_bus.p_wind_mw > 0],
                                                                column_to_plot="p_wind_mw", marker_scaling=2,
                                                                trace_name='P_wind')

        # create bus trace with color map for wind generation
        bt = pp.plotting.create_bus_trace(net, buses=net.res_bus.index[net.res_bus.p_wind_mw > 0], cmap='ylorrd',
                                          colormap_column='p_wind_mw', cbar_title='Wind generation [MW]',
                                          cmax=max_wind_gen)

        # add marker size information to bus trace
        bt[0]['marker']['size'] = marker_trace['marker']['size']
        bt[0]['marker']['sizemode'] = marker_trace['marker']['sizemode']
        bt[0]['text'] = marker_trace['text']

        pp.plotting.simple_plotly(net, bus_size=3, bus_color='black', additional_traces=bt, showlegend=False)


def _create_sgen_load_trace(net):
    load_trace = pp.plotting.create_weighted_marker_trace(net, "load", column_to_plot="p_mw", marker_scaling=2,
                                                          color='green')

    if 'wind_hc' in net.sgen:
        sgen_ids = net.sgen.index[net.sgen.wind_hc == False]
        sgen_trace = pp.plotting.create_weighted_marker_trace(net, "sgen", elm_ids=sgen_ids, column_to_plot="p_mw",
                                                              marker_scaling=2, color="blue")
    else:
        sgen_trace = pp.plotting.create_weighted_marker_trace(net, "sgen", column_to_plot="p_mw",
                                                              marker_scaling=2, color='blue')

    return sgen_trace, load_trace


def plot_wind_hc_sgens_loads(net):
    sgen_trace, load_trace = _create_sgen_load_trace(net)

    sgen_wind = net.sgen[net.sgen.wind_hc]
    net.res_bus.loc[sgen_wind.bus, 'p_wind_mw'] = sgen_wind['p_mw'].values
    wind_trace = pp.plotting.create_weighted_marker_trace(net, "res_bus",
                                                          elm_ids=net.res_bus.index[net.res_bus.p_wind_mw > 0],
                                                          column_to_plot="p_wind_mw", marker_scaling=2,
                                                          trace_name='P_wind')

    pp.plotting.simple_plotly(net, bus_size=3, additional_traces=[wind_trace, sgen_trace, load_trace])


def plot_sgen_load(net):
    sgen_trace, load_trace = _create_sgen_load_trace(net)
    pp.plotting.simple_plotly(net, bus_size=3, additional_traces=[sgen_trace, load_trace])


def plot_pq_gridcodes():
    plt.style.use('rwth-word')

    clrs = ['#00549F', '#000000', '#E30066', '#FFED00', '#006165',
            '#0098A1', '#57AB27', '#BDCD00', '#F6A800', '#CC071E',
            '#A11035', '#612158', '#7A6FAC']

    fig, ax = plt.subplots()

    # Variant 1
    ax.vlines(1., -0.23, 0.48, color=clrs[8])  # Pmax
    ax.vlines(0.1, -0.1, 0.1, color=clrs[8])  # Pmin
    ax.plot([0.1, 0.2, 1], [-0.1, -0.23, -0.23], color=clrs[8], label='Variante 1')  # Qmin(P)
    ax.plot([0.1, 0.2, 1], [0.1, 0.48, 0.48], color=clrs[8])  # Qmax(P)

    # Variant 2
    ax.vlines(1., -0.33, 0.41, color=clrs[7])  # Pmax
    ax.vlines(0.1, -0.1, 0.1, color=clrs[7])  # Pmin
    ax.plot([0.1, 0.2, 1], [-0.1, -0.33, -0.33], color=clrs[7], label='Variante 2')  # Qmin(P)
    ax.plot([0.1, 0.2, 1], [0.1, 0.41, 0.41], color=clrs[7])  # Qmax(P)

    # Variant 3
    ax.vlines(1., -0.41, 0.33, color=clrs[0])  # Pmax
    ax.vlines(0.1, -0.1, 0.1, color=clrs[0])  # Pmin
    ax.plot([0.1, 0.2, 1], [-0.1, -0.41, -0.41], color=clrs[0], label='Variante 3')  # Qmin(P)
    ax.plot([0.1, 0.2, 1], [0.1, 0.33, 0.33], color=clrs[0])  # Qmax(P)

    plt.xlabel('$P_{mom} / P_{inst}$')
    plt.ylabel('$Q / P_{inst}$')
    plt.yticks(np.arange(-0.5, 0.6, 0.1))
    plt.ylim(-0.5, 0.5)
    plt.legend(loc='center left', bbox_to_anchor=(0.6, 0.5))
    plt.show()


def plot_qu_gridcodes():
    plt.style.use('rwth-word')

    clrs = ['#00549F', '#000000', '#E30066', '#FFED00', '#006165',
            '#0098A1', '#57AB27', '#BDCD00', '#F6A800', '#CC071E',
            '#A11035', '#612158', '#7A6FAC']

    fig, ax = plt.subplots()

    # Variant 1
    ax.plot([96, 120, 127], [0.48, 0.48, -0.23], color=clrs[8], label='Variante 1')  # Qmin(P)
    ax.plot([96, 103, 127], [0.48, -0.23, -0.23], color=clrs[8])  # Qmax(P)

    # Variant 2
    ax.plot([96, 120, 127], [0.41, 0.41, -0.33], color=clrs[7], label='Variante 2')  # Qmin(P)
    ax.plot([96, 103, 127], [0.41, -0.33, -0.33], color=clrs[7])  # Qmax(P)

    # Variant 3
    ax.plot([96, 120, 127], [0.33, 0.33, -0.41], color=clrs[0], label='Variante 3')  # Qmin(P)
    ax.plot([96, 103, 127], [0.33, -0.41, -0.41], color=clrs[0])  # Qmax(P)

    plt.xlabel('$U_{b} [kV]$')
    plt.ylabel('$Q / P_{inst}$')
    plt.yticks(np.arange(-0.5, 0.6, 0.1))
    plt.ylim(-0.5, 0.5)
    plt.legend(loc='center')
    plt.show()


def plot_qu_res(nets, labels=None):
    # values for variant 1
    m_qu = (0.48 + 0.23) / (96 - 103) * 110
    qu_min = -m_qu * 96 / 110 + 0.48
    qu_max = -m_qu * 120 / 110 + 0.48
    ymax = m_qu * 1.1 + qu_max
    ymin = m_qu * 0.9 + qu_min

    plt.style.use('rwth-word')
    clrs = ['#00549F', '#E30066', '#BDCD00', '#000000', '#FFED00', '#006165',
            '#0098A1', '#57AB27', '#F6A800', '#CC071E',
            '#A11035', '#612158', '#7A6FAC']
    marker = ['o', '*', '+', '.', ',']

    fig, ax = plt.subplots()
    ax.fill_between([0.9, 103 / 110, 120 / 110, 1.1], [0.48, 0.48, 0.48, ymax], [ymin, -0.23, -0.23, -0.23],
                    color=clrs[8], alpha=0.2)

    ax.plot([0.9, 120 / 110, 1.1, 1.1], [0.48, 0.48, ymax, -0.23], color=clrs[8])  # Qmax(P)
    ax.plot([0.9, 0.9, 103 / 110, 1.1], [0.48, ymin, -0.23, -0.23], color=clrs[8])  # Qmin(P)

    plt.xticks(np.arange(0.9, 1.1, 0.05))
    plt.xlabel('$U_{b}$ [p.u.]')
    plt.ylabel('$Q_{w} / P_{w}$')

    for i, net in enumerate(nets):
        line = ax.scatter(net.res_bus.vm_pu[net.sgen.bus[net.sgen.wind_hc]],
                          net.res_sgen.q_mvar[net.sgen.wind_hc].values / net.res_sgen.p_mw[net.sgen.wind_hc].values,
                          color=clrs[i], marker=marker[i])
        if labels:
            line.set_label(labels[i])

    if labels:
        plt.legend()


def plot_all_pG(hcs):
    # Initialize an empty dictionary to store the values
    values = {}
    # Exclude list
    exclude_w = [105, 106, 110]
    # Iterate over all hcs
    for hc in hcs:
        # Iterate over all w in hc.model.WIND
        for w in hc.model.WIND:
            # Skip if w is in the exclude list
            if w in exclude_w:
                continue
            # If w is not in the dictionary, add it with an empty list as value
            if w not in values:
                values[w] = []
            # Append hc.model.pG[w] to the list corresponding to w
            values[w].append(hc.model.pG[w].value)
    # Create a bar plot for each w
    plot_index = 0
    plot_w = []
    for i, (w, vals) in enumerate(values.items()):
        # Skip if all values are smaller than 1e-3
        if all(val < 1e-3 for val in vals):
            continue
        # Plot a light-colored bar as background
        for j in range(len(vals)):
            plt.bar(plot_index + j / len(vals), max(vals), width=1 / len(vals), color='lightgray')
        # Plot the actual bars
        plt.bar([plot_index + j / len(vals) for j in range(len(vals))], vals, width=1 / len(vals))
        plot_w.append(w)
        plot_index += 1

    plt.xlabel('Generator')
    plt.xticks(ticks=[])
    plt.ylabel('$P_{wind}$ [MW]')
    plt.show()



def plot_all_pG(hcs):
    # Initialize an empty dictionary to store the values
    values = [[] for _ in range(len(hcs[0].model.WIND))]
    ws = list(hcs[0].model.WIND)

    # Exclude list
    exclude_w = [105, 106, 110]
    exclude_i = [ws.index(i) for i in exclude_w]

    # Iterate over all hcs
    for hc in hcs:
        # Iterate over all w in hc.model.WINDhcs
        for i, w in enumerate(ws):
            if hc.model.y[w].value:
                # Append hc.model.pG[w] to the list corresponding to w
                values[i].append(hc.model.pG[w].value)
            else:
                values[i].append(0)
    # Create a bar plot for each w
    for i, (i, vals) in enumerate(values):
        # Skip if all values are smaller than 1e-3
        if all(val < 1e-3 for val in vals):
            continue
        # Plot a light-colored bar as background
        for j in range(len(vals)):
            plt.bar(i + j / len(vals), max(vals), width=1 / len(vals), color='lightgray')
        # Plot the actual bars
        plt.bar([i + j / len(vals) for j in range(len(vals))], vals, width=1 / len(vals))
    plt.xlabel('Generator Number')
    plt.ylabel('pG Value')
    plt.title('pG values for all w')
    plt.show()