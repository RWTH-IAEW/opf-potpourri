import matplotlib.pyplot as plt
import pyomo.environ as pe
import pandapower as pp
import numpy as np
import math

from potpourri.models.class_based.HC_ACOPF import HC_ACOPF


def init_pGW_qGW(hc_model, pGW=1, qGW=0):
    for w in hc_model.WIND:
        hc_model.pG[w] = pGW / hc_model.baseMVA
        hc_model.qG[w] = qGW / hc_model.baseMVA


def set_SLmax(self, max_loading, unit: ['percent', 'MW']):
    if unit == 'percent':
        SLmax = max_loading / 100. * self.line_data['SLmax_data']
    elif unit == 'MW':
        self.line_data['SLmax_data'] = max_loading * np.ones(len(self.model.L)) / self.baseMVA
        SLmax = self.line_data['SLmax_data']
    for l in self.model.L:
        self.model.SLmax[l] = SLmax[l]


def vary_SLmax(hc):
    res = {'p_wind_mw': {}, 'q_wind_mvar': {}, 'v_bus_pu': {}, 'va_bus_degree': {}, 'line_loading_%': {},
           'p_line_mw': {}, 'q_line_mvar': {}, 'p_slack_mw': {}, 'q_slack_mvar': {}}
    obj = []
    for i in range(10, 110, 10):
        init_pGW_qGW(hc.model, pGW=10, qGW=0)
        #        hc = HC_ACOPF(net)
        set_SLmax(hc, i, 'percent')
        print("max loading: " + str(i))
        hc.solve()

        store_values(hc, res)
        obj.append(pe.value(hc.model.obj_hc))

    return res


def diff_SLmax(SLmax, hc=None, SWmax=None, peGmax=None):
    if hc:
        hc_act = hc
        init_pGW_qGW(hc_act.model, pGW=10, qGW=0)
        if SWmax:
            for w in hc_act.model.WIND:
                hc_act.model.SWmax[w] = SWmax
    else:
        if SWmax:
            hc_act = HC_ACOPF(net, SWmax=SWmax)
        else:
            hc_act = HC_ACOPF(net)
    if peGmax:
        hc_act.model.peGmax[0] = peGmax
    set_SLmax(hc_act, SLmax, 'MW')
    return hc_act


def diff_peGmax(peGmax, hc=None, SWmax=None, SLmax=None):
    if hc:
        hc_act = hc
        init_pGW_qGW(hc_act.model, pGW=10, qGW=0)
        if SWmax:
            for w in hc_act.model.WIND:
                hc_act.model.SWmax[w] = SWmax
        hc_act.model.peGmax[0] = peGmax
    else:
        if SWmax:
            hc_act = HC_ACOPF(net, SWmax=SWmax, peGmax=peGmax)
        else:
            hc_act = HC_ACOPF(net, peGmax=peGmax)
    if SLmax:
        set_SLmax(hc_act, SLmax, 'MW')
    return hc_act


def diff_SWmax(SWmax, hc=None, SLmax=None, peGmax=None):
    if hc:
        hc_act = hc
        init_pGW_qGW(hc_act.model, pGW=10, qGW=0)
        for w in hc_act.model.WIND:
            hc_act.model.SWmax[w] = SWmax
    else:
        hc_act = HC_ACOPF(net, SWmax=SWmax)

    if SLmax:
        set_SLmax(hc_act, SLmax)
    if peGmax:
        hc_act.model.peGmax[0] = peGmax

    return hc_act


def vary_slack_and_SLmax(values=None, key=None, hc=None, SL_percent=[100], SWmax=None, SLmax=None, peGmax=None,
                         save_net: bool = False,
                         folder=None, solver='ipopt', mip_solver=None):
    results = []
    if not values and not key:
        values = [0]
        key = 'SL_%'

    for val in values:
        if key == 'SL_%':
            val = []
        res = {'p_wind_mw': {}, 'q_wind_mvar': {}, 'v_bus_pu': {}, 'va_bus_degree': {}, 'line_loading_%': {},
               'p_line_mw': {}, 'q_line_mvar': {}, 'p_slack_mw': {}, 'q_slack_mvar': {}}
        obj = []
        x = []
        model = []
        limits = {'peGmax': [], 'SWmax': [], 'SLmax': []}

        for i in SL_percent:
            if key == 'peGmax':
                hc_act = diff_peGmax(val, hc, SWmax, SLmax)
            elif key == 'SLmax':
                hc_act = diff_SLmax(val, hc, SWmax, peGmax)
            elif key == 'SWmax':
                hc_act = diff_SWmax(val, hc, SLmax, peGmax)
            else:
                if hc:
                    hc_act = hc
                else:
                    hc_act = HC_ACOPF(net)
                if key == 'SL_%':
                    val.append(i)

            set_SLmax(hc_act, i, 'percent')

            print("max loading: " + str(i))
            if mip_solver:
                hc_act.solve(solver=solver, mip_solver=mip_solver)
            else:
                hc_act.solve(solver=solver)

            store_values(hc_act, res)
            obj.append(pe.value(hc_act.model.obj_hc))
            x.append(i)
            model.append(hc_act)
            baseMVA = hc_act.model.baseMVA.value
            limits['peGmax'].append(hc_act.model.peGmax[0].value * baseMVA)
            limits['SWmax'].append(hc_act.model.SWmax[3].value * baseMVA)
            limits['SLmax'].append(hc_act.model.SLmax[0].value * baseMVA)

            if save_net:
                pp.to_excel(hc_act.net, folder + key + '_' + str(val) + '_SWmax_' + str(SWmax) + '_SLmax_' + str(
                    SLmax) + '_peGmax_' + str(peGmax) + '.xlsx')

        results.append(
            {'varied': {'key': key, 'val': val}, 'measurements': res, 'obj': obj, 'x': x, 'model': model,
             'limits': limits})

    return results


def store_values(hc, res):
    for w in hc.model.WIND:
        res['p_wind_mw'].setdefault(w, []).append(hc.net.res_sgen.p_mw[w])
        res['q_wind_mvar'].setdefault(w, []).append(hc.net.res_sgen.q_mvar[w])

    for b in hc.model.B:
        res['v_bus_pu'].setdefault(b, []).append(hc.net.res_bus.vm_pu[b])
        res['va_bus_degree'].setdefault(b, []).append(hc.net.res_bus.va_degree[b])

    for l in hc.model.L:
        res['line_loading_%'].setdefault(l, []).append(hc.net.res_line.loading_percent[l])
        res['p_line_mw'].setdefault(l, []).append(hc.net.res_line.pl_mw[l])
        res['q_line_mvar'].setdefault(l, []).append(hc.net.res_line.ql_mvar[l])

    for g in hc.model.eG:
        res['p_slack_mw'].setdefault(g, []).append(hc.net.res_ext_grid.p_mw[g])
        res['q_slack_mvar'].setdefault(g, []).append(hc.net.res_ext_grid.q_mvar[g])

    return res


def plt_diff_SLmax(dicts, keys_to_plt=None, save: bool = False, folder='../hc_test_results/figures/'):
    marker = ['o', '*', '+', '.', ',']
    # x = np.linspace(100, 10, 10)
    res = dicts[0]
    x = res['x']

    if not keys_to_plt:
        keys_to_plt = list(res['measurements'].keys())

    plt.style.use('rwth-word')
    for (key, entries) in res['measurements'].items():
        if key in keys_to_plt:
            for (b, values) in entries.items():
                plt.plot(x, values, label=str(b) + ' - ' + res['varied']['key'] + ' ' + str(res['varied']['val']),
                         marker=marker[0], markersize=10)
                for i in range(len(dicts) - 1):
                    res_act = dicts[i + 1]
                    val = res_act['measurements'][key][b]
                    plt.plot(x, val,
                             label=str(b) + ' - ' + res_act['varied']['key'] + ' ' + str(res_act['varied']['val']),
                             marker=marker[i + 1], markersize=10)
            plt.title(str(key))
            # plt.suptitle('p_max slack = ' + str(hc.model.peGmax[0].value) + " MW")
            plt.xlabel("max line loading [%]")
            plt.grid()
            if save:
                plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
                plt.savefig(folder + str(key) + '_' + res['varied']['key'] + '.png', bbox_inches="tight")
            else:
                plt.legend()
            plt.show()



def plt_pq_wind_from_res(res, save: bool = False, folder=None):
    plt.style.use('rwth-word')

    markers = ['.', '*', '+', 'o', ',']
    clrs = ['#00549F', '#000000', '#E30066', '#FFED00', '#006165',
            '#0098A1', '#57AB27', '#BDCD00', '#F6A800', '#CC071E',
            '#A11035', '#612158', '#7A6FAC']

    fig, ax = plt.subplots()

    # --- PQ region from power factor ---
    # x = np.arange(0, 1.1, 0.1)
    # y = math.tan(math.acos(0.95)) * x
    #
    # ax.plot(x, y, x, -y, color=clrs[0])
    # ax.fill_between(x, y, -y, alpha=0.2)

    # --- PQ according to VDE Variante ---
    ax.vlines(1., -0.23, 0.48, color=clrs[0])  # Pmax
    ax.vlines(0.1, -0.1, 0.1, color=clrs[0])  # Pmin
    ax.plot([0.1, 0.2, 1], [-0.1, -0.23, -0.23], color=clrs[0])  # Qmin(P)
    ax.plot([0.1, 0.2, 1], [0.1, 0.48, 0.48], color=clrs[0])  # Qmax(P)

    ax.fill_between([0.1, 0.2, 1], [-0.1, -0.23, -0.23], [0.1, 0.48, 0.48], color=clrs[0], alpha=0.2)

    for key, dict in enumerate(res):
        for bus in dict['measurements']['p_wind_mw']:
            for i in range(len(dict['measurements']['p_wind_mw'][bus])):
                p = dict['measurements']['p_wind_mw'][bus][i]
                q = dict['measurements']['q_wind_mvar'][bus][i]
                SWmax = dict['limits']['SWmax'][i]

                if dict['varied']['key'] == 'SL_%':
                    value = dict['varied']['val'][i]
                    color = clrs[i + 1]
                else:
                    value = dict['varied']['val']
                    color = clrs[key + 1]

                ax.plot(p / p, q / p, label=str(bus) + ' - ' + dict['varied']['key'] + ' ' + str(value),
                        marker='.', color=color, markersize=10)

    plt.xlabel('P_wind / P_wind_inst')
    plt.ylabel('Q_wind / P_wind_inst')
    if save:
        plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
        plt.savefig(folder + 'pq_' + dict['varied']['key'] + '.png', bbox_inches="tight")
    else:
        plt.legend()

    plt.show()


# def plt_diff_SLmax(dicts):
#     marker = [',', '+', '.', 'o', '*']
#     x = np.linspace(100, 10, 10)
#     res = dicts[0]
#     for (key, entries) in res.items():
#         for (b, values) in entries.items():
#             plt.plot(x, values, label="bus " + str(b), marker=marker[0])
#             for i in range(len(dicts) - 1):
#                 res_act = dicts[i + 1]
#                 val = res_act[key][b]
#                 plt.plot(x, val, label=str(b) + ' ' + str(i+1))
#         plt.title(str(key))
#         # plt.suptitle('p_max slack = ' + str(hc.model.peGmax[0].value) + " MW")
#         plt.xlabel("max line loading [%]")
#         plt.grid()
#         plt.legend()
#         plt.show()


def plt_results(res, x, xlabel, suptitle):
    for (meas, entries) in res.items():
        for (b, values) in entries.items():
            plt.plot(x, values, label=str(b), marker='.')
        plt.title(str(meas))
        plt.suptitle(suptitle)
        plt.xlabel(xlabel)
        plt.grid()
        plt.legend()
        plt.show()


def plt_pq_wind_from_hc(hc):
    plt.style.use('rwth-word')

    x = np.arange(0, 1.1, 0.1)
    fig, ax = plt.subplots()

    y = math.tan(math.acos(0.95)) * x

    ax.plot(x, y, x, -y)
    ax.fill_between(x, y, -y, alpha=0.2)

    for w in hc.model.WIND:
        ax.plot(hc.model.pG[w].value / hc.model.SWmax[w].value, hc.model.qG[w].value / hc.model.SWmax[w].value,
                label=str(w), marker='o')
    plt.xlabel('P_wind / S_wind_max')
    plt.ylabel('Q_wind / S_wind_max')
    plt.legend()
    plt.show()



if __name__ == '__main__':
    net = pp.networks.simple_four_bus_system()