import pickle
import pyomo.environ as pe

from potpourri.models.HC_ACOPF import HC_ACOPF


def create_hc_cc_basics(net, n_scenarios, eps, Delta):
    hc = HC_ACOPF(net)
    # hc.solve()
    hc.add_OPF(SWmin=0)
    # hc.solve(solver='mindtpy')
    hc.fix_vars('y', 1.)

    hc.model.eps = pe.Param(mutable=True, initialize=eps)
    hc.model.Delta = pe.Param(mutable=True, initialize=Delta)

    hc.model.S = pe.Set(initialize=range(n_scenarios))

    hc.model.u = pe.Var(hc.model.S, within=pe.Binary, initialize=0)

    @hc.model.Constraint()
    def success_probability(model):
        return sum(model.u[sc] for sc in model.S) <= model.eps * len(model.S)

    return hc


def create_model_cc_node_branch(net, scenarios, case, eps, Delta):
    n_scenarios = len(scenarios['load'][0][0])
    scenario_busses = scenarios['load_sgen_busses']

    hc = create_hc_cc_basics(net, n_scenarios, eps, Delta)

    hc.model.Bsc = pe.Set(within=hc.model.B, initialize=hc.bus_lookup[scenario_busses])

    sgen_p_dict = {(b, sc): scenarios['sgen'][case][i, sc] for i, b in enumerate(hc.model.Bsc) for sc in hc.model.S}
    sgen_q_dict = {(b, sc): scenarios['sgen_q'][case][i, sc] for i, b in enumerate(hc.model.Bsc) for sc in hc.model.S}
    load_p_dict = {(b, sc): scenarios['load'][case][i, sc] for i, b in enumerate(hc.model.Bsc) for sc in hc.model.S}
    load_q_dict = {(b, sc): scenarios['load_q'][case][i, sc] for i, b in enumerate(hc.model.Bsc) for sc in hc.model.S}

    hc.model.PG_sc = pe.Param(hc.model.Bsc, hc.model.S, within=pe.NonNegativeReals, initialize=sgen_p_dict)
    hc.model.QG_sc = pe.Param(hc.model.Bsc, hc.model.S, within=pe.NonNegativeReals, initialize=sgen_q_dict)
    hc.model.PD_sc = pe.Param(hc.model.Bsc, hc.model.S, within=pe.NonNegativeReals, initialize=load_p_dict)
    hc.model.QD_sc = pe.Param(hc.model.Bsc, hc.model.S, within=pe.NonNegativeReals, initialize=load_q_dict)

    hc.model.v = pe.Var(hc.model.B, hc.model.S, within=pe.NonNegativeReals, initialize=1.)
    hc.model.delta = pe.Var(hc.model.B, hc.model.S, within=pe.Reals, initialize=0.)

    hc.model.pLfrom = pe.Var(hc.model.L, hc.model.S, within=pe.Reals)
    hc.model.pLto = pe.Var(hc.model.L, hc.model.S, within=pe.Reals)
    hc.model.qLfrom = pe.Var(hc.model.L, hc.model.S, within=pe.Reals)
    hc.model.qLto = pe.Var(hc.model.L, hc.model.S, within=pe.Reals)

    hc.model.pThv = pe.Var(hc.model.TRANSF, hc.model.S, within=pe.Reals)
    hc.model.pTlv = pe.Var(hc.model.TRANSF, hc.model.S, within=pe.Reals)
    hc.model.qThv = pe.Var(hc.model.TRANSF, hc.model.S, within=pe.Reals)
    hc.model.qTlv = pe.Var(hc.model.TRANSF, hc.model.S, within=pe.Reals)

    def M_P(model, b, sc):
        return (model.PG_sc[b, sc] - model.PD_sc[b, sc] +
                sum(model.psG[g] for g in (model.sG - model.sGc - model.WIND_HC) if (g, b) in model.sGbs) +
                sum(model.sPGmax[g] for g in model.sGc if (g, b) in model.sGbs) +
                sum(model.pWmax[w] for w in model.WIND_HC if (w, b) in model.sGbs) +
                sum(model.PGmax[g] for g in model.G if (g, b) in model.Gbs) +
                sum(model.SLmax[l] for l in model.L if model.A[l, 1] == b) +
                sum(model.SLmax[l] for l in model.L if model.A[l, 2] == b) +
                sum(model.SLmaxT[l] for l in model.TRANSF if model.AT[l, 1] == b) +
                sum(model.SLmaxT[l] for l in model.TRANSF if model.AT[l, 2] == b) +
                - (sum(model.pD[d] for d in (model.D - model.Dc) if (b, d) in model.Dbs) +
                   sum(model.PDmin[d] for d in model.Dc if (b, d) in model.Dbs) +
                   sum(model.GB[s] * model.Vmin[b] ** 2 for s in model.SHUNT if
                       (b, s) in model.SHUNTbs and model.GB[s] != 0))
                - model.Delta)

    hc.model.M_P = pe.Expression(hc.model.Bsc, hc.model.S, rule=M_P)

    def M_Q(model, b, sc):
        return (model.QG_sc[b, sc] - model.QD_sc[b, sc] +
                sum(model.qsG[g] for g in (model.sG - model.sGc - model.WIND_HC) if (g, b) in model.sGbs) +
                sum(model.sQGmax[g] for g in model.sGc if (g, b) in model.sGbs) +
                sum(model.pWmax[w] for w in model.WIND_HC if (w, b) in model.sGbs) +
                sum(model.QGmax[g] for g in model.G if (g, b) in model.Gbs) +
                sum(model.SLmax[l] for l in model.L if model.A[l, 1] == b) +
                sum(model.SLmax[l] for l in model.L if model.A[l, 2] == b) +
                sum(model.SLmaxT[l] for l in model.TRANSF if model.AT[l, 1] == b) +
                sum(model.SLmaxT[l] for l in model.TRANSF if model.AT[l, 2] == b) +
                - (sum(model.qD[d] for d in (model.D - model.Dc) if (b, d) in model.Dbs) +
                   sum(model.QDmin[d] for d in model.Dc if (b, d) in model.Dbs) +
                   sum(model.BB[s] * model.Vmin[b] ** 2 for s in model.SHUNT if
                       (b, s) in model.SHUNTbs and model.GB[s] != 0))
                - model.Delta)

    hc.model.M_Q = pe.Expression(hc.model.Bsc, hc.model.S, rule=M_Q)

    def KCL_real_def(model, b, sc):
        kcl = ((model.PG_sc[b, sc] - model.PD_sc[b, sc]) if b in model.Bsc else 0) + \
              sum(model.psG[g] for g in model.sG if (g, b) in model.sGbs) + \
              sum(model.pG[g] for g in model.G if (g, b) in model.Gbs) - \
              (sum(model.pD[d] for d in model.D if (b, d) in model.Dbs) +
               sum(model.pLfrom[l, sc] for l in model.L if model.A[l, 1] == b) +
               sum(model.pLto[l, sc] for l in model.L if model.A[l, 2] == b) +
               sum(model.pThv[l, sc] for l in model.TRANSF if model.AT[l, 1] == b) +
               sum(model.pTlv[l, sc] for l in model.TRANSF if model.AT[l, 2] == b) +
               sum(model.GB[s] * model.v[b, sc] ** 2 for s in model.SHUNT if
                   (b, s) in model.SHUNTbs and model.GB[s] != 0)) <= \
              model.Delta + (model.u[sc] * model.M_P[b, sc] if b in model.Bsc else 0)
        if isinstance(kcl, bool):
            return pe.Constraint.Skip
        return kcl

    def KCL_real_def_neg(model, b, sc):
        kcl = ((model.PG_sc[b, sc] - model.PD_sc[b, sc]) if b in model.Bsc else 0) + \
              sum(model.psG[g] for g in model.sG if (g, b) in model.sGbs) + \
              sum(model.pG[g] for g in model.G if (g, b) in model.Gbs) - \
              (sum(model.pD[d] for d in model.D if (b, d) in model.Dbs) +
               sum(model.pLfrom[l, sc] for l in model.L if model.A[l, 1] == b) +
               sum(model.pLto[l, sc] for l in model.L if model.A[l, 2] == b) +
               sum(model.pThv[l, sc] for l in model.TRANSF if model.AT[l, 1] == b) +
               sum(model.pTlv[l, sc] for l in model.TRANSF if model.AT[l, 2] == b) +
               sum(model.GB[s] * model.v[b, sc] ** 2 for s in model.SHUNT if
                   (b, s) in model.SHUNTbs and model.GB[s] != 0)) >= \
              - model.Delta - (model.u[sc] * model.M_P[b, sc] if b in model.Bsc else 0)
        if isinstance(kcl, bool):
            return pe.Constraint.Skip
        return kcl

    def KCL_reactive_def(model, b, sc):
        kcl = ((model.QG_sc[b, sc] - model.QD_sc[b, sc]) if b in model.Bsc else 0) + \
              sum(model.qsG[g] for g in model.sG if (g, b) in model.sGbs) + \
              sum(model.qG[g] for g in model.G if (g, b) in model.Gbs) - \
              (sum(model.qD[d] for d in model.D if (b, d) in model.Dbs) +
               sum(model.qLfrom[l, sc] for l in model.L if model.A[l, 1] == b) +
               sum(model.qLto[l, sc] for l in model.L if model.A[l, 2] == b) +
               sum(model.qThv[l, sc] for l in model.TRANSF if model.AT[l, 1] == b) +
               sum(model.qTlv[l, sc] for l in model.TRANSF if model.AT[l, 2] == b) -
               sum(model.BB[s] * model.v[b, sc] ** 2 for s in model.SHUNT if
                   (b, s) in model.SHUNTbs and model.BB[s] != 0)) <= \
              model.Delta + (model.u[sc] * model.M_Q[b, sc] if b in model.Bsc else 0)
        if isinstance(kcl, bool):
            return pe.Constraint.Skip
        return kcl

    def KCL_reactive_def_neg(model, b, sc):
        kcl = ((model.QG_sc[b, sc] - model.QD_sc[b, sc]) if b in model.Bsc else 0) + \
              sum(model.qsG[g] for g in model.sG if (g, b) in model.sGbs) + \
              sum(model.qG[g] for g in model.G if (g, b) in model.Gbs) - \
              (sum(model.qD[d] for d in model.D if (b, d) in model.Dbs) +
               sum(model.qLfrom[l, sc] for l in model.L if model.A[l, 1] == b) +
               sum(model.qLto[l, sc] for l in model.L if model.A[l, 2] == b) +
               sum(model.qThv[l, sc] for l in model.TRANSF if model.AT[l, 1] == b) +
               sum(model.qTlv[l, sc] for l in model.TRANSF if model.AT[l, 2] == b) -
               sum(model.BB[s] * model.v[b, sc] ** 2 for s in model.SHUNT if
                   (b, s) in model.SHUNTbs and model.BB[s] != 0)) >= \
              - model.Delta - (model.u[sc] * model.M_Q[b, sc] if b in model.Bsc else 0)
        if isinstance(kcl, bool):
            return pe.Constraint.Skip
        return kcl

    hc.model.KCL_real = pe.Constraint(hc.model.B, hc.model.S, rule=KCL_real_def)
    hc.model.KCL_real_neg = pe.Constraint(hc.model.B, hc.model.S, rule=KCL_real_def_neg)
    hc.model.KCL_reactive = pe.Constraint(hc.model.B, hc.model.S, rule=KCL_reactive_def)
    hc.model.KCL_reactive_neg = pe.Constraint(hc.model.B, hc.model.S, rule=KCL_reactive_def_neg)

    # --- Kirchoff's voltage law on each line ---
    def KVL_real_fromend(model, l, sc):
        return model.pLfrom[l, sc] == model.Gii[l] * (model.v[model.A[l, 1], sc] ** 2) + \
            model.v[model.A[l, 1], sc] * model.v[model.A[l, 2], sc] * (
                    model.Bik[l] * pe.sin(model.delta[model.A[l, 1], sc] - model.delta[model.A[l, 2], sc]) +
                    model.Gik[l] * pe.cos(model.delta[model.A[l, 1], sc] - model.delta[model.A[l, 2], sc]))

    def KVL_real_toend(model, l, sc):
        return model.pLto[l, sc] == model.Gii[l] * (model.v[model.A[l, 2], sc] ** 2) + \
            model.v[model.A[l, 1], sc] * model.v[model.A[l, 2], sc] * (
                    model.Bik[l] * pe.sin(model.delta[model.A[l, 2], sc] - model.delta[model.A[l, 1], sc]) +
                    model.Gik[l] * pe.cos(model.delta[model.A[l, 2], sc] - model.delta[model.A[l, 1], sc]))

    def KVL_reactive_fromend(model, l, sc):
        return model.qLfrom[l, sc] == -model.Bii[l] * (model.v[model.A[l, 1], sc] ** 2) + \
            model.v[model.A[l, 1], sc] * model.v[model.A[l, 2], sc] * (
                    model.Gik[l] * pe.sin(model.delta[model.A[l, 1], sc] - model.delta[model.A[l, 2], sc]) -
                    model.Bik[l] * pe.cos(model.delta[model.A[l, 1], sc] - model.delta[model.A[l, 2], sc]))

    def KVL_reactive_toend(model, l, sc):
        return model.qLto[l, sc] == -model.Bii[l] * (model.v[model.A[l, 2], sc] ** 2) + \
            model.v[model.A[l, 1], sc] * model.v[model.A[l, 2], sc] * (
                    model.Gik[l] * pe.sin(model.delta[model.A[l, 2], sc] - model.delta[model.A[l, 1], sc]) -
                    model.Bik[l] * pe.cos(model.delta[model.A[l, 2], sc] - model.delta[model.A[l, 1], sc]))

    hc.model.KVL_real_from = pe.Constraint(hc.model.L, hc.model.S, rule=KVL_real_fromend)
    hc.model.KVL_real_to = pe.Constraint(hc.model.L, hc.model.S, rule=KVL_real_toend)
    hc.model.KVL_reactive_from = pe.Constraint(hc.model.L, hc.model.S, rule=KVL_reactive_fromend)
    hc.model.KVL_reactive_to = pe.Constraint(hc.model.L, hc.model.S, rule=KVL_reactive_toend)

    # --- Kirchoff's voltage law on each transformer line ---
    def KVL_real_fromendTransf(model, l, sc):
        if model.shift[l]:
            return model.pThv[l, sc] == model.GiiT[l] / model.Tap[l] ** 2 * (model.v[model.AT[l, 1], sc] ** 2) + \
                model.v[model.AT[l, 1], sc] * model.v[model.AT[l, 2], sc] / model.Tap[l] * (
                        model.GikT[l] * pe.cos(
                    model.delta[model.AT[l, 1], sc] - model.delta[model.AT[l, 2], sc] - model.shift[l]) +
                        model.BikT[l] * pe.sin(
                    model.delta[model.AT[l, 1], sc] - model.delta[model.AT[l, 2], sc] - model.shift[l]))

        return model.pThv[l, sc] == model.GiiT[l] / model.Tap[l] ** 2 * (model.v[model.AT[l, 1], sc] ** 2) + \
            model.v[model.AT[l, 1], sc] * model.v[model.AT[l, 2], sc] / model.Tap[l] * (
                    model.GikT[l] * pe.cos(model.delta[model.AT[l, 1], sc] - model.delta[model.AT[l, 2], sc]) +
                    model.BikT[l] * pe.sin(model.delta[model.AT[l, 1], sc] - model.delta[model.AT[l, 2], sc]))

    def KVL_real_toendTransf(model, l, sc):
        if model.shift[l]:
            return model.pTlv[l, sc] == model.GiiT[l] * (model.v[model.AT[l, 2], sc] ** 2) + \
                model.v[model.AT[l, 1], sc] * model.v[model.AT[l, 2], sc] / model.Tap[l] * (
                        model.BikT[l] * pe.sin(
                    model.delta[model.AT[l, 2], sc] - model.delta[model.AT[l, 1], sc] + model.shift[l]) +
                        model.GikT[l] * pe.cos(
                    model.delta[model.AT[l, 2], sc] - model.delta[model.AT[l, 1], sc] + model.shift[l]))

        return model.pTlv[l, sc] == model.GiiT[l] * (model.v[model.AT[l, 2], sc] ** 2) + \
            model.v[model.AT[l, 1], sc] * model.v[model.AT[l, 2], sc] / model.Tap[l] * (
                    model.BikT[l] * pe.sin(model.delta[model.AT[l, 2], sc] - model.delta[model.AT[l, 1], sc]) +
                    model.GikT[l] * pe.cos(model.delta[model.AT[l, 2], sc] - model.delta[model.AT[l, 1], sc]))

    def KVL_reactive_fromendTransf(model, l, sc):
        if model.shift[l]:
            return model.qThv[l, sc] == -model.BiiT[l] / model.Tap[l] ** 2 * (model.v[model.AT[l, 1], sc] ** 2) + \
                model.v[model.AT[l, 1], sc] * model.v[model.AT[l, 2], sc] / model.Tap[l] * (
                        - model.BikT[l] * pe.cos(
                    model.delta[model.AT[l, 1], sc] - model.delta[model.AT[l, 2], sc] - model.shift[l]) +
                        model.GikT[l] * pe.sin(
                    model.delta[model.AT[l, 1], sc] - model.delta[model.AT[l, 2], sc] - model.shift[l]))

        return model.qThv[l, sc] == -model.BiiT[l] / model.Tap[l] ** 2 * (model.v[model.AT[l, 1], sc] ** 2) + \
            model.v[model.AT[l, 1], sc] * model.v[model.AT[l, 2], sc] / model.Tap[l] * (
                    - model.BikT[l] * pe.cos(model.delta[model.AT[l, 1], sc] - model.delta[model.AT[l, 2], sc]) +
                    model.GikT[l] * pe.sin(model.delta[model.AT[l, 1], sc] - model.delta[model.AT[l, 2], sc]))

    def KVL_reactive_toendTransf(model, l, sc):
        if model.shift[l]:
            return model.qTlv[l, sc] == -model.BiiT[l] * (model.v[model.AT[l, 2], sc] ** 2) + \
                model.v[model.AT[l, 1], sc] * model.v[model.AT[l, 2]] / model.Tap[l] * (
                        - model.BikT[l] * pe.cos(
                    model.delta[model.AT[l, 2], sc] - model.delta[model.AT[l, 1], sc] + model.shift[l]) +
                        model.GikT[l] * pe.sin(
                    model.delta[model.AT[l, 2], sc] - model.delta[model.AT[l, 1], sc] + model.shift[l]))

        return model.qTlv[l, sc] == -model.BiiT[l] * (model.v[model.AT[l, 2], sc] ** 2) + \
            model.v[model.AT[l, 1], sc] * model.v[model.AT[l, 2], sc] / model.Tap[l] * (
                    - model.BikT[l] * pe.cos(model.delta[model.AT[l, 2], sc] - model.delta[model.AT[l, 1], sc]) +
                    model.GikT[l] * pe.sin(model.delta[model.AT[l, 2], sc] - model.delta[model.AT[l, 1], sc]))

    hc.model.KVL_real_fromTransf = pe.Constraint(hc.model.TRANSF, hc.model.S, rule=KVL_real_fromendTransf)
    hc.model.KVL_real_toTransf = pe.Constraint(hc.model.TRANSF, hc.model.S, rule=KVL_real_toendTransf)
    hc.model.KVL_reactive_fromTransf = pe.Constraint(hc.model.TRANSF, hc.model.S, rule=KVL_reactive_fromendTransf)
    hc.model.KVL_reactive_toTransf = pe.Constraint(hc.model.TRANSF, hc.model.S, rule=KVL_reactive_toendTransf)

    # --- wind ---
    def QU_min_hc(model, w, sc):
        for (g, b) in model.sGbs:
            if g == w:
                return model.qsG[w] >= (hc.m_qu_min * model.v[b, sc] + hc.qu_min) * model.psG[w]

    hc.model.QU_min_hc_constraint = pe.Constraint(hc.model.WIND_HC, hc.model.S, rule=QU_min_hc)

    def QU_max_hc(model, w, sc):
        for (g, b) in model.sGbs:
            if g == w:
                return model.qsG[w] <= (hc.m_qu_max * model.v[b, sc] + hc.qu_max) * model.psG[w]

    hc.model.QU_max_hc_constraint = pe.Constraint(hc.model.WIND_HC, hc.model.S, rule=QU_max_hc)

    # --- voltage bounds ---
    def v_bounds(model, b, sc):
        return model.Vmin[b], model.v[b, sc], model.Vmax[b]

    hc.model.v_constraint = pe.Constraint(hc.model.B, hc.model.S, rule=v_bounds)

    # --- line power limits ---
    def line_lim_from_def(model, l, sc):
        return model.pLfrom[l, sc] ** 2 + model.qLfrom[l, sc] ** 2 <= model.SLmax[l] ** 2 * model.v[model.A[l, 1], sc] ** 2

    def line_lim_to_def(model, l, sc):
        return model.pLto[l, sc] ** 2 + model.qLto[l, sc] ** 2 <= model.SLmax[l] ** 2 * model.v[model.A[l, 2], sc] ** 2

    hc.model.line_lim_from = pe.Constraint(hc.model.L, hc.model.S, rule=line_lim_from_def)
    hc.model.line_lim_to = pe.Constraint(hc.model.L, hc.model.S, rule=line_lim_to_def)

    # --- power flow limits on transformer lines---
    def transf_lim1_def(model, l, sc):
        return model.pThv[l, sc] ** 2 + model.qThv[l, sc] ** 2 <= model.SLmaxT[l] ** 2 * model.v[model.AT[l, 1], sc] ** 2

    def transf_lim2_def(model, l, sc):
        return model.pTlv[l, sc] ** 2 + model.qTlv[l, sc] ** 2 <= model.SLmaxT[l] ** 2 * model.v[model.AT[l, 2], sc] ** 2

    hc.model.transf_lim1 = pe.Constraint(hc.model.TRANSF, hc.model.S, rule=transf_lim1_def)
    hc.model.transf_lim2 = pe.Constraint(hc.model.TRANSF, hc.model.S, rule=transf_lim2_def)

    # --- objective ---
    def obj_wind_loss_rule(model):
        return sum(model.psG[w] for w in model.WIND_HC)
            #     - sum(model.pLfrom[l, sc] + model.pLto[l, sc] for l in model.L) - sum(
            # model.pThv[t, sc] + model.pTlv[t, sc] for t in model.TRANSF))

    hc.model.obj = pe.Objective(rule=obj_wind_loss_rule, sense=pe.maximize)

    return hc



def create_model_cc_node(net, scenarios, case, eps, Delta):
    n_scenarios = len(scenarios['load'][0][0])
    scenario_busses = scenarios['load_sgen_busses']

    hc = create_hc_cc_basics(net, n_scenarios, eps, Delta)

    hc.model.Bsc = pe.Set(within=hc.model.B, initialize=hc.bus_lookup[scenario_busses])

    sgen_p_dict = {(b, sc): scenarios['sgen'][case][i, sc] for i, b in enumerate(hc.model.Bsc) for sc in hc.model.S}
    sgen_q_dict = {(b, sc): scenarios['sgen_q'][case][i, sc] for i, b in enumerate(hc.model.Bsc) for sc in hc.model.S}
    load_p_dict = {(b, sc): scenarios['load'][case][i, sc] for i, b in enumerate(hc.model.Bsc) for sc in hc.model.S}
    load_q_dict = {(b, sc): scenarios['load_q'][case][i, sc] for i, b in enumerate(hc.model.Bsc) for sc in hc.model.S}

    hc.model.PG_sc = pe.Param(hc.model.Bsc, hc.model.S, within=pe.NonNegativeReals, initialize=sgen_p_dict)
    hc.model.QG_sc = pe.Param(hc.model.Bsc, hc.model.S, within=pe.NonNegativeReals, initialize=sgen_q_dict)
    hc.model.PD_sc = pe.Param(hc.model.Bsc, hc.model.S, within=pe.NonNegativeReals, initialize=load_p_dict)
    hc.model.QD_sc = pe.Param(hc.model.Bsc, hc.model.S, within=pe.NonNegativeReals, initialize=load_q_dict)

    def M_P(model, b, sc):
        return (model.PG_sc[b, sc] - model.PD_sc[b, sc] +
                sum(model.psG[g] for g in (model.sG - model.sGc - model.WIND_HC) if (g, b) in model.sGbs) +
                sum(model.sPGmax[g] for g in model.sGc if (g, b) in model.sGbs) +
                sum(model.pWmax[w] for w in model.WIND_HC if (w, b) in model.sGbs) +
                sum(model.PGmax[g] for g in model.G if (g, b) in model.Gbs) +
                sum(model.SLmax[l] for l in model.L if model.A[l, 1] == b) +
                sum(model.SLmax[l] for l in model.L if model.A[l, 2] == b) +
                sum(model.SLmaxT[l] for l in model.TRANSF if model.AT[l, 1] == b) +
                sum(model.SLmaxT[l] for l in model.TRANSF if model.AT[l, 2] == b) +
                - (sum(model.pD[d] for d in (model.D - model.Dc) if (b, d) in model.Dbs) +
                   sum(model.PDmin[d] for d in model.Dc if (b, d) in model.Dbs) +
                   sum(model.GB[s] * model.Vmin[b] ** 2 for s in model.SHUNT if
                       (b, s) in model.SHUNTbs and model.GB[s] != 0))
                - model.Delta)

    hc.model.M_P = pe.Expression(hc.model.Bsc, hc.model.S, rule=M_P)

    def M_Q(model, b, sc):
        return (model.QG_sc[b, sc] - model.QD_sc[b, sc] +
                sum(model.qsG[g] for g in (model.sG - model.sGc - model.WIND_HC) if (g, b) in model.sGbs) +
                sum(model.sQGmax[g] for g in model.sGc if (g, b) in model.sGbs) +
                sum(model.pWmax[w] for w in model.WIND_HC if (w, b) in model.sGbs) +
                sum(model.QGmax[g] for g in model.G if (g, b) in model.Gbs) +
                sum(model.SLmax[l] for l in model.L if model.A[l, 1] == b) +
                sum(model.SLmax[l] for l in model.L if model.A[l, 2] == b) +
                sum(model.SLmaxT[l] for l in model.TRANSF if model.AT[l, 1] == b) +
                sum(model.SLmaxT[l] for l in model.TRANSF if model.AT[l, 2] == b) +
                - (sum(model.qD[d] for d in (model.D - model.Dc) if (b, d) in model.Dbs) +
                   sum(model.QDmin[d] for d in model.Dc if (b, d) in model.Dbs) +
                   sum(model.BB[s] * model.Vmin[b] ** 2 for s in model.SHUNT if
                       (b, s) in model.SHUNTbs and model.GB[s] != 0))
                - model.Delta)

    hc.model.M_Q = pe.Expression(hc.model.Bsc, hc.model.S, rule=M_Q)

    def KCL_real_def(model, b, sc):
        kcl = ((model.PG_sc[b, sc] - model.PD_sc[b, sc]) if b in model.Bsc else 0) + \
              sum(model.psG[g] for g in model.sG if (g, b) in model.sGbs) + \
              sum(model.pG[g] for g in model.G if (g, b) in model.Gbs) - \
              (sum(model.pD[d] for d in model.D if (b, d) in model.Dbs) +
               sum(model.pLfrom[l] for l in model.L if model.A[l, 1] == b) +
               sum(model.pLto[l] for l in model.L if model.A[l, 2] == b) +
               sum(model.pThv[l] for l in model.TRANSF if model.AT[l, 1] == b) +
               sum(model.pTlv[l] for l in model.TRANSF if model.AT[l, 2] == b) +
               sum(model.GB[s] * model.v[b] ** 2 for s in model.SHUNT if
                   (b, s) in model.SHUNTbs and model.GB[s] != 0)) <= \
              model.Delta + (model.u[sc] * model.M_P[b, sc] if b in model.Bsc else 0)
        if isinstance(kcl, bool):
            return pe.Constraint.Skip
        return kcl

    def KCL_real_def_neg(model, b, sc):
        kcl = ((model.PG_sc[b, sc] - model.PD_sc[b, sc]) if b in model.Bsc else 0) + \
              sum(model.psG[g] for g in model.sG if (g, b) in model.sGbs) + \
              sum(model.pG[g] for g in model.G if (g, b) in model.Gbs) - \
              (sum(model.pD[d] for d in model.D if (b, d) in model.Dbs) +
               sum(model.pLfrom[l] for l in model.L if model.A[l, 1] == b) +
               sum(model.pLto[l] for l in model.L if model.A[l, 2] == b) +
               sum(model.pThv[l] for l in model.TRANSF if model.AT[l, 1] == b) +
               sum(model.pTlv[l] for l in model.TRANSF if model.AT[l, 2] == b) +
               sum(model.GB[s] * model.v[b] ** 2 for s in model.SHUNT if
                   (b, s) in model.SHUNTbs and model.GB[s] != 0)) >= \
              - model.Delta - (model.u[sc] * model.M_P[b, sc] if b in model.Bsc else 0)
        if isinstance(kcl, bool):
            return pe.Constraint.Skip
        return kcl

    def KCL_reactive_def(model, b, sc):
        kcl = ((model.QG_sc[b, sc] - model.QD_sc[b, sc]) if b in model.Bsc else 0) + \
              sum(model.qsG[g] for g in model.sG if (g, b) in model.sGbs) + \
              sum(model.qG[g] for g in model.G if (g, b) in model.Gbs) - \
              (sum(model.qD[d] for d in model.D if (b, d) in model.Dbs) +
               sum(model.qLfrom[l] for l in model.L if model.A[l, 1] == b) +
               sum(model.qLto[l] for l in model.L if model.A[l, 2] == b) +
               sum(model.qThv[l] for l in model.TRANSF if model.AT[l, 1] == b) +
               sum(model.qTlv[l] for l in model.TRANSF if model.AT[l, 2] == b) -
               sum(model.BB[s] * model.v[b] ** 2 for s in model.SHUNT if
                   (b, s) in model.SHUNTbs and model.BB[s] != 0)) <= \
              model.Delta + (model.u[sc] * model.M_Q[b, sc] if b in model.Bsc else 0)
        if isinstance(kcl, bool):
            return pe.Constraint.Skip
        return kcl

    def KCL_reactive_def_neg(model, b, sc):
        kcl = ((model.QG_sc[b, sc] - model.QD_sc[b, sc]) if b in model.Bsc else 0) + \
              sum(model.qsG[g] for g in model.sG if (g, b) in model.sGbs) + \
              sum(model.qG[g] for g in model.G if (g, b) in model.Gbs) - \
              (sum(model.qD[d] for d in model.D if (b, d) in model.Dbs) +
               sum(model.qLfrom[l] for l in model.L if model.A[l, 1] == b) +
               sum(model.qLto[l] for l in model.L if model.A[l, 2] == b) +
               sum(model.qThv[l] for l in model.TRANSF if model.AT[l, 1] == b) +
               sum(model.qTlv[l] for l in model.TRANSF if model.AT[l, 2] == b) -
               sum(model.BB[s] * model.v[b] ** 2 for s in model.SHUNT if
                   (b, s) in model.SHUNTbs and model.BB[s] != 0)) >= \
              - model.Delta - (model.u[sc] * model.M_Q[b, sc] if b in model.Bsc else 0)
        if isinstance(kcl, bool):
            return pe.Constraint.Skip
        return kcl

    hc.model.KCL_real = pe.Constraint(hc.model.B, hc.model.S, rule=KCL_real_def)
    hc.model.KCL_real_neg = pe.Constraint(hc.model.B, hc.model.S, rule=KCL_real_def_neg)
    hc.model.KCL_reactive = pe.Constraint(hc.model.B, hc.model.S, rule=KCL_reactive_def)
    hc.model.KCL_reactive_neg = pe.Constraint(hc.model.B, hc.model.S, rule=KCL_reactive_def_neg)

    return hc


def create_model_cc_sum(net, scenarios, case, eps, Delta):
    n_scenarios = len(scenarios['load'][0][0])

    hc = create_hc_cc_basics(net, n_scenarios, eps, Delta)

    sc = 10
    p_base_sc = sum(scenarios['sgen'][case][:, sc]) - sum(scenarios['load'][case][:, sc])
    q_base_sc = sum(scenarios['sgen_q'][case][:, sc]) - sum(scenarios['load_q'][case][:, sc])
    p_sc = []
    q_sc = []
    for sc in range(len(scenarios['load'][0][0])):
        p_sc.append(p_base_sc - (sum(scenarios['sgen'][case][:, sc]) - sum(scenarios['load'][case][:, sc])))
        q_sc.append(q_base_sc - (sum(scenarios['sgen_q'][case][:, sc]) - sum(scenarios['load_q'][case][:, sc])))

    hc.model.P_sc = pe.Param(hc.model.S, within=pe.Reals, initialize=p_sc)
    hc.model.Q_sc = pe.Param(hc.model.S, within=pe.Reals, initialize=q_sc)

    def M_P(model, sc):
        return (model.P_sc[sc] +
                sum(model.psG[g] for g in (model.sG - model.sGc - model.WIND_HC)) +
                sum(model.sPGmax[g] for g in model.sGc) +
                sum(model.pWmax[w] for w in model.WIND_HC) +
                sum(model.PGmax[g] for g in model.G) +
                sum(model.SLmax[l] for l in model.L) +
                sum(model.SLmax[l] for l in model.L) +
                sum(model.SLmaxT[l] for l in model.TRANSF) +
                sum(model.SLmaxT[l] for l in model.TRANSF) +
                - (sum(model.pD[d] for d in (model.D - model.Dc)) +
                   sum(model.PDmin[d] for d in model.Dc) +
                   sum(model.GB[s] * model.Vmin[b] ** 2 for s in model.SHUNT for b in model.B if
                       (b, s) in model.SHUNTbs and model.GB[s] != 0))
                - model.Delta)

    def M_Q(model, sc):
        return (model.Q_sc[sc] +
                sum(model.qsG[g] for g in (model.sG - model.sGc - model.WIND_HC)) +
                sum(model.sQGmax[g] for g in model.sGc) +
                sum(model.pWmax[w] for w in model.WIND_HC) +
                sum(model.QGmax[g] for g in model.G) +
                sum(model.SLmax[l] for l in model.L) +
                sum(model.SLmax[l] for l in model.L) +
                sum(model.SLmaxT[l] for l in model.TRANSF) +
                sum(model.SLmaxT[l] for l in model.TRANSF) +
                - (sum(model.qD[d] for d in (model.D - model.Dc)) +
                      sum(model.QDmin[d] for d in model.Dc) +
                        sum(model.BB[s] * model.Vmin[b] ** 2 for s in model.SHUNT for b in model.B if
                            (b, s) in model.SHUNTbs and model.BB[s] != 0))
                - model.Delta)

    hc.model.M_P = pe.Param(hc.model.S, within=pe.Reals, initialize=M_P)
    hc.model.M_Q = pe.Param(hc.model.S, within=pe.Reals, initialize=M_Q)

    def p_leq(model, sc):
        return (model.P_sc[sc] +
                sum(model.pG[g] for g in model.G) +
                sum(model.psG[g] for g in model.sG) -
                (sum(model.pD[d] for d in model.D) +
                 sum(model.pLfrom[l] + model.pLto[l] for l in model.L) +
                 sum(model.pThv[t] + model.pTlv[t] for t in model.TRANSF) +
                 sum(model.GB[s] * model.v[b] ** 2 for s in model.SHUNT for b in model.B if
                     (b, s) in model.SHUNTbs and model.GB[s] != 0))) <= \
            model.Delta + (model.u[sc] * model.M_P[sc])

    def p_geq(model, sc):
        return (model.P_sc[sc] +
                sum(model.pG[g] for g in model.G) +
                sum(model.psG[g] for g in model.sG) -
                (sum(model.pD[d] for d in model.D) +
                 sum(model.pLfrom[l] + model.pLto[l] for l in model.L) +
                 sum(model.pThv[t] + model.pTlv[t] for t in model.TRANSF) +
                 sum(model.GB[s] * model.v[b] ** 2 for s in model.SHUNT for b in model.B if
                     (b, s) in model.SHUNTbs and model.GB[s] != 0))) >= \
            - model.Delta - (model.u[sc] * model.M_P[sc])

    def q_leq(model, sc):
        return (model.Q_sc[sc] +
                sum(model.qG[g] for g in model.G) +
                sum(model.qsG[g] for g in model.sG) -
                (sum(model.qD[d] for d in model.D) +
                 sum(model.qLfrom[l] + model.qLto[l] for l in model.L) +
                 sum(model.qThv[t] + model.qTlv[t] for t in model.TRANSF) -
                 sum(model.BB[s] * model.v[b] ** 2 for s in model.SHUNT for b in model.B if
                     (b, s) in model.SHUNTbs and model.BB[s] != 0))) <= \
            model.Delta + (model.u[sc] * model.M_Q[sc])

    def q_geq(model, sc):
        return (model.Q_sc[sc] +
                sum(model.qG[g] for g in model.G) +
                sum(model.qsG[g] for g in model.sG) -
                (sum(model.qD[d] for d in model.D) +
                 sum(model.qLfrom[l] + model.qLto[l] for l in model.L) +
                 sum(model.qThv[t] + model.qTlv[t] for t in model.TRANSF) -
                 sum(model.BB[s] * model.v[b] ** 2 for s in model.SHUNT for b in model.B if
                     (b, s) in model.SHUNTbs and model.BB[s] != 0))) >= \
            - model.Delta - (model.u[sc] * model.M_Q[sc])

    hc.model.p_leq = pe.Constraint(hc.model.S, rule=p_leq)
    hc.model.p_geq = pe.Constraint(hc.model.S, rule=p_geq)
    hc.model.q_leq = pe.Constraint(hc.model.S, rule=q_leq)
    hc.model.q_geq = pe.Constraint(hc.model.S, rule=q_geq)

    return hc


if __name__ == '__main__':
    net = pickle.load(
        open(
            'C:\\Users\\f.lohse\PycharmProjects\potpourri\potpourri\data\windpot\simbench_hv_grid_with_potential_2.pkl',
            'rb'))
    scenarios = pickle.load(open(
        '/potpourri/data/scenarios/sb_hv_grid_loadcases_assumption_scenarios_with_q.pkl',
        'rb'))

    # net.sgen['in_service'] = False
    # net.load['in_service'] = False
    # busses = scenarios['load_sgen_busses']
    #
    # hc = create_model_cc_node(net, scenarios, 3, 0.2, 10)
    # load_ind, sgen_ind = add_scenario_load_sgen_to_net(net, busses)

    loadcase = 'lW'
    factors = net.loadcases.loc[loadcase]
    net.load.p_mw *= factors['pload']
    net.load.q_mvar *= factors['qload']
    net.sgen.scaling[net.sgen.type == 'Wind'] = factors['Wind_p']
    net.sgen.scaling[net.sgen.type == 'PV'] = factors['PV_p']
    net.sgen.scaling[(net.sgen.type != 'Wind') & (net.sgen.type != 'Solar')] = factors['RES_p']
    net.ext_grid.vm_pu = factors['Slack_vm']

    hc = create_model_cc_sum(net, scenarios, 3, 0.2, 10)
