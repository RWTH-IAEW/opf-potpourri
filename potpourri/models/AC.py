from pyomo.environ import *
from potpourri.models.basemodel import Basemodel


class AC(Basemodel):
    def __init__(self, net):
        super().__init__(net)

        self.BB_data = - self.net.shunt.q_mvar * self.net.shunt.step / self.baseMVA

        # line and transformer addmittances
        r = self.net._ppc["branch"][:, 2].real
        x = self.net._ppc["branch"][:, 3].real
        y = (self.net._ppc["branch"][:, 4] * 1j)
        gt_ik = r / (r ** 2 + x ** 2)
        bt_ik = -x / (r ** 2 + x ** 2)
        BiiT = bt_ik + y.imag / 2
        BikT = -bt_ik
        GiiT = gt_ik + y.real / 2
        GikT = -gt_ik
        trafo_start = len(self.net.line)
        trafo_end = trafo_start + len(self.net.trafo)

        self.line_data = self.line_data.assign(
            **{'Bii_data': BiiT[:trafo_start], 'Bik_data': BikT[:trafo_start], 'Gii_data': GiiT[:trafo_start],
               'Gik_data': GikT[:trafo_start]})
        self.trafo_data = self.trafo_data.assign(
            **{'BiiT_data': BiiT[trafo_start:trafo_end], 'BikT_data': BikT[trafo_start:trafo_end],
               'GiiT_data': GiiT[trafo_start:trafo_end], 'GikT_data': GikT[trafo_start:trafo_end]})

        # generator and external grids voltage set points
        self.generation_data['v'] = self.net._ppc['gen'][:, 5]

        qsg = (self.net.sgen.q_mvar.fillna(0) * self.net.sgen.scaling / self.baseMVA).values
        self.static_generation_data['q'] = qsg

        self.QD_data = self.net.load.q_mvar * self.net.load.scaling / self.baseMVA

        self.create_model()

    def create_model(self):
        super().create_model()

        self.model.name = "AC"

        # shunt
        self.model.BB = Param(self.model.SHUNT, within=Reals,
                              initialize=self.BB_data[self.model.SHUNT])  # shunt susceptance

        # derived line parameters
        self.model.Bii = Param(self.model.L, within=Reals, initialize=self.line_data.Bii_data[self.model.L])
        self.model.Bik = Param(self.model.L, within=Reals, initialize=self.line_data.Bik_data[self.model.L])
        self.model.Gii = Param(self.model.L, within=Reals, initialize=self.line_data.Gii_data[self.model.L])
        self.model.Gik = Param(self.model.L, within=Reals, initialize=self.line_data.Gik_data[self.model.L])

        ## derived transformer parameters
        self.model.BiiT = Param(self.model.TRANSF, within=Reals,
                                initialize=self.trafo_data.BiiT_data[self.model.TRANSF])
        self.model.BikT = Param(self.model.TRANSF, within=Reals,
                                initialize=self.trafo_data.BikT_data[self.model.TRANSF])
        self.model.GiiT = Param(self.model.TRANSF, within=Reals,
                                initialize=self.trafo_data.GiiT_data[self.model.TRANSF])
        self.model.GikT = Param(self.model.TRANSF, within=Reals,
                                initialize=self.trafo_data.GikT_data[self.model.TRANSF])

        # reactive generation
        self.model.QsG = Param(self.model.sG, initialize=self.static_generation_data['q'][self.model.sG])

        self.model.v_bPV = Param(self.model.bPV, within=NonNegativeReals, initialize=self.bus_data.v_m[self.model.bPV])

        # reactive demand
        self.model.QD = Param(self.model.D, initialize=self.QD_data[self.model.D])

        # external grid voltage
        self.model.v_b0 = Param(self.model.b0, within=NonNegativeReals, initialize=self.bus_data.v_m[self.model.b0])

        # --- control variables ---
        self.model.qsG = Var(self.model.sG, domain=Reals)  # reactive power of static generators

        self.model.qD = Var(self.model.D, domain=Reals)  # reactive power absorbed by demand

        self.model.qLfrom = Var(self.model.L, domain=Reals)  # reactive power injected at b onto line
        self.model.qLto = Var(self.model.L, domain=Reals)  # reactive power injected at b' onto line
        self.model.qThv = Var(self.model.TRANSF, domain=Reals)  # reactive power injected at b onto transformer
        self.model.qTlv = Var(self.model.TRANSF, domain=Reals)  # reactive power injected at b' onto transformer

        self.model.v = Var(self.model.B, domain=NonNegativeReals, initialize=1.0)  # voltage magnitude at bus b, rad

        self.model.qG = Var(self.model.G, domain=Reals)

        # --- Kirchoff's current law at each bus b ---
        def KCL_real_def(model, b):
            kcl = sum(model.psG[g] for g in model.sG if (g, b) in model.sGbs) + \
                  sum(model.pG[g] for g in model.G if (g, b) in model.Gbs) == \
                  sum(model.pD[d] for d in model.D if (b, d) in model.Dbs) + \
                  sum(model.pLfrom[l] for l in model.L if model.A[l, 1] == b) + \
                  sum(model.pLto[l] for l in model.L if model.A[l, 2] == b) + \
                  sum(model.pThv[l] for l in model.TRANSF if model.AT[l, 1] == b) + \
                  sum(model.pTlv[l] for l in model.TRANSF if model.AT[l, 2] == b) + \
                  sum(model.GB[s] * model.v[b] ** 2 for s in model.SHUNT if
                      (b, s) in model.SHUNTbs and model.GB[s] != 0)
            if isinstance(kcl, bool):
                return Constraint.Skip
            return kcl

        def KCL_reactive_def(model, b):
            kcl = sum(model.qsG[g] for g in model.sG if (g, b) in model.sGbs) + \
                  sum(model.qG[g] for g in model.G if (g, b) in model.Gbs) == \
                  sum(model.qD[d] for d in model.D if (b, d) in model.Dbs) + \
                  sum(model.qLfrom[l] for l in model.L if model.A[l, 1] == b) + \
                  sum(model.qLto[l] for l in model.L if model.A[l, 2] == b) + \
                  sum(model.qThv[l] for l in model.TRANSF if model.AT[l, 1] == b) + \
                  sum(model.qTlv[l] for l in model.TRANSF if model.AT[l, 2] == b) - \
                  sum(model.BB[s] * model.v[b] ** 2 for s in model.SHUNT if
                      (b, s) in model.SHUNTbs and model.BB[s] != 0)
            if isinstance(kcl, bool):
                return Constraint.Skip
            return kcl

        self.model.KCL_real = Constraint(self.model.B, rule=KCL_real_def)
        self.model.KCL_reactive = Constraint(self.model.B, rule=KCL_reactive_def)

        # --- Kirchoff's voltage law on each line ---
        def KVL_real_fromend(model, l):
            return model.pLfrom[l] == model.Gii[l] * (model.v[model.A[l, 1]] ** 2) + \
                model.v[model.A[l, 1]] * model.v[model.A[l, 2]] * (
                        model.Bik[l] * sin(model.delta[model.A[l, 1]] - model.delta[model.A[l, 2]]) +
                        model.Gik[l] * cos(model.delta[model.A[l, 1]] - model.delta[model.A[l, 2]]))

        def KVL_real_toend(model, l):
            return model.pLto[l] == model.Gii[l] * (model.v[model.A[l, 2]] ** 2) + \
                model.v[model.A[l, 1]] * model.v[model.A[l, 2]] * (
                        model.Bik[l] * sin(model.delta[model.A[l, 2]] - model.delta[model.A[l, 1]]) +
                        model.Gik[l] * cos(model.delta[model.A[l, 2]] - model.delta[model.A[l, 1]]))

        def KVL_reactive_fromend(model, l):
            return model.qLfrom[l] == -model.Bii[l] * (model.v[model.A[l, 1]] ** 2) + \
                model.v[model.A[l, 1]] * model.v[model.A[l, 2]] * (
                        model.Gik[l] * sin(model.delta[model.A[l, 1]] - model.delta[model.A[l, 2]]) -
                        model.Bik[l] * cos(model.delta[model.A[l, 1]] - model.delta[model.A[l, 2]]))

        def KVL_reactive_toend(model, l):
            return model.qLto[l] == -model.Bii[l] * (model.v[model.A[l, 2]] ** 2) + \
                model.v[model.A[l, 1]] * model.v[model.A[l, 2]] * (
                        model.Gik[l] * sin(model.delta[model.A[l, 2]] - model.delta[model.A[l, 1]]) -
                        model.Bik[l] * cos(model.delta[model.A[l, 2]] - model.delta[model.A[l, 1]]))

        self.model.KVL_real_from = Constraint(self.model.L, rule=KVL_real_fromend)
        self.model.KVL_real_to = Constraint(self.model.L, rule=KVL_real_toend)
        self.model.KVL_reactive_from = Constraint(self.model.L, rule=KVL_reactive_fromend)
        self.model.KVL_reactive_to = Constraint(self.model.L, rule=KVL_reactive_toend)

        # --- Kirchoff's voltage law on each transformer line ---
        def KVL_real_fromendTransf(model, l):
            if model.shift[l]:
                return model.pThv[l] == model.GiiT[l] / model.Tap[l] ** 2 * (model.v[model.AT[l, 1]] ** 2) + \
                    model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] / model.Tap[l] * (
                            model.GikT[l] * cos(
                        model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]] - model.shift[l]) +
                            model.BikT[l] * sin(
                        model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]] - model.shift[l]))

            return model.pThv[l] == model.GiiT[l] / model.Tap[l] ** 2 * (model.v[model.AT[l, 1]] ** 2) + \
                model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] / model.Tap[l] * (
                        model.GikT[l] * cos(model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]]) +
                        model.BikT[l] * sin(model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]]))

        def KVL_real_toendTransf(model, l):
            if model.shift[l]:
                return model.pTlv[l] == model.GiiT[l] * (model.v[model.AT[l, 2]] ** 2) + \
                    model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] / model.Tap[l] * (
                            model.BikT[l] * sin(
                        model.delta[model.AT[l, 2]] - model.delta[model.AT[l, 1]] + model.shift[l]) +
                            model.GikT[l] * cos(
                        model.delta[model.AT[l, 2]] - model.delta[model.AT[l, 1]] + model.shift[l]))

            return model.pTlv[l] == model.GiiT[l] * (model.v[model.AT[l, 2]] ** 2) + \
                model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] / model.Tap[l] * (
                        model.BikT[l] * sin(model.delta[model.AT[l, 2]] - model.delta[model.AT[l, 1]]) +
                        model.GikT[l] * cos(model.delta[model.AT[l, 2]] - model.delta[model.AT[l, 1]]))

        def KVL_reactive_fromendTransf(model, l):
            if model.shift[l]:
                return model.qThv[l] == -model.BiiT[l] / model.Tap[l] ** 2 * (model.v[model.AT[l, 1]] ** 2) + \
                    model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] / model.Tap[l] * (
                            - model.BikT[l] * cos(
                        model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]] - model.shift[l]) +
                            model.GikT[l] * sin(
                        model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]] - model.shift[l]))

            return model.qThv[l] == -model.BiiT[l] / model.Tap[l] ** 2 * (model.v[model.AT[l, 1]] ** 2) + \
                model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] / model.Tap[l] * (
                        - model.BikT[l] * cos(model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]]) +
                        model.GikT[l] * sin(model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]]))

        def KVL_reactive_toendTransf(model, l):
            if model.shift[l]:
                return model.qTlv[l] == -model.BiiT[l] * (model.v[model.AT[l, 2]] ** 2) + \
                    model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] / model.Tap[l] * (
                            - model.BikT[l] * cos(
                        model.delta[model.AT[l, 2]] - model.delta[model.AT[l, 1]] + model.shift[l]) +
                            model.GikT[l] * sin(
                        model.delta[model.AT[l, 2]] - model.delta[model.AT[l, 1]] + model.shift[l]))

            return model.qTlv[l] == -model.BiiT[l] * (model.v[model.AT[l, 2]] ** 2) + \
                model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] / model.Tap[l] * (
                        - model.BikT[l] * cos(model.delta[model.AT[l, 2]] - model.delta[model.AT[l, 1]]) +
                        model.GikT[l] * sin(model.delta[model.AT[l, 2]] - model.delta[model.AT[l, 1]]))

        self.model.KVL_real_fromTransf = Constraint(self.model.TRANSF, rule=KVL_real_fromendTransf)
        self.model.KVL_real_toTransf = Constraint(self.model.TRANSF, rule=KVL_real_toendTransf)
        self.model.KVL_reactive_fromTransf = Constraint(self.model.TRANSF, rule=KVL_reactive_fromendTransf)
        self.model.KVL_reactive_toTransf = Constraint(self.model.TRANSF, rule=KVL_reactive_toendTransf)

        # --- reactive generator power limits ---
        # for g in self.model.sG:
        #    self.model.qsG[g].fix(self.model.QsG[g])  # reactive power of static generators fixed

        # --- reactive demand limits ---
        for d in self.model.D:
            self.model.qD[d].fix(self.model.QD[d])

        # --- generator voltage operating point ---
        def v_bPV_setpoint_rule(model, b):
            return model.v[b] == model.v_bPV[b]
        self.model.v_bPV_setpoint = Constraint(self.model.bPV, rule=v_bPV_setpoint_rule)

        # --- reference bus voltage constraint ---
        for b0 in self.model.b0:
            self.model.v[b0].fix(self.model.v_b0[b0])
