import pandas as pd
from pyomo.environ import *
from potpourri.models.class_based.basemodel import Basemodel
import numpy as np


class AC(Basemodel):
    def __init__(self, net):
        super().__init__(net)

        self.BB_data = - self.net.shunt.q_mvar * net.shunt.step / self.baseMVA

        ZN = self.net.bus.vn_kv ** 2 / self.baseMVA

        b_ik_data = -self.net.line.x_ohm_per_km / (
                (self.net.line.x_ohm_per_km ** 2 + self.net.line.r_ohm_per_km ** 2) * self.net.line.length_km) * ZN[
                        self.net.line.from_bus].values
        b_ik_C_data = 2 * np.pi * self.net.f_hz * self.net.line.c_nf_per_km * 1e-9 * self.net.line.length_km * ZN[
            self.net.line.from_bus].values
        g_ik_data = self.net.line.r_ohm_per_km / (
                (self.net.line.x_ohm_per_km ** 2 + self.net.line.r_ohm_per_km ** 2) * self.net.line.length_km) * ZN[
                        self.net.line.from_bus].values
        self.Bii_data = b_ik_data + 0.5 * b_ik_C_data
        self.Bik_data = -b_ik_data
        self.Gii_data = g_ik_data
        self.Gik_data = -g_ik_data

        gt_ik = self.trafo_parameters["r"] / (self.trafo_parameters["r"] ** 2 + self.trafo_parameters["x"] ** 2)
        bt_ik = -self.trafo_parameters["x"] / (self.trafo_parameters["r"] ** 2 + self.trafo_parameters["x"] ** 2)

        self.BiiT_data = pd.Series(bt_ik + self.trafo_parameters["y"].imag / 2, self.trafo_set)
        self.BikT_data = pd.Series(-bt_ik, self.trafo_set)
        self.GiiT_data = pd.Series(gt_ik + self.trafo_parameters["y"].real / 2, self.trafo_set)
        self.GikT_data = pd.Series(-gt_ik, self.trafo_set)

        self.v_gG_data = self.generators.vm_pu[self.gen_set]  # voltage set points for generators (not static)

        self.QsG_data = self.generators.q_mvar[self.sgen_set].fillna(0) * self.generators.scaling[
            self.sgen_set] / self.baseMVA  # reactive generation set points for static generators

        self.QD_data = self.net.load.q_mvar * self.net.load.scaling / self.baseMVA

        self.v_eG_data = self.net.ext_grid.vm_pu

        self.create_model()

    def create_model(self):
        super().create_model()

        self.model.name = "AC"

        # shunt
        self.model.BB = Param(self.model.SHUNT, within=Reals, initialize=self.BB_data[self.model.SHUNT])  # shunt susceptance

        # derived line parameters
        self.model.Bii = Param(self.model.L, within=Reals, initialize=self.Bii_data[self.model.L])
        self.model.Bik = Param(self.model.L, within=Reals, initialize=self.Bik_data[self.model.L])
        self.model.Gii = Param(self.model.L, within=Reals, initialize=self.Gii_data[self.model.L])
        self.model.Gik = Param(self.model.L, within=Reals, initialize=self.Gik_data[self.model.L])

        ## derived transformer parameters
        self.model.BiiT = Param(self.model.TRANSF, within=Reals, initialize=self.BiiT_data[self.model.TRANSF])
        self.model.BikT = Param(self.model.TRANSF, within=Reals, initialize=self.BikT_data[self.model.TRANSF])
        self.model.GiiT = Param(self.model.TRANSF, within=Reals, initialize=self.GiiT_data[self.model.TRANSF])
        self.model.GikT = Param(self.model.TRANSF, within=Reals, initialize=self.GikT_data[self.model.TRANSF])

        # reactive generation
        self.model.QsG = Param(self.model.sG, initialize=self.QsG_data[self.model.sG])

        self.model.v_gG = Param(self.model.gG, initialize=self.v_gG_data[self.model.gG])

        # reactive demand
        self.model.QD = Param(self.model.D, initialize=self.QD_data[self.model.D])

        # external grid voltage
        self.model.v_eG = Param(self.model.eG, within=NonNegativeReals, initialize=self.v_eG_data[self.model.eG])

        # --- control variables ---
        self.model.qG = Var(self.model.G, domain=Reals)  # reactive generator power

        self.model.qD = Var(self.model.D, domain=Reals)  # reactive power absorbed by demand

        self.model.qLfrom = Var(self.model.L, domain=Reals)  # reactive power injected at b onto line
        self.model.qLto = Var(self.model.L, domain=Reals)  # reactive power injected at b' onto line
        self.model.qThv = Var(self.model.TRANSF, domain=Reals)  # reactive power injected at b onto transformer
        self.model.qTlv = Var(self.model.TRANSF, domain=Reals)  # reactive power injected at b' onto transformer

        self.model.v = Var(self.model.B, domain=NonNegativeReals, initialize=1.0)  # voltage magnitude at bus b, rad

        self.model.qeG = Var(self.model.eG, domain=Reals)

        # --- Kirchoff's current law at each bus b ---
        def KCL_real_def(model, b):
            kcl = sum(model.pG[g] for g in model.G if (b, g) in model.Gbs) + \
                  sum(model.peG[g] for g in model.eG if (b, g) in model.eGbs) == \
                  sum(model.pD[d] for d in model.D if (b, d) in model.Dbs) + \
                  sum(model.pLfrom[l] for l in model.L if model.A[l, 1] == b) + \
                  sum(model.pLto[l] for l in model.L if model.A[l, 2] == b) + \
                  sum(model.pThv[l] for l in model.TRANSF if model.AT[l, 1] == b) + \
                  sum(model.pTlv[l] for l in model.TRANSF if model.AT[l, 2] == b) + \
                  sum(model.GB[s] * model.v[b] ** 2 for s in model.SHUNT if (b, s) in model.SHUNTbs and model.GB[s] != 0)
            if isinstance(kcl, bool):
                return Constraint.Skip
            return kcl

        def KCL_reactive_def(model, b):
            kcl = sum(model.qG[g] for g in model.G if (b, g) in model.Gbs) + \
                  sum(model.qeG[g] for g in model.eG if (b, g) in model.eGbs) == \
                  sum(model.qD[d] for d in model.D if (b, d) in model.Dbs) + \
                  sum(model.qLfrom[l] for l in model.L if model.A[l, 1] == b) + \
                  sum(model.qLto[l] for l in model.L if model.A[l, 2] == b) + \
                  sum(model.qThv[l] for l in model.TRANSF if model.AT[l, 1] == b) + \
                  sum(model.qTlv[l] for l in model.TRANSF if model.AT[l, 2] == b) - \
                  sum(model.BB[s] * model.v[b] ** 2 for s in model.SHUNT if (b, s) in model.SHUNTbs and model.BB[s] != 0)
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
            # if model.shift[l]:
            #     return model.pThv[l] == model.GiiT[l] / model.Tap[l] ** 2 * model.v[model.AT[l, 1]] ** 2 + \
            #         model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] / model.Tap[l] * (
            #             (model.GikT[l] * cos(model.shift[l]) - model.BikT[l] * sin(model.shift[l])) * cos(model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]]) +
            #             (model.GikT[l] * sin(model.shift[l]) + model.BikT[l] * cos(model.shift[l])) * sin(model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]]))

            if model.shift[l]:
                return model.pThv[l] == model.GiiT[l] / model.Tap[l] ** 2 * (model.v[model.AT[l, 1]] ** 2) + \
                model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] / model.Tap[l] * (
                        model.GikT[l] * cos(model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]] - model.shift[l]) +
                        model.BikT[l] * sin(model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]] - model.shift[l]))

            return model.pThv[l] == model.GiiT[l] / model.Tap[l] ** 2 * (model.v[model.AT[l, 1]] ** 2) + \
                model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] / model.Tap[l] * (
                        model.GikT[l] * cos(model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]]) +
                        model.BikT[l] * sin(model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]]))

        def KVL_real_toendTransf(model, l):
            # if model.shift[l]:
            #     return model.pTlv[l] == model.GiiT[l] * model.v[model.AT[l, 2]] ** 2 + \
            #         model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] / model.Tap[l] * (
            #             (model.GikT[l] * cos(model.shift[l]) + model.BikT[l] * sin(model.shift[l])) * cos(model.delta[model.AT[l, 2]] - model.delta[model.AT[l, 1]]) +
            #             (-model.GikT[l] * sin(model.shift[l]) + model.BikT[l] * cos(model.shift[l])) * sin(model.delta[model.AT[l, 2]] - model.delta[model.AT[l, 1]]))

            if model.shift[l]:
                return model.pTlv[l] == model.GiiT[l] * (model.v[model.AT[l, 2]] ** 2) + \
                model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] / model.Tap[l] * (
                        model.BikT[l] * sin(model.delta[model.AT[l, 2]] - model.delta[model.AT[l, 1]] + model.shift[l]) +
                        model.GikT[l] * cos(model.delta[model.AT[l, 2]] - model.delta[model.AT[l, 1]] + model.shift[l]))

            return model.pTlv[l] == model.GiiT[l] * (model.v[model.AT[l, 2]] ** 2) + \
                model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] / model.Tap[l] * (
                        model.BikT[l] * sin(model.delta[model.AT[l, 2]] - model.delta[model.AT[l, 1]]) +
                        model.GikT[l] * cos(model.delta[model.AT[l, 2]] - model.delta[model.AT[l, 1]]))

        def KVL_reactive_fromendTransf(model, l):
            # if model.shift[l]:
            #     return model.qThv[l] == -model.BiiT[l] / model.Tap[l] ** 2 * (model.v[model.AT[l, 1]] ** 2) + \
            #             model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] / model.Tap[l] * (
            #             -(model.GikT[l] * sin(model.shift[l]) + model.BikT[l] * cos(model.shift[l])) * cos(model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]]) +
            #             (model.GikT[l] * cos(model.shift[l]) - model.BikT[l] * sin(model.shift[l])) * sin(model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]]))

            if model.shift[l]:
                return model.qThv[l] == -model.BiiT[l] / model.Tap[l] ** 2 * (model.v[model.AT[l, 1]] ** 2) + \
                model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] / model.Tap[l] * (
                        - model.BikT[l] * cos(model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]] - model.shift[l]) +
                        model.GikT[l] * sin(model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]] - model.shift[l]))


            return model.qThv[l] == -model.BiiT[l] / model.Tap[l] ** 2 * (model.v[model.AT[l, 1]] ** 2) + \
                model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] / model.Tap[l] * (
                        - model.BikT[l] * cos(model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]]) +
                        model.GikT[l] * sin(model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]]))

        def KVL_reactive_toendTransf(model, l):
            # if model.shift[l]:
            #     return model.qTlv[l] == -model.BiiT[l] * (model.v[model.AT[l, 2]] ** 2) + \
            #         model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] / model.Tap[l] * (
            #                 (model.GikT[l] * sin(model.shift[l]) - model.BikT[l] * cos(model.shift[l])) * cos(
            #             model.delta[model.AT[l, 2]] - model.delta[model.AT[l, 1]]) +
            #                 (model.GikT[l] * cos(model.shift[l]) + model.BikT[l] * sin(model.shift[l])) * sin(
            #             model.delta[model.AT[l, 2]] - model.delta[model.AT[l, 1]]))

            if model.shift[l]:
                return model.qTlv[l] == -model.BiiT[l] * (model.v[model.AT[l, 2]] ** 2) + \
                model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] / model.Tap[l] * (
                        - model.BikT[l] * cos(model.delta[model.AT[l, 2]] - model.delta[model.AT[l, 1]] + model.shift[l]) +
                        model.GikT[l] * sin(model.delta[model.AT[l, 2]] - model.delta[model.AT[l, 1]] + model.shift[l]))

            return model.qTlv[l] == -model.BiiT[l] * (model.v[model.AT[l, 2]] ** 2) + \
                model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] / model.Tap[l] * (
                        - model.BikT[l] * cos(model.delta[model.AT[l, 2]] - model.delta[model.AT[l, 1]]) +
                        model.GikT[l] * sin(model.delta[model.AT[l, 2]] - model.delta[model.AT[l, 1]]))

        self.model.KVL_real_fromTransf = Constraint(self.model.TRANSF, rule=KVL_real_fromendTransf)
        self.model.KVL_real_toTransf = Constraint(self.model.TRANSF, rule=KVL_real_toendTransf)
        self.model.KVL_reactive_fromTransf = Constraint(self.model.TRANSF, rule=KVL_reactive_fromendTransf)
        self.model.KVL_reactive_toTransf = Constraint(self.model.TRANSF, rule=KVL_reactive_toendTransf)

        # --- reactive generator power limits ---
        for g in self.model.sG:
            self.model.qG[g].fix(self.model.QsG[g])  # reactive power of static generators fixed

        # --- reactive demand limits ---
        for d in self.model.D:
            self.model.qD[d].fix(self.model.QD[d])

        # --- generator voltage operation point ---
        for (b, g) in self.model.Gbs:
            if g in self.model.gG:
                self.model.v[b].fix(self.model.v_gG[g])

        # --- reference bus voltage constraint ---
        for (b, g) in self.model.eGbs:
            self.model.v[b].fix(self.model.v_eG[g])

        # def ref_bus_def(model, b, g):
        #     return model.v[b] == model.v_eG[g]
        #
        # self.model.refbus_v = Constraint(self.model.eGbs, rule=ref_bus_def)
