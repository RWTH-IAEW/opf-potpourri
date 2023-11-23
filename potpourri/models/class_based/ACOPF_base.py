from pyomo.environ import *
from potpourri.models.class_based.AC import AC
from potpourri.models.class_based.OPF import OPF


class ACOPF(AC, OPF):
    def __init__(self, net):
        super().__init__(net)

        self.Vmax_data = net.bus.max_vm_pu
        self.Vmin_data = net.bus.min_vm_pu

    def create_model(self):
        super().create_model()

        self.model.name = "ACOPF"

        self.model.Vmax = Param(self.model.B, within=NonNegativeReals, initialize=self.Vmax_data)  # max voltage angle
        self.model.Vmin = Param(self.model.B, within=NonNegativeReals, initialize=self.Vmin_data)  # min voltage angle

        # --- cost function ---
        def objective(model):
            obj = sum(
                model.c2[g] * (model.baseMVA * model.pG[g]) ** 2 + model.c1[g] * model.baseMVA * model.pG[g] +
                model.c0[g] for g in model.G) + \
                  sum(model.VOLL[d] * (model.PD[d] - model.pD[d]) * model.baseMVA for d in model.D)
            return obj

        self.model.OBJ = Objective(rule=objective, sense=minimize)
        # --- Kirchoff's voltage law on each transformer line ---
        def KVL_real_fromendTransf(model, l):
            return model.pLfromT[l] == model.GiiT[l] * (model.v[model.AT[l, 1]] ** 2) + \
                model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] * (
                        model.BikT[l] * sin(model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]]) +
                        model.GikT[l] * cos(model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]]))

        def KVL_real_toendTransf(model, l):
            return model.pLtoT[l] == model.GiiT[l] * (model.v[model.AT[l, 2]] ** 2) + \
                model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] * (
                        model.BikT[l] * sin(model.delta[model.AT[l, 2]] - model.delta[model.AT[l, 1]]) +
                        model.GikT[l] * cos(model.delta[model.AT[l, 2]] - model.delta[model.AT[l, 1]]))

        def KVL_reactive_fromendTransf(model, l):
            return model.qLfromT[l] == -model.BiiT[l] * (model.v[model.AT[l, 1]] ** 2) + \
                model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] * (
                        model.GikT[l] * sin(model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]]) -
                        model.BikT[l] * cos(model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]]))

        def KVL_reactive_toendTransf(model, l):
            return model.qLtoT[l] == -model.BiiT[l] * (model.v[model.AT[l, 2]] ** 2) + \
                model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] * (
                        model.GikT[l] * sin(model.delta[model.AT[l, 2]] - model.delta[model.AT[l, 1]]) -
                        model.BikT[l] * cos(model.delta[model.AT[l, 2]] - model.delta[model.AT[l, 1]]))

        self.model.KVL_real_fromTransf = Constraint(self.model.TRANSF, rule=KVL_real_fromendTransf)
        self.model.KVL_real_toTransf = Constraint(self.model.TRANSF, rule=KVL_real_toendTransf)
        self.model.KVL_reactive_fromTransf = Constraint(self.model.TRANSF, rule=KVL_reactive_fromendTransf)
        self.model.KVL_reactive_toTransf = Constraint(self.model.TRANSF, rule=KVL_reactive_toendTransf)

        # --- line power limits ---
        def line_lim1_def(model, l):
            return model.pLfrom[l] ** 2 + model.qLfrom[l] ** 2 <= model.SLmax[l] ** 2

        def line_lim2_def(model, l):
            return model.pLto[l] ** 2 + model.qLto[l] ** 2 <= model.SLmax[l] ** 2

        self.model.line_lim1 = Constraint(self.model.L, rule=line_lim1_def)
        self.model.line_lim2 = Constraint(self.model.L, rule=line_lim2_def)

        # --- power flow limits on transformer lines---
        def transf_lim1_def(model, l):
            return model.pLfromT[l] ** 2 + model.qLfromT[l] ** 2 <= model.SLmaxT[l] ** 2

        def transf_lim2_def(model, l):
            return model.pLtoT[l] ** 2 + model.qLtoT[l] ** 2 <= model.SLmaxT[l] ** 2

        self.model.transf_lim1 = Constraint(self.model.TRANSF, rule=transf_lim1_def)
        self.model.transf_lim2 = Constraint(self.model.TRANSF, rule=transf_lim2_def)

        # --- voltage constraints ---
        def v_bounds(model, b):
            return model.Vmin[b], model.v[b], model.Vmax[b]

        self.model.v_constraint = Constraint(self.model.B, rule=v_bounds)
