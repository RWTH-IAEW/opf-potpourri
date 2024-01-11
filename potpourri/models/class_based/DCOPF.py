from pyomo.environ import *
from potpourri.models.class_based.DC import DC
from potpourri.models.class_based.OPF import OPF


class DCOPF(DC, OPF):
    def __init__(self, net):
        super().__init__(net)

        self.create_model()


    def create_model(self):
        super().create_model()
        self.model.name = "DCOPF"

        # --- cost function ---
        def objective(model):
            obj = sum(model.c1[g] * (model.baseMVA * model.pG[g]) + model.c0[g] for g in model.G) + \
                  sum(model.VOLL[d] * (model.PD[d] - model.pD[d]) * model.baseMVA for d in model.D)
            return obj

        self.model.OBJ = Objective(rule=objective, sense=minimize)

        # --- line power limits ---
        def line_lim1_def(model, l):
            return model.pL[l] <= model.SLmax[l]

        def line_lim2_def(model, l):
            return model.pL[l] >= -model.SLmax[l]

        self.model.line_lim_from = Constraint(self.model.L, rule=line_lim1_def)
        self.model.line_lim_to = Constraint(self.model.L, rule=line_lim2_def)

        # --- power flow limits on transformer lines---
        def transf_lim1_def(model, l):
            return model.pLT[l] <= model.SLmaxT[l]

        def transf_lim2_def(model, l):
            return model.pLT[l] >= -model.SLmaxT[l]

        self.model.transf_lim1 = Constraint(self.model.TRANSF, rule=transf_lim1_def)
        self.model.transf_lim2 = Constraint(self.model.TRANSF, rule=transf_lim2_def)


