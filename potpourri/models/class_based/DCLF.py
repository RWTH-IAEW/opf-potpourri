from pyomo.environ import *
from potpourri.models.class_based.DC import DC


class DCLF(DC):
    def __init__(self, net):
        super().__init__(net)

        self.create_model()

    def create_model(self):
        super().create_model()

        self.model.name = "DCLF"

        # --- variables that model small violations on PG values. this is to help with the convergence
        self.model.epsG_up = Var(self.model.G, domain=NonNegativeReals)
        self.model.epsG_down = Var(self.model.G, domain=NonNegativeReals)

        # --- cost function ---
        def objective(model):
            obj = sum(model.baseMVA * (model.epsG_down[g] + model.epsG_up[g]) for (b, g) in model.Gbs) + \
                  sum(model.baseMVA * model.VOLL[d] * (model.PD[d] - model.pD[d]) for d in model.D)
            return obj

        self.model.OBJ = Objective(rule=objective, sense=minimize)
