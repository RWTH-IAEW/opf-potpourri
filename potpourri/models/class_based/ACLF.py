from pyomo.environ import *
from potpourri.models.class_based.AC import AC


class ACLF(AC):
    def __init__(self, net):
        super().__init__(net)

        self.create_model()

    def create_model(self):
        super().create_model()
        self.model.name = "ACLF"

        # voltage targets
        self.model.VS = Param(self.model.G, within=NonNegativeReals, initialize=1.)  # voltage target at generator buses

        # variables that model the violations (to help with the convergence)
        self.model.epsV_up = Var(self.model.G, initialize=0.0, domain=NonNegativeReals)
        self.model.epsV_down = Var(self.model.G, initialize=0.0, domain=NonNegativeReals)
        self.model.epsPG_up = Var(self.model.eG, initialize=0.0, domain=NonNegativeReals)
        self.model.epsPG_down = Var(self.model.eG, initialize=0.0, domain=NonNegativeReals)

        # --- cost function ---
        def objective(model):
            obj = 0.8 * sum(model.epsV_up[g] + model.epsV_down[g] for g in model.G) + \
                  0.1 * sum((model.pLto[l] + model.pLfrom[l]) for l in model.L) + \
                  0.1 * sum((model.epsPG_up[g] + model.epsPG_down[g]) for g in model.eG)
            return obj

        self.model.OBJ = Objective(rule=objective, sense=minimize)

        # --- voltage target constraints ---
        def bus_max_voltage(model, b, g):
            return model.v[b] <= model.VS[g] + model.epsV_up[g]

        def bus_min_voltage(model, b, g):
            return model.v[b] >= model.VS[g] - model.epsV_down[g]

        self.model.Vmaxc = Constraint(self.model.Gbs, rule=bus_max_voltage)
        self.model.Vminc = Constraint(self.model.Gbs, rule=bus_min_voltage)





