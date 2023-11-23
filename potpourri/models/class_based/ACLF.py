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

        # --- Kirchoff's voltage law on each transformer line ---
        def KVL_real_fromendTransf(model, l):
            if (model.AT[l, 1] in model.Bvolt) or (model.AT[l, 1] in model.Bvolt):
                return model.pLfromT[l] == (model.g[l] / (model.tap[l] ** 2)) * (model.v[model.AT[l, 1]] ** 2) - \
                    (1 / model.tap[l]) * model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] * (
                            model.b[l] * sin(model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]]) + model.g[
                        l] * cos(
                        model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]]))
            else:
                return model.pLfromT[l] == (model.g[l] / (model.Tap[l] ** 2)) * (model.v[model.AT[l, 1]] ** 2) - \
                    (1 / model.Tap[l]) * model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] * (
                            model.b[l] * sin(model.delta[model.AT[l, 1]] -
                                             model.delta[model.AT[l, 2]]) + model.g[l] * cos(
                        model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]]))

        def KVL_real_toendTransf(model, l):
            if (model.AT[l, 1] in model.Bvolt) or (model.AT[l, 1] in model.Bvolt):
                return model.pLtoT[l] == model.g[l] * (model.v[model.AT[l, 2]] ** 2) - \
                    (1 / model.tap[l]) * model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] * (
                            model.b[l] * sin(model.delta[model.AT[l, 2]] -
                                             model.delta[model.AT[l, 1]]) + model.g[l] * cos(
                        model.delta[model.AT[l, 2]] - model.delta[model.AT[l, 1]]))
            else:
                return model.pLtoT[l] == model.g[l] * (model.v[model.AT[l, 2]] ** 2) - \
                    (1 / model.Tap[l]) * model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] * (
                            model.b[l] * sin(model.delta[model.AT[l, 2]] - model.delta[model.AT[l, 1]]) + model.g[
                        l] * cos(
                        model.delta[model.AT[l, 2]] - model.delta[model.AT[l, 1]]))

        def KVL_reactive_fromendTransf(model, l):
            if (model.AT[l, 1] in model.Bvolt) or (model.AT[l, 1] in model.Bvolt):
                return model.qLfromT[l] == -(model.b[l] + 0.5 * model.bC[l]) * (model.v[model.AT[l, 1]] ** 2) / (
                        model.tap[l] ** 2) - \
                    (1 / model.tap[l]) * model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] * (
                            model.g[l] * sin(model.delta[model.AT[l, 1]] -
                                             model.delta[model.AT[l, 2]]) - model.b[l] * cos(
                        model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]]))
            else:
                return model.qLfromT[l] == -(model.b[l] + 0.5 * model.bC[l]) * (model.v[model.AT[l, 1]] ** 2) / (
                        model.Tap[l] ** 2) - \
                    (1 / model.Tap[l]) * model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] * (
                            model.g[l] * sin(model.delta[model.AT[l, 1]] -
                                             model.delta[model.AT[l, 2]]) - model.b[l] * cos(
                        model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]]))

        def KVL_reactive_toendTransf(model, l):
            if (model.AT[l, 1] in model.Bvolt) or (model.AT[l, 1] in model.Bvolt):
                return model.qLtoT[l] == -(model.b[l] + 0.5 * model.bC[l]) * (model.v[model.AT[l, 2]] ** 2) - \
                    (1 / model.tap[l]) * model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] * (
                            model.g[l] * sin(model.delta[model.AT[l, 2]] -
                                             model.delta[model.AT[l, 1]]) - model.b[l] * cos(
                        model.delta[model.AT[l, 2]] - model.delta[model.AT[l, 1]]))
            else:
                return model.qLtoT[l] == -(model.b[l]) * (model.v[model.AT[l, 2]] ** 2) - \
                    (1 / model.Tap[l]) * model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] * (
                            model.g[l] * sin(model.delta[model.AT[l, 2]] -
                                             model.delta[model.AT[l, 1]]) - model.b[l] * cos(
                        model.delta[model.AT[l, 2]] - model.delta[model.AT[l, 1]]))

        self.model.KVL_real_fromTransf = Constraint(self.model.TRANSF, rule=KVL_real_fromendTransf)
        self.model.KVL_real_toTransf = Constraint(self.model.TRANSF, rule=KVL_real_toendTransf)
        self.model.KVL_reactive_fromTransf = Constraint(self.model.TRANSF, rule=KVL_reactive_fromendTransf)
        self.model.KVL_reactive_toTransf = Constraint(self.model.TRANSF, rule=KVL_reactive_toendTransf)

        # --- voltage target constraints ---
        def bus_max_voltage(model, b, g):
            return model.v[b] <= model.VS[g] + model.epsV_up[g]

        def bus_min_voltage(model, b, g):
            return model.v[b] >= model.VS[g] - model.epsV_down[g]

        self.model.Vmaxc = Constraint(self.model.Gbs, rule=bus_max_voltage)
        self.model.Vminc = Constraint(self.model.Gbs, rule=bus_min_voltage)





