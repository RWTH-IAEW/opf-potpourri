from pyomo.environ import *
from potpourri.models.class_based.ACOPF_base import ACOPF
from potpourri.models.class_based.basemodel import Basemodel


class HCOPF(ACOPF, Basemodel):
    def __init__(self, net):
        super().__init__(net)

        self.SWmax_data = 150 / self.baseMVA
        self.SWmin_data = 10 / self.baseMVA

        self.create_model()

    def create_model(self):
        super().create_model()

        self.model.name = "HCOPF"

        self.model.W = Set(initialize=self.model.B)

        self.model.pW = Var(self.model.W, domain=NonNegativeReals)
        self.model.qW = Var(self.model.W, domain=Reals)

        self.model.SWmax = Param(self.model.W, initialize=self.SWmax_data)
        self.model.SWmin = Param(self.model.W, initialize=self.SWmin_data)
        def objective(model):
            return sum(model.pW[w] for w in self.model.W)

        self.model.OBJ = Objective(rule=objective, sense=maximize)

        # --- Kirchoff's current law at each bus b ---
        def KCL_real_def(model, b):
            return (sum(model.pG[g] for g in model.G if (b, g) in model.Gbs) + model.pW[b] +
                    sum(model.peG[g] for g in model.eG if (b, g) in model.eGbs) ==
                    sum(model.pD[d] for d in model.D if (b, d) in model.Dbs) +
                    sum(model.pLfrom[l] for l in model.L if model.A[l, 1] == b) +
                    sum(model.pLto[l] for l in model.L if model.A[l, 2] == b) +
                    sum(model.pLfromT[l] for l in model.TRANSF if model.AT[l, 1] == b) +
                    sum(model.pLtoT[l] for l in model.TRANSF if model.AT[l, 2] == b) +
                    sum(model.GB[s] * model.v[b] ** 2 for s in model.SHUNT if (b, s) in model.SHUNTbs))

        def KCL_reactive_def(model, b):
            return (sum(model.qG[g] for g in model.G if (b, g) in model.Gbs) + model.qW[b] +
                    sum(model.qb0[s] for s in model.eG if s == b) + \
                    sum(model.qW[w] for w in model.WIND if (b, w) in model.Wbs) == \
                    sum(model.qD[d] for d in model.D if (b, d) in model.Dbs) + \
                    sum(model.qLfrom[l] for l in model.L if model.A[l, 1] == b) + \
                    sum(model.qLto[l] for l in model.L if model.A[l, 2] == b) + \
                    sum(model.qLfromT[l] for l in model.TRANSF if model.AT[l, 1] == b) + \
                    sum(model.qLtoT[l] for l in model.TRANSF if model.AT[l, 2] == b) - \
                    sum(model.BB[s] * model.v[b] ** 2 for s in model.SHUNT if (b, s) in model.SHUNTbs))

        self.model.KCL_real = Constraint(self.model.B, rule=KCL_real_def)
        self.model.KCL_reactive = Constraint(self.model.B, rule=KCL_reactive_def)

        def slack_bounds(model, b):
            return -200, model.peG[b], 200
        self.model.slack_constraint = Constraint(self.model.eG, rule=slack_bounds)

        def SW_bounds(model, w):
            return model.SWmin[w] ** 2, model.pW[w] ** 2 + model.qW[w] ** 2, model.SWmax[w] ** 2

        self.model.SW_constraint = Constraint(self.model.W, rule=SW_bounds)

        def power_factor(model, w):
            return -sin(acos(0.95)), (model.qW[w] / sqrt(model.pW[w] ** 2 + model.qW[w] ** 2)), sin(acos(0.95))

        self.model.power_factor_constraint = Constraint(self.model.W, rule=power_factor)

    def add_objectice(self):
        def objective(model):
            return sum(model.pW[w] for w in self.model.W)

        self.model.OBJ = Objective(rule=objective, sense=maximize)
