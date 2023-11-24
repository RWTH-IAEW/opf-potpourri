from pyomo.environ import *
from potpourri.models.class_based.ACOPF_base import ACOPF
from potpourri.models.class_based.basemodel import Basemodel


class HC_ACOPF(ACOPF):
    def __init__(self, net):
        super().__init__(net)

        net.sgen.wind_hc.fillna(False, inplace=True)

        self.SWmax_data = 150 / self.baseMVA
        self.SWmin_data = 10 / self.baseMVA

        self.wind_hc_set = net.sgen.index[net.sgen.wind_hc]
        self.bus_whc_set = list(zip(net.sgen.bus[net.sgen.wind_hc], self.wind_hc_set))

        self.create_model()

    def create_model(self):
        super().create_model()

        self.model.name = "HC_ACOPF"

        self.model.W = Set(within=self.model.G, initialize=self.wind_hc_set)
        self.model.Wbs = Set(within=self.model.B * self.model.G, initialize=self.bus_whc_set)

        self.model.SWmax = Param(self.model.W, initialize=self.SWmax_data)
        self.model.SWmin = Param(self.model.W, initialize=self.SWmin_data)

        def objective(model):
            return sum(model.pG[w] for w in self.model.W)
        self.model.OBJ = Objective(rule=objective, sense=maximize)

        def slack_bounds(model, b):
            return -200, model.peG[b], 200
        self.model.slack_constraint = Constraint(self.model.eG, rule=slack_bounds)

        def SW_bounds(model, w):
            return model.SWmin[w] ** 2, model.pG[w] ** 2 + model.qG[w] ** 2, model.SWmax[w] ** 2

        self.model.SW_constraint = Constraint(self.model.W, rule=SW_bounds)

        def power_factor(model, w):
            return -sin(acos(0.95)), (model.qG[w] / sqrt(model.pG[w] ** 2 + model.qG[w] ** 2)), sin(acos(0.95))

        self.model.power_factor_constraint = Constraint(self.model.W, rule=power_factor)

        # --- generator power limits ---
        def real_power_bounds(model, g):
            if g in model.Gc:
                return model.PGmin[g], model.pG[g], model.PGmax[g]
            elif g in model.W:
                model.pG[g].unfix()
                model.pG[g] = 1.
                return Constraint.Skip
            else:
                model.pG[g].fix(model.PG[g])
                return Constraint.Skip

        self.model.PG_Constraint = Constraint(self.model.G, rule=real_power_bounds)

        # --- reactive generator power limits ---
        def reactive_power_bounds(model, g):
            if g in model.Gc:
                return model.QGmin[g], model.qG[g], model.QGmax[g]
            elif g in model.W:
                model.qG[g].unfix()
                return Constraint.Skip
            else:
                model.qG[g].fix(model.QG[g])
                return Constraint.Skip

        self.model.QG_Constraint = Constraint(self.model.G, rule=reactive_power_bounds)
