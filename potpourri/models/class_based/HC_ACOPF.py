from pyomo.environ import *
from potpourri.models.class_based.ACOPF_base import ACOPF
import pandapower as pp

class HC_ACOPF(ACOPF):
    def __init__(self, net, SWmax=150, SWmin=10, peGmax=10000):
        if 'wind_hc' not in net.sgen:
            buses_excl_extGrids = net.bus.loc[~net.bus.index.isin(net.ext_grid.bus)].index

            pp.create_sgens(net, buses_excl_extGrids, p_mw=0, wind_hc=True)
        net.sgen.wind_hc.fillna(False, inplace=True)

        super().__init__(net)

        net.sgen.wind_hc.fillna(False, inplace=True)

        self.SWmax_data = SWmax / self.baseMVA
        self.SWmin_data = SWmin / self.baseMVA

        self.peGmax_data = peGmax / self.baseMVA

        self.wind_hc_set = net.sgen.index[net.sgen.wind_hc]
        self.bus_whc_set = list(zip(net.sgen.bus[net.sgen.wind_hc], self.wind_hc_set))

        self.create_model()

    def add_OPF(self):
        super().add_OPF()

        self.model.name = "HC_ACOPF"

        self.model.WIND = Set(within=self.model.G, initialize=self.wind_hc_set)
        self.model.Wbs = Set(within=self.model.B * self.model.G, initialize=self.bus_whc_set)

        self.model.SWmax = Param(self.model.WIND, initialize=self.SWmax_data, mutable=True)
        self.model.SWmin = Param(self.model.WIND, initialize=self.SWmin_data, mutable=True)

        self.model.peGmax = Param(self.model.eG, initialize=self.peGmax_data, mutable=True)

        self.model.y = Var(self.model.WIND, within=Binary)

        def objective(model):
            return sum(model.pG[w] for w in model.WIND)

        self.model.OBJ = Objective(rule=objective, sense=maximize)

        def slack_p_bounds(model, b):
            return -model.peGmax[b], model.peG[b], model.peGmax[b]

        self.model.slack_p_min_max_constr = Constraint(self.model.eG, rule=slack_p_bounds)

        # def SW_bounds(model, w):
        #     return model.SWmin[w] ** 2, model.pG[w] ** 2 + model.qG[w] ** 2, model.SWmax[w] ** 2
        #
        # self.model.SW_constraint = Constraint(self.model.WIND, rule=SW_bounds)

        def SW_max(model, w):
            return model.pG[w] ** 2 + model.qG[w] ** 2 <= model.SWmax[w] ** 2 * model.y[w]

        def SW_min(model, w):
            return model.pG[w] ** 2 + model.qG[w] ** 2 >= model.SWmin[w] ** 2 * model.y[w]

        self.model.SW_max_constraint = Constraint(self.model.WIND, rule=SW_max)
        self.model.SW_min_constraint = Constraint(self.model.WIND, rule=SW_min)

        # def power_factor(model, w):
        #     return -sin(acos(0.95)), (model.qG[w] / (sqrt(model.pG[w]**2 + model.qG[w]**2) + 1e-3)), sin(acos(0.95))
        #
        # self.model.power_factor_constraint = Constraint(self.model.WIND, rule=power_factor)

        # --- generator power ---
        for w in self.model.WIND:
            self.model.pG[w].unfix()
            self.model.pG[w] = 10.
            self.model.qG[w].unfix()

        # --- PQ VED Variante 1 ---
        def PW_min_max(model, w):
            return model.SWmax[w] * 0.1, model.pG[w], model.SWmax[w]

        def QW_pos(model, w):
            return model.qG[w] <= -0.3 * model.SWmax[w] + 3.8 * model.pG[w]

        def QW_neg(model, w):
            return model.qG[w] >= -1.3 * model.pG[w]

        def QW_min_max(model, w):
            return -model.SWmax[w] * 0.23, model.qG[w], model.SWmax[w] * 0.48

        def QW_min(model, w):
            return model.qG[w] >= -0.23 * model.pG[w]

        def QW_max(model, w):
            return model.qG[w] <= 0.48 * model.pG[w]

#        self.model.PW_min_max_constraint = Constraint(self.model.WIND, rule=PW_min_max)
#        self.model.QW_pos_constraint = Constraint(self.model.WIND, rule=QW_pos)
#        self.model.QW_neg_constraint = Constraint(self.model.WIND, rule=QW_neg)
#        self.model.QW_min_max_constraint = Constraint(self.model.WIND, rule=QW_min_max)
        self.model.QW_min_constraint = Constraint(self.model.WIND, rule=QW_min)
        self.model.QW_max_constraint = Constraint(self.model.WIND, rule=QW_max)

    def fix_y(self, value=1.):
        for w in self.model.WIND:
            self.model.y[w].fix(value)

    def unfix_y(self, value=1.):
        for w in self.model.WIND:
            self.model.y[w].unfix()
            self.model.y[w] = value

