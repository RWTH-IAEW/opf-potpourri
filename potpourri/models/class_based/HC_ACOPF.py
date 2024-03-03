import numpy as np
import pandas as pd
from pyomo.environ import *
from potpourri.models.class_based.ACOPF_base import ACOPF
import pandapower as pp

class HC_ACOPF(ACOPF):
    def __init__(self, net):
        if 'wind_hc' not in net.sgen:
            buses_excl_extGrids = net.bus.loc[~net.bus.index.isin(net.ext_grid.bus)].index

            pp.create_sgens(net, buses_excl_extGrids, p_mw=0, wind_hc=True)
        net.sgen.wind_hc.fillna(False, inplace=True)

        super().__init__(net)

    def _calc_opf_parameters(self, SWmax=10000, SWmin=0):
        super()._calc_opf_parameters()

        wind_sgen_index = self.net.sgen.index[self.net.sgen.wind_hc & self.net.sgen.in_service]
        self.wind_hc_set = np.arange(len(self.net.sgen))[self.net.sgen.wind_hc & self.net.sgen.in_service]
        self.static_generation_data['wind_hc'] = self.net.sgen.wind_hc

        if 'windpot_p_mw' in self.net.bus:
            self.static_generation_data['windpot'] = self.net.bus.windpot_p_mw[self.net.sgen.bus.values].values
            # self.windpot = pd.Series(self.net.bus.windpot_p_mw[self.net.sgen.bus[self.net.sgen.wind_hc & self.net.sgen.in_service].values].values, self.wind_hc_set)

        self.SWmax_data = pd.Series(SWmax / self.baseMVA, self.wind_hc_set)
        self.SWmin_data = pd.Series(SWmin / self.baseMVA, self.wind_hc_set)

        self.m_qu_max = (0.48 + 0.23) / (96 - 103) * 110    # Variante 1
        self.qu_max = -self.m_qu_max * 120 / 110 + 0.48
        qu_min_1 = (0.33, 96/110)
        qu_min_2 = (-0.41, 103/110)
        m_qu_min1 = (qu_min_1[0] - qu_min_2[0]) / (qu_min_1[1] - qu_min_2[1])
        b_qu_min1 = qu_min_1[0] - m_qu_min1 * qu_min_1[1]
        self.m_qu_min = (0.33 + 0.41) / (96 - 103) * 110    # Variante 3
        self.qu_min = -self.m_qu_min * 96 / 110 + 0.33

        x = pd.DataFrame([[96, 103], [120, 127]], index=['x1', 'x2'], columns=['min', 'max'])/110
        y = pd.DataFrame([[0.48, 0.41, 0.33], [-0.23, -0.33, -0.41]], index=['y1', 'y2'], columns=['v1', 'v2', 'v3'])


        def get_lin_constr(p1, p2):
            m = (p1[0] - p2[0]) / (p1[1] - p2[1])
            b = p1[0] - m * p1[1]
            return m, b


    def add_OPF(self, **kwargs):
        super().add_OPF(**kwargs)

        self.model.name = "HC_ACOPF"

        self.model.WIND = Set(within=self.model.sG, initialize=self.static_generation_data.index[self.static_generation_data['wind_hc'] & self.static_generation_data.in_service])

        self.model.SWmax = Param(self.model.WIND, initialize=self.SWmax_data[self.model.WIND], mutable=True)
        self.model.SWmin = Param(self.model.WIND, initialize=self.SWmin_data[self.model.WIND], mutable=True)

        if 'windpot_p_mw' in self.net.bus:
            self.model.pWmax = Param(self.model.WIND, initialize=self.static_generation_data['windpot'][self.model.WIND], mutable=True)

        self.model.y = Var(self.model.WIND, within=Binary, initialize=1.)

        def objective(model):
            return sum(model.psG[w] for w in model.WIND)

        self.model.OBJ = Objective(rule=objective, sense=maximize)

        # def SW_bounds(model, w):
        #     return model.SWmin[w] ** 2, model.pG[w] ** 2 + model.qG[w] ** 2, model.SWmax[w] ** 2
        #
        # self.model.SW_constraint = Constraint(self.model.WIND, rule=SW_bounds)

        def SW_max(model, w):
            return model.psG[w] ** 2 + model.qsG[w] ** 2 <= model.SWmax[w] ** 2 * model.y[w]

        def SW_min(model, w):
            return model.psG[w] ** 2 + model.qsG[w] ** 2 >= model.SWmin[w] ** 2 * model.y[w]

        self.model.SW_max_constraint = Constraint(self.model.WIND, rule=SW_max)
        self.model.SW_min_constraint = Constraint(self.model.WIND, rule=SW_min)

        # def power_factor(model, w):
        #     return -sin(acos(0.95)), (model.qG[w] / (sqrt(model.pG[w]**2 + model.qG[w]**2) + 1e-3)), sin(acos(0.95))
        #
        # self.model.power_factor_constraint = Constraint(self.model.WIND, rule=power_factor)

        # --- generator power ---
        for w in self.model.WIND:
            self.model.psG[w].unfix()
            self.model.qsG[w].unfix()

        # --- PQ VED Variante 1 ---
        def PW_min_max(model, w):
            return model.SWmax[w] * 0.1, model.psG[w], model.SWmax[w]

        def QW_pos(model, w):
            return model.qsG[w] <= -0.3 * model.SWmax[w] + 3.8 * model.psG[w]

        def QW_neg(model, w):
            return model.qsG[w] >= -1.3 * model.psG[w]

        def QW_min_max(model, w):
            return -model.SWmax[w] * 0.23, model.qsG[w], model.SWmax[w] * 0.48

        def QW_min(model, w):
            return model.qsG[w] >= -0.41 * model.psG[w]

        def QW_max(model, w):
            return model.qsG[w] <= 0.48 * model.psG[w]

#        self.model.PW_min_max_constraint = Constraint(self.model.WIND, rule=PW_min_max)
#        self.model.QW_pos_constraint = Constraint(self.model.WIND, rule=QW_pos)
#        self.model.QW_neg_constraint = Constraint(self.model.WIND, rule=QW_neg)
#        self.model.QW_min_max_constraint = Constraint(self.model.WIND, rule=QW_min_max)
        self.model.QW_min_constraint = Constraint(self.model.WIND, rule=QW_min)
        self.model.QW_max_constraint = Constraint(self.model.WIND, rule=QW_max)

        def QU_min(model, w):
            for (g, b) in model.sGbs:
                if g == w:
                    return model.qsG[w] >= (self.m_qu_min * model.v[b] + self.qu_min) * model.psG[w]
        self.model.QU_min_constraint = Constraint(self.model.WIND, rule=QU_min)

        def QU_max(model, w):
            for (g, b) in model.sGbs:
                if g == w:
                    return model.qsG[w] <= (self.m_qu_max * model.v[b] + self.qu_max) * model.psG[w]
        self.model.QU_max_constraint = Constraint(self.model.WIND, rule=QU_max)

        if 'windpot_p_mw' in self.net.bus:
            def PW_max(model, w):
                return model.psG[w] <= model.pWmax[w]
            self.model.PW_max_constraint = Constraint(self.model.WIND, rule=PW_max)

    def add_loss_obj(self):
        self.model.eps = Param(domain=Reals, initialize=1., mutable=True)

        def objective_pwind_loss(model):
            return model.eps * sum(model.psG[w] for w in model.WIND) + (1-model.eps)*(- sum(model.pLfrom[l] + model.pLto[l] for l in model.L))

        self.model.OBJ.deactivate()
        self.model.OBJ_with_loss = Objective(rule=objective_pwind_loss, sense=maximize)

