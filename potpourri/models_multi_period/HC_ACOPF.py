import copy

import numpy as np
import pandas as pd
from pyomo.environ import *
from potpourri.models.ACOPF_base import ACOPF
import pandapower as pp


class HC_ACOPF(ACOPF):
    def __init__(self, net):
        if 'wind_hc' not in net.sgen:
            net = copy.deepcopy(net)
            buses_excl_extGrids = net.bus.loc[~net.bus.index.isin(net.ext_grid.bus)].index

            pp.create_sgens(net, buses_excl_extGrids, p_mw=0, wind_hc=True)
        net.sgen.wind_hc.fillna(False, inplace=True)

        super().__init__(net)

    def _calc_opf_parameters(self, SWmax=10000, SWmin=0):
        super()._calc_opf_parameters()

        if 'windpot_p_mw' in self.net.bus:
            self.static_generation_data['windpot'] = self.net.bus.windpot_p_mw[self.net.sgen.bus.values].values

        wind_hc_set = np.arange(len(self.net.sgen))[self.net.sgen.wind_hc & self.net.sgen.in_service]
        self.SWmax_data = pd.Series(SWmax / self.baseMVA, wind_hc_set)
        self.SWmin_data = pd.Series(SWmin / self.baseMVA, wind_hc_set)

        self.m_qu_max = (0.48 + 0.23) / (96 - 103) * 110  # Variante 1
        self.qu_max = -self.m_qu_max * 120 / 110 + 0.48
        self.m_qu_min = (0.33 + 0.41) / (96 - 103) * 110  # Variante 3
        self.qu_min = -self.m_qu_min * 96 / 110 + 0.33

    def add_OPF(self, **kwargs):
        super().add_OPF(**kwargs)

        self.model.name = "HC_ACOPF"

        self.model.SWmax = Param(self.model.WIND_HC, initialize=self.SWmax_data[self.model.WIND_HC], mutable=True)
        self.model.SWmin = Param(self.model.WIND_HC, initialize=self.SWmin_data[self.model.WIND_HC], mutable=True)

        if 'windpot_p_mw' in self.net.bus:
            self.model.pWmax = Param(self.model.WIND_HC,
                                     initialize=self.static_generation_data['windpot'][self.model.WIND_HC],
                                     mutable=True)

        self.model.y = Var(self.model.WIND_HC, within=Binary, initialize=1.)

        # def objective(model):
        #     return sum(model.psG[w] for w in model.WIND_HC)
        # self.model.obj_hc = Objective(rule=objective, sense=maximize)
        # self.model.obj_hc.deactivate()

        def obj_wind_loss_rule(model):
            return sum(model.psG[w] for w in model.WIND_HC) - sum(model.pLfrom[l] + model.pLto[l] for l in model.L) - sum(model.pThv[t] + model.pTlv[t] for t in model.TRANSF)

        self.model.obj = Objective(rule=obj_wind_loss_rule, sense=maximize)

        def SW_max(model, w):
            return model.psG[w] ** 2 + model.qsG[w] ** 2 <= model.SWmax[w] ** 2 * model.y[w]

        def SW_min(model, w):
            return model.psG[w] ** 2 + model.qsG[w] ** 2 >= model.SWmin[w] ** 2 * model.y[w]

        self.model.SW_max_constraint = Constraint(self.model.WIND_HC, rule=SW_max)
        self.model.SW_min_constraint = Constraint(self.model.WIND_HC, rule=SW_min)

        # def power_factor(model, w):
        #     return -sin(acos(0.95)), (model.qG[w] / (sqrt(model.pG[w]**2 + model.qG[w]**2) + 1e-3)), sin(acos(0.95))
        #
        # self.model.power_factor_constraint = Constraint(self.model.WIND_HC, rule=power_factor)

        # --- generator power ---
        for w in self.model.WIND_HC:
            self.model.psG[w].unfix()
            self.model.qsG[w].unfix()

        # --- QU Variante 1 ---
        def QW_min(model, w):
            return model.qsG[w] >= -0.41 * model.psG[w]

        def QW_max(model, w):
            return model.qsG[w] <= 0.48 * model.psG[w]

        self.model.QW_min_constraint = Constraint(self.model.WIND_HC, rule=QW_min)
        self.model.QW_max_constraint = Constraint(self.model.WIND_HC, rule=QW_max)

        def QU_min_hc(model, w):
            for (g, b) in model.sGbs:
                if g == w:
                    return model.qsG[w] >= (self.m_qu_min * model.v[b] + self.qu_min) * model.psG[w]

        self.model.QU_min_hc_constraint = Constraint(self.model.WIND_HC, rule=QU_min_hc)

        def QU_max_hc(model, w):
            for (g, b) in model.sGbs:
                if g == w:
                    return model.qsG[w] <= (self.m_qu_max * model.v[b] + self.qu_max) * model.psG[w]

        self.model.QU_max_hc_constraint = Constraint(self.model.WIND_HC, rule=QU_max_hc)

        if 'windpot_p_mw' in self.net.bus:
            def PW_max(model, w):
                return model.psG[w] <= model.pWmax[w]

            self.model.PW_max_constraint = Constraint(self.model.WIND_HC, rule=PW_max)

    def add_loss_obj(self):
        self.model.eps = Param(domain=Reals, initialize=1., mutable=True)

        def objective_pwind_loss(model):
            return model.eps * sum(model.psG[w] for w in model.WIND_HC) + (1 - model.eps) * (
                - sum(model.pLfrom[l] + model.pLto[l] for l in model.L))

        self.model.obj_hc.deactivate()
        self.model.OBJ_with_loss = Objective(rule=objective_pwind_loss, sense=maximize)
