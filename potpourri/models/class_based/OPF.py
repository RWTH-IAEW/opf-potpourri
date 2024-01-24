from pyomo.environ import *
from potpourri.models.class_based.basemodel import Basemodel
import numpy as np
import pandas as pd


class OPF(Basemodel):
    def __init__(self, net):
        super().__init__(net)

        max_load = net.line.max_loading_percent.values if "max_loading_percent" in net.line else 100.
        self.SLmax_data = self.__calc_SLmax(max_load)

        max_load_T = self.net.trafo.max_loading_percent.fillna(
            100.) / 100. if "max_loading_percent" in net.trafo else 1.
        sn_mva = self.net.trafo.sn_mva
        df_T = self.net.trafo.df
        self.SLmaxT_data = max_load_T * sn_mva * df_T * self.net.trafo.parallel / self.baseMVA

        self.c0_data = [1] * len(self.gen_all_set)
        self.c1_data = [2] * len(self.gen_all_set)
        self.c2_data = [3] * len(self.gen_all_set)

        self.get_generator_real_power_data()
        self.get_demand_real_power_data()

        self.tap_min_data = pd.Series(self.trafo_parameters["ratio"][1].round(15), self.trafo_set)
        self.tap_max_data = pd.Series(self.trafo_parameters["ratio"][2].round(15), self.trafo_set)

        self.tap_neutral = self.net.trafo.tap_neutral
        self.tap_step = self.net.trafo.tap_step_percent / 100.
        self.tap_pos_max = self.net.trafo.tap_max
        self.tap_pos_min = self.net.trafo.tap_min
        self.tap_side_data = pd.Series(np.where(self.net.trafo.tap_side == 'lv', 1, 0),
                                       self.trafo_set)  # tap side data = 1 if tap side is lv, 0 if tap side is hv

    def get_generator_real_power_data(self):
        if 'controllable' not in self.generators:
            self.gen_controllable_set = pd.Index([])  # create empty Set if no controllable generators exist
        else:
            self.gen_controllable_set = self.generators.index[self.generators.controllable == True]

        # add rows with active generation limits if not existing
        if 'max_p_mw' not in self.generators:
            self.generators['max_p_mw'] = self.generators.p_mw

        if 'min_p_mw' not in self.net.sgen:
            self.generators['min_p_mw'] = [0] * len(self.generators.index)

        # generation limits for sgens
        self.PGmax_data = self.generators.max_p_mw.fillna(self.generators.p_mw) / self.baseMVA
        self.PGmin_data = self.generators.min_p_mw.fillna(0) / self.baseMVA

    def get_demand_real_power_data(self):
        if 'controllable' not in self.net.load:
            self.demand_controllable_set = None  # create empty Set if no controllable load exist
        else:
            self.demand_controllable_set = self.net.load.index[self.net.load.controllable == True]

        # add rows with active demand limits if not existing
        if 'max_p_mw' not in self.net.load:
            self.net.load['max_p_mw'] = self.net.load.p_mw

        if 'min_p_mw' not in self.net.load:
            self.net.load['min_p_mw'] = [0] * len(self.net.load.index)

        # demand limits for loads
        self.PDmax_data = self.net.load.max_p_mw.fillna(self.net.load.p_mw) / self.baseMVA
        self.PDmin_data = self.net.load.min_p_mw.fillna(0) / self.baseMVA

    def __calc_SLmax(self, max_loading_percent=100):
        vr = self.net.bus.loc[self.net.line["from_bus"].values, "vn_kv"].values * np.sqrt(3.)
        max_i_ka = self.net.line.max_i_ka.values
        df = self.net.line.df.values
        return max_loading_percent / 100. * max_i_ka * df * self.net.line.parallel.values * vr / self.baseMVA

    # def _calc_tap_range(self):
    #     # only longitudinal regulators considered regarding control, no phase shifters, no cross regulators
    #
    #     tap_pos_min = self.net.trafo.tap_min
    #     tap_pos_max = self.net.trafo.tap_max
    #     tap_neutral = self.net.trafo.tap_neutral
    #     tap_step = self.net.trafo.tap_step_percent / 100.
    #
    #     n_tap_min = 1 + (tap_pos_min - tap_neutral) * tap_step
    #     n_tap_max = 1 + (tap_pos_max - tap_neutral) * tap_step
    #
    #     tap_range = []
    #     for i in self.net.trafo.index:
    #         tap_range.append(np.arange(n_tap_min[i], n_tap_max[i], tap_step[i]))
    #
    #     trafos_lv_index = self.net.trafo.index[self.net.trafo.tap_side == 'lv']
    #     for i in trafos_lv_index:
    #         tap_range[i] = 1 / tap_range[i]
    #
    #     return tap_range
    #
    #     # trafos_step_degree = self.net.trafo.index[(self.net.trafo.tap_step_degree == 0) & (self.net.trafo.tap_side == 'hv')]
    #     #
    #     # trafos_lv_degree = self.net.trafo.index[~trafos_hv_nodegree]
    #     # for i in trafos_lv_degree:
    #     #     vn_trafo_hv, vn_trafo_lv, shift = self._calc_tap_shift(tap_pos=i)
    #     #     ratio = self._calc_nominal_ratio_from_dataframe(vn_trafo_hv, vn_trafo_lv)
    #     #     tap_range[i]= ([i, ratio])
    #     #     shift_range.append([i, shift])

    def set_SLmax(self, max_loading, unit: ['percent', 'MW']):
        if unit == 'percent':
            SLmax = max_loading / 100. * self.SLmax_data
        elif unit == 'MW':
            self.SLmax_data = max_loading * np.ones(len(self.model.L)) / self.baseMVA
            SLmax = self.SLmax_data
        for l in self.model.L:
            self.model.SLmax[l] = SLmax[l]

    def add_OPF(self):
        # controllable generation
        self.model.Gc = Set(within=self.model.G, initialize=self.gen_controllable_set)  # controllable generators
        self.model.PGmax = Param(self.model.G, initialize=self.PGmax_data[self.model.G])
        self.model.PGmin = Param(self.model.G, initialize=self.PGmin_data[self.model.G])

        # controllable loads
        self.model.Dc = Set(within=self.model.D, initialize=self.demand_controllable_set)  # controllable loads
        self.model.PDmax = Param(self.model.D, initialize=self.PDmax_data[self.model.D])
        self.model.PDmin = Param(self.model.D, initialize=self.PDmin_data[self.model.D])

        # lines and transformer chracteristics and ratings
        self.model.SLmax = Param(self.model.L, within=NonNegativeReals,
                                 initialize=self.SLmax_data[self.model.L], mutable=True)  # real power line limit
        self.model.SLmaxT = Param(self.model.TRANSF, within=NonNegativeReals,
                                  initialize=self.SLmaxT_data[self.model.TRANSF],
                                  mutable=True)  # real power transformer limit

        # cost data
        self.model.c2 = Param(self.model.G, within=NonNegativeReals,
                              initialize=self.c2_data)  # generator cost coefficient c2 (*pG^2)
        self.model.c1 = Param(self.model.G, within=NonNegativeReals,
                              initialize=self.c1_data)  # generator cost coefficient c1 (*pG)
        self.model.c0 = Param(self.model.G, within=NonNegativeReals,
                              initialize=self.c0_data)  # generator cost coefficient c0
        self.model.VOLL = Param(self.model.D, within=Reals, initialize=10000)  # value of lost load

        # --- generator power limits ---
        def real_power_bounds(model, g):
            model.pG[g].unfix()
            return model.PGmin[g], model.pG[g], model.PGmax[g]

        self.model.PGc_Constraint = Constraint(self.model.Gc, rule=real_power_bounds)

        # --- demand limits ---
        def real_demand_bounds(model, d):
            return model.PDmin[d], model.pD[d], model.PDmax[d]

        self.model.PD_Constraint = Constraint(self.model.Dc, rule=real_demand_bounds)

        # --- transformer tap ratio limits ---


    def add_tap_changer_linear(self):
        self.model.Tap_min = Param(self.model.TRANSF, within=Reals, initialize=self.tap_min_data[self.model.TRANSF])
        self.model.Tap_max = Param(self.model.TRANSF, within=Reals, initialize=self.tap_max_data[self.model.TRANSF])

        def trafo_tap_linear_bounds(model, t):
            return model.Tap_min[t], model.Tap[t], model.Tap_max[t]

        self.model.Tap_linear_constr = Constraint(self.model.TRANSF, rule=trafo_tap_linear_bounds)

        self.unfix_vars('Tap')

    def add_tap_changer_discrete(self):
        self.model.Tap_pos = Var(self.model.TRANSF, within=Integers, initialize=0.)  # transformer tap position
        self.model.Tap_pos_min = Param(self.model.TRANSF, within=Integers, initialize=self.tap_pos_min)
        self.model.Tap_pos_max = Param(self.model.TRANSF, within=Integers, initialize=self.tap_pos_max)
        self.model.Tap_neutral = Param(self.model.TRANSF, within=Integers,
                                       initialize=self.tap_neutral)  # transformer tap neutral position
        self.model.Tap_step = Param(self.model.TRANSF, within=Reals,
                                    initialize=self.tap_step)  # transformer tap step size
        self.model.Tap_side = Param(self.model.TRANSF, initialize=self.tap_side_data[self.model.TRANSF])

        def trafo_tap_pos_min_max(model, t):
            return model.Tap_pos_min[t], model.Tap_pos[t], model.Tap_pos_max[t]

        self.model.Tap_pos_constr = Constraint(self.model.TRANSF, rule=trafo_tap_pos_min_max)

        def trafo_tap_discrete(model, t):
            if model.Tap_side[t]:
                # tap side: lv
                return model.Tap[t] == 1 / (1 + (model.Tap_pos[t] - model.Tap_neutral[t]) * model.Tap_step[t])

            return model.Tap[t] == 1. + (model.Tap_pos[t] - model.Tap_neutral[t]) * model.Tap_step[t]

        self.model.Tap_discrete_constr = Constraint(self.model.TRANSF, rule=trafo_tap_discrete)

        self.unfix_vars('Tap')
