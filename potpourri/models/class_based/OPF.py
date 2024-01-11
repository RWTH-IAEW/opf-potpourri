from pyomo.environ import *
from potpourri.models.class_based.basemodel import Basemodel
import numpy as np
import pandas as pd


class OPF(Basemodel):
    def __init__(self, net):
        super().__init__(net)

        max_load = net.line.max_loading_percent.values if "max_loading_percent" in net.line else 100.
        self.SLmax_data = self.__calc_SLmax(max_load)

        max_load_T = net.trafo.max_loading_percent.values / 100. if "max_loading_percent" in net.trafo else 1.
        sn_mva = net.trafo.sn_mva.values
        df_T = net.trafo.df.values
        self.SLmaxT_data = max_load_T * sn_mva * df_T * net.trafo.parallel.values / self.baseMVA

        self.c0_data = [1] * len(self.gen_all_set)
        self.c1_data = [2] * len(self.gen_all_set)
        self.c2_data = [3] * len(self.gen_all_set)

        self.get_generator_real_power_data()
        self.get_demand_real_power_data()

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
                                  initialize=self.SLmaxT_data[self.model.TRANSF])  # real power transformer limit

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

