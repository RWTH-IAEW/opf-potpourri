from pyomo.environ import *
from potpourri.models.class_based.basemodel import Basemodel
import numpy as np


class OPF(Basemodel):
    def __init__(self, net):
        super().__init__(net)

        max_load = net.line.max_loading_percent.values / 100. if "max_loading_percent" in net.line else 1.
        vr = net.bus.loc[net.line["from_bus"].values, "vn_kv"].values * np.sqrt(3.)
        max_i_ka = net.line.max_i_ka.values
        df = net.line.df.values
        self.SLmax_data = max_load * max_i_ka * df * net.line.parallel.values * vr / self.baseMVA

        max_load_T = net.trafo.max_loading_percent.values / 100. if "max_loading_percent" in net.trafo else 1.
        sn_mva = net.trafo.sn_mva.values
        df_T = net.trafo.df.values
        self.SLmaxT_data = max_load_T * sn_mva * df_T * net.trafo.parallel.values / self.baseMVA

        self.c0_data = [1] * len(self.sgen_set)
        self.c1_data = [2] * len(self.sgen_set)
        self.c2_data = [3] * len(self.sgen_set)

    def create_model(self):
        super().create_model()

        # lines and transformer chracteristics and ratings
        self.model.SLmax = Param(self.model.L, within=NonNegativeReals,
                                 initialize=self.SLmax_data)  # real power line limit
        self.model.SLmaxT = Param(self.model.TRANSF, within=NonNegativeReals,
                                  initialize=self.SLmaxT_data)  # real power transformer limit

        # cost data
        self.model.c2 = Param(self.model.G, within=NonNegativeReals,
                              initialize=self.c2_data)  # generator cost coefficient c2 (*pG^2)
        self.model.c1 = Param(self.model.G, within=NonNegativeReals,
                              initialize=self.c1_data)  # generator cost coefficient c1 (*pG)
        self.model.c0 = Param(self.model.G, within=NonNegativeReals,
                              initialize=self.c0_data)  # generator cost coefficient c0
