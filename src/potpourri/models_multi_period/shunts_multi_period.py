import pandas as pd
from pyomo.environ import *
from math import pi
import copy
import numpy as np
import pandapower as pp
from src.potpourri.models_multi_period.flexibility_multi_period import Flexibility_multi_period
from src.potpourri.models_multi_period.pyo_to_net_multi_period import pyo_sol_to_net_res
import simbench as sb
import os


class Shunts_multi_period(Flexibility_multi_period):

    def __init__(self,net, T=None, scenario=None):
        super().__init__(net, T, scenario)

        self.shunt_set = self.net.shunt.index[self.net.shunt.in_service]
        self.bus_shunt_set = list(zip(self.bus_lookup[self.net.shunt.bus[self.shunt_set].values], self.shunt_set))

    def get_all(self, model):
        """"
        **get_all** \n
        This function gets the sets, parameters and variables of the class
        """
        self.get_sets(model)
        self.get_parameters(model)

    def get_all_opf(self, model):
        """"
        **get_all** \n
        This function gets the sets, parameters and variables of the class
        """
    def get_all_acopf(self, model):
        pass

    def get_sets(self, model):
        super().get_sets(model)
        model.SHUNT = Set(initialize=self.shunt_set) # set of shunts
        # loads linked to each bus b
        model.SHUNTbs = Set(within=model.B * model.SHUNT,
                                 initialize=self.bus_shunt_set)  # set of shunt-bus mapping
        return True

    def get_parameters(self, model):
        # --- Parameters ---
        model.GB = Param(model.SHUNT, within=Reals,
                              initialize=self.GB_data[model.SHUNT])  # shunt conductance
        return True



