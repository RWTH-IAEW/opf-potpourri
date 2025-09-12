import pandas as pd
from pyomo.environ import *
from math import pi
import copy
import numpy as np
import pandapower as pp
from potpourri.models_multi_period.flexibility_multi_period import Flexibility_multi_period
from potpourri.models_multi_period.pyo_to_net_multi_period import pyo_sol_to_net_res
import simbench as sb

class PV_multi_period(Flexibility_multi_period):
    """
    **PV** \n
    This class creates the sets, parameters, variables and constraints for the PV generation units
    percentages taken from paper: Anforderungen an aktuelle Verteilnetze und deren zukuenftige Versorgungsaufgabe
    """

    def __init__(self,net, T=None, scenario=None):
        super().__init__(net, T, scenario)
        self.net = net

        if scenario == 0:
            self.pv_percentage = 13.4
        elif scenario == 1:
            self.pv_percentage = 22.4
        elif scenario == 2:
            self.pv_percentage = 24.4
        elif scenario == 3:
            self.pv_percentage = 25.4

        # take a random columnn of input sb net renewables as pv load profile
        self.pv_load_profiles = (-1)*(self.net.pv_load_profiles)
        self.pv_load_profile = self.pv_load_profiles["PV5"]

        #self.pv_pmax = -0.01447 #MW
        self.pv_pmax = self.pv_load_profile
        self.pv_pmin = 0.0
        num_indexes = round(len(self.buses_excl_extGrids) * self.pv_percentage / 100)
        self.random_indexes = np.random.choice(self.buses_excl_extGrids, num_indexes, replace=False)

       # pp.create_sgens(net, self.random_indexes, p_mw=self.pv_pmax, pv=True, type='PV_own', name=num_indexes)
        #net.sgen.pv.fillna(False, inplace=True)

    def get_all(self, model):
        """"
        **get_all** \n
        This function gets the sets, parameters and variables of the class, unfixes variables
        """
        self.get_sets(model)
        self.get_parameters(model)
        self.get_variables(model)
        self.get_all_constraints(model)
        self.unfix_variables(model)

    def unfix_variables(self, model):
        for t in model.T:
            for pv in model.PV:
                model.pPV[pv, t].unfix()

    def get_sets(self, model):
        super().get_sets(model)
        model.PV = Set(initialize=list(range(len(self.random_indexes))))
        model.PV_bus = Set(initialize=list(enumerate(self.random_indexes)))
        return True

    def get_parameters(self, model):
        # PV Parameters
        self.PV_Pmax_dict, self.PV_Pmax_tuple = self.make_to_dict(model.PV, model.T, self.pv_pmax)
        self.PV_Pmin_dict, self.PV_Pmin_tuple = self.make_to_dict(model.PV, model.T, self.pv_pmin, False)

        model.PV_Pmax = Param(self.PV_Pmax_tuple, within=Reals, initialize=self.PV_Pmax_dict)
        model.PV_Pmin = Param(self.PV_Pmin_tuple, within=Reals, initialize=self.PV_Pmin_dict)
        # model.PV_Qmax = Param(model.PV, within=NegativeReals, initialize=0)
        # model.PV_Qmin = Param(model.PV, within=NegativeReals, initialize=0)

    def get_variables(self, model):
        # PV Variables
        self.pPV_data_dict, self.pPV_tuple = self.make_to_dict(model.PV, model.T, self.pv_load_profile)
        model.pPV = Var(self.pPV_tuple, within=Reals, initialize=self.pPV_data_dict)

    def get_all_constraints(self, model):
        @model.Constraint(model.PV, model.T)
        def PV_real_power_bounds(model, pv, t):
            return model.PV_Pmax[pv,t], model.pPV[pv, t], model.PV_Pmin[pv,t]

    #  empty placeholders, because they are called in get_all
    def get_all_acopf(self, model):
        pass
    def get_all_ac(self, model):
        pass
    def get_all_opf(self, model):
        pass
