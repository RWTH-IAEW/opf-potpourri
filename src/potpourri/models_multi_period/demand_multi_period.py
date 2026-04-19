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


class Demand_multi_period(Flexibility_multi_period):

    def __init__(self,net, T=None, scenario=None):
        super().__init__(net, T, scenario)

        self.demand_set = self.net.load.index[self.net.load.in_service]
        self.bus_demand_set = list(zip(self.bus_lookup[self.net.load.bus[self.demand_set].values], self.demand_set))


    def get_all(self, model):
        """"
        **get_all** \n
        This function gets the sets, parameters and variables of the class
        """
        self.get_sets(model)
        self.get_parameters(model)
        self.get_variables(model)
        self.fix_variables(model)


    def get_all_opf(self, model):
        """"
        **get_all** \n
        This function gets the sets, parameters and variables of the class
        """
        self.get_opf_sets(model)
        self.get_opf_parameters(model)
        self.get_all_constraints_opf(model)

    def get_all_acopf(self, model):
        self.get_acopf_parameters(model)
        self.get_all_constraints_acopf(model)

    def get_all_ac(self, model):
        self.get_ac_parameters(model)
        self.get_ac_variables(model)

    def get_sets(self, model):
        super().get_sets(model)
        model.D = Set(initialize=self.demand_set)  #set of demands
        # loads linked to each bus b
        model.Dbs = Set(within=model.B * model.D,
                             initialize=self.bus_demand_set)  # set of demand-bus mapping
        return True

    def get_parameters(self, model):
        # make PD_data_dict index-able through the tuple list
        self.PD_data_dict, self.PD_tuple = self.make_to_dict(model.D, model.T, self.PD_data)
        # demand at each bus
        model.PD = Param(self.PD_tuple, initialize=self.PD_data_dict)
        return True

    def get_ac_parameters(self, model):
        # reactive demand
        self.QD_data_dict, self.QD_tuple = self.make_to_dict(model.D, model.T, self.QD_data)
        model.QD = Param(self.QD_tuple, initialize=self.QD_data_dict)
        return True

    def get_ac_variables(self, model):
        # reactive demand
        model.qD = Var(self.QD_tuple, domain=Reals)
        return True

    def get_variables(self, model):
        # --- Variables ---
        # demand at each bus
        model.pD = Var(self.PD_tuple, domain=Reals)  # real power demand delivered
        return True

    def fix_variables(self, model):
        # fix the demand values over all d in D and t in T
        for d_t in self.PD_tuple:
            model.pD[d_t].fix(model.PD[d_t])
        # for d in model.D for t in model.T:
        #     model.pD[(d,t)].fix(model.PD[(d,t)])

    # def unfix_variables(self, model):
    #     # unfix the static generation values
    #     for d in model.D:
    #         for t in model.T:
    #             model.pD[(d, t)].unfix()


    def get_opf_sets(self, model):
        model.Dc = Set(within=model.D, initialize=self.demand_controllable_set)  # controllable loads

    def get_opf_parameters(self, model):
        # real demand
        model.PDmax = Param(self.PDmax_tuple, initialize=self.PDmax_data_dict)
        model.PDmin = Param(self.PDmin_tuple, initialize=self.PDmin_data_dict)

    def get_acopf_parameters(self, model):
        # reactive demand
        model.QDmax = Param( self.QDmax_tuple, initialize=self.QDmax_data_dict)
        model.QDmin = Param(self.QDmin_tuple, initialize=self.QDmin_data_dict)

    def get_all_constraints_opf(self, model):
        @model.Constraint(model.Dc, model.T)
        def real_demand_bounds(model, d, t):
            model.pD[(d,t)].unfix()
            return model.PDmin[(d,t)], model.pD[(d,t)], model.PDmax[(d,t)]

    def get_all_constraints_acopf(self, model):
        @model.Constraint(model.Dc, model.T)
        def reactive_demand_bounds(model, d, t):
            model.qD[(d,t)].unfix()
            return model.QDmin[(d,t)], model.qD[(d,t)], model.QDmax[(d,t)]

    def get_demand_real_power_data(self, model, max_p_mw=None, min_p_mw=None):

        if 'controllable' not in self.net.load:
             self.demand_controllable_set = None  # create empty Set if no controllable load exist
        else:
             self.demand_controllable_set = self.net.load.index[self.net.load.controllable == True]

        # Aus dem Netzwerk die maximalen und minimalen Lasten holen
        if max_p_mw is not None:
            self.PDmax_data_dict, self.PDmax_tuple = self.make_to_dict(model.D, model.T, max_p_mw)
        elif max_p_mw is None:
            self.PDmax_data_dict, self.PDmax_tuple = self.make_to_dict(model.D, model.T, self.PD_data)

         # add rows with active demand limits if not existing
        if min_p_mw is not None:
            self.PDmin_data_dict, self.PDmin_tuple = self.make_to_dict(model.D, model.T, min_p_mw, False)
        elif min_p_mw is None:
            self.PDmin_data_dict, self.PDmin_tuple = self.make_to_dict(model.D, model.T, 0, False)

    def get_demand_reactive_data(self, model, max_q_mvar=None, min_q_mvar=None):
        # reactive power demand
        if max_q_mvar is not None:
            self.QDmax_data_dict, self.QDmax_tuple = self.make_to_dict(model.D, model.T, max_q_mvar)
        elif max_q_mvar is None:
            self.QDmax_data_dict, self.QDmax_tuple = self.make_to_dict(model.D, model.T, self.QD_data.abs())

        if min_q_mvar is not None:
            self.QDmin_data_dict, self.QDmin_tuple = self.make_to_dict(model.D, model.T, min_q_mvar)
        elif min_q_mvar is None:
            self.QDmin_data_dict, self.QDmin_tuple = self.make_to_dict(model.D, model.T,-(self.QD_data.abs()))
        # demand limits for loads