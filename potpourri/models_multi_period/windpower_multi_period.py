import pandas as pd
from pyomo.environ import *
from math import pi
import copy
import numpy as np
import pandapower as pp
from potpourri.models_multi_period.sgens_multi_period import Sgens_multi_period
from potpourri.models_multi_period.pyo_to_net_multi_period import pyo_sol_to_net_res
import simbench as sb
import os
import logging
import pandapower as pp

#TODO: check if properly made multiperiod
class Windpower_multi_period(Sgens_multi_period):
    def __init__(self, net, T=None, scenario=None):
        super().__init__(net, T, scenario)
        #taken from HC_ACOPF, formerly used in _calc_opf_parameters
        if 'windpot_p_mw' in net.bus:
            self.static_generation_data['windpot'] = net.bus.windpot_p_mw[net.sgen.bus.values].values
            #taken from ACOPF_base, formerly used in get static_generation_reactive_power_limits()
            self.static_generation_data['type'] = net.sgen.type.values

    def get_all(self, model):
        """"
        **get_all** \n
        This function gets the sets, parameters and variables of the class
        """

    def get_all_opf(self, model):
        """"
        **get_all** \n
        This function gets the sets, parameters and variables of the class
        """
        self.get_opf_sets(model)
        self.get_opf_parameters(model)

    def get_opf_sets(self, model):
        # --- sets ---
        # generators for hc calculation
        model.WIND_HC = Set(within=model.sG, initialize=self.static_generation_data.index[
            self.static_generation_data['wind_hc'] & self.static_generation_data.in_service])
        # all wind generators
        model.WIND = model.WIND_HC | Set( within=model.sG, initialize=self.static_generation_data.index[(self.static_generation_data['type'] == 'Wind') & self.static_generation_data.in_service])
        # controllable wind generators, not for hc calculation
        model.WINDc = model.WIND & model.sGc & Set(initialize=self.static_generation_data.index[self.static_generation_data['var_q'].values != None])
        return True

    def _calc_wind_opf_parameters(self, model, SWmax = 10000, SWmin = 0):
        if 'windpot_p_mw' in self.net.bus:
            self.static_generation_data['windpot'] = self.net.bus.windpot_p_mw[self.net.sgen.bus.values].values

        wind_hc_set = np.arange(len(self.net.sgen))[self.net.sgen.wind_hc & self.net.sgen.in_service]
        self.SWmax_data = pd.Series(SWmax / self.baseMVA, wind_hc_set)
        self.SWmax_data_dict, self.SWmax_tuple = self.make_to_dict(model.WIND_HC, model.T, self.SWmax_data)
        self.SWmin_data = pd.Series(SWmin / self.baseMVA, wind_hc_set)
        self.SWmin_data_dict, self.SWmin_tuple = self.make_to_dict(model.WIND_HC, model.T, self.SWmin_data)

        self.m_qu_max = (0.48 + 0.23) / (96 - 103) * 110  # Variante 1
        self.qu_max = -self.m_qu_max * 120 / 110 + 0.48
        self.m_qu_min = (0.33 + 0.41) / (96 - 103) * 110  # Variante 3
        self.qu_min = -self.m_qu_min * 96 / 110 + 0.33
        return True

    def get_hc_acopf_parameters(self, model, net):
    # --- Parameters ---
        model.SWmax = Param(self.SWmax_tuple, initialize=self.SWmax_data_dict, mutable=True)
        model.SWmin = Param(self.SWmin_data_dict, initialize=self.SWmin_data_dict, mutable=True)

        if 'windpot_p_mw' in self.net.bus:
            self.Windpot_data_dict, self.Windpot_tuple = self.make_to_dict(model.WIND_HC, model.T, self.static_generation_data['windpot'])
            model.pWmax = Param(self.Windpot_tuple,
                                     initialize=self.Windpot_data_dict,
                                     mutable=True)
        return True

    def get_hc_acopf_variables(self, model):
        # --- Variables ---
        model.y = Var(self.Windpot_tuple, within=Binary, initialize=1.)
        return True

    def get_opf_parameters(self, model):
        # --- Parameters ---
        #super().get_opf_parameters(model)

        model.var_q = Param(model.WINDc, model.T, initialize=self.static_generation_data['var_q'][model.WINDc])
        model.PsG_inst = Param(model.WINDc, model.T, initialize=self.static_generation_data['p_inst'][model.WINDc])
        return True

    def static_generation_wind_var_q(self, net):
        x = np.array([[96, 103], [120, 127]]) / 110
        y = np.array([[0.48, 0.41, 0.33], [-0.23, -0.33, -0.41]])
        m = (y[1] - y[0]) / (x[0, 1] - x[0, 0])
        b = np.array([y[0] - m * x[i, 0] for i in range(len(x))]).T
        # self.m_v = m
        # self.b_min = b[:, 0]
        # self.b_max = b[:, 1]

        m_qp_max = (0.1 - y[0]) / (0.1 - 0.2)
        m_qp_min = (-0.1 - y[1]) / (0.1 - 0.2)
        b_qp_max = 0.1 - m_qp_max * 0.1
        b_qp_min = -0.1 - m_qp_min * 0.1

        self.q_limit_parameter = pd.DataFrame(
            {'m_qv': m, 'b_qv_min': b[:, 0], 'b_qv_max': b[:, 1], 'm_qp_max': m_qp_max, 'm_qp_min': m_qp_min,
             'b_qp_max': b_qp_max, 'b_qp_min': b_qp_min})

        if 'var_q' in self.net.sgen:
            self.static_generation_data['var_q'] = self.net.sgen.var_q.values
            sgens_var_q = self.static_generation_data.index[self.static_generation_data.var_q.notna()]

            try:
                p_inst = self.net.sgen.p_inst_mw.values / self.baseMVA
            except AttributeError:
                logging.warning(
                    "No p_inst_mw attribute found in self.net.sgen. Using p as p_inst for wind generators power limits.")
                p_inst = self.static_generation_data['p']

            self.static_generation_data['p_inst'] = p_inst

            self.static_generation_data['max_q'][sgens_var_q] = [
                y[0, int(self.static_generation_data.var_q[g])] * self.static_generation_data['p_inst'][g] for g in
                sgens_var_q]
            self.static_generation_data['min_q'][sgens_var_q] = [
                y[1, int(self.static_generation_data.var_q[g])] * self.static_generation_data['p_inst'][g] for g in
                sgens_var_q]

            self.static_generation_data['max_p'][sgens_var_q] = p_inst[sgens_var_q]
            self.static_generation_data['min_p'][sgens_var_q] = p_inst[sgens_var_q] * 0.1

        else:
            self.static_generation_data['var_q'] = None
            self.static_generation_data['p_inst'] = None

        if 'wind_hc' in self.net.sgen:
            self.static_generation_data['wind_hc'] = self.net.sgen.wind_hc.values
        else:
            self.static_generation_data['wind_hc'] = False


         # --- Constraints and objective ---

    def get_objective(self, model):
        # --- Objective ---
        # @model.Objective(model.WIND_HC, sense=maximize)
        # def obj(model, w):
        #     return sum(model.psG[w] for w in model.WIND_HC) - sum(
        #         model.pLfrom[l] + model.pLto[l] for l in model.L) - sum(
        #         model.pThv[t] + model.pTlv[t] for t in model.TRANSF)


        # old objective function because solver has problems with the new one in decorators.
        def obj_wind_loss_rule(model):
             return sum(model.psG[w] for w in model.WIND_HC) - sum(model.pLfrom[l] + model.pLto[l] for l in model.L) - sum(model.pThv[t] + model.pTlv[t] for t in model.TRANSF)

        model.obj = Objective(rule=obj_wind_loss_rule, sense=maximize)

    def get_constraints(self, model, net):
        # --- Constraints ---

        # --- constraints for wind genrators from ACOPF_base ---

        @model.Constraint(model.WINDc)
        def QW_pos(model, w):
            return model.qsG[w] <= self.q_limit_parameter.b_qp_max[model.var_q[w]] * model.PsG_inst[w] + \
                self.q_limit_parameter.m_qp_max[model.var_q[w]] * model.psG[w]

        @model.Constraint(model.WINDc)
        def QW_neg(model, w):
            return model.qsG[w] >= self.q_limit_parameter.b_qp_min[model.var_q[w]] * model.PsG_inst[w] + \
                self.q_limit_parameter.m_qp_min[model.var_q[w]] * model.psG[w]

        @model.Constraint(model.WINDc)
        def QV_min(model, w):
            for(g,b) in model.sGbs:
                if g == w:
                    return model.qsG[w] >= (self.q_limit_parameter.m_qv[model.var_q[w]] * model.v[b] + self.q_limit_parameter.b_qv_min[model.var_q[w]]) * model.PsG_inst[w]

        @model.Constraint(model.WINDc)
        def QV_max(model, w):
            for(g,b) in model.sGbs:
                if g == w:
                    return model.qsG[w] <= (self.q_limit_parameter.m_qv[model.var_q[w]] * model.v[b] + self.q_limit_parameter.b_qv_max[model.var_q[w]]) * model.PsG_inst[w]

        # --- constraints for wind genrators from HC_ACOPF_base ---
        @model.Constraint(model.WIND_HC)
        def SW_max(model, w):
            return model.psG[w] ** 2 + model.qsG[w] ** 2 <= model.SWmax[w] ** 2 * model.y[w]

        @model.Constraint(model.WIND_HC)
        def SW_min(model, w):
            return model.psG[w] ** 2 + model.qsG[w] ** 2 >= model.SWmin[w] ** 2 * model.y[w]

        # --- QU Variante 1 ---
        @model.Constraint(model.WIND_HC)
        def QW_min(model, w):
            return model.qsG[w] >= -0.41 * model.psG[w]

        @model.Constraint(model.WIND_HC)
        def QW_max(model, w):
            return model.qsG[w] <= 0.48 * model.psG[w]

        @model.Constraint(model.WIND_HC)
        def QU_min_hc(model, w):
            for (g, b) in model.sGbs:
                if g == w:
                    return model.qsG[w] >= (self.m_qu_min * model.v[b] + self.qu_min) * model.psG[w]

        @model.Constraint(model.WIND_HC)
        def QU_max_hc(model, w):
            for (g, b) in model.sGbs:
                if g == w:
                    return model.qsG[w] <= (self.m_qu_max * model.v[b] + self.qu_max) * model.psG[w]

        if 'windpot_p_mw' in net.bus:
            @model.Constraint(model.WIND_HC)
            def PW_max(model, w):
                return model.psG[w] <= model.pWmax[w]

    def unfix_variables(self, model):
        # --- generator power ---
        for w in model.WIND_HC:
            model.psG[w].unfix()
            model.qsG[w].unfix()

    def get_all_acopf(self,model):
        '''
        getter in windpower created so no error occurs when running get_all_acopf for flexibilities
        '''
        pass