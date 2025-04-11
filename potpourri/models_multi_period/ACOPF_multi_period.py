import pandas as pd
from pyomo.environ import *
from potpourri.models_multi_period.AC_multi_period import AC_multi_period
from potpourri.models_multi_period.OPF_multi_period import OPF_multi_period
from potpourri.models_multi_period.generator_multi_period import Generator_multi_period
from potpourri.models_multi_period.demand_multi_period import Demand_multi_period
from potpourri.models_multi_period.windpower_multi_period import Windpower_multi_period
from potpourri.models_multi_period.sgens_multi_period import Sgens_multi_period
from potpourri.models_multi_period.EVs_multi_period import EV_multi_period


import numpy as np
import logging


class ACOPF_multi_period(AC_multi_period, OPF_multi_period):
    """
    This class implements the AC Optimal Power Flow (ACOPF) model.
    """
    def __init__(self, net, toT,  fromT=None, pf=1, vehicles=None, locations=None, chargingpoints=None, scenario=None):
        super().__init__(net, toT,  fromT, pf, vehicles, locations, chargingpoints, scenario)

    def _calc_opf_parameters(self):
        super()._calc_opf_parameters()

        max_vm_pu, min_vm_pu = self.get_v_limits()
        self.v_limits = (max_vm_pu, min_vm_pu)

        # create sgen instance
        sgens_object = next((obj for obj in self.flexibilities if isinstance(obj, Sgens_multi_period)), None)
        sgens_object.static_generation_reactive_power_limits(self.model) # gives the model now instead of the net

        # get the object of class 'Windpower' from the 'flexibilities' list
        if 'windpot_p_mw' in self.net.bus:
            windpower_object = next((obj for obj in self.flexibilities if isinstance(obj, Windpower_multi_period)), None)
            windpower_object.static_generation_wind_var_q(self.net)

        # create generator instance
        generator_object = next((obj for obj in self.flexibilities if isinstance(obj, Generator_multi_period)), None)
        generator_object.generation_reactive_power_limits_acopf()

        # create demand instance
        demand_object = next((obj for obj in self.flexibilities if isinstance(obj, Demand_multi_period)), None)
        demand_object.get_demand_reactive_data(self.model)

    def get_v_limits(self):
        if 'max_vm_pu' in self.net.bus:
            max_vm_pu = self.net.bus.max_vm_pu.values
        else:
            max_vm_pu = np.full(len(self.net.bus.index), 1.1)

        if 'min_vm_pu' in self.net.bus:
            min_vm_pu = self.net.bus.min_vm_pu.values
        else:
            min_vm_pu = np.full(len(self.net.bus.index), 0.9)

        if any(self.net.gen.index):
            self.add_generator_v_limits(max_vm_pu, min_vm_pu)

        return max_vm_pu, min_vm_pu

    def add_generator_v_limits(self, max_vm_pu, min_vm_pu):
        # check max_vm_pu / min_vm_pu bus limit violation by gens
        gen_buses = self.bus_lookup[self.net.gen.bus.values]
        if "max_vm_pu" in self.net["gen"].columns:
            v_max_bound = max_vm_pu[gen_buses] < self.net["gen"]["max_vm_pu"].values
            if np.any(v_max_bound):
                bound_gens = self.net["gen"].index.values[v_max_bound]
                print("gen max_vm_pu > bus max_vm_pu for gens {}. "
                      "Setting bus limit for these gens.".format(bound_gens))
                # set only vm of gens which do not violate the limits
                max_vm_pu[gen_buses[~v_max_bound]] = self.net["gen"]["max_vm_pu"].values[~v_max_bound]
            else:
                # set vm of all gens
                max_vm_pu[gen_buses] = self.net["gen"]["max_vm_pu"].values

        if "min_vm_pu" in self.net["gen"].columns:
            v_min_bound = self.net["gen"]["min_vm_pu"].values < min_vm_pu[gen_buses]
            if np.any(v_min_bound):
                bound_gens = self.net["gen"].index.values[v_min_bound]
                print("gen min_vm_pu < bus min_vm_pu for gens {}. "
                      "Setting bus limit for these gens.".format(bound_gens))
                # set only vm of gens which do not violate the limits
                min_vm_pu[gen_buses[~v_min_bound]] = self.net["gen"]["min_vm_pu"].values[~v_min_bound]
            else:
                # set vm of all gens
                min_vm_pu[gen_buses] = self.net["gen"]["min_vm_pu"].values

        if 'controllable' in self.net.gen:
            controllable = self.net["gen"]["controllable"].values
            not_controllable = ~controllable.astype(bool)

            # get voltage setpoints for not controllable generators
            if np.any(not_controllable):
                bus = self.net["gen"]["bus"].values[not_controllable]
                vm_pu = self.net["gen"]["vm_pu"].values[not_controllable]

                not_controllable_buses = self.bus_lookup[bus]
                max_vm_pu[not_controllable_buses] = vm_pu
                min_vm_pu[not_controllable_buses] = vm_pu

        return max_vm_pu, min_vm_pu

    def add_OPF(self, **kwargs):
        super().add_OPF(**kwargs)

        self.model.name = "ACOPF"

        # run get_all_opf from flexibility instances
        for flex in self.flexibilities:
            flex.get_all_acopf(self.model)

        # voltage limits DONE: make non time dependent
        self.model.Vmax = Param(self.model.B, within=NonNegativeReals, initialize=self.v_limits[0][self.model.B], mutable=True)  # max voltage (p.u.)
        self.model.Vmin = Param(self.model.B, within=NonNegativeReals,initialize=self.v_limits[1][self.model.B], mutable=True)  # min voltage (p.u.)

        # --- cost function ---
        # def objective(model):
        #     obj = sum(
        #         model.c2[g] * (model.baseMVA * model.pG[g]) ** 2 + model.c1[g] * model.baseMVA * model.pG[g] +
        #         model.c0[g] for g in model.G) + \
        #           sum(model.VOLL[d] * (model.PD[d] - model.pD[d]) * model.baseMVA for d in model.D)
        #     return obj
        #
        # self.model.OBJ = Objective(rule=objective, sense=minimize)

        # --- line power limits --- DONE non time dependent
        def line_lim_from_def(model, l, t):
            return model.pLfrom[l,t] ** 2 + model.qLfrom[l,t] ** 2 <= model.SLmax[l] ** 2 * model.v[model.A[l, 1],t] ** 2

        def line_lim_to_def(model, l,t):
            return model.pLto[l,t] ** 2 + model.qLto[l,t] ** 2 <= model.SLmax[l] ** 2 * model.v[model.A[l, 2],t] ** 2

        self.model.line_lim_from = Constraint(self.model.L, self.model.T, rule=line_lim_from_def)
        self.model.line_lim_to = Constraint(self.model.L, self.model.T, rule=line_lim_to_def)

        # --- power flow limits on transformer lines--- DONE non time dependent
        def transf_lim1_def(model, l, t):
            return model.pThv[l,t] ** 2 + model.qThv[l,t] ** 2 <= model.SLmaxT[l] ** 2 * model.v[model.AT[l, 1], t] ** 2

        def transf_lim2_def(model, l,t):
            return model.pTlv[l,t] ** 2 + model.qTlv[l,t] ** 2 <= model.SLmaxT[l] ** 2 * model.v[model.AT[l, 2], t] ** 2

        self.model.transf_lim1 = Constraint(self.model.TRANSF, self.model.T, rule=transf_lim1_def)
        self.model.transf_lim2 = Constraint(self.model.TRANSF, self.model.T, rule=transf_lim2_def)

        # --- voltage constraints --- can be removed
        # self.model.v_bPV_setpoint.deactivate()

        #should be time dependent
        def v_bounds(model, b, t):
            return model.Vmin[b], model.v[b, t], model.Vmax[b]

        self.model.v_constraint = Constraint(self.model.B, self.model.T, rule=v_bounds)

        # # --- wind generation q requirements variant 3---
        # def QW_pos(model, w):
        #     return model.qsG[w] <= self.q_limit_parameter.b_qp_max[model.var_q[w]] * model.PsG_inst[w] + self.q_limit_parameter.m_qp_max[model.var_q[w]] * model.psG[w]
        #
        # def QW_neg(model, w):
        #     return model.qsG[w] >= self.q_limit_parameter.b_qp_min[model.var_q[w]] * model.PsG_inst[w] + self.q_limit_parameter.m_qp_min[model.var_q[w]] * model.psG[w]
        #
        # self.model.QW_pos_constraint = Constraint(self.model.WINDc, rule=QW_pos)
        # self.model.QW_neg_constraint = Constraint(self.model.WINDc, rule=QW_neg)
        #
        # #
        # def QV_min(model, w):
        #     for (g, b) in model.sGbs:
        #         if g == w:
        #             return model.qsG[w] >= (self.q_limit_parameter.m_qv[model.var_q[w]] * model.v[b] + self.q_limit_parameter.b_qv_min[model.var_q[w]]) * model.PsG_inst[w]
        # self.model.QU_min_constraint = Constraint(self.model.WINDc, rule=QV_min)
        #
        # def QV_max(model, w):
        #     for (g, b) in model.sGbs:
        #         if g == w:
        #             return model.qsG[w] <= (self.q_limit_parameter.m_qv[model.var_q[w]] * model.v[b] + self.q_limit_parameter.b_qv_max[model.var_q[w]]) * model.PsG_inst[w]
        # self.model.QU_max_constraint = Constraint(self.model.WINDc, rule=QV_max)
        #

    def add_voltage_deviation_objective(self):
        self.model.vm = Param(self.model.B, initialize=self.bus_data['v_m'][self.model.B])

        def voltage_deviation_objective(model, t):
            return sum((model.v[b, t] - 1.) ** 2 for b in model.B - model.b0 for t in model.T) + \
                sum((model.v[b, t] - model.v_b0[b]) ** 2 for b in model.b0 for t in model.T)

        self.model.obj_v_deviation = Objective(rule=voltage_deviation_objective, sense=minimize)

    # Assuming model, G, T, and pf are already defined

    #test objective function
    def add_minimize_power_objective(self):
        def power_minimization_objective(model):
            return sum(model.pD[d, t] for d in model.D for t in model.T)

        self.model.Objective = Objective(rule=power_minimization_objective, sense=minimize)

    def add_generation_objective(self):
        def minimize_generation(model):
            return sum(model.pG[(g, t)]**2 for g in model.G for t in model.T)

        self.model.obj = Objective(rule=minimize_generation, sense=minimize)

    # possible objective: weighted generation and discharging
    def add_weighted_generation_objective(self):
        def weighted_generation_objective(model):
            c1 = 4
            c2 = 1
            c3 = 1
            return ((c1 * sum(model.pG[(g, t)] for g in model.G for t in model.T)
                     + c2 * sum(model.p_discharging[(v, t)] for v in model.veh for t in model.T))
                    + c3 * sum(model.psG[(g, t)] for g in model.sG for t in model.T))
        self.model.obj = Objective(rule=weighted_generation_objective, sense=minimize)

    def add_charging_power_obj(self):
        def maximize_charging_power(model):
            return sum(model.p_opf[v, t] for v in model.veh for t in model.T)

        self.model.charging_power_obj = Objective(rule=maximize_charging_power, sense=maximize)

    def add_discharging_power_obj(self):
        def maximize_discharging_power(model):
            return sum(model.p_opf[v, t] for v in model.veh for t in model.T)

        self.model.discharging_power_obj = Objective(rule=maximize_discharging_power, sense=minimize)

    def add_arbitrage_objective(self):
        ev_object = next((obj for obj in self.flexibilities if isinstance(obj, EV_multi_period)), None)
        ev_object.get_market_constraints(self.model)

        # day ahead market arbitrage
        def arbitrage_objective(model):
            return sum(model.p_opf[v, t] * model.p_da[t] * model.timestep_size/60 + model.buffer_soc[v, t] * model.bigM for v in model.veh for t in model.T)  # p_da in [€/MWh], timestep_size in [min], + model.buffer_soc[v, t] * model.bigM

        self.model.arbitrage_obj = Objective(rule=arbitrage_objective, sense=minimize)

