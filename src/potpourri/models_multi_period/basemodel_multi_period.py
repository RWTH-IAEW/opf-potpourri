import pandas as pd
from pyomo.environ import *
from math import pi
import copy
import os
import numpy as np
import pandapower as pp
import simbench as sb
import time as ctime
from src.potpourri.models_multi_period.generator_multi_period import Generator_multi_period
from src.potpourri.models_multi_period.shunts_multi_period import Shunts_multi_period
from src.potpourri.models_multi_period.sgens_multi_period import Sgens_multi_period
from src.potpourri.models_multi_period.demand_multi_period import Demand_multi_period
from src.potpourri.models_multi_period.windpower_multi_period import Windpower_multi_period
from src.potpourri.models_multi_period.EVs_multi_period import EV_multi_period
from src.potpourri.models_multi_period.flexibility_multi_period import Flexibility_multi_period
from src.potpourri.models_multi_period.pyo_to_net_multi_period import pyo_sol_to_net_res


class Basemodel_multi_period:
    def __init__(self, net, toT,  fromT=None, pf=1, num_vehicles=None):
        """
        **Constructor Basemodel_multi_period** \n
        *net:* simbench/ T.B.D. pandapower ... network \n
        *toT:* Model runs till this Time-step, beware: Index 3 equals time step 4, int \n
        *fromT:* Model runs from this Time-step beware: Index 0 equals time step 1, standard value=None, int \n
        *pf:* Power Factor (cos between active and reactive power, equals 1 if reactive power is zero), standard value: 1,float \n
        *scenario:* Scenario for the model, standard value=None, 0-2020, 1-2030, 2-2040, 3-2050, int \n
        """
        self.net = copy.deepcopy(net)
        pp.runpp(self.net)

        # --- initialize flexibility list ---
        self.flexibilities = []

        # --- Sets ---
        bus_set = self.net._ppc['bus'][:, [0, 1, 7, 8]]
        bus_set[:, -1] *= pi / 180
        self.bus_data = pd.DataFrame(bus_set[:, 1:], index=bus_set[:, 0].astype(int), columns=['type', 'v_m', 'v_a_rad'])
        self.bus_lookup = self.net._pd2ppc_lookups["bus"]

        # --- Param Data ---
        self.baseMVA = self.net.sn_mva

        # --- line ---
        hv_bus = self.net._ppc['branch'][:, 0].real
        lv_bus = self.net._ppc['branch'][:, 1].real
        trafo_start = len(self.net.line.index)
        trafo_end = trafo_start + len(self.net.trafo.index)

        hv_bus_line = hv_bus[:trafo_start]
        lv_bus_line = lv_bus[:trafo_start]
        self.line_data = pd.DataFrame({'in_service': self.net.line.in_service.values})
        line_ind = self.line_data.index[self.line_data.in_service]
        self.bus_line_dict = dict(zip(list(zip(line_ind, [1] * len(line_ind))) + list(
            zip(line_ind, [2] * len(line_ind))), np.concatenate([hv_bus_line[line_ind], lv_bus_line[line_ind]])))

        # --- transformer ---
        shift = self.net._ppc['branch'][trafo_start:trafo_end, 9].real * pi / 180
        tap = self.net._ppc['branch'][trafo_start:trafo_end, 8].real
        self.trafo_data = pd.DataFrame({'in_service': self.net.trafo.in_service.values, 'shift_rad': shift, 'tap': tap})

        hv_bus_trafo = hv_bus[trafo_start:trafo_end]
        lv_bus_trafo = lv_bus[trafo_start:trafo_end]
        trafo_ind = self.trafo_data.index[self.trafo_data.in_service]
        self.bus_trafo_dict = dict(zip(list(zip(trafo_ind, [1] * len(trafo_ind))) + list(
            zip(trafo_ind, [2] * len(trafo_ind))), np.concatenate([hv_bus_trafo[trafo_ind], lv_bus_trafo[trafo_ind]])))

        self.toT = toT
        if fromT is not None:
            self.fromT = fromT
        else:
            self.fromT = 0
        self.T = toT - fromT if fromT else toT
        self.pf = pf # powerfactor for reactive power calculation
        if toT != None and not hasattr(self.net, "profiles"):
            raise ValueError("The net object does not have profiles. Please provide a net object with profiles.")
        self.profiles = sb.get_absolute_values(self.net, profiles_instead_of_study_cases=True)
        # TODO: q static generation data with power factor ( check for power factor in pp net)

        # get the load profiles out of the sb net for pv before they are deleted
        # in the future profiles might be coming from another source
        self.net.pv_load_profiles = self.net.profiles['renewables'].iloc[self.fromT:self.toT]

        self.net.profiles.clear() #  TODO: needed here?
        # cut off profiles for T, depending on fromT and toT
        if self.fromT == None:
            for profile in self.profiles.keys():
                self.net.profiles[profile] = self.profiles[profile].iloc[:self.toT]  # cut off profiles to T
        elif self.fromT < 0:
            raise ValueError("fromT must be positive integer")
        elif self.toT < 0:
            raise ValueError("toT must be positive integer")
        elif self.toT > len(self.profiles[list(self.profiles.keys())[0]]):
            raise ValueError("toT must be smaller than the length of the profiles")
        elif self.toT < self.fromT:
            raise ValueError("toT must be greater than fromT")
        elif self.fromT != None:
            for profile in self.profiles.keys():
                self.net.profiles[profile] = self.profiles[profile].iloc[self.fromT:self.toT]  # cuts out the profile from fromT to toT

        # print("Model runs from time step " + str(self.fromT+1) + " to time step " + str(self.toT))

        # calculation of reactive power for static generators with power factor
        if not ('sgen', 'q_mvar') in self.net.profiles.keys():
            self.calc_reactive_sgen_power(self.pf)

        # --- create flexibility object and append to list ---

        self.flexibilities.append(Demand_multi_period(self.net))
        self.flexibilities.append(Shunts_multi_period(self.net))
        self.flexibilities.append(Sgens_multi_period(self.net))
        self.flexibilities.append(Generator_multi_period(self.net))
        # only create instance of EV if necessary parameters are provided
        if num_vehicles is not None :
            self.flexibilities.append(EV_multi_period(self.net, num_vehicles))

        if 'windpot_p_mw' in self.net.bus:
            self.flexibilities.append(Windpower_multi_period(self.net))

    # calculates reactive power for static generators with given power factor
    def calc_reactive_sgen_power(self, pf=1):
        self.net.profiles[('sgen', 'q_mvar')] = self.net.profiles[('sgen', 'p_mw')] * np.tan(np.arccos(pf))

    def create_model(self):
        print("Creating Model at: ", ctime.ctime())
        self.model = ConcreteModel()

        # time dependency to model

        self.model.T = Set(initialize=range(self.fromT, self.toT), ordered=True)  # time periods
        self.deltaT = (1/4) # time step length in hours = 15 minutes
        self.model.deltaT = Param(initialize=self.deltaT, within=PositiveReals)  # time step length

        # --- iterate through flexibility list ---
        for flex in self.flexibilities:
            flex.get_all(self.model)

        # --- SETS ---
        self.model.b0 = Set(initialize=self.bus_data.index[self.bus_data.type == 3],
                            within=self.model.B)  # reference buses
        self.model.bPV = Set(initialize=self.bus_data.index[self.bus_data.type == 2], within=self.model.B)  # PV buses

        self.model.L = Set(initialize=self.line_data.index[self.line_data.in_service])
        self.model.LE = Set(initialize=[1, 2])
        self.model.TRANSF = Set(initialize=self.trafo_data.index[self.trafo_data.in_service])

        # --- parameters ---
        # line and trafo matrix
        self.model.A = Param(self.model.L*self.model.LE, initialize=self.bus_line_dict)  # bus-line matrix
        self.model.AT = Param(self.model.TRANSF*self.model.LE, initialize=self.bus_trafo_dict)  # bus-transformer matrix

        # trafo
        self.model.shift = Param(self.model.TRANSF, within=Reals, initialize=self.trafo_data.shift_rad[self.model.TRANSF])  # transformer phase shift in rad DONE:remove Time dependency

        # external grid voltage angle DONE: remove time depdendency
        self.model.delta_b0 = Param(self.model.b0, within=Reals, initialize=self.bus_data.v_a_rad[self.model.b0])

        # baseMVA of the net
        self.model.baseMVA = Param(within=NonNegativeReals, initialize=self.baseMVA)

        # --- variables ---
        #Done stay multiperiod
        self.delta_data_dict, self.delta_tuple = self.make_to_dict(self.model.B, self.model.T, 0.0, False) #False or true?
        self.pLfrom_tuple = self.make_to_tuple(self.model.L, self.model.T)
        self.pLto_tuple = self.make_to_tuple(self.model.L, self.model.T)
        self.pThv_tuple = self.make_to_tuple(self.model.TRANSF, self.model.T)
        self.pTlv_tuple = self.make_to_tuple(self.model.TRANSF, self.model.T)
        self.Tap_data_dict, self.Tap_tuple = self.make_to_dict(self.model.TRANSF, self.model.T, self.trafo_data.tap[self.model.TRANSF], False) #False or true on time dpendency?


        self.model.delta = Var(self.delta_tuple, domain=Reals, initialize=self.delta_data_dict,
                               bounds=(-pi, pi))  # voltage phase angle at bus b, rad
        self.model.pLfrom = Var(self.pLfrom_tuple, domain=Reals)  # real power injected at b onto line
        self.model.pLto = Var(self.pLto_tuple, domain=Reals)  # real power injected at b' onto line
        self.model.pThv = Var(self.pThv_tuple, domain=Reals)  # real power injected at b onto transformer
        self.model.pTlv = Var(self.pTlv_tuple, domain=Reals)  # real power injected at b' onto transformer
        self.model.Tap = Var(self.Tap_tuple, domain=Reals, initialize=self.Tap_data_dict)  # transformer tap ratio


        # transformer tap ratio
        for tr in self.model.TRANSF:
            for t in self.model.T:
                self.model.Tap[tr,t].fix()

        # --- reference bus constraint ---
        for b in self.model.b0:
            self.model.delta[b, t].fix(self.model.delta_b0[b])

    def solve(self, to_net: bool = True, print_solver_output: bool = True, solver='ipopt', load_solutions: bool = True,
              mip_solver='gurobi', max_iter=None, time_limit=600, init_strategy='rNLP', neos_opt='bonmin'):
        optimizer = SolverFactory(solver)

        if solver == 'mindtpy':
            if not max_iter:
                max_iter = 50

            if mip_solver == 'gurobi':
                mip_solver = 'gurobi_persistent'

            try:
                self.results = optimizer.solve(self.model, mip_solver=mip_solver, nlp_solver='ipopt',
                                               tee=print_solver_output, iteration_limit=max_iter, time_limit=time_limit, init_strategy=init_strategy)
            except ValueError as err:
                print(err)

            # if self.results.solver.termination_condition == TerminationCondition.feasible:
            #     print('Model is feasible but not optimal. Trying to solve a second time.')
            #     self.results = optimizer.solve(self.model, mip_solver=mip_solver, nlp_solver='ipopt',
            #                                    tee=print_solver_output, iteration_limit=max_iter, time_limit=time_limit)
        elif solver == 'neos':
            os.environ['NEOS_EMAIL'] = "xyz@rwth-aachen.de"
            solver_manager = SolverManagerFactory('neos')
            self.results = solver_manager.solve(self.model, opt=neos_opt, tee=True)

        else:
            if max_iter:
                optimizer.options['max_iter'] = max_iter

            self.results = optimizer.solve(self.model, load_solutions=load_solutions, tee=print_solver_output)

        try:
            if check_optimal_termination(self.results) & to_net:
                #pyo_sol_to_net_res(self.net, self.model)
                print("Solved successfully")
        except AttributeError as err:
            print(err)

    def change_vals(self, key, value):
        component = self.model.component(key)
        if not component:
            print("Model " + self.model.name + " has no component " + key)
            return
        try:
            for index in component:
                component[index] = value
        except TypeError as err:
            print(err)

    def fix_vars(self, key, value=None):
        component = self.model.component(key)
        if not component:
            print("Model " + self.model.name + " has no component " + key)
            return
        try:
            for index in component:
                if value is not None:
                    component[index].fix(value)
                else:
                    component[index].fix()
        except AttributeError as err:
            print(err)

    def unfix_vars(self, key, value=None):
        component = self.model.component(key)
        if not component:
            print("Model " + self.model.name + " has no component " + key)
            return
        try:
            for index in component:
                component[index].unfix()
                if value:
                    component[index] = value
        except AttributeError as err:
            print(err)

    def make_to_dict(self, model_obj, model_time, data, time_dependent=True):
        """
        **make_to_dict** \n
        Args:
            model_obj: Object Indices
            model_time: Time Indices
            data: Data to be converted to dictionary
            time_dependent: should the data be time dependent True or constant for all time steps False

        Returns:
            data_dict: Dictionary with the data
            tuple_list: List of tuples with the object and time index

        """

        # if isinstance(data, float):
        #     data_dict = {(o, t): data for o in model_obj for t in model_time}
        #     tuple_list = list([(o, t) for o in model_obj for t in model_time])
        #     return data_dict, tuple_list

        if data is 0:
            data_dict = {(o, t): 0 for o in model_obj for t in model_time}
            tuple_list = list([(o, t) for o in model_obj for t in model_time])
            return data_dict, tuple_list

        if isinstance(data, np.ndarray):
            if time_dependent:
                data_dict = {(o, t): data[o][t] for o in model_obj for t in model_time}
                tuple_list = [(o, t) for o in model_obj for t in model_time]
            else:
                data_dict = {(o, t): data[o] for o in model_obj for t in model_time}
                tuple_list = [(o, t) for o in model_obj for t in model_time]
            return data_dict, tuple_list

        # when data is float put it in dict, correct?
        if isinstance(data, float):
            data_dict = {(o, t): data for o in model_obj for t in model_time}
            tuple_list = list([(o, t) for o in model_obj for t in model_time])
            return data_dict, tuple_list

        #if data is already a dict just put it in data_dict and make tuple_list with time
        if isinstance(data, dict):
            if time_dependent:
                data_dict = {(o, t): data[o][t] for o in model_obj for t in model_time}
                tuple_list = list([(o, t) for o in model_obj for t in model_time])
                return data_dict, tuple_list
            else:
                data_dict = {(o, t): data[o] for o in model_obj for t in model_time}
                tuple_list = list([(o, t) for o in model_obj for t in model_time])
                return data_dict, tuple_list

        # if isinstance(data, dict):
        #     return data, list(data.keys())
        # make data_dict with constant values over time if not time dependent
        if time_dependent:
            data_dict = data.to_dict()
            data_dict = {(o, t): data_dict[o][t] for o in model_obj for t in model_time}
            tuple_list = list([(o, t) for o in model_obj for t in model_time])
        else:
            data_dict = data.to_dict()
            data_dict = {(o, t): data_dict[o] for o in model_obj for t in model_time}
            tuple_list = list([(o, t) for o in model_obj for t in model_time])

        return data_dict, tuple_list

    def make_to_tuple(self, model_obj, model_time):
        """
        **make_to_tuple** \n
        Args:
            model_obj: Object Indices
            model_time: Time Indices

        Returns:
            tuple_list: List of tuples with the object and time index

        """
        tuple_list = list([(o, t) for o in model_obj for t in model_time])
        return tuple_list