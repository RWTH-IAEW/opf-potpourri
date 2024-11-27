import pandas as pd
from pyomo.environ import *
from math import pi
import copy
import numpy as np
import pandapower as pp
import os

from potpourri.models.pyo_to_net import pyo_sol_to_net_res


class Basemodel:
    def __init__(self, net):
        self.net = copy.deepcopy(net)
        pp.runpp(self.net)

        # --- Sets ---
        bus_set = self.net._ppc['bus'][:, [0, 1, 7, 8]]
        bus_set[:, -1] *= pi / 180
        self.bus_data = pd.DataFrame(bus_set[:, 1:], index=bus_set[:, 0].astype(int), columns=['type', 'v_m', 'v_a_rad'])

        self.bus_lookup = self.net._pd2ppc_lookups["bus"]
        self.demand_set = self.net.load.index[self.net.load.in_service]
        self.shunt_set = self.net.shunt.index[self.net.shunt.in_service]

        self.bus_demand_set = list(zip(self.bus_lookup[self.net.load.bus[self.demand_set].values], self.demand_set))
        self.bus_shunt_set = list(zip(self.bus_lookup[self.net.shunt.bus[self.shunt_set].values], self.shunt_set))

        # --- Param Data ---
        self.baseMVA = self.net.sn_mva

        self.PD_data = self.net.load.p_mw * self.net.load.scaling / self.baseMVA

        self.GB_data = self.net.shunt.p_mw * self.net.shunt.step / self.baseMVA

        # --- generation ---
        pg = self.net._ppc['gen'][:, 1] / self.baseMVA
        ref_gens = self.net._ppc['internal']['ref_gens']
        in_service_gens = self.net._ppc['gen'][:, 7].astype(bool)
        gen_bus = self.net._ppc['gen'][:, 0].astype(int)
        self.generation_data = pd.DataFrame({'pg': pg, 'ref': False, 'in_service': in_service_gens, 'bus': gen_bus})
        self.generation_data.loc[ref_gens, 'ref'] = True
        gen_bus_tuples = list(enumerate(self.generation_data['bus']))
        self.generation_data['gen_bus'] = gen_bus_tuples

        # --- static generation ---
        psg = self.net.sgen.p_mw * self.net.sgen.scaling / self.baseMVA
        sgen_bus = self.bus_lookup[self.net.sgen.bus.values]
        self.static_generation_data = pd.DataFrame(
            {'p': psg.values, 'in_service': self.net.sgen.in_service.values, 'bus': sgen_bus})
        self.static_generation_data['gen_bus'] = list(enumerate(self.static_generation_data['bus']))

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

    def create_model(self):
        self.model = ConcreteModel()

        # --- SETS ---
        self.model.B = Set(initialize=self.bus_data.index)  # buses
        self.model.b0 = Set(initialize=self.bus_data.index[self.bus_data.type == 3],
                            within=self.model.B)  # reference buses
        self.model.bPV = Set(initialize=self.bus_data.index[self.bus_data.type == 2], within=self.model.B)  # PV buses
        self.model.sG = Set(
            initialize=self.static_generation_data.index[self.static_generation_data.in_service])  # static generators
        self.model.G = Set(
            initialize=self.generation_data.index[self.generation_data.in_service])  # external grids and generators
        self.model.eG = Set(initialize=self.generation_data.index[self.generation_data.ref],
                            within=self.model.G)  # external grids and slack generators
        self.model.gG = Set(initialize=self.generation_data.index[self.generation_data.ref == False],
                            within=self.model.G)  # generators (not static) not slack
        self.model.D = Set(initialize=self.demand_set)
        self.model.L = Set(initialize=self.line_data.index[self.line_data.in_service])
        self.model.SHUNT = Set(initialize=self.shunt_set)
        self.model.LE = Set(initialize=[1, 2])
        self.model.TRANSF = Set(initialize=self.trafo_data.index[self.trafo_data.in_service])

        # generators, buses, loads linked to each bus b
        self.model.Dbs = Set(within=self.model.B * self.model.D,
                             initialize=self.bus_demand_set)  # set of demand-bus mapping
        self.model.SHUNTbs = Set(within=self.model.B * self.model.SHUNT,
                                 initialize=self.bus_shunt_set)  # set of shunt-bus mapping
        self.model.Gbs = Set(within=self.model.G * self.model.B,
                             initialize=self.generation_data['gen_bus'][self.model.G])
        self.model.sGbs = Set(within=self.model.sG * self.model.B,
                              initialize=self.static_generation_data['gen_bus'][self.model.sG])

        # --- parameters ---
        # line and trafo matrix
        self.model.A = Param(self.model.L * self.model.LE, initialize=self.bus_line_dict)  # bus-line matrix
        self.model.AT = Param(self.model.TRANSF * self.model.LE,
                              initialize=self.bus_trafo_dict)  # bus-transformer matrix

        # generation
        self.model.PsG = Param(self.model.sG, initialize=self.static_generation_data.p[self.model.sG])
        self.model.PG = Param(self.model.G, initialize=self.generation_data.pg[self.model.G])

        # demand
        self.model.PD = Param(self.model.D, initialize=self.PD_data[self.model.D])

        # shunt
        self.model.GB = Param(self.model.SHUNT, within=Reals,
                              initialize=self.GB_data[self.model.SHUNT])  # shunt conductance

        # trafo
        self.model.shift = Param(self.model.TRANSF, within=Reals, initialize=self.trafo_data.shift_rad[
            self.model.TRANSF])  # transformer phase shift in rad

        # external grid voltage angle
        self.model.delta_b0 = Param(self.model.b0, within=Reals, initialize=self.bus_data.v_a_rad[self.model.b0])

        # baseMVA of the net
        self.model.baseMVA = Param(within=NonNegativeReals, initialize=self.baseMVA)

        # --- variables ---
        self.model.delta = Var(self.model.B, domain=Reals, initialize=0.0,
                               bounds=(-pi, pi))  # voltage phase angle at bus b, rad
        self.model.pD = Var(self.model.D, domain=Reals)  # real power demand delivered
        self.model.psG = Var(self.model.sG, domain=NonNegativeReals)  # real static generator power
        self.model.pG = Var(self.model.G, domain=Reals,
                            initialize=self.model.PG)  # real power injection from static generators
        self.model.pLfrom = Var(self.model.L, domain=Reals)  # real power injected at b onto line
        self.model.pLto = Var(self.model.L, domain=Reals)  # real power injected at b' onto line
        self.model.pThv = Var(self.model.TRANSF, domain=Reals)  # real power injected at b onto transformer
        self.model.pTlv = Var(self.model.TRANSF, domain=Reals)  # real power injected at b' onto transformer
        self.model.Tap = Var(self.model.TRANSF, domain=Reals,
                             initialize=self.trafo_data.tap[self.model.TRANSF])  # transformer tap ratio

        # transformer tap ratio
        for t in self.model.TRANSF:
            self.model.Tap[t].fix()

        # --- generator power ---
        for g in self.model.sG:
            self.model.psG[g].fix(self.model.PsG[g])
        for g in self.model.gG:
            self.model.pG[g].fix(self.model.PG[g])

        # --- demand ---
        for d in self.model.D:
            self.model.pD[d].fix(self.model.PD[d])

        # --- reference bus constraint ---
        for b in self.model.b0:
            self.model.delta[b].fix(self.model.delta_b0[b])

    def solve(self, to_net: bool = True, print_solver_output: bool = False, solver='ipopt', load_solutions: bool = True,
              mip_solver='gurobi', max_iter=None, time_limit=600, init_strategy='rNLP', neos_opt='ipopt'):

        if solver == 'mindtpy':
            optimizer = SolverFactory(solver)
            if not max_iter:
                max_iter = 50

            if mip_solver == 'gurobi':
                mip_solver = 'gurobi_persistent'

            try:
                self.results = optimizer.solve(self.model, mip_solver=mip_solver, nlp_solver='ipopt',
                                               tee=print_solver_output, iteration_limit=max_iter, time_limit=time_limit,
                                               init_strategy=init_strategy)
            except ValueError as err:
                print(err)

            # if self.results.solver.termination_condition == TerminationCondition.feasible:
            #     print('Model is feasible but not optimal. Trying to solve a second time.')
            #     self.results = optimizer.solve(self.model, mip_solver=mip_solver, nlp_solver='ipopt',
            #                                    tee=print_solver_output, iteration_limit=max_iter, time_limit=time_limit)
        elif solver == 'neos':
            os.environ['NEOS_EMAIL'] = "ben.jamin@bluem-chen.de"
            solver_manager = SolverManagerFactory('neos')
            self.results = solver_manager.solve(self.model, opt=neos_opt, tee=True)
        else:
            if max_iter:
                optimizer.options['max_iter'] = max_iter

            self.results = optimizer.solve(self.model, load_solutions=load_solutions, tee=print_solver_output)

        try:
            if check_optimal_termination(self.results) & to_net:
                pyo_sol_to_net_res(self.net, self.model)
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
