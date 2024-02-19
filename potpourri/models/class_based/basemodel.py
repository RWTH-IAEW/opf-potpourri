import pandas as pd
from pyomo.environ import *
from math import pi
import copy
import numpy as np
import pandapower as pp

from potpourri.models.class_based.pyo_to_net import pyo_sol_to_net_res


class Basemodel:
    def __init__(self, net):
        self.net = copy.deepcopy(net)
        pp.runpp(self.net)
        self.net.gen.index += len(self.net.sgen.index)

        self.generators = pd.concat([self.net.sgen, self.net.gen])

        # --- Sets ---
        self.gen_set = self.net.gen.index[self.net.gen.in_service]
        self.sgen_set = self.net.sgen.index[self.net.sgen.in_service]
        self.gen_all_set = self.generators.index[self.generators.in_service]

        # self.bus_set = self.net.bus.index[self.net.bus.in_service]
        self.bus_lookup = self.net._pd2ppc_lookups["bus"]
        self.bus_set = self.net._ppc['bus'][:, 0].astype(int)
        self.demand_set = self.net.load.index[self.net.load.in_service]
        self.line_set = self.net.line.index[self.net.line.in_service]
        self.shunt_set = self.net.shunt.index[self.net.shunt.in_service]
        self.ext_grid_set = self.net.ext_grid.index[self.net.ext_grid.in_service]

        self.bus_gen_set = list(zip(self.bus_lookup[self.generators.bus[self.gen_all_set].values], self.gen_all_set))
        self.bus_demand_set = list(zip(self.bus_lookup[self.net.load.bus[self.demand_set].values], self.demand_set))
        self.bus_shunt_set = list(zip(self.bus_lookup[self.net.shunt.bus[self.shunt_set].values], self.shunt_set))
        self.bus_ext_grid_set = list(zip(self.bus_lookup[self.net.ext_grid.bus[self.ext_grid_set]], self.ext_grid_set))

        self.bus_line_dict = dict(
            zip(list(zip(self.line_set, [1] * len(self.line_set))) + list(zip(self.line_set, [2] * len(self.line_set))),
                np.concatenate([self.bus_lookup[self.net.line.from_bus[self.line_set].values],
                                self.bus_lookup[self.net.line.to_bus[self.line_set].values]])))

        # --- Param Data ---
        self.baseMVA = self.net.sn_mva

        self.PG_data = self.generators.p_mw * self.generators.scaling / self.baseMVA  # active power generation set points

        self.PD_data = self.net.load.p_mw * self.net.load.scaling / self.baseMVA

        self.GB_data = self.net.shunt.p_mw * self.net.shunt.step / self.baseMVA

        self.delta_eG_data = self.net.ext_grid.va_degree * pi / 180

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
        shift = self.net._ppc['branch'][trafo_start:trafo_end, 9].real * pi/180
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
        self.model.B = Set(initialize=self.bus_set)
        self.model.eG = Set(initialize=self.ext_grid_set)  # external grids
        self.model.G = Set(initialize=self.gen_all_set)  # static generators and generators
        self.model.sG = Set(within=self.model.G, initialize=self.sgen_set)  # static generators
        self.model.gG = Set(within=self.model.G, initialize=self.gen_set)  # generators (not static)
        self.model.D = Set(initialize=self.demand_set)
        # self.model.L = Set(initialize=self.line_set)
        self.model.L = Set(initialize=self.line_data.index[self.line_data.in_service])
        self.model.SHUNT = Set(initialize=self.shunt_set)
        self.model.LE = Set(initialize=[1, 2])
        self.model.TRANSF = Set(initialize=self.trafo_data.index[self.trafo_data.in_service])

        # generators, buses, loads linked to each bus b
        self.model.Gbs = Set(within=self.model.B * self.model.G,
                             initialize=self.bus_gen_set)  # set of generator-bus mapping
        self.model.Dbs = Set(within=self.model.B * self.model.D,
                             initialize=self.bus_demand_set)  # set of demand-bus mapping
        self.model.SHUNTbs = Set(within=self.model.B * self.model.SHUNT,
                                 initialize=self.bus_shunt_set)  # set of shunt-bus mapping
        self.model.eGbs = Set(within=self.model.B * self.model.eG,
                              initialize=self.bus_ext_grid_set)  # set of external grid-bus mapping

        # --- parameters ---
        # line and trafo matrix
        self.model.A = Param(self.model.L * self.model.LE, initialize=self.bus_line_dict)  # bus-line matrix
        self.model.AT = Param(self.model.TRANSF * self.model.LE,
                              initialize=self.bus_trafo_dict)  # bus-transformer matrix

        # generation
        self.model.PG = Param(self.model.G, initialize=self.PG_data[self.model.G])

        # demand
        self.model.PD = Param(self.model.D, initialize=self.PD_data[self.model.D])

        # shunt
        self.model.GB = Param(self.model.SHUNT, within=Reals,
                              initialize=self.GB_data[self.model.SHUNT])  # shunt conductance

        # trafo
        # self.model.shift = Param(self.model.TRANSF, within=Reals,
        #                          initialize=self.shift_rad_data[self.model.TRANSF])  # transformer phase shift in rad
        self.model.shift = Param(self.model.TRANSF, within=Reals,
                                 initialize=self.trafo_data.shift_rad[self.model.TRANSF])  # transformer phase shift in rad

        # external grid voltage angle
        self.model.delta_eG = Param(self.model.eG, within=Reals, initialize=self.delta_eG_data[self.model.eG])

        # baseMVA of the net
        self.model.baseMVA = Param(within=NonNegativeReals, initialize=self.baseMVA)

        # --- variables ---
        self.model.delta = Var(self.model.B, domain=Reals, initialize=0.0,
                               bounds=(-pi, pi))  # voltage phase angle at bus b, rad
        self.model.pD = Var(self.model.D, domain=Reals)  # real power demand delivered
        self.model.pG = Var(self.model.G, domain=NonNegativeReals)  # real generator power
        self.model.peG = Var(self.model.eG, domain=Reals)  # real power injection from external grids
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
        for g in self.model.G:
            self.model.pG[g].fix(self.model.PG[g])

        # --- demand ---
        for d in self.model.D:
            self.model.pD[d].fix(self.model.PD[d])

        # --- reference bus constraint ---
        for (b, g) in self.model.eGbs:
            self.model.delta[b].fix(self.model.delta_eG[g])
        # def ref_bus_def(model, b, g):
        #     return model.delta[b] == self.model.delta_eG[g]
        #
        # self.model.refbus_delta = Constraint(self.model.eGbs, rule=ref_bus_def)

    def solve(self, to_net: bool = True, print_solver_output: bool = False, solver='ipopt', load_solutions: bool = True,
              mip_solver='gurobi', max_iter=None, time_limit=600):
        optimizer = SolverFactory(solver)

        if solver == 'mindtpy':
            if not max_iter:
                max_iter = 50

            if mip_solver == 'gurobi':
                mip_solver = 'gurobi_persistent'

            try:
                self.results = optimizer.solve(self.model, mip_solver=mip_solver, nlp_solver='ipopt',
                                               tee=print_solver_output, iteration_limit=max_iter, time_limit=time_limit)
            except ValueError as err:
                print(err)

            # if self.results.solver.termination_condition == TerminationCondition.feasible:
            #     print('Model is feasible but not optimal. Trying to solve a second time.')
            #     self.results = optimizer.solve(self.model, mip_solver=mip_solver, nlp_solver='ipopt',
            #                                    tee=print_solver_output, iteration_limit=max_iter, time_limit=time_limit)
        else:
            if max_iter:
                optimizer.options['max_iter'] = max_iter

            self.results = optimizer.solve(self.model, load_solutions=load_solutions, tee=print_solver_output)

        if check_optimal_termination(self.results) & to_net:
            pyo_sol_to_net_res(self.net, self.model)

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
                if value:
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
