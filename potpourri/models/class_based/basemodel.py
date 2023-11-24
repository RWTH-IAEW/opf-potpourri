import pandas as pd
from pyomo.environ import *


class Basemodel:
    def __init__(self, net, controllable_sgens: bool = False):
        self.net = net

        # --- Set indices ---
        self.bus_set = net.bus.index
        self.sgen_set = net.sgen.index
        self.demand_set = net.load.index
        self.line_set = net.line.index
        self.shunt_set = net.shunt.index
        self.trafo_set = net.trafo.index
        self.ext_grid_set = net.ext_grid.index
        self.slack_set = net.ext_grid.bus

        self.bus_gen_set = list(zip(net.sgen.bus, self.sgen_set))
        self.bus_demand_set = list(zip(net.load.bus, self.demand_set))
        self.bus_shunt_set = list(zip(net.shunt.bus, self.shunt_set))
        self.bus_ext_grid_set = list(zip(net.ext_grid.bus, self.ext_grid_set))

        self.bus_line_dict = dict(
            zip(list(zip(self.line_set, [1] * len(self.line_set))) + list(zip(self.line_set, [2] * len(self.line_set))),
                pd.concat([net.line.from_bus, net.line.to_bus])))
        self.bus_trafo_dict = dict(zip(list(zip(self.trafo_set, [1] * len(self.trafo_set))) + list(
            zip(self.trafo_set, [2] * len(self.trafo_set))), pd.concat([net.trafo.hv_bus, net.trafo.lv_bus])))

        # --- Param Data ---
        self.baseMVA = net.sn_mva

        self.demand_data = net.load.p_mw / self.baseMVA

        self.GB_data = net.shunt.p_mw * net.shunt.step / self.baseMVA

        self.get_generator_active_data()
        self.get_demand_active_data()

    def get_generator_active_data(self):
        # generation set points
        self.PG_data = self.net.sgen.p_mw / self.baseMVA

        if 'controllable' not in self.net.sgen:
            self.sgen_controllable_set = None  # create empty Set if no controllable generators exist
        else:
            self.sgen_controllable_set = self.net.sgen.index[self.net.sgen.controllable]

        # add rows with active generation limits if not existing
        if 'max_p_mw' not in self.net.sgen:
            self.net.sgen['max_p_mw'] = self.net.sgen.p_mw

        if 'min_p_mw' not in self.net.sgen:
            self.net.sgen['min_p_mw'] = [0] * len(self.net.sgen.index)

        # generation limits for sgens
        self.PGmax_data = self.net.sgen.max_p_mw.fillna(self.net.sgen.p_mw) / self.baseMVA
        self.PGmin_data = self.net.sgen.min_p_mw.fillna(0) / self.baseMVA



    def get_demand_active_data(self):
        # active power demand
        self.PD_data = self.net.load.p_mw / self.baseMVA

        if 'controllable' not in self.net.load:
            self.demand_controllable_set = None  # create empty Set if no controllable load exist
        else:
            self.demand_controllable_set = self.net.load.index[self.net.load.controllable]

        # add rows with active demand limits if not existing
        if 'max_p_mw' not in self.net.load:
            self.net.load['max_p_mw'] = self.net.load.p_mw

        if 'min_p_mw' not in self.net.load:
            self.net.load['min_p_mw'] = [0] * len(self.net.load.index)

        # demand limits for loads
        self.PDmax_data = self.net.load.max_p_mw.fillna(self.net.load.p_mw) / self.baseMVA
        self.PDmin_data = self.net.load.min_p_mw.fillna(0) / self.baseMVA

    def create_model(self):
        self.model = ConcreteModel()

        # --- SETS ---
        self.model.B = Set(initialize=self.bus_set)
        self.model.eG = Set(within=self.model.B, initialize=self.ext_grid_set)     # external grids
        self.model.G = Set(initialize=self.sgen_set)
        self.model.Gc = Set(within=self.model.G, initialize=self.sgen_controllable_set)     # controllable static generators
        self.model.D = Set(initialize=self.demand_set)
        self.model.Dc = Set(within=self.model.D, initialize=self.demand_controllable_set)     # controllable loads
        self.model.L = Set(initialize=self.line_set)
        self.model.SHUNT = Set(initialize=self.shunt_set)
        self.model.LE = Set(initialize=[1, 2])
        self.model.TRANSF = Set(initialize=self.trafo_set)
        self.model.WIND = Set()  # set of wind generators

        # generators, buses, loads linked to each bus b
        self.model.Gbs = Set(within=self.model.B * self.model.G,
                             initialize=self.bus_gen_set)  # set of generator-bus mapping
        self.model.Dbs = Set(within=self.model.B * self.model.D,
                             initialize=self.bus_demand_set)  # set of demand-bus mapping
        self.model.SHUNTbs = Set(within=self.model.B * self.model.SHUNT,
                                 initialize=self.bus_shunt_set)  # set of shunt-bus mapping
        self.model.eGbs = Set(within=self.model.B * self.model.eG, initialize=self.bus_ext_grid_set) # set of external grid-bus mapping

        # --- parameters ---
        # line matrix
        self.model.A = Param(self.model.L * self.model.LE, initialize=self.bus_line_dict)  # bus-line matrix
        self.model.AT = Param(self.model.TRANSF * self.model.LE,
                              initialize=self.bus_trafo_dict)  # bus-transformer matrix

        # generation
        self.model.PG = Param(self.model.G, initialize=self.PG_data)
        self.model.PGmax = Param(self.model.G, initialize=self.PGmax_data)
        self.model.PGmin = Param(self.model.G, initialize=self.PGmin_data)

        # demand
        self.model.PD = Param(self.model.D, initialize=self.PD_data)
        self.model.PDmax = Param(self.model.D, initialize=self.PDmax_data)
        self.model.PDmin = Param(self.model.D, initialize=self.PDmin_data)

        self.model.VOLL = Param(self.model.D, within=Reals, initialize=10000)  # value of lost load

        # shunt
        self.model.GB = Param(self.model.SHUNT, within=Reals, initialize=self.GB_data)  # shunt conductance

        # baseMVA of the net
        self.model.baseMVA = Param(within=NonNegativeReals, initialize=self.baseMVA)

        # --- variables ---
        self.model.delta = Var(self.model.B, domain=Reals, initialize=0.0)  # voltage phase angle at bus b, rad
        self.model.pD = Var(self.model.D, domain=Reals)  # real power demand delivered
        self.model.pG = Var(self.model.G, domain=Reals)  # real generator power
        self.model.peG = Var(self.model.eG, domain=Reals)  # real power injection from external grids
        self.model.pLfrom = Var(self.model.L, domain=Reals)  # real power injected at b onto line
        self.model.pLto = Var(self.model.L, domain=Reals)  # real power injected at b' onto line
        self.model.pLfromT = Var(self.model.TRANSF, domain=Reals)  # real power injected at b onto transformer
        self.model.pLtoT = Var(self.model.TRANSF, domain=Reals)  # real power injected at b' onto transformer

        # --- generator power limits ---
        def real_power_bounds(model, g):
            if g in model.Gc:
                return model.PGmin[g], model.pG[g], model.PGmax[g]
            else:
                model.pG[g].fix(model.PG[g])
                return Constraint.Skip

        self.model.PG_Constraint = Constraint(self.model.G, rule=real_power_bounds)

        # --- demand limits ---
        def real_demand_bounds(model, d):
            if d in model.Dc:
                return model.PDmin[d], model.pD[d], model.PDmax[d]
            else:
                model.pD[d].fix(model.PD[d])
                return Constraint.Skip

        self.model.PD_Constraint = Constraint(self.model.D, rule=real_demand_bounds)

        # --- reference bus constraint ---
        def ref_bus_def(model, b):
            return model.delta[b] == 0

        self.model.refbus_delta = Constraint(self.model.eG, rule=ref_bus_def)

    def solve(self):
        solver = SolverFactory('ipopt')
        self.results = solver.solve(self.model, load_solutions=True)


