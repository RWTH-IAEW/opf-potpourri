"""Generator mix-in: attaches external grid and generator variables and OPF limits to a multi-period model."""

import pandas as pd
from pyomo.environ import *
import numpy as np

from src.potpourri.technologies.flexibility import Flexibility_multi_period


class Generator_multi_period(Flexibility_multi_period):
    """Multi-period generator device module for external grids and PV generators."""

    def __init__(self, net, T=None, scenario=None):
        super().__init__(net, T, scenario)

        # --- generation ---
        pg = self.net._ppc['gen'][:, 1] / self.baseMVA  # active power generation
        ref_gens = self.net._ppc['internal']['ref_gens']  # reference generators indices
        in_service_gens = self.net._ppc['gen'][:, 7].astype(bool)  # in service generators
        gen_bus = self.net._ppc['gen'][:, 0].astype(int)  # generator buses
        self.generation_data = pd.DataFrame({'pg': pg, 'ref': False,
                                             'in_service': in_service_gens, 'bus': gen_bus})  # generation data
        self.generation_data.loc[ref_gens, 'ref'] = True  # reference generators are marked as such
        gen_bus_tuples = list(enumerate(self.generation_data['bus']))  # generator bus tuples for indexing
        self.generation_data['gen_bus'] = gen_bus_tuples  # generator bus tuples for indexing

        # generator and external grids voltage set points
        self.generation_data['v'] = self.net._ppc['gen'][:, 5]

    def get_all(self, model):
        """Attach generator sets, parameters, variables and fix non-slack generators."""
        self.get_sets(model)
        self.get_parameters(model)
        self.get_variables(model)
        self.fix_variables(model)

    def get_all_opf(self, model):
        self.get_opf_parameters(model)
        self.get_all_Constraints_opf(model)

    def get_all_acopf(self, model):
        self.get_acopf_parameters(model)
        self.get_all_Constraints_acopf(model)
        pass

    def get_sets(self, model):

        # external grids and generators
        model.G = Set(initialize=self.generation_data.index[self.generation_data.in_service])

        # external grids and slack generators
        model.eG = Set(initialize=self.generation_data.index[self.generation_data.ref], within=model.G)

        # generators (not static) not slack generators
        model.gG = Set(initialize=self.generation_data.index[self.generation_data.ref == False], within=model.G)

        # generators linked to each bus b
        model.Gbs = Set(within=model.G * model.B, initialize=self.generation_data['gen_bus'][model.G])

    def get_parameters(self, model):
        """Create time-indexed PG parameter from generation_data['pg']."""
        self.PG_data_dict, self.PG_tuple = self.make_to_dict(model.G, model.T, self.generation_data['pg'], False)
        model.PG = Param(self.PG_tuple, initialize=self.PG_data_dict)

    def get_opf_parameters(self, model):
        """Create PGmax and PGmin parameters for OPF real-power limits."""

        self.Pgmax_data_dict, self.Pgmax_tuple = self.make_to_dict(model.G, model.T,
                                                                   self.generation_data['max_p'], False)
        self.PGmin_data_dict, self.PGmin_tuple = self.make_to_dict(model.G, model.T,
                                                                   self.generation_data['min_p'], False)

        # generation real power limits multi period, TODO not necessarily time dependent
        model.PGmax = Param(self.Pgmax_tuple, initialize=self.Pgmax_data_dict)
        model.PGmin = Param(self.PGmin_tuple, initialize=self.PGmin_data_dict)

    def get_acopf_parameters(self, model):
        """Create QGmax and QGmin parameters for AC OPF reactive-power limits."""
        # create dict and tuple for reactive power limits multi period
        self.QGmax_data_dict, self.QGmax_tuple = self.make_to_dict(model.G, model.T, self.generation_data['max_q'], False)
        self.QGmin_data_dict, self.QGmin_tuple = self.make_to_dict(model.G, model.T, self.generation_data['min_q'], False)

        # generation reactive power limits multi period
        model.QGmax = Param(self.QGmax_tuple, initialize=self.QGmax_data_dict)
        model.QGmin = Param(self.QGmin_tuple, initialize=self.QGmin_data_dict)

    def get_variables(self, model):
        """Create pG variable for real power injection over all generators and time steps."""
        model.pG = Var(self.PG_tuple, domain=Reals)  # real power injection from static generators

    def fix_variables(self, model):
        """Fix non-slack (non-reference) generators to their profile values."""
        for g in model.gG:
            for t in model.T:
                model.pG[(g, t)].fix(model.PG[(g, t)])


    def generation_real_power_limits_opf(self, model):
        """Compute per-generator real-power limits and store in generation_data."""
        max_p = np.full(len(self.generation_data), 1e9) / self.baseMVA
        min_p = np.full(len(self.generation_data), -1e9) / self.baseMVA

        for element, (f, t) in self.net._gen_order.items():
            if 'max_p_mw' in self.net[element]:
                max_p[f:t] = self.net[element].max_p_mw.fillna(1e9).values / self.baseMVA
            if 'min_p_mw' in self.net[element]:
                min_p[f:t] = self.net[element].min_p_mw.fillna(-1e9).values / self.baseMVA

        if 'controllable' in self.net.gen:
            controllable = self.net["gen"]["controllable"].values
            not_controllable = ~controllable.astype(bool)

            if np.any(not_controllable):
                f, t = self.net._gen_order['gen']

                p_mw = self.net["gen"]["p_mw"].values[not_controllable]

                not_controllable_gens = np.arange(f, t)[not_controllable]
                max_p[not_controllable_gens] = p_mw / self.baseMVA
                min_p[not_controllable_gens] = p_mw / self.baseMVA

        self.generation_data['max_p'] = max_p
        self.generation_data['min_p'] = min_p

    def generation_reactive_power_limits_acopf(self):
        """Compute per-generator reactive-power limits and store in generation_data."""
        max_q = np.full(len(self.generation_data), 1e9) / self.baseMVA
        min_q = np.full(len(self.generation_data), -1e9) / self.baseMVA

        for element, (f, t) in self.net._gen_order.items():
            if 'max_q_mvar' in self.net[element]:
                max_q[f:t] = self.net[element].max_q_mvar.fillna(1e9) / self.baseMVA
            if 'min_q_mvar' in self.net[element]:
                min_q[f:t] = self.net[element].min_q_mvar.fillna(-1e9) / self.baseMVA

        self.generation_data['max_q'] = max_q
        self.generation_data['min_q'] = min_q

    # TODO make non time dependent
    def get_all_Constraints_opf(self, model):
        """Add real-power bound constraints for all generators over all time steps."""
        # --- generation real power limits ---
        #PG Constraint
        @model.Constraint(model.G, model.T)
        def real_power_bounds(model, g, t):
            model.pG[(g, t)].unfix()
            return model.PGmin[(g, t)], model.pG[(g, t)], model.PGmax[(g, t)]

    def get_all_Constraints_acopf(self, model):
        """Add reactive-power bound constraints for all generators over all time steps."""
        # --- reactive generator power limits ---
        @model.Constraint(model.G, model.T)
        def reactive_power_bounds(model, g, t):
            model.qG[(g, t)].unfix()
            return model.QGmin[(g, t)], model.qG[(g, t)], model.QGmax[(g, t)]