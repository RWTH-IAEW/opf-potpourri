"""PV mix-in: attaches PV generation sets, parameters, variables, and power-bound constraints to a multi-period model."""

import numpy as np
import pyomo.environ as pyo
from src.potpourri.technologies.flexibility import Flexibility_multi_period


class PV_multi_period(Flexibility_multi_period):
    """Multi-period PV device module with scenario-based penetration levels.

    Penetration percentages are taken from:
    "Anforderungen an aktuelle Verteilnetze und deren zukuenftige Versorgungsaufgabe".
    """

    def __init__(self, net, T=None, scenario=None):
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

        # take a random column of input sb net renewables as pv load profile
        self.pv_load_profiles = (-1) * (self.net.pv_load_profiles)
        self.pv_load_profile = self.pv_load_profiles["PV5"]

        self.pv_pmax = self.pv_load_profile
        self.pv_pmin = 0.0
        num_indexes = round(len(self.buses_excl_extGrids) * self.pv_percentage / 100)
        self.random_indexes = np.random.choice(self.buses_excl_extGrids, num_indexes, replace=False)

    def get_all(self, model):
        """Attach PV sets, parameters, variables, constraints, and unfix variables."""
        self.get_sets(model)
        self.get_parameters(model)
        self.get_variables(model)
        self.get_all_constraints(model)
        self.unfix_variables(model)

    def unfix_variables(self, model):
        """Unfix pPV variables for all PV units over the full time horizon."""
        for t in model.T:
            for pv in model.PV:
                model.pPV[pv, t].unfix()

    def get_sets(self, model):
        """Define PV and PV_bus sets from randomly placed PV units."""
        super().get_sets(model)
        model.PV = pyo.Set(initialize=list(range(len(self.random_indexes))))
        model.PV_bus = pyo.Set(initialize=list(enumerate(self.random_indexes)))
        return True

    def get_parameters(self, model):
        """Attach time-indexed PV_Pmax and PV_Pmin parameters."""
        self.PV_Pmax_dict, self.PV_Pmax_tuple = self.make_to_dict(model.PV, model.T, self.pv_pmax)
        self.PV_Pmin_dict, self.PV_Pmin_tuple = self.make_to_dict(model.PV, model.T, self.pv_pmin, False)

        model.PV_Pmax = pyo.Param(self.PV_Pmax_tuple, within=pyo.Reals, initialize=self.PV_Pmax_dict)
        model.PV_Pmin = pyo.Param(self.PV_Pmin_tuple, within=pyo.Reals, initialize=self.PV_Pmin_dict)

    def get_variables(self, model):
        """Create pPV variable initialised from the load profile."""
        self.pPV_data_dict, self.pPV_tuple = self.make_to_dict(model.PV, model.T, self.pv_load_profile)
        model.pPV = pyo.Var(self.pPV_tuple, within=pyo.Reals, initialize=self.pPV_data_dict)

    def get_all_constraints(self, model):
        """Add real-power bound constraints for all PV units over all time steps."""
        @model.Constraint(model.PV, model.T)
        def PV_real_power_bounds(model, pv, t):
            return model.PV_Pmax[pv, t], model.pPV[pv, t], model.PV_Pmin[pv, t]

    # empty placeholders, because they are called in get_all
    def get_all_acopf(self, model):
        pass

    def get_all_ac(self, model):
        pass

    def get_all_opf(self, model):
        pass
