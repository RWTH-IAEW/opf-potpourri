"""Battery mix-in: attaches battery storage sets, parameters, variables, and constraints to a multi-period model."""

import numpy as np
from pyomo.environ import *
from src.potpourri.technologies.flexibility import Flexibility_multi_period


class Battery_multi_period(Flexibility_multi_period):
    """Multi-period battery storage device module with scenario-based penetration levels."""

    def __init__(self, net, T=None, scenario=None):
        super().__init__(net, T, scenario)

        if scenario == 0:
            self.bat_percentage = 1
        elif scenario == 1:
            self.bat_percentage = 7.9
        elif scenario == 2:
            self.bat_percentage = 9.9
        elif scenario == 3:
            self.bat_percentage = 10.6

        # randomly select buses for battery placement based on the percentage of scenario
        num_indexes = round(len(self.buses_excl_extGrids) * self.bat_percentage / 100)
        self.random_indexes = np.random.choice(self.buses_excl_extGrids, num_indexes, replace=False)

        # --- battery parameter values ---
        self.bat_power = 0.006
        self.bat_soc_max = 1
        self.bat_soc_min = 0.2
        self.bat_cap = 0.015

    def get_all(self, model):
        """Attach battery sets, parameters, variables, and constraints to the model."""
        self.get_sets(model)
        self.get_parameters(model)
        self.get_variables(model)
        self.get_all_constraints(model)

    def get_sets(self, model):
        """Define BAT and BAT_bus sets from randomly placed batteries."""
        super().get_sets(model)
        model.BAT = Set(initialize=list(range(len(self.random_indexes))))  # set of batteries
        model.BAT_bus = Set(initialize=list(enumerate(self.random_indexes)))  # set of battery-bus mapping
        return True

    def get_parameters(self, model):
        """Attach battery power, SOC, capacity, and efficiency parameters."""
        model.BAT_Pmax = Param(model.BAT, within=Reals, initialize=self.bat_power)
        model.BAT_Pmin = Param(model.BAT, within=Reals, initialize=-self.bat_power)
        model.BAT_SOCmax = Param(model.BAT, within=Reals, initialize=self.bat_soc_max)
        model.BAT_SOCmin = Param(model.BAT, within=Reals, initialize=self.bat_soc_min)
        model.BAT_Cap = Param(model.BAT, within=Reals, initialize=self.bat_cap)
        model.BAT_Eff = Param(model.BAT, within=Reals, initialize=0.9)
        return True

    def get_variables(self, model):
        """Create BAT_P (power) and BAT_SOC (state of charge) variables."""
        model.BAT_P = Var(model.BAT, model.T, within=Reals)
        model.BAT_SOC = Var(model.BAT, model.T, within=Reals)
        return True

    def get_all_constraints(self, model):
        """Add power limits, SOC limits, and SOC update constraints for all batteries."""
        def bat_power_rule(model, b, t):
            return model.BAT_Pmin[b], model.BAT_P[b, t], model.BAT_Pmax[b]
        model.bat_power_con = Constraint(model.BAT, model.T, rule=bat_power_rule)

        def bat_soc_rule(model, b, t):
            if t == model.T.at(1):
                return model.BAT_SOC[b, t] == 0.5 * model.BAT_SOCmax[b]
            return model.BAT_SOCmin[b], model.BAT_SOC[b, t], model.BAT_SOCmax[b]
        model.bat_soc_con = Constraint(model.BAT, model.T, rule=bat_soc_rule)

        def bat_soc_update_rule(model, b, t):
            if t == model.T.at(1):
                return Constraint.Skip
            return (model.BAT_SOC[b, t] == model.BAT_SOC[b, t - 1]
                    + model.deltaT * (model.BAT_P[b, t] * model.BAT_Eff[b]) / model.BAT_Cap[b])
        model.bat_soc_update_con = Constraint(model.BAT, model.T, rule=bat_soc_update_rule)
        return True

    def get_all_acopf(self, model):
        pass

    def get_all_ac(self, model):
        pass

    def get_all_opf(self, model):
        pass
