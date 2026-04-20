"""Electric vehicle mix-in: attaches simplified EV charging sets, parameters, and constraints to a multi-period model."""

import numpy as np
from pyomo.environ import *
from loguru import logger
from src.potpourri.technologies.flexibility import Flexibility_multi_period


# noinspection PyAttributeOutsideInit
class Electric_vehicle_multi_period(Flexibility_multi_period):
    """Multi-period electric vehicle device module with scenario-based penetration and simple charging model."""

    def __init__(self, net, T=None, scenario=None):
        super().__init__(net, T, scenario)

        if scenario == 0:
            self.ev_percentage = 0
        elif scenario == 1:
            self.ev_percentage = 52.8
        elif scenario == 2:
            self.ev_percentage = 91
        elif scenario == 3:
            self.ev_percentage = 134.3

        # randomly select buses for electric vehicle placement based on the percentage of scenario
        num_indexes = round(len(self.buses_excl_extGrids) * self.ev_percentage / 100)
        self.random_indexes = np.random.choice(self.buses_excl_extGrids, num_indexes, replace=False)

        self.ev_power_to_bat = 0.011  # power of electric vehicle in MW (AC, 11 kW)
        self.ev_power_to_grid = -0.011
        self.ev_soc_max = 1
        self.ev_soc_min = 0.2
        self.ev_batt_eff = 0.95
        self.ev_cap = 0.08  # capacity of electric vehicle in MWh

        # night steps per year calculation:
        days_per_year = 365
        steps_per_day = 96
        steps_in_night = list(range(0, 24)) + list(range(88, 96))  # steps corresponding to night hours (0-6 and 22-24)
        self.load_set = []
        for day in range(days_per_year):
            for step in range(steps_per_day):
                if step in steps_in_night:
                    self.load_set.append(True)
                else:
                    self.load_set.append(False)

        self.not_load_set = [not i for i in self.load_set]
        self.drive_away_list = []

        logger.debug("night_indicator array constructed")

    def get_all(self, model):
        """Attach EV sets, parameters, variables, and constraints to the model."""
        self.create_drive_profile(model)
        self.get_sets(model)
        self.get_parameters(model)
        self.get_variables(model)
        self.get_all_constraints(model)
        return True

    def create_drive_profile(self, model):
        """Build plugged-in and drive-away binary profiles for the optimisation horizon."""
        self.plugged_in_set = []
        for t in model.T:
            self.plugged_in_set.append(self.load_set[t])

        self.not_plugged_in_set = [not i for i in self.plugged_in_set]
        self.drive_away_set = []

        for t in range(len(self.plugged_in_set)):
            # end of set
            if t == len(self.plugged_in_set) - 1:
                self.drive_away_set.append(1)
            elif (self.plugged_in_set[t], self.plugged_in_set[t + 1]) == (True, False):
                self.drive_away_set.append(1)
            else:
                self.drive_away_set.append(0)

    def get_sets(self, model):
        """Define EV and EV_bus sets from randomly placed electric vehicles."""
        super().get_sets(model)
        model.EV = Set(initialize=list(range(len(self.random_indexes))))
        model.EV_bus = Set(initialize=list(enumerate(self.random_indexes)))
        return True

    def get_parameters(self, model):
        """Attach EV power limits, SOC bounds, capacity, efficiency, and drive profiles."""
        model.EV_P_to_grid = Param(model.EV, within=Reals, initialize=self.ev_power_to_grid)
        model.EV_P_to_bat = Param(model.EV, within=Reals, initialize=-self.ev_power_to_bat)
        model.EV_SOCmax = Param(model.EV, within=Reals, initialize=self.ev_soc_max)
        model.EV_SOCmin = Param(model.EV, within=Reals, initialize=self.ev_soc_min)
        model.EV_Cap = Param(model.EV, within=Reals, initialize=self.ev_cap)
        model.EV_batt_eff = Param(model.EV, within=Reals, initialize=self.ev_batt_eff)

        self.plugged_in_dict = {t: value for t, value in zip(model.T, self.plugged_in_set)}
        model.Plugged_in_profile = Param(model.T, initialize=self.plugged_in_dict)
        self.not_plugged_in_dict = {t: value for t, value in zip(model.T, self.not_plugged_in_set)}
        model.Not_plugged_in_profile = Param(model.T, initialize=self.not_plugged_in_dict)
        self.drive_away_dict = {t: value for t, value in zip(model.T, self.drive_away_set)}
        model.Drive_profile = Param(model.T, initialize=self.drive_away_dict)
        return True

    def get_variables(self, model):
        """Create EV_P (power) and EV_SOC (state of charge) variables."""
        model.EV_P = Var(model.EV, model.T, within=Reals)
        model.EV_SOC = Var(model.EV, model.T, within=Reals)
        return True

    def get_all_constraints(self, model):
        """Add power limits, SOC limits, and SOC update constraints for all EVs."""
        def ev_power_rule(model, e, t):
            return model.EV_P_to_grid[e], model.EV_P[e, t], model.EV_P_to_bat[e]
        model.ev_power_constr = Constraint(model.EV, model.T, rule=ev_power_rule)

        def ev_soc_rule(model, e, t):
            if t == model.T.at(1):
                return model.EV_SOC[e, t] == 0.5
            return model.EV_SOCmin[e], model.EV_SOC[e, t], model.EV_SOCmax[e]
        model.ev_soc_constr = Constraint(model.EV, model.T, rule=ev_soc_rule)

        def ev_soc_update_rule(model, e, t):
            if t == model.T.at(1):
                return Constraint.Skip
            return (model.EV_SOC[e, t] == model.EV_SOC[e, t - 1]
                    + model.deltaT * (model.EV_P[e, t] * model.EV_batt_eff[e]) / model.EV_Cap[e])
        model.ev_soc_update_constr = Constraint(model.EV, model.T, rule=ev_soc_update_rule)

    def get_all_opf(self, model):
        pass

    def get_all_ac(self, model):
        pass

    def get_all_acopf(self, model):
        pass
