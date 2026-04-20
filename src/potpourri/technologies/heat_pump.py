"""Heat pump mix-in: attaches heat pump storage sets, parameters, variables, and constraints to a multi-period model."""

import numpy as np
from pyomo.environ import *
from src.potpourri.technologies.flexibility import Flexibility_multi_period


class Heatpump_multi_period(Flexibility_multi_period):
    """Multi-period heat pump device module with thermal building model and scenario-based penetration."""

    def __init__(self, net, T=None, scenario=None):
        super().__init__(net, T, scenario)

        if scenario == 0:
            self.hp_percentage = 6.3
        elif scenario == 1:
            self.hp_percentage = 26.4
        elif scenario == 2:
            self.hp_percentage = 44.6
        elif scenario == 3:
            self.hp_percentage = 59.9

        # randomly select buses for heat pump placement based on the percentage of scenario
        num_indexes = round(len(self.buses_excl_extGrids) * self.hp_percentage / 100)
        self.random_indexes = np.random.choice(self.buses_excl_extGrids, num_indexes, replace=False)

        # --- house calculation ---
        self.temp_max = 20  # average max temperature in K
        self.temp_min = 15  # average min temperature in K
        self.temp_outside = 10  # average outside temperature in K
        self.avg_house_size = 500  # average house size in m^3
        self.air_density = 1.2  # air density in kg/m^3
        self.air_capacity = 1.005  # air capacity in J/(kg*K)
        self.concrete_density = 2400  # concrete density in kg/m^3
        self.concrete_capacity = 0.84
        self.wall_area = self.avg_house_size * (1 / 3)  # wall area in m^2 (house is a cube)
        self.wall_thickness = 0.2  # wall thickness in m

        # --- mass and heat_cap calculation ---
        self.air_mass = self.avg_house_size * self.air_density  # air mass in kg
        self.concrete_mass = 4 * self.wall_area * self.wall_thickness * self.concrete_density
        # heat capacity of the house in MWh/K
        self.heat_cap = (self.air_mass * self.air_capacity + self.concrete_mass * self.concrete_capacity) / 3600000

        # TODO: Translate to temperature parameters for heat pump

        # --- heat pump parameter values ---
        self.hp_power_max = 0.005  # in MW_heat
        self.hp_cop = 4  # coefficient of performance in MW_el/MW_heat

        # --- heat loss calculation ---
        # max load out of profiles
        self.max_load = np.max(net.profiles[('load', 'p_mw')][0] / self.baseMVA)
        self.Qloss_max = 0.01
        self.heat_scaling_fac = self.Qloss_max / self.max_load
        self.heat_load = self.heat_scaling_fac * self.net.profiles[('load', 'p_mw')][0]

    def get_all(self, model):
        """Attach heat pump sets, parameters, variables, and constraints to the model."""
        self.get_sets(model)
        self.get_parameters(model)
        self.get_variables(model)
        self.get_all_constraints(model)

    def get_sets(self, model):
        """Define HP and HP_bus sets from randomly placed heat pumps."""
        super().get_sets(model)
        model.HP = Set(initialize=list(range(len(self.random_indexes))))  # set of heat pumps
        model.HP_bus = Set(initialize=list(enumerate(self.random_indexes)))  # set of heat pump-bus mapping
        return True

    def get_parameters(self, model):
        """Attach heat pump power limits, temperature bounds, and heat loss profile."""
        model.HP_Pmax = Param(model.HP, within=Reals, initialize=self.hp_power_max)
        model.HP_Pmin = Param(model.HP, within=Reals, initialize=0)
        model.TempMax = Param(model.HP, within=Reals, initialize=self.temp_max)
        model.TempMin = Param(model.HP, within=Reals, initialize=self.temp_min)

        self.Qloss_data_dict, self.Qloss_tuple = self.make_to_dict(model.HP, model.T, self.heat_load, True)
        model.Qloss = Param(self.Qloss_tuple, initialize=self.Qloss_data_dict)
        return True

    def get_variables(self, model):
        """Create hp_p (electrical power) and temp (indoor temperature) variables."""
        model.hp_p = Var(model.HP, model.T, within=Reals)
        model.temp = Var(model.HP, model.T, within=Reals)
        return True

    def get_all_constraints(self, model):
        """Add power limits, temperature bounds, and thermal state update constraints."""
        def hp_power_rule(model, h, t):
            return model.HP_Pmin[h], model.hp_p[h, t], model.HP_Pmax[h]
        model.temperature_rule = Constraint(model.HP, model.T, rule=hp_power_rule)

        def hp_temp_rule(model, h, t):
            return model.TempMin[h], model.temp[h, t], model.TempMax[h]
        model.hp_temp_con = Constraint(model.HP, model.T, rule=hp_temp_rule)

        def hp_temp_update_rule(model, h, t):
            if t == model.T.at(1):
                return Constraint.Skip
            return (model.temp[h, t] == model.temp[h, t - 1]
                    + model.deltaT * ((model.hp_p[h, t] * self.hp_cop) - model.Qloss[h, t]) / self.heat_cap)
        model.hp_temp_update_con = Constraint(model.HP, model.T, rule=hp_temp_update_rule)
        return True

    def get_all_opf(self, model):
        pass

    def get_all_ac(self, model):
        pass

    def get_all_acopf(self, model):
        pass
