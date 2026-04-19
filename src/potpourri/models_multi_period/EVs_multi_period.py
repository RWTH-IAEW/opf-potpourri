import copy
import random
import numpy as np
import pandas as pd
import pickle
from pathlib import Path
from pyomo.environ import *
from src.potpourri.models_multi_period.flexibility_multi_period import Flexibility_multi_period
from pyomo.common.errors import ApplicationError

_DATA_DIR = Path(__file__).parent.parent / "data"

# linearized model --> weakness: efficiencies(if not market oriented) , calculation of p for soc > soclim

class EV_multi_period(Flexibility_multi_period):
    def __init__(self, net, num_vehicles):
        """
        Args:
            vehicles: dictionary of the vehicles created in the charging profile generator simulation
            locations: dictionary of the locations created in the charging profile generator simulation
            chargingpoints: dictionary of the chargingpoints created in the charging profile generator simulation
            scenario: scenario storing the simulation configuration
        """
        super().__init__(net)

        # Extract the vehicles, locations, and charging points from the pkl file
        file_path = _DATA_DIR / "emob_profiles.pkl"

        with open(file_path, 'rb') as f:
            vehicles, locations, chargingpoints, scenario = pickle.load(f)

        def extract_vehicles(original_dict, num_vehicles, start_index=0):

            # Extract a subset of vehicles from a dictionary of arrays.

            new_dict = {}
            for key, array in original_dict.items():

                if isinstance(array, list):
                    new_dict[key] = array[start_index:start_index + num_vehicles]
                elif isinstance(array, np.ndarray):
                    if array.ndim == 2:
                        new_dict[key] = array[:num_vehicles, :]
                    elif array.ndim == 1:
                        new_dict[key] = array[:num_vehicles]
                    else:
                        raise ValueError(f"Array for key '{key}' is not 2D.")
                elif isinstance(array, pd.DataFrame):
                    new_dict[key] = array.iloc[:num_vehicles]
                else:
                    new_dict[key] = array
            return new_dict
        # Extract a subset of vehicles from the original dictionary
        vehicles = extract_vehicles(vehicles, num_vehicles, start_index=0)

        # Sets
        self.vehicles = copy.deepcopy(vehicles)
        self.locations = copy.deepcopy(locations)
        self.chargingpoints = copy.deepcopy(chargingpoints)
        self.scenario = copy.deepcopy(scenario)

        # VARIANTE 1: bus-cp mapping, connect cps to random buses in the given net
        random.seed(88)
        buses_excl_ext_grid = list(net.bus.loc[~net.bus.index.isin(net.ext_grid.bus)].index)
        self.chargingpoints['bus'] = np.full(self.chargingpoints['n'], 0, dtype=int)
        for i in range(self.chargingpoints['n']):
            self.chargingpoints['bus'][i] = random.choice(buses_excl_ext_grid)

        self.bus_cp_set = list(zip(self.bus_lookup[self.chargingpoints['bus']], self.chargingpoints['id']))

        # VARIANTE 2: bus-cp mapping, look up ev buses in simbench
        # buses_ev_home = net.load.bus[net.load[net.load.profile.str.startswith("H0")].index].values
        # buses_ev_work = net.load.bus[net.load[net.load.profile.str.startswith("APLS")].index].values
        # cps_home = self.chargingpoints['loc_ids'] == 3
        # cps_work = self.chargingpoints['loc_ids'] == 0
        # # assign home charging points to buses with household loads but don't repeat buses
        # for i in range(self.chargingpoints['n']):
        #     if i in cps_home:
        #         self.chargingpoints['bus'][i] = random.choice(buses_ev_home)
        #     elif i in cps_work:
        #         self.chargingpoints['bus'][i] = random.choice(buses_ev_work)


        # driving profile parameters:
        self.distance = self.vehicles['distance']  # in km
        self.consumption = self.vehicles['consumption'] / 100000  # in MWh/km
        self.t_arr = self.vehicles['t_arr']
        self.t_dep = self.vehicles['t_dep']
        self.t_connect = self.vehicles['t_connect']
        self.t_disconnect = self.vehicles['t_disconnect']
        self.veh_cps = self.vehicles['cps']
        self.param = self.vehicles['param1']*1000  # convert from kWh to MWh
        self.cp_locids = self.chargingpoints['loc_ids']  # location id of the chargingpoint
        self.timestep_size = self.scenario['config']['timestep_size']  # in minutes
        self.loctypes = self.locations['type']
        self.cptype = self.chargingpoints['type']


        # charging profile parameters
        self.soc = {(i, t): self.vehicles['soc'][i, t] for i in vehicles['id'] for t in range(self.scenario['config']['n_t_steps'])}
        self.soclim = self.scenario['config']['soclim']
        self.e_max = self.vehicles['e_max'] / 1000  # battery capacity of the vehicle (in MWh)
        self.p = self.vehicles['p'] / 1000  # current power demand vehicle
        self.p_req_ev = self.vehicles['p_req'] / 1000  # required power for the vehicle
        self.p_dischargable = self.vehicles['p_dischargable'] / 1000  # dischargable power of the vehicle

        self.p_max_veh = self.vehicles['p_max'] / 1000  # max. charging power vehicle
        self.p_max_cp = self.chargingpoints['p_max'] / 1000  # max. charging power chargingpoint
        self.p_sk = self.vehicles['p_sk'] / 1000
        self.p_ac = self.vehicles['p_ac'] / 1000
        self.p_dc = self.vehicles['p_dc'] / 1000
        self.p_req_cp = self.chargingpoints['p_req'] / 1000
        self.p_cp = self.chargingpoints['p'] / 1000
        self.eta = self.chargingpoints['eta']  # charging efficiency of the chargingpoints
        self.soc_lim = self.vehicles['soc_lim']  # soc limit for the vehicle
        self.tau = self.vehicles['tau']  # time constant for the exponential decrease of charging power

        # derived parameters
        ev_connected = self.veh_cps != -1  # vehicle is not connected to any charging points at timestep t
        self.ev_connected = np.where(ev_connected, 1, 0)  # binary indicating if vehicle is connected to a chargingpoint at timestep t
        self.v2g = self.vehicles['v2g']  # bool indicating if vehicle has V2G capability
        self.p_drive = np.where(self.p < 0, self.p, 0)

        # self.controllable = self.vehicles['controllable']  # bool indicating if vehicle is controllable


        self.v2g_vehicles = np.where(np.array(self.v2g))[0]
        self.v1g_vehicles = np.where(~(np.array(self.v2g)))[0]
        self.nc_vehicles = []  # non controllable vehicles

        # market parameters
        file_path = _DATA_DIR / "da_prices_hourly_2022.xlsx"
        # self.p_id = pd.read_excel(file_path, sheet_name='Intraday')['Price (€/MWh)'].values
        # self.p_da_hourly = pd.read_excel(file_path, sheet_name='Day Ahead')['Price (€/MWh)'].values
        self.p_da_hourly = pd.read_excel(file_path)['Deutschland/Luxemburg [€/MWh]'].values
        self.p_da = np.repeat(self.p_da_hourly, 4)  # repeat the hourly prices 4 times to match the timesteps

    def get_opf_sets(self, model):
        # --SETS--
        model.veh = Set(initialize=self.vehicles['id'])  # vehicles
        model.cp = Set(initialize=self.chargingpoints['id'])  # charging points
        model.Bcp = Set(within=model.B * model.cp, initialize=self.bus_cp_set)  # buses of the charging points

        model.veh_v2g = Set(within=model.veh, initialize=self.v2g_vehicles)  # bidirectional capable vehicles
        model.veh_v1g = Set(within=model.veh, initialize=self.v1g_vehicles)  # unidirectional capable vehicles
        model.veh_nc = Set(within=model.veh, initialize=self.nc_vehicles)  # vehicles with no controlled charging capability (p fixed to cpg reference charging)

    def get_opf_parameters(self, model):
        # --- parameters ---
        # non-time dependent parameters
        model.timestep_size = self.timestep_size
        model.p_max_veh = Param(model.veh, mutable=True, initialize=self.p_max_veh[model.veh])
        model.p_max_cp = Param(model.cp, mutable=True, initialize=self.p_max_cp)
        model.p_sk = Param(model.veh, initialize=self.p_sk[model.veh])
        model.p_ac = Param(model.veh, mutable=True, initialize=self.p_ac[model.veh])
        model.p_dc = Param(model.veh, initialize=self.p_dc[model.veh])
        model.e_max = Param(model.veh, initialize=self.e_max[model.veh])
        model.cp_locids = Param(model.cp, initialize=self.cp_locids)
        model.cp_type = Param(model.cp, initialize=self.cptype)
        model.eta = Param(initialize=self.eta)
        model.consumption = Param(model.veh, initialize=self.consumption[model.veh])
        model.param = Param(model.veh, mutable=True, initialize=self.param[model.veh])
        model.soc_lim = Param(initialize=self.soc_lim)
        model.tau = Param(model.veh, initialize=self.tau[model.veh])
        model.min_soc_parking = Param(initialize=0.3)
        model.min_soc_departing = Param(initialize=0.8)  #mutable=True,
        model.bigM = Param(initialize=1e6)

        # time dependent parameters
        veh_cps_dict = {(i, t): self.veh_cps[i, t] for i in model.veh for t in model.T}  # flatten array for easier mapping
        model.veh_cps = Param(model.veh, model.T, initialize=veh_cps_dict)
        p_dict = {(i, t): self.p[i, t] for i in model.veh for t in model.T}
        model.p = Param(model.veh, model.T, initialize=p_dict)
        p_dischargable_dict = {(i, t): self.p_dischargable[i, t] for i in model.veh for t in model.T}
        model.p_dischargable_ev = Param(model.veh, model.T, initialize=p_dischargable_dict)
        ev_connected_dict = {(i, t): self.ev_connected[i, t] for i in model.veh for t in model.T}
        model.ev_connected = Param(model.veh, model.T, initialize=ev_connected_dict)
        soc_init_dict = {(i, t): self.soc[i, t] for i in model.veh for t in model.T}
        model.soc_init = Param(model.veh, model.T, initialize=soc_init_dict)
        p_req_ev_dict = {(i, t): self.p_req_ev[i, t] for i in model.veh for t in model.T}
        model.p_req_ev = Param(model.veh, model.T, initialize=p_req_ev_dict)
        p_drive_dict = {(i, t): self.p_drive[i, t] for i in model.veh for t in model.T}
        model.p_drive = Param(model.veh, model.T, initialize=p_drive_dict)

        # trip no. dependent parameters
        model.t_dep = Param(model.veh, initialize=self.t_dep)
        model.t_arr = Param(model.veh, initialize=self.t_arr)
        model.distance = Param(model.veh, initialize=self.distance)

        # prices on day ahead market
        model.p_da = Param(model.T, initialize=lambda model, t: self.p_da[t])

    def get_opf_variables(self, model):
        # --- variables ---
        model.soc = Var(model.veh, model.T, bounds=(0, 1))  # vehicle state of charge
        model.p_opf = Var(model.veh, model.T)  # EV controlled charging power
        model.p_charging = Var(model.veh, model.T, within=NonNegativeReals)  # charging power
        model.p_discharging = Var(model.veh, model.T, within=NonPositiveReals)  # discharging power
        model.buffer_soc = Var(model.veh, model.T, within=NonNegativeReals)  # buffer soc for the min soc constraints

    def get_opf_constraints(self, model):
        # ---mobility constraints ---

        #  charging and discharging power limitations
        def lower_power_bounds_v2g(model, v, t):
            return model.p_dischargable_ev[v, t] * -1 <= model.p_discharging[v, t]
        model.discharging_power_lb = Constraint(model.veh_v2g, model.T, rule=lower_power_bounds_v2g)

        def lower_power_bounds_v1g(model, v, t):
            return model.p_discharging[v, t] == 0
        model.charging_power_lb_v1g = Constraint(model.veh_v1g, model.T, rule=lower_power_bounds_v1g)

        def upper_power_bounds_cp(model, v, t):
            # 1.1 maximum power of vehicle depending on soc
            # p_max_veh = model.p_max_veh[v]

            # 2. maximum power of vehicle depending on connection type
            if model.veh_cps[v, t] != -1:
                veh_cp_type = model.cp_type[model.veh_cps[v, t]]
            else:
                veh_cp_type = 'None'

            if veh_cp_type == 'SK':
                p_max_veh_cptype = model.p_sk[v]
            elif veh_cp_type == 'AC':
                p_max_veh_cptype = model.p_ac[v]
            elif veh_cp_type == 'DC':
                p_max_veh_cptype = model.p_dc[v]
            else:
                p_max_veh_cptype = 0

            # Resulting maximum power is the minimum of all max powers
            p_max_charging = p_max_veh_cptype
            return model.p_charging[v, t] <= p_max_charging * model.ev_connected[v, t]
        model.charging_power_ub = Constraint(model.veh, model.T, rule=upper_power_bounds_cp)

        def upper_power_bounds_cp2(model, v, t):
            # 3. maximum cp power
            if model.veh_cps[v, t] == -1:
                p_max_cp = 0
            else:
                p_max_cp = model.p_max_cp[model.veh_cps[v, t]]
            return model.p_charging[v, t] <= p_max_cp * model.ev_connected[v, t]
        model.charging_power_ub3 = Constraint(model.veh, model.T, rule=upper_power_bounds_cp2)
        # temp for mutable param
        def upper_power_bounds_veh(model, v, t):
            # 1.1 maximum power of vehicle depending on soc
            return model.p_charging[v, t] <= model.p_max_veh[v] * model.ev_connected[v, t]
        model.charging_power_ub2 = Constraint(model.veh, model.T, rule=upper_power_bounds_veh)


        def power_balance(model, v, t):
            return model.p_opf[v, t] == model.p_charging[v, t] + model.p_discharging[v, t]  # remove p_opf eventually
        model.power_balance_constraint = Constraint(model.veh, model.T, rule=power_balance)

        # initial SOC defined as reference SOC at the beginning of the optimization horizon
        def initial_soc_rule(model, v):
            return model.soc[v, model.T.first()] == model.soc_init[v, model.T.first()]
        model.initial_soc_constraint = Constraint(model.veh, rule=initial_soc_rule)

        # final SOC must be equal to initial SOC to have comparable net energy charging
        def final_soc_rule(model, v):
            return model.soc[v, model.T.last()] == model.soc_init[v, model.T.last()]
        model.final_soc_constraint = Constraint(model.veh, rule=final_soc_rule)


        # SOC after discharging the battery must be above preferred by owner
        def soc_limits(model, v, t):
            if t == model.T.first():
                is_driving = self.p[v, t-1] < 0
            else:
                is_driving = model.p[v, t-1] < 0
            # vehicle is not connected at t-1 or t
            ev_not_connected = self.veh_cps[v, t-1] == -1 or self.veh_cps[v, t] == -1

            # Find the next departure and last arrival time for the vehicle
            dep_times = model.t_dep[v]
            next_deps = dep_times[dep_times >= t]
            next_dep_time = next_deps.min() if next_deps.size > 0 else None
            arr_times = model.t_arr[v]
            past_arrs = arr_times[arr_times <= t]
            last_arr_time = past_arrs.max() if past_arrs.size > 0 else 0

            recovery_period = 6
            recovery_time_end = last_arr_time + recovery_period if last_arr_time is not None else t

            # if the last 4 timesteps prior to arrival time the vehicle has been driving create a parameter 'long Trip' and set it to 1
            long_trip = 1 if all(self.p[v, last_arr_time - i] < 0 for i in range(1, 5)) else 0

            # Skipping step if discharging is due to driving -> then EV is allowed to go below these limits
            if ev_not_connected or t == last_arr_time:
                return Constraint.Skip

            # ensure feasibility of model by temporarily allowing soc to be below min soc for parking
            elif last_arr_time is not None and t < recovery_time_end:
                min_soc_dynamic = model.min_soc_parking * (t - last_arr_time) / (recovery_period * 2)
                return model.soc[v, t] >= min_soc_dynamic * long_trip + model.min_soc_parking * (1 - long_trip)

            # if vehicle is departing in the next timestep,the soc must be above the min soc for departing
            elif t == next_dep_time:
                return model.soc[v, t] >= model.min_soc_departing - model.buffer_soc[v, t]

            # if vehicle is still parking in the next time step, the soc must be above the min soc for parking
            else:
                return model.soc[v, t] >= model.min_soc_parking - model.buffer_soc[v, t]
        model.soc_limits = Constraint(model.veh, model.T, rule=soc_limits)


        # SOC update
        def update_soc(model, v, t):
            if t == model.T.first():
                # Initial state of charge (e.g., based on initial scenario conditions)
                return Constraint.Skip
            else:
                # soc is updated based on the ch./dch. power if ev is connected or driving power if ev is on the road
                return (model.soc[v, t] == model.soc[v, t-1] + (model.p_charging[v, t-1] * model.eta
                                                                + model.p_discharging[v, t-1] * (1/model.eta))
                        * model.param[v] + model.p_drive[v, t-1] * model.param[v])

        model.update_soc = Constraint(model.veh, model.T, rule=update_soc)

        # take soc after last timestep into account
        def limit_last_power(model, v):
            return (0, model.soc[v, model.T.last()] + (model.p_charging[v, model.T.last()] * model.eta
                                                       + model.p_discharging[v, model.T.last()] * (1/model.eta)) * model.param[v]
                    + model.p_drive[v, model.T.last()] * model.param[v], 1)
        model.limit_last_power = Constraint(model.veh, rule=limit_last_power)


        def non_controllable_power(model, v, t):
            return model.p_charging[v, t] == model.p_req_ev[v, t]
        model.non_controllable_power_constraint = Constraint(model.veh_nc, model.T, rule=non_controllable_power)

    def get_market_constraints(self, model):
        # --- market constraints ---
        # aggregated power has to be constant over 4 time steps (1 hour = 4 x 15 min)
        def constant_aggregated_power_per_hour(model, t):
            if (t - model.T.first()) % 4 != 3:
                return sum(model.p_opf[v, t] for v in model.veh) - sum(model.p_opf[v, t + 1] for v in model.veh) == 0
            else:
                return Constraint.Skip
        model.constant_aggregated_power_per_hour = Constraint(model.T, rule=constant_aggregated_power_per_hour)

    def get_opf_objective(self, model):
        pass

    def get_all(self, model):
        self.get_opf_sets(model)
        self.get_opf_parameters(model)
        self.get_opf_variables(model)

    def get_all_opf(self, model):
        self.get_opf_constraints(model)
        self.get_opf_objective(model)
        self.set_simbench_ev_to_zero(model)
        # self.active_soc_constraints(model)

    def get_all_acopf(self, model):
        pass

    def fix_variables(self, model):
        # fix the ev charging values as in cpg without optimizing , no discharging possibilty
        for v in model.veh:
            for t in model.T:
                model.p_opf[v, t].fix(model.p_req_ev[v, t])
        return True

    def set_simbench_ev_to_zero(self, model):
        hls_loads = self.net.load[self.net.load.profile.str.startswith("HLS")].index
        apls_loads = self.net.load[self.net.load.profile.str.startswith("APLS")].index
        for d in model.D:
            for t in model.T:
                if d in hls_loads:
                    model.pD[d, t].fix(0)
                elif d in apls_loads:
                    model.pD[d, t].fix(0)
        return True
