"""Base mix-in class for all multi-period flexibility device modules
(batteries, EVs, demand, etc.)."""

import pandas as pd
import pyomo.environ as pyo
from math import pi
import numpy as np


class Flexibility_multi_period:
    """Base class for technology mix-in objects that attach Pyomo components
    to a multi-period model.

    Reads network topology and profile data from net in __init__. Subclasses
    implement get_all(), get_sets(), get_parameters(), get_variables(), and
    constraint methods.
    """

    # TODO see if correct time dependency implementation
    def __init__(self, net, T=None, scenario=None):
        self.net = net
        self.t_range = (
            range(T) if T else [0]
        )  # If T is None, default to a single time period [0]

        # buses that are not ext_grids
        self.buses_excl_extGrids = self.net.bus.loc[
            ~self.net.bus.index.isin(net.ext_grid.bus)
        ].index

        # --- Sets ---
        bus_set = self.net._ppc["bus"][:, [0, 1, 7, 8]]
        bus_set[:, -1] *= pi / 180
        self.bus_data = pd.DataFrame(
            bus_set[:, 1:],
            index=bus_set[:, 0].astype(int),
            columns=["type", "v_m", "v_a_rad"],
        )
        self.bus_lookup = self.net._pd2ppc_lookups["bus"]

        # --- Param Data ---
        self.baseMVA = self.net.sn_mva

        self.PD_data = self.net.profiles[("load", "p_mw")] / self.baseMVA
        self.QD_data = self.net.profiles[("load", "q_mvar")] / self.baseMVA
        # self.QD_data = self.net.load.q_mvar / self.baseMVA

        self.GB_data = self.net.shunt.p_mw * self.net.shunt.step / self.baseMVA

    def get_sets(self, model):
        """Initialise (or re-initialise) the bus set B on the Pyomo model
        from network topology."""
        if hasattr(model, "B"):
            model.del_component(model.B)
        # Make B time-dependent # why?
        model.B = pyo.Set(initialize=self.bus_data.index)
        return True

    def make_to_dict(self, model_obj, model_time, data, time_dependent=True):
        """
        **make_to_dict** \n
        Args:
            model_obj: Object Indices
            model_time: Time Indices
            data: Data to be converted to dictionary
            time_dependent: should the data be time dependent True or
                constant for all time steps False

        Returns:
            data_dict: Dictionary with the data
            tuple_list: List of tuples with the object and time index

        """
        # TODO Check if it is dict already and then just take dict and pack
        # indexes in tuple_list
        if time_dependent:
            if isinstance(data, np.ndarray):
                data_dict = {
                    (o, t): data[o] for o in model_obj for t in model_time
                }
                tuple_list = [(o, t) for o in model_obj for t in model_time]
                return data_dict, tuple_list

            if isinstance(data, dict):
                data_dict = {
                    (o, t): data[o][t] for o in model_obj for t in model_time
                }
                tuple_list = list(
                    [(o, t) for o in model_obj for t in model_time]
                )
                return data_dict, tuple_list

            if isinstance(data, pd.Series):
                data_dict = {
                    (o, t): data[t] for o in model_obj for t in model_time
                }
                tuple_list = list(
                    [(o, t) for o in model_obj for t in model_time]
                )
                return data_dict, tuple_list
            if isinstance(data, list):
                data_dict = {
                    (o, t): data[t] for o in model_obj for t in model_time
                }
                tuple_list = list(
                    [(o, t) for o in model_obj for t in model_time]
                )
                return data_dict, tuple_list

            else:
                data_dict = data.to_dict()
                data_dict = {
                    (o, t): data_dict[o][t]
                    for o in model_obj
                    for t in model_time
                }
                tuple_list = list(
                    [(o, t) for o in model_obj for t in model_time]
                )
                return data_dict, tuple_list

        else:
            if data == 0:
                data_dict = {(o, t): 0 for o in model_obj for t in model_time}
                tuple_list = list(
                    [(o, t) for o in model_obj for t in model_time]
                )
                return data_dict, tuple_list

            if isinstance(data, pd.Series):
                data_dict = {
                    (o, t): data[o] for o in model_obj for t in model_time
                }
                tuple_list = list(
                    [(o, t) for o in model_obj for t in model_time]
                )
                return data_dict, tuple_list

            if isinstance(data, np.ndarray):
                data_dict = {
                    (o, t): data[o] for o in model_obj for t in model_time
                }
                tuple_list = [(o, t) for o in model_obj for t in model_time]
                return data_dict, tuple_list

            # when data is float put it in dict, correct?
            if isinstance(data, float):
                data_dict = {
                    (o, t): data for o in model_obj for t in model_time
                }
                tuple_list = list(
                    [(o, t) for o in model_obj for t in model_time]
                )
                return data_dict, tuple_list

            # if data is already a dict just put it in data_dict and make
            # tuple_list with time
            if isinstance(data, dict):
                data_dict = {
                    (o, t): data[o] for o in model_obj for t in model_time
                }
                tuple_list = list(
                    [(o, t) for o in model_obj for t in model_time]
                )
                return data_dict, tuple_list

            if isinstance(data, list):
                data_dict = {
                    (o, t): data[o] for o in model_obj for t in model_time
                }
                tuple_list = list(
                    [(o, t) for o in model_obj for t in model_time]
                )
                return data_dict, tuple_list
            else:
                data_dict = data.to_dict()
                data_dict = {
                    (o, t): data_dict[o] for o in model_obj for t in model_time
                }
                tuple_list = list(
                    [(o, t) for o in model_obj for t in model_time]
                )
                return data_dict, tuple_list
