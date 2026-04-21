"""Base mix-in class for all multi-period flexibility device modules
(batteries, heat pumps, PV, demand, etc.)."""

import pandas as pd
import pyomo.environ as pyo
from math import pi
import numpy as np


class Flexibility_multi_period:
    """Base class for technology mix-in objects that attach Pyomo components
    to a multi-period model.

    Reads network topology and profile data from *net* in ``__init__``.
    Subclasses implement ``get_all()``, ``get_sets()``, ``get_parameters()``,
    ``get_variables()``, and constraint methods.
    """

    def __init__(self, net, T=None, scenario=None):
        self.net = net

        # buses that are not ext_grids (eligible for technology placement)
        self.buses_excl_extGrids = self.net.bus.loc[
            ~self.net.bus.index.isin(net.ext_grid.bus)
        ].index

        bus_set = self.net._ppc["bus"][:, [0, 1, 7, 8]]
        bus_set[:, -1] *= pi / 180
        self.bus_data = pd.DataFrame(
            bus_set[:, 1:],
            index=bus_set[:, 0].astype(int),
            columns=["type", "v_m", "v_a_rad"],
        )
        self.bus_lookup = self.net._pd2ppc_lookups["bus"]

        self.baseMVA = self.net.sn_mva

        self.PD_data = self.net.profiles[("load", "p_mw")] / self.baseMVA
        self.QD_data = self.net.profiles[("load", "q_mvar")] / self.baseMVA

        self.GB_data = self.net.shunt.p_mw * self.net.shunt.step / self.baseMVA

    def get_sets(self, model):
        """Initialise (or re-initialise) the bus set B on the Pyomo model
        from network topology."""
        if hasattr(model, "B"):
            model.del_component(model.B)
        model.B = pyo.Set(initialize=self.bus_data.index)
        return True

    def make_to_dict(self, model_obj, model_time, data, time_dependent=True):
        """Convert data into a ``{(object_index, time_index): value}`` dict
        and a matching list of index tuples for Pyomo Param initialisation.

        Args:
            model_obj: Iterable of object indices (e.g. a Pyomo Set).
            model_time: Iterable of time indices (e.g. ``model.T``).
            data: Data values.  Accepted types:

                * ``np.ndarray`` — indexed by object index.
                * ``pd.Series`` — indexed by time (``time_dependent=True``)
                  or by object index (``time_dependent=False``).
                * ``dict`` — ``{obj: value}`` (time-independent) or
                  ``{obj: {t: value}}`` (time-dependent).
                * ``list`` — same layout as ndarray.
                * ``float`` / ``int`` — scalar, broadcast to all (obj, t).
                * ``pd.DataFrame`` — calls ``.to_dict()`` internally.

            time_dependent: When ``True`` (default) the same value is
                repeated for every time step for each object.  When
                ``False`` the data varies only over objects (not time).

        Returns:
            tuple: ``(data_dict, tuple_list)`` where *data_dict* maps
            ``(obj, t)`` → value and *tuple_list* is the ordered list of
            index pairs.
        """
        tuple_list = [(o, t) for o in model_obj for t in model_time]

        if time_dependent:
            if isinstance(data, np.ndarray):
                data_dict = {(o, t): data[o] for o, t in tuple_list}
            elif isinstance(data, dict):
                # Support both flat {obj: val} and nested {obj: {t: val}}
                first_val = next(iter(data.values())) if data else None
                if isinstance(first_val, dict):
                    data_dict = {(o, t): data[o][t] for o, t in tuple_list}
                else:
                    data_dict = {(o, t): data[o] for o, t in tuple_list}
            elif isinstance(data, pd.Series):
                data_dict = {(o, t): data[t] for o, t in tuple_list}
            elif isinstance(data, list):
                data_dict = {(o, t): data[t] for o, t in tuple_list}
            elif isinstance(data, (int, float)):
                data_dict = {(o, t): data for o, t in tuple_list}
            else:
                raw = data.to_dict()
                data_dict = {(o, t): raw[o][t] for o, t in tuple_list}
        else:
            if isinstance(data, (int, float)) and data == 0:
                data_dict = {(o, t): 0 for o, t in tuple_list}
            elif isinstance(data, np.ndarray):
                data_dict = {(o, t): data[o] for o, t in tuple_list}
            elif isinstance(data, pd.Series):
                data_dict = {(o, t): data[o] for o, t in tuple_list}
            elif isinstance(data, dict):
                data_dict = {(o, t): data[o] for o, t in tuple_list}
            elif isinstance(data, list):
                data_dict = {(o, t): data[o] for o, t in tuple_list}
            elif isinstance(data, (int, float)):
                data_dict = {(o, t): data for o, t in tuple_list}
            else:
                raw = data.to_dict()
                data_dict = {(o, t): raw[o] for o, t in tuple_list}

        return data_dict, tuple_list
