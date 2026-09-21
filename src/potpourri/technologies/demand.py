# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Demand mix-in: load profiles and demand bounds over time."""

import pyomo.environ as pyo
from potpourri.technologies.flexibility import Flexibility_multi_period


class Demand_multi_period(Flexibility_multi_period):
    """Multi-period demand device, driven by load profiles.

    Multi-period demand device module; reads load profiles from
    net.profiles.
    """

    def __init__(self, net, T=None, scenario=None):
        super().__init__(net, T, scenario)

        self.demand_set = self.net.load.index[self.net.load.in_service]
        self.bus_demand_set = list(
            zip(
                self.bus_lookup[self.net.load.bus[self.demand_set].values],
                self.demand_set,
            )
        )

    def get_all(self, model):
        """Attach the demand and fix it to its profile.

        Attach demand sets, parameters, variables and fix them to profile
        values.
        """
        self.get_sets(model)
        self.get_parameters(model)
        self.get_variables(model)
        self.fix_variables(model)

    def get_all_opf(self, model):
        """Attach the demand OPF bounds.

        Attach OPF demand sets, parameters, and real-power bound
        constraints.
        """
        self.get_opf_sets(model)
        self.get_opf_parameters(model)
        self.get_all_constraints_opf(model)

    def get_all_acopf(self, model):
        """Attach the AC OPF layer for this device.

        Adds the reactive-power bounds on top of the active-power ones.

        Returns:
                None. The components are added to `model` in place.
        """
        self.get_acopf_parameters(model)
        self.get_all_constraints_acopf(model)

    def get_all_ac(self, model):
        """Attach the AC layer for this device.

        Adds the reactive-power parameters and variables.

        Returns:
                None. The components are added to `model` in place.
        """
        self.get_ac_parameters(model)
        self.get_ac_variables(model)

    def get_sets(self, model):
        """Add this device's index sets to `model`.

        Extends the base sets with the device's own element set and its
        element-to-bus mapping.

        Returns:
                True, so a caller can chain the lifecycle steps. The real
                result is the components added to `model`.
        """
        super().get_sets(model)
        model.D = pyo.Set(initialize=self.demand_set)  # set of demands
        # loads linked to each bus b
        model.Dbs = pyo.Set(
            within=model.B * model.D, initialize=self.bus_demand_set
        )  # set of demand-bus mapping
        return True

    def get_parameters(self, model):
        """Add the device's active-power parameters to `model`.

        Values are spread over the time index by `make_to_dict`, so each
        parameter is keyed by `(element, time)` and carried in per unit.

        Returns:
                True, so a caller can chain the lifecycle steps. The real
                result is the components added to `model`.
        """
        # make PD_data_dict index-able through the tuple list
        self.PD_data_dict, self.PD_tuple = self.make_to_dict(
            model.D, model.T, self.PD_data
        )
        # demand at each bus
        model.PD = pyo.Param(self.PD_tuple, initialize=self.PD_data_dict)
        return True

    def get_ac_parameters(self, model):
        """Add the device's reactive-power parameters to `model`.

        Only meaningful on an AC model; a DC model has no reactive power to
        bound.

        Returns:
                True, so a caller can chain the lifecycle steps. The real
                result is the components added to `model`.
        """
        # reactive demand
        self.QD_data_dict, self.QD_tuple = self.make_to_dict(
            model.D, model.T, self.QD_data
        )
        model.QD = pyo.Param(self.QD_tuple, initialize=self.QD_data_dict)
        return True

    def get_ac_variables(self, model):
        """Add the device's reactive-power variables to `model`.

        Indexed by `(element, time)`, in per unit. AC models only.

        Returns:
                True, so a caller can chain the lifecycle steps. The real
                result is the components added to `model`.
        """
        # reactive demand
        model.qD = pyo.Var(self.QD_tuple, domain=pyo.Reals)
        return True

    def get_variables(self, model):
        """Add the device's active-power variables to `model`.

        Indexed by `(element, time)`, in per unit.

        Returns:
                True, so a caller can chain the lifecycle steps. The real
                result is the components added to `model`.
        """
        # --- Variables ---
        # demand at each bus
        model.pD = pyo.Var(
            self.PD_tuple, domain=pyo.Reals
        )  # real power demand delivered
        return True

    def fix_variables(self, model):
        """Pin the device's power to its profile.

        Fixes each variable at the corresponding parameter value, so the device
        behaves as a fixed injection until an OPF layer frees it again.

        Returns:
                None. The components are added to `model` in place.
        """
        # fix the demand values over all d in D and t in T
        for d_t in self.PD_tuple:
            model.pD[d_t].fix(model.PD[d_t])
        # for d in model.D for t in model.T:
        #     model.pD[(d,t)].fix(model.PD[(d,t)])

    def get_opf_sets(self, model):
        """Add the controllable subset of this device to `model`.

        Only the elements flagged controllable in the pandapower table get
        operating bounds; the rest stay at their profile.

        Returns:
                None. The components are added to `model` in place.
        """
        model.Dc = pyo.Set(
            within=model.D, initialize=self.demand_controllable_set
        )  # controllable loads

    def get_opf_parameters(self, model):
        """Add the active-power operating bounds to `model`.

        Keyed by `(element, time)` and in per unit.

        Returns:
                None. The components are added to `model` in place.
        """
        # real demand
        model.PDmax = pyo.Param(
            self.PDmax_tuple, initialize=self.PDmax_data_dict
        )
        model.PDmin = pyo.Param(
            self.PDmin_tuple, initialize=self.PDmin_data_dict
        )

    def get_acopf_parameters(self, model):
        """Add the reactive-power operating bounds to `model`.

        Keyed by `(element, time)` and in per unit. AC models only.

        Returns:
                None. The components are added to `model` in place.
        """
        # reactive demand
        model.QDmax = pyo.Param(
            self.QDmax_tuple, initialize=self.QDmax_data_dict
        )
        model.QDmin = pyo.Param(
            self.QDmin_tuple, initialize=self.QDmin_data_dict
        )

    def get_all_constraints_opf(self, model):
        """Add the active-power bound constraints.

        Each rule **unfixes** its variable before returning the bound, so the
        profile value becomes a starting point rather than a fixed value, and
        returns Pyomo's `(lower, expr, upper)` ranged form.

        Returns:
                None. The components are added to `model` in place.
        """

        @model.Constraint(model.Dc, model.T)
        def real_demand_bounds(model, d, t):
            """Bound a controllable load's active power, freeing it first.

            Args:
                model: The Pyomo model being built.
                d: Load index.
                t: Time index.

            Returns:
                A Pyomo expression.
            """
            model.pD[(d, t)].unfix()
            return model.PDmin[(d, t)], model.pD[(d, t)], model.PDmax[(d, t)]

    def get_all_constraints_acopf(self, model):
        """Add the reactive-power bound constraints.

        As `get_all_constraints_opf`, for reactive power.

        Returns:
                None. The components are added to `model` in place.
        """

        @model.Constraint(model.Dc, model.T)
        def reactive_demand_bounds(model, d, t):
            """Bound a controllable load's reactive power, freeing it first.

            Args:
                model: The Pyomo model being built.
                d: Load index.
                t: Time index.

            Returns:
                A Pyomo expression.
            """
            model.qD[(d, t)].unfix()
            return model.QDmin[(d, t)], model.qD[(d, t)], model.QDmax[(d, t)]

    def get_demand_real_power_data(self, model, max_p_mw=None, min_p_mw=None):
        """Work out the active-power bounds for controllable loads.

        Loads flagged `controllable` in `net.load` become adjustable between a
        minimum and a maximum; everything else stays at its profile. Bounds
        come from the explicit arguments when given, otherwise from the
        `max_p_mw` / `min_p_mw` columns.

        Sets `self.demand_controllable_set` and the `PDmax`/`PDmin`
        dictionaries used by `get_opf_parameters`. When the network has no
        `controllable` column the set is empty and no load is adjustable.

        Args:
            model: The multi-period model being extended.
            max_p_mw: Upper active-power bound in MW, applied to every
                controllable load. `None` reads the network column.
            min_p_mw: Lower active-power bound in MW. `None` reads the network
                column.

        Returns:
            None. Stores its results on `self`.
        """
        if "controllable" not in self.net.load:
            self.demand_controllable_set = (
                None  # create empty Set if no controllable load exist
            )
        else:
            self.demand_controllable_set = self.net.load.index[
                self.net.load.controllable.astype(bool)
            ]

        # Aus dem Netzwerk die maximalen und minimalen Lasten holen
        if max_p_mw is not None:
            self.PDmax_data_dict, self.PDmax_tuple = self.make_to_dict(
                model.D, model.T, max_p_mw
            )
        elif max_p_mw is None:
            self.PDmax_data_dict, self.PDmax_tuple = self.make_to_dict(
                model.D, model.T, self.PD_data
            )

        # add rows with active demand limits if not existing
        if min_p_mw is not None:
            self.PDmin_data_dict, self.PDmin_tuple = self.make_to_dict(
                model.D, model.T, min_p_mw, False
            )
        elif min_p_mw is None:
            self.PDmin_data_dict, self.PDmin_tuple = self.make_to_dict(
                model.D, model.T, 0, False
            )

    def get_demand_reactive_data(
        self, model, max_q_mvar=None, min_q_mvar=None
    ):
        """Work out the reactive-power bounds for controllable loads.

        The reactive counterpart of `get_demand_real_power_data`, used only by
        AC models.

        Args:
            model: The multi-period model being extended.
            max_q_mvar: Upper reactive bound in MVAr, or `None` to read the
                network column.
            min_q_mvar: Lower reactive bound in MVAr, or `None`.

        Returns:
            None. Stores its results on `self`.
        """
        # reactive power demand
        if max_q_mvar is not None:
            self.QDmax_data_dict, self.QDmax_tuple = self.make_to_dict(
                model.D, model.T, max_q_mvar
            )
        elif max_q_mvar is None:
            self.QDmax_data_dict, self.QDmax_tuple = self.make_to_dict(
                model.D, model.T, self.QD_data.abs()
            )

        if min_q_mvar is not None:
            self.QDmin_data_dict, self.QDmin_tuple = self.make_to_dict(
                model.D, model.T, min_q_mvar
            )
        elif min_q_mvar is None:
            self.QDmin_data_dict, self.QDmin_tuple = self.make_to_dict(
                model.D, model.T, -(self.QD_data.abs())
            )
        # demand limits for loads
