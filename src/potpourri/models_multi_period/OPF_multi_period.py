# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Multi-period OPF mix-in: operating limits over a horizon.

Adds generator, load and branch limits, and thermal ratings,
indexed over the model's time set.
"""

import numpy as np
import pyomo.environ as pyo

from potpourri.models.oltc import OLTCControlMixin
from potpourri.models.shunt_control import ShuntControlMixin
from potpourri.models_multi_period.basemodel_multi_period import (
    Basemodel_multi_period,
)
from potpourri.technologies.demand import Demand_multi_period
from potpourri.technologies.generator import Generator_multi_period
from potpourri.technologies.sgens import Sgens_multi_period


class OPF_multi_period(
    OLTCControlMixin, ShuntControlMixin, Basemodel_multi_period
):
    """OPF mixin for multi-period models.

    Provides line/transformer ratings and generator/demand limits.
    """

    def __init__(self, net, toT, fromT=None, pf=1):
        super().__init__(net, toT, fromT, pf)

    def __calc_SLmax(self, max_loading_percent=100):
        """Apparent-power rating of every branch, in p.u.

        A line's rating is a current limit in pandapower, converted with the
        from-bus nominal voltage; a transformer already carries an MVA rating.

        Returns:
            A Pyomo expression.
        """
        # Native lines: max_i_ka @ from-bus vn_kv → MVA limit, p.u.
        vr = self.net.bus.loc[
            self.net.line["from_bus"].values, "vn_kv"
        ].values * np.sqrt(3.0)
        max_i_ka = self.net.line.max_i_ka.values
        df = self.net.line.df.values
        line_lim = (
            max_loading_percent
            / 100.0
            * max_i_ka
            * df
            * self.net.line.parallel.values
            * vr
            / self.baseMVA
        )

        # Impedance branches: per-unit S limit from net.impedance.sn_mva
        # (0 → unrated). They appear in model.L after the native lines (see
        # Basemodel_multi_period.__init__). max_loading_percent may be an
        # array over net.line only — apply 100% by default for impedance
        # entries (pandapower has no max_loading column on net.impedance).
        imp = self.net.get("impedance")
        if imp is not None and not imp.empty:
            sn_imp = imp["sn_mva"].astype(float).values.copy()
            unrated = sn_imp <= 0
            sn_imp[unrated] = 1e6
            imp_lim = 1.0 * sn_imp / self.baseMVA
            return np.concatenate([line_lim, imp_lim])
        return line_lim

    def _calc_opf_parameters(self, **kwargs):
        """Compute branch ratings and the device operating limits.

        Compute line/transformer ratings and call generator/demand limit
        methods on flexibility objects.

        Args:
            **kwargs: No options are consumed here. Anything left over has
                been forwarded from ``add_OPF`` without a consumer, so it is
                rejected rather than ignored — silently swallowing, e.g.,
                ``thermal_limit`` on the DC path made an option look
                supported when it changed nothing.

        Raises:
            TypeError: If any keyword argument reaches this point.
        """
        if kwargs:
            unsupported = ", ".join(sorted(kwargs))
            raise TypeError(
                f"{type(self).__name__}.add_OPF() got unsupported option(s): "
                f"{unsupported}. This model kind does not implement them; "
                f"see the single-period vs multi-period table in "
                f"docs/user-guide/multi-period.md."
            )
        max_load = (
            self.net.line.max_loading_percent.values
            if "max_loading_percent" in self.net.line
            else 100.0
        )
        self.line_data["SLmax_data"] = self.__calc_SLmax(max_load)

        # maximum transformer loading
        max_load_T = (
            self.net.trafo.max_loading_percent.fillna(100.0) / 100.0
            if "max_loading_percent" in self.net.trafo
            else 1.0
        )
        sn_mva = self.net.trafo.sn_mva
        df_T = self.net.trafo.df
        SLmaxT_data = (
            max_load_T * sn_mva * df_T * self.net.trafo.parallel / self.baseMVA
        )
        self.trafo_data["SLmaxT_data"] = SLmaxT_data.values

        # create generator instance and call method
        # generation_real_power_limits_opf
        generator_object = next(
            (
                obj
                for obj in self.flexibilities
                if isinstance(obj, Generator_multi_period)
            ),
            None,
        )
        generator_object.generation_real_power_limits_opf(self.model)

        # create instance of class 'Sgens' from the 'flexibilities' list and
        # call method 'static_generation_real_power_limits'
        sgens_object = next(
            (
                obj
                for obj in self.flexibilities
                if isinstance(obj, Sgens_multi_period)
            ),
            None,
        )
        sgens_object.static_generation_real_power_limits(self.model)

        # Get the object of class 'Demand' from the 'flexibilities' list and
        # call method 'get_demand_real_power_data'
        demand_object = next(
            (
                obj
                for obj in self.flexibilities
                if isinstance(obj, Demand_multi_period)
            ),
            None,
        )
        # Change: gives model instead of net, because model is needed and net
        # is already in flexibilities
        demand_object.get_demand_real_power_data(self.model)

    def add_OPF(self, **kwargs):
        """Attach OPF parameters and constraints.

        Ratings, generator limits, demand limits.

        Args:
            **kwargs: Forwarded to _calc_opf_parameters.
        """
        self._calc_opf_parameters(**kwargs)

        # get all opf parameters from flexibility objects
        for flex in (
            self.flexibilities
        ):  # Gets all Opf parmas ands sets from flexibility objects
            flex.get_all_opf(self.model)

        # Devices attach to the model after the power-flow equations were
        # built, so the balance is rebuilt here to pick up their injections.
        # A no-op when nothing registered a coupling term.
        self.rebuild_kcl()

        # lines and transformer chracteristics and ratings
        self.model.SLmax = pyo.Param(
            self.model.L,
            within=pyo.NonNegativeReals,
            initialize=self.line_data["SLmax_data"][self.model.L],
            mutable=True,
        )  # real power line limit
        self.model.SLmaxT = pyo.Param(
            self.model.TRANSF,
            within=pyo.NonNegativeReals,
            initialize=self.trafo_data.SLmaxT_data[self.model.TRANSF],
            mutable=True,
        )  # real power transformer limit

    def _branch_angle_bounds(self, table, idx_set, hv_col, lv_col):
        """Read finite per-branch angle bounds, in radians, keyed by index.

        Mirrors the reader in ``ACOPF._add_branch_angle_limits``: branches
        with no angle columns, non-finite bounds, or the ±360° placeholder
        MATPOWER uses for "unconstrained" are left out, as are the synthetic
        impedance indices in ``model.L`` that have no row in ``net.line``.
        """
        angmin_col, angmax_col = "angmin_degree", "angmax_degree"
        if angmin_col not in table.columns or angmax_col not in table.columns:
            return {}
        valid = set(table.index)
        out = {}
        for ix in idx_set:
            if ix not in valid:
                continue
            amin = float(table.at[ix, angmin_col])
            amax = float(table.at[ix, angmax_col])
            if (
                not np.isfinite(amin)
                or not np.isfinite(amax)
                or abs(amin) >= 359.0
                or abs(amax) >= 359.0
            ):
                continue
            out[ix] = (
                self.bus_lookup[int(table.at[ix, hv_col])],
                self.bus_lookup[int(table.at[ix, lv_col])],
                np.deg2rad(amin),
                np.deg2rad(amax),
            )
        return out

    def _add_branch_angle_limits(self):
        """Attach time-indexed branch phase-angle-difference constraints.

        PowerModels.jl convention: ``angmin ≤ δ_from − δ_to ≤ angmax``, held
        at every time step. The single-period equivalent is
        ``ACOPF._add_branch_angle_limits``.
        """
        line_bounds = self._branch_angle_bounds(
            self.net.line, list(self.model.L), "from_bus", "to_bus"
        )
        trafo_bounds = self._branch_angle_bounds(
            self.net.trafo, list(self.model.TRANSF), "hv_bus", "lv_bus"
        )

        if line_bounds:
            self.model.LineAngleSet = pyo.Set(initialize=list(line_bounds))

            @self.model.Constraint(self.model.LineAngleSet, self.model.T)
            def line_angle_diff(model, l, t):
                r"""Phase-angle-difference limit on line `l`.

                $\alpha_{min} \le \theta_f - \theta_t \le \alpha_{max}$, in
                radians.

                Args:
                    model: The Pyomo model being built.
                    l: Branch index.
                    t: Time index.

                Returns:
                    A Pyomo expression.
                """
                f, to, amin, amax = line_bounds[l]
                return amin, model.delta[f, t] - model.delta[to, t], amax

        if trafo_bounds:
            self.model.TrafoAngleSet = pyo.Set(initialize=list(trafo_bounds))

            @self.model.Constraint(self.model.TrafoAngleSet, self.model.T)
            def trafo_angle_diff(model, l, t):
                """Phase-angle-difference limit on transformer `l`.

                Args:
                    model: The Pyomo model being built.
                    l: Branch index.
                    t: Time index.

                Returns:
                    A Pyomo expression.
                """
                f, to, amin, amax = trafo_bounds[l]
                return amin, model.delta[f, t] - model.delta[to, t], amax
