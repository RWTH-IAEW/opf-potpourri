# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Multi-period OPF mixin: adds operational limit constraints and thermal
ratings over time."""

import copy

import numpy as np
from pyomo.environ import *

from potpourri.models_multi_period.basemodel_multi_period import (
    Basemodel_multi_period,
)
from potpourri.technologies.demand import Demand_multi_period
from potpourri.technologies.generator import Generator_multi_period
from potpourri.technologies.sgens import Sgens_multi_period


def _calc_tap_min_max_mp(obj):
    """Compute off-nominal tap-ratio bounds from net.trafo tap_min/max."""
    vnh_min, vnl_min, _ = _calc_tap_shift_mp(
        obj, tap_pos=obj.net.trafo.tap_min
    )
    vnh_max, vnl_max, _ = _calc_tap_shift_mp(
        obj, tap_pos=obj.net.trafo.tap_max
    )
    ratio_min = _calc_nominal_ratio_mp(obj, vnh_min, vnl_min)
    ratio_max = _calc_nominal_ratio_mp(obj, vnh_max, vnl_max)
    return ratio_min, ratio_max


def _calc_nominal_ratio_mp(obj, vn_hv_kv, vn_lv_kv):
    tap_rat = vn_hv_kv / vn_lv_kv
    hv_bus = obj.net.trafo.hv_bus
    lv_bus = obj.net.trafo.lv_bus
    nom_rat = (
        obj.net.bus.vn_kv[hv_bus].values / obj.net.bus.vn_kv[lv_bus].values
    )
    return tap_rat / nom_rat


def _calc_tap_shift_mp(obj, tap_pos=None):
    """Compute adjusted vn_hv_kv, vn_lv_kv, and phase shift for tap_pos."""
    vnh = copy.deepcopy(obj.net.trafo.vn_hv_kv.values)
    vnl = copy.deepcopy(obj.net.trafo.vn_lv_kv.values)
    trafo_shift = obj.net.trafo.shift_degree.values.copy()

    if tap_pos is None:
        tap_pos = obj.net.trafo.tap_pos
    tap_neutral = obj.net.trafo.tap_neutral
    tap_diff = tap_pos - tap_neutral
    tap_phase_shifter = obj.net.trafo.get(
        "tap_phase_shifter",
        __import__("pandas").Series(False, index=obj.net.trafo.index),
    )
    tap_side = obj.net.trafo.tap_side
    tap_step_percent = obj.net.trafo.tap_step_percent
    tap_step_degree = obj.net.trafo.tap_step_degree

    for side, vn, direction in [("hv", vnh, 1), ("lv", vnl, -1)]:
        phase_shifters = tap_phase_shifter & (tap_side == side)
        tap_complex = (
            np.isfinite(tap_step_percent)
            & np.isfinite(tap_pos)
            & (tap_side == side)
            & ~phase_shifters
        )
        if tap_complex.any():
            tap_steps = (
                tap_step_percent[tap_complex] * tap_diff[tap_complex] / 100
            )
            tap_angles = tap_step_degree[tap_complex].fillna(0)
            u1 = vn[tap_complex]
            du = u1 * tap_steps.fillna(0)
            cos_a = np.cos(np.deg2rad(tap_angles))
            sin_a = np.sin(np.deg2rad(tap_angles))
            vn[tap_complex] = np.sqrt(
                (u1 + du * cos_a) ** 2 + (du * sin_a) ** 2
            )
            trafo_shift[tap_complex] += np.rad2deg(
                np.arctan(direction * du * sin_a / (u1 + du * cos_a))
            )
        if phase_shifters.any():
            deg_set = tap_step_degree[phase_shifters].fillna(0) != 0
            pct_set = tap_step_percent[phase_shifters].fillna(0) != 0
            if (deg_set & pct_set).any():
                raise ValueError(
                    "Both tap_step_degree and tap_step_percent set for "
                    "ideal phase shifter."
                )
            trafo_shift[phase_shifters] += np.where(
                deg_set,
                direction
                * tap_diff[phase_shifters]
                * tap_step_degree[phase_shifters],
                direction
                * 2
                * np.rad2deg(
                    np.arcsin(
                        tap_diff[phase_shifters]
                        * tap_step_percent[phase_shifters]
                        / 100
                        / 2
                    )
                ),
            )
    return vnh, vnl, trafo_shift


class OPF_multi_period(Basemodel_multi_period):
    """OPF mixin for multi-period models: provides line/transformer ratings
    and generator/demand limits."""

    def __init__(self, net, toT, fromT=None, pf=1):
        super().__init__(net, toT, fromT, pf)

    def __calc_SLmax(self, max_loading_percent=100):
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
        """Compute line/transformer ratings and call generator/demand limit
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
        """Attach OPF parameters and constraints: ratings, generator limits,
        demand limits.

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
        self.model.SLmax = Param(
            self.model.L,
            within=NonNegativeReals,
            initialize=self.line_data["SLmax_data"][self.model.L],
            mutable=True,
        )  # real power line limit
        self.model.SLmaxT = Param(
            self.model.TRANSF,
            within=NonNegativeReals,
            initialize=self.trafo_data.SLmaxT_data[self.model.TRANSF],
            mutable=True,
        )  # real power transformer limit

        # --- transformer tap ratio limits ---

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
            self.model.LineAngleSet = Set(initialize=list(line_bounds))

            def _line_angle_rule(model, l, t):
                f, to, amin, amax = line_bounds[l]
                return amin, model.delta[f, t] - model.delta[to, t], amax

            self.model.line_angle_diff = Constraint(
                self.model.LineAngleSet, self.model.T, rule=_line_angle_rule
            )

        if trafo_bounds:
            self.model.TrafoAngleSet = Set(initialize=list(trafo_bounds))

            def _trafo_angle_rule(model, l, t):
                f, to, amin, amax = trafo_bounds[l]
                return amin, model.delta[f, t] - model.delta[to, t], amax

            self.model.trafo_angle_diff = Constraint(
                self.model.TrafoAngleSet, self.model.T, rule=_trafo_angle_rule
            )

    def add_tap_changer_linear(self, max_tap_change_per_step=None):
        """Enable continuous (linear) OLTC tap-ratio optimisation.

        Unfixes the time-indexed ``Tap[tr, t]`` variables and adds per-step
        tap-ratio bounds derived from ``net.trafo.tap_min`` /
        ``net.trafo.tap_max``.

        Args:
            max_tap_change_per_step: Maximum allowed change in tap ratio
                between consecutive time steps.  ``None`` (default) means
                unconstrained.  Adds ``model.tap_rate_up`` /
                ``model.tap_rate_down`` constraints when set.
        """
        ratio_min, ratio_max = _calc_tap_min_max_mp(self)
        self.trafo_data = self.trafo_data.assign(
            **{"tap_min_data": ratio_min, "tap_max_data": ratio_max}
        )

        self.model.Tap_min = Param(
            self.model.TRANSF,
            within=Reals,
            initialize=self.trafo_data.tap_min_data[self.model.TRANSF],
        )
        self.model.Tap_max = Param(
            self.model.TRANSF,
            within=Reals,
            initialize=self.trafo_data.tap_max_data[self.model.TRANSF],
        )

        def trafo_tap_linear_bounds(model, tr, t):
            return model.Tap_min[tr], model.Tap[tr, t], model.Tap_max[tr]

        self.model.Tap_linear_constr = Constraint(
            self.model.TRANSF, self.model.T, rule=trafo_tap_linear_bounds
        )
        self.unfix_vars("Tap")

        if max_tap_change_per_step is not None:
            t_list = sorted(self.model.T)
            t_pairs = list(zip(t_list[:-1], t_list[1:]))

            def tap_rate_up(model, tr, t_prev, t_next):
                return (
                    model.Tap[tr, t_next] - model.Tap[tr, t_prev]
                    <= max_tap_change_per_step
                )

            def tap_rate_down(model, tr, t_prev, t_next):
                return (
                    model.Tap[tr, t_next] - model.Tap[tr, t_prev]
                    >= -max_tap_change_per_step
                )

            self.model.tap_rate_up = Constraint(
                self.model.TRANSF, t_pairs, rule=tap_rate_up
            )
            self.model.tap_rate_down = Constraint(
                self.model.TRANSF, t_pairs, rule=tap_rate_down
            )

    def add_tap_changer_discrete(self):
        """Enable discrete OLTC tap-position optimisation.

        Introduces an integer variable ``Tap_pos[tr, t]`` per transformer and
        time step, and links it to the continuous ``Tap[tr, t]`` via an
        equality constraint.  Requires a MIP/MINLP solver (MindtPy, Gurobi).
        """
        tap_neutral = self.net.trafo.tap_neutral
        tap_step = self.net.trafo.tap_step_percent / 100.0
        tap_pos_max = self.net.trafo.tap_max
        tap_pos_min = self.net.trafo.tap_min
        tap_side_data = np.where(self.net.trafo.tap_side == "lv", 1, 0)

        self.trafo_data = self.trafo_data.assign(
            **{
                "tap_neutral": tap_neutral,
                "tap_step": tap_step,
                "tap_pos_max": tap_pos_max,
                "tap_pos_min": tap_pos_min,
                "tap_side_data": tap_side_data,
            }
        )

        self.model.Tap_pos = Var(
            self.model.TRANSF,
            self.model.T,
            within=Integers,
            initialize=0,
        )
        self.model.Tap_pos_min = Param(
            self.model.TRANSF,
            within=Integers,
            initialize=self.trafo_data.tap_pos_min[self.model.TRANSF],
        )
        self.model.Tap_pos_max = Param(
            self.model.TRANSF,
            within=Integers,
            initialize=self.trafo_data.tap_pos_max[self.model.TRANSF],
        )
        self.model.Tap_neutral = Param(
            self.model.TRANSF,
            within=Integers,
            initialize=self.trafo_data.tap_neutral[self.model.TRANSF],
        )
        self.model.Tap_step = Param(
            self.model.TRANSF,
            within=Reals,
            initialize=self.trafo_data.tap_step[self.model.TRANSF],
        )
        self.model.Tap_side = Param(
            self.model.TRANSF,
            initialize=self.trafo_data.tap_side_data[self.model.TRANSF],
        )

        def trafo_tap_pos_min_max(model, tr, t):
            return (
                model.Tap_pos_min[tr],
                model.Tap_pos[tr, t],
                model.Tap_pos_max[tr],
            )

        self.model.Tap_pos_constr = Constraint(
            self.model.TRANSF, self.model.T, rule=trafo_tap_pos_min_max
        )

        def trafo_tap_discrete(model, tr, t):
            if model.Tap_side[tr]:
                return model.Tap[tr, t] == 1 / (
                    1
                    + (model.Tap_pos[tr, t] - model.Tap_neutral[tr])
                    * model.Tap_step[tr]
                )
            return (
                model.Tap[tr, t]
                == 1.0
                + (model.Tap_pos[tr, t] - model.Tap_neutral[tr])
                * model.Tap_step[tr]
            )

        self.model.Tap_discrete_constr = Constraint(
            self.model.TRANSF, self.model.T, rule=trafo_tap_discrete
        )
        self.unfix_vars("Tap")
