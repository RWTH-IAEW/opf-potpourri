"""Multi-period OPF mixin: adds operational limit constraints and thermal
ratings over time."""

import copy

from pyomo.environ import *
from potpourri.models_multi_period.basemodel_multi_period import (
    Basemodel_multi_period,
)
from potpourri.technologies.generator import Generator_multi_period
from potpourri.technologies.sgens import Sgens_multi_period
from potpourri.technologies.demand import Demand_multi_period
import numpy as np


class OPF_multi_period(Basemodel_multi_period):
    """OPF mixin for multi-period models: provides line/transformer ratings
    and generator/demand limits."""

    def __init__(self, net, toT, fromT=None, pf=1):
        super().__init__(net, toT, fromT, pf)

    def __calc_SLmax(self, max_loading_percent=100):
        vr = self.net.bus.loc[
            self.net.line["from_bus"].values, "vn_kv"
        ].values * np.sqrt(3.0)
        max_i_ka = self.net.line.max_i_ka.values
        df = self.net.line.df.values
        return (
            max_loading_percent
            / 100.0
            * max_i_ka
            * df
            * self.net.line.parallel.values
            * vr
            / self.baseMVA
        )

    def _calc_opf_parameters(self, **kwargs):
        """Compute line/transformer ratings and call generator/demand limit
        methods on flexibility objects."""
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

    def add_tap_changer_linear(self):
        def _calc_tap_min_max(self):
            vn_trafo_hv_min, vn_trafo_lv_min, shift_min = _calc_tap_shift(
                self, tap_pos=self.net.trafo.tap_min
            )
            vn_trafo_hv_max, vn_trafo_lv_max, shift_max = _calc_tap_shift(
                self, tap_pos=self.net.trafo.tap_max
            )
            ratio_min = _calc_nominal_ratio_from_dataframe(
                self, vn_trafo_hv_min, vn_trafo_lv_min
            )
            ratio_max = _calc_nominal_ratio_from_dataframe(
                self, vn_trafo_hv_max, vn_trafo_lv_max
            )
            return ratio_min, ratio_max

        def _calc_nominal_ratio_from_dataframe(self, vn_hv_kv, vn_lv_kv):
            """
            Calculates (Vectorized) the off nominal tap ratio::

                          (vn_hv_kv / vn_lv_kv) / (ub1_in_kv / ub2_in_kv)

            INPUT:
                **net** (Dataframe) - The net for which to calc the tap ratio.

                **vn_hv_kv** (1d array, float) - The adjusted nominal high
                voltages

                **vn_lv_kv** (1d array, float) - The adjusted nominal low
                voltages

            OUTPUT:
                **tab** (1d array, float) - The off-nominal tap ratio
            """
            # Calculating tab (transformer off nominal turns ratio)
            tap_rat = vn_hv_kv / vn_lv_kv
            hv_bus = self.net.trafo.hv_bus
            lv_bus = self.net.trafo.lv_bus
            nom_rat = (
                self.net.bus.vn_kv[hv_bus].values
                / self.net.bus.vn_kv[lv_bus].values
            )
            return tap_rat / nom_rat

        def _calc_tap_shift(self, tap_pos=None):
            """
            Adjust the nominal voltage vnh and vnl to the active tab position
            "tap_pos". If "side" is 1 (high-voltage side) the high voltage
            vnh is adjusted. If "side" is 2 (low-voltage side) the low
            voltage vnl is adjusted

            INPUT:
                **net** - The pandapower format network

                **trafo** (Dataframe) - The dataframe in
                pd_net["structure"]["trafo"]
                which contains transformer calculation values.

            OUTPUT:
                **vn_hv_kv** (1d array, float) - The adusted high voltages

                **vn_lv_kv** (1d array, float) - The adjusted low voltages

                **trafo_shift** (1d array, float) - phase shift angle

            """
            vnh = copy.deepcopy(self.net.trafo.vn_hv_kv.values)
            vnl = copy.deepcopy(self.net.trafo.vn_lv_kv.values)
            trafo_shift = self.net.trafo.shift_degree.values

            if tap_pos is None:
                tap_pos = self.net.trafo.tap_pos
            tap_neutral = self.net.trafo.tap_neutral
            tap_diff = tap_pos - tap_neutral
            tap_phase_shifter = self.net.trafo.tap_phase_shifter
            tap_side = self.net.trafo.tap_side
            tap_step_percent = self.net.trafo.tap_step_percent
            tap_step_degree = self.net.trafo.tap_step_degree

            def cos(x):
                return np.cos(np.deg2rad(x))

            def sin(x):
                return np.sin(np.deg2rad(x))

            def arctan(x):
                return np.rad2deg(np.arctan(x))

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
                        tap_step_percent[tap_complex]
                        * tap_diff[tap_complex]
                        / 100
                    )
                    tap_angles = (tap_step_degree[tap_complex]).fillna(0)
                    u1 = vn[tap_complex]
                    du = u1 * tap_steps.fillna(0)
                    vn[tap_complex] = np.sqrt(
                        (u1 + du * cos(tap_angles)) ** 2
                        + (du * sin(tap_angles)) ** 2
                    )
                    trafo_shift[tap_complex] += arctan(
                        direction
                        * du
                        * sin(tap_angles)
                        / (u1 + du * cos(tap_angles))
                    )
                if phase_shifters.any():
                    degree_is_set = (
                        tap_step_degree[phase_shifters].fillna(0) != 0
                    )
                    percent_is_set = (
                        tap_step_percent[phase_shifters].fillna(0) != 0
                    )
                    if (degree_is_set & percent_is_set).any():
                        raise UserWarning(
                            "Both tap_step_degree and tap_step_percent set"
                            " for ideal phase shifter"
                        )
                    trafo_shift[phase_shifters] += np.where(
                        (degree_is_set),
                        (
                            direction
                            * tap_diff[phase_shifters]
                            * tap_step_degree[phase_shifters]
                        ),
                        (
                            direction
                            * 2
                            * np.rad2deg(
                                np.arcsin(
                                    tap_diff[phase_shifters]
                                    * tap_step_percent[phase_shifters]
                                    / 100
                                    / 2
                                )
                            )
                        ),
                    )

            return vnh, vnl, trafo_shift

        ratio_min, ratio_max = _calc_tap_min_max(self)
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

        def trafo_tap_linear_bounds(model, tr):
            return model.Tap_min[tr], model.Tap[tr], model.Tap_max[tr]

        self.model.Tap_linear_constr = Constraint(
            self.model.TRANSF, rule=trafo_tap_linear_bounds
        )

        self.unfix_vars("Tap")

    def add_tap_changer_discrete(self):
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
            self.model.TRANSF, within=Integers, initialize=0.0
        )  # transformer tap position
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
        )  # transformer tap neutral position
        self.model.Tap_step = Param(
            self.model.TRANSF,
            within=Reals,
            initialize=self.trafo_data.tap_step[self.model.TRANSF],
        )  # transformer tap step size
        self.model.Tap_side = Param(
            self.model.TRANSF,
            initialize=self.trafo_data.tap_side_data[self.model.TRANSF],
        )  # transformer tap side; 0: hv, 1: lv

        def trafo_tap_pos_min_max(model, tr):
            return (
                model.Tap_pos_min[tr],
                model.Tap_pos[tr],
                model.Tap_pos_max[tr],
            )

        self.model.Tap_pos_constr = Constraint(
            self.model.TRANSF, rule=trafo_tap_pos_min_max
        )

        def trafo_tap_discrete(model, tr):
            if model.Tap_side[tr]:
                # tap side: lv
                return model.Tap[tr] == 1 / (
                    1
                    + (model.Tap_pos[tr] - model.Tap_neutral[tr])
                    * model.Tap_step[tr]
                )

            return (
                model.Tap[tr]
                == 1.0
                + (model.Tap_pos[tr] - model.Tap_neutral[tr])
                * model.Tap_step[tr]
            )

        self.model.Tap_discrete_constr = Constraint(
            self.model.TRANSF, rule=trafo_tap_discrete
        )

        self.unfix_vars("Tap")
