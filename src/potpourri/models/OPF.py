import copy

import pyomo.environ as pyo
from src.potpourri.models.basemodel import Basemodel
import numpy as np


class OPF(Basemodel):
    """OPF mixin that adds power and thermal limit constraints to a power flow model.

    Intended for use via multiple inheritance alongside AC or DC:

        class ACOPF(AC, OPF): ...

    Provides methods to read operational limits from a pandapower network and
    attach the corresponding Pyomo sets, parameters, and constraints.
    """

    def __init__(self, net):
        super().__init__(net)

    def generation_real_power_limits(self):
        """Read generator real power limits from net into generation_data.

        Populates generation_data['max_p'] and ['min_p'] (per-unit on baseMVA).
        Non-controllable generators are pinned to their current p_mw setpoint.
        """
        max_p = np.full(len(self.generation_data), 1e9) / self.baseMVA
        min_p = np.full(len(self.generation_data), -1e9) / self.baseMVA

        for element, (f, t) in self.net._gen_order.items():
            if 'max_p_mw' in self.net[element]:
                max_p[f:t] = self.net[element].max_p_mw.fillna(1e9).values / self.baseMVA
            if 'min_p_mw' in self.net[element]:
                min_p[f:t] = self.net[element].min_p_mw.fillna(-1e9).values / self.baseMVA

        if 'controllable' in self.net.gen:
            controllable = self.net["gen"]["controllable"].values
            not_controllable = ~controllable.astype(bool)

            if np.any(not_controllable):
                f, t = self.net._gen_order['gen']

                p_mw = self.net["gen"]["p_mw"].values[not_controllable]

                not_controllable_gens = np.arange(f, t)[not_controllable]
                max_p[not_controllable_gens] = p_mw / self.baseMVA
                min_p[not_controllable_gens] = p_mw / self.baseMVA

        self.generation_data['max_p'] = max_p
        self.generation_data['min_p'] = min_p

    def static_generation_real_power_limits(self):
        """Read static generator real power limits from net.sgen.

        Populates static_generation_data['max_p'], ['min_p'], and
        ['controllable'] (per-unit). Defaults: max = p_mw, min = 0.
        """
        if 'controllable' in self.net.sgen:
            self.static_generation_data['controllable'] = self.net.sgen.controllable.values
        else:
            self.static_generation_data['controllable'] = False

        if 'max_p_mw' in self.net.sgen:
            self.static_generation_data['max_p'] = self.net.sgen.max_p_mw.fillna(
                self.net.sgen.p_mw).values / self.baseMVA
        else:
            self.static_generation_data['max_p'] = self.net.sgen.p_mw.values / self.baseMVA
        if 'min_p_mw' in self.net.sgen:
            self.static_generation_data['min_p'] = self.net.sgen.min_p_mw.fillna(0).values / self.baseMVA
        else:
            self.static_generation_data['min_p'] = np.zeros(len(self.net.sgen.index))

    def get_demand_real_power_data(self):
        """Read load real power bounds from net.load.

        Populates self.PDmax_data and self.PDmin_data (per-unit).
        Falls back to p_mw and 0 respectively if columns are absent.
        """
        if 'controllable' not in self.net.load:
            self.demand_controllable_set = None  # create empty  pyo.Set if no controllable load exist
        else:
            self.demand_controllable_set = self.net.load.index[self.net.load.controllable.fillna(False).astype(bool)]

        # add rows with active demand limits if not existing
        if 'max_p_mw' not in self.net.load:
            self.net.load['max_p_mw'] = self.net.load.p_mw

        if 'min_p_mw' not in self.net.load:
            self.net.load['min_p_mw'] = np.zeros(len(self.net.load.index))

        # demand limits for loads
        self.PDmax_data = self.net.load.max_p_mw.fillna(self.net.load.p_mw) / self.baseMVA
        self.PDmin_data = self.net.load.min_p_mw.fillna(0) / self.baseMVA

    def __calc_SLmax(self, max_loading_percent=100):
        vr = self.net.bus.loc[self.net.line["from_bus"].values, "vn_kv"].values * np.sqrt(3.)
        max_i_ka = self.net.line.max_i_ka.values
        df = self.net.line.df.values
        return max_loading_percent / 100. * max_i_ka * df * self.net.line.parallel.values * vr / self.baseMVA

    def _calc_opf_parameters(self, **kwargs):
        """Compute all OPF limit data from the network before model construction.

        Calculates line apparent power limits (SLmax) and transformer limits
        (SLmaxT), then reads generator, static generator, and demand limits.
        """
        max_load = self.net.line.max_loading_percent.values if "max_loading_percent" in self.net.line else 100.
        self.line_data['SLmax_data'] = self.__calc_SLmax(max_load)

        # maximum transformer loading
        max_load_T = self.net.trafo.max_loading_percent.fillna(
            100.) / 100. if "max_loading_percent" in self.net.trafo else 1.
        sn_mva = self.net.trafo.sn_mva
        df_T = self.net.trafo.df
        SLmaxT_data = max_load_T * sn_mva * df_T * self.net.trafo.parallel / self.baseMVA
        self.trafo_data['SLmaxT_data'] = SLmaxT_data.values

        # self.c0_data = pd.Series([1] * len(self.gen_all_set), self.gen_all_set)
        # self.c1_data = pd.Series([2] * len(self.gen_all_set), self.gen_all_set)
        # self.c2_data = pd.Series([3] * len(self.gen_all_set), self.gen_all_set)

        # self.get_generator_real_power_data()
        self.static_generation_real_power_limits()
        self.generation_real_power_limits()
        self.get_demand_real_power_data()

    def add_OPF(self, **kwargs):
        """Attach OPF sets, parameters, and constraints to self.model.

        Calls _calc_opf_parameters(), then adds sets sGc and Dc, parameters
        PGmax/min, sPGmax/min, PDmax/min, SLmax, SLmaxT, and real-power bound
        constraints for generators, static generators, and controllable loads.

        Args:
            **kwargs: Forwarded to _calc_opf_parameters (e.g. max_loading_percent).
        """
        self._calc_opf_parameters(**kwargs)

        # controllable generation
        self.model.sGc =  pyo.Set(within=self.model.sG,
                             initialize=self.static_generation_data.index[self.static_generation_data.controllable & self.static_generation_data.in_service])
        self.model.sPGmax =  pyo.Param(self.model.sGc, initialize=self.static_generation_data.max_p[self.model.sGc])
        self.model.sPGmin =  pyo.Param(self.model.sGc, initialize=self.static_generation_data.min_p[self.model.sGc])
        self.model.PGmax =  pyo.Param(self.model.G, initialize=self.generation_data['max_p'][self.model.G])
        self.model.PGmin =  pyo.Param(self.model.G, initialize=self.generation_data['min_p'][self.model.G])

        # controllable loads
        self.model.Dc =  pyo.Set(within=self.model.D, initialize=self.demand_controllable_set)  # controllable loads
        self.model.PDmax =  pyo.Param(self.model.D, initialize=self.PDmax_data[self.model.D])
        self.model.PDmin =  pyo.Param(self.model.D, initialize=self.PDmin_data[self.model.D])

        # lines and transformer chracteristics and ratings
        self.model.SLmax =  pyo.Param(self.model.L, within=pyo.NonNegativeReals,
                                 initialize=self.line_data['SLmax_data'][self.model.L],
                                 mutable=True)  # real power line limit
        self.model.SLmaxT =  pyo.Param(self.model.TRANSF, within=pyo.NonNegativeReals,
                                  initialize=self.trafo_data.SLmaxT_data[self.model.TRANSF],
                                  mutable=True)  # real power transformer limit

        # cost data
        # self.model.c2 =  pyo.Param(self.model.G, within=pyo.NonNegativeReals,
        #                       initialize=self.c2_data)  # generator cost coefficient c2 (*pG^2)
        # self.model.c1 =  pyo.Param(self.model.G, within=pyo.NonNegativeReals,
        #                       initialize=self.c1_data)  # generator cost coefficient c1 (*pG)
        # self.model.c0 =  pyo.Param(self.model.G, within=pyo.NonNegativeReals,
        #                       initialize=self.c0_data)  # generator cost coefficient c0
        # self.model.VOLL =  pyo.Param(self.model.D, within=pyo.Reals, initialize=10000)  # value of lost load

        # --- static generator power limits ---
        def static_generation_real_power_bounds(model, g):
            model.psG[g].unfix()
            return model.sPGmin[g], model.psG[g], model.sPGmax[g]

        self.model.PsG_Constraint = pyo.Constraint(self.model.sGc, rule=static_generation_real_power_bounds)

        # --- generation real power limits ---
        def real_power_bounds(model, g):
            model.pG[g].unfix()
            return model.PGmin[g], model.pG[g], model.PGmax[g]

        self.model.PG_Constraint = pyo.Constraint(self.model.G, rule=real_power_bounds)

        # --- demand limits ---
        def real_demand_bounds(model, d):
            return model.PDmin[d], model.pD[d], model.PDmax[d]

        self.model.PD_Constraint = pyo.Constraint(self.model.Dc, rule=real_demand_bounds)

        # --- transformer tap ratio limits ---

    def add_tap_changer_linear(self):
        """Enable continuous (linear) transformer tap ratio optimisation.

        Unfixes Tap variables and adds [Tap_min, Tap_max] bounds derived from
        net.trafo tap_min / tap_max columns.
        """
        def _calc_tap_min_max(self):
            vn_trafo_hv_min, vn_trafo_lv_min, shift_min = _calc_tap_shift(self, tap_pos=self.net.trafo.tap_min)
            vn_trafo_hv_max, vn_trafo_lv_max, shift_max = _calc_tap_shift(self, tap_pos=self.net.trafo.tap_max)
            ratio_min = _calc_nominal_ratio_from_dataframe(self, vn_trafo_hv_min, vn_trafo_lv_min)
            ratio_max = _calc_nominal_ratio_from_dataframe(self, vn_trafo_hv_max, vn_trafo_lv_max)
            return ratio_min, ratio_max

        def _calc_nominal_ratio_from_dataframe(self, vn_hv_kv, vn_lv_kv):
            """
            Calculates (Vectorized) the off nominal tap ratio::

                          (vn_hv_kv / vn_lv_kv) / (ub1_in_kv / ub2_in_kv)

            INPUT:
                **net** (Dataframe) - The net for which to calc the tap ratio.

                **vn_hv_kv** (1d array, float) - The adjusted nominal high voltages

                **vn_lv_kv** (1d array, float) - The adjusted nominal low voltages

            OUTPUT:
                **tab** (1d array, float) - The off-nominal tap ratio
            """
            # Calculating tab (transformer off nominal turns ratio)
            tap_rat = vn_hv_kv / vn_lv_kv
            hv_bus = self.net.trafo.hv_bus
            lv_bus = self.net.trafo.lv_bus
            nom_rat = self.net.bus.vn_kv[hv_bus].values / self.net.bus.vn_kv[lv_bus].values
            return tap_rat / nom_rat

        def _calc_tap_shift(self, tap_pos=None):
            """
            Adjust the nominal voltage vnh and vnl to the active tab position "tap_pos".
            If "side" is 1 (high-voltage side) the high voltage vnh is adjusted.
            If "side" is 2 (low-voltage side) the low voltage vnl is adjusted

            INPUT:
                **net** - The pandapower format network

                **trafo** (Dataframe) - The dataframe in pd_net["structure"]["trafo"]
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

            cos = lambda x: np.cos(np.deg2rad(x))
            sin = lambda x: np.sin(np.deg2rad(x))
            arctan = lambda x: np.rad2deg(np.arctan(x))

            for side, vn, direction in [("hv", vnh, 1), ("lv", vnl, -1)]:
                phase_shifters = tap_phase_shifter & (tap_side == side)
                tap_complex = np.isfinite(tap_step_percent) & np.isfinite(tap_pos) & (tap_side == side) & \
                              ~phase_shifters
                if tap_complex.any():
                    tap_steps = tap_step_percent[tap_complex] * tap_diff[tap_complex] / 100
                    tap_angles = (tap_step_degree[tap_complex]).fillna(0)
                    u1 = vn[tap_complex]
                    du = u1 * tap_steps.fillna(0)
                    vn[tap_complex] = np.sqrt((u1 + du * cos(tap_angles)) ** 2 + (du * sin(tap_angles)) ** 2)
                    trafo_shift[tap_complex] += (arctan(direction * du * sin(tap_angles) /
                                                        (u1 + du * cos(tap_angles))))
                if phase_shifters.any():
                    degree_is_set = tap_step_degree[phase_shifters].fillna(0) != 0
                    percent_is_set = tap_step_percent[phase_shifters].fillna(0) != 0
                    if (degree_is_set & percent_is_set).any():
                        raise ValueError(
                            "Both tap_step_degree and tap_step_percent set for ideal phase shifter.")
                    trafo_shift[phase_shifters] += np.where(
                        (degree_is_set),
                        (direction * tap_diff[phase_shifters] * tap_step_degree[phase_shifters]),
                        (direction * 2 * np.rad2deg(np.arcsin(tap_diff[phase_shifters] *
                                                              tap_step_percent[phase_shifters] / 100 / 2)))
                    )

            return vnh, vnl, trafo_shift

        ratio_min, ratio_max = _calc_tap_min_max(self)
        self.trafo_data = self.trafo_data.assign(**{"tap_min_data": ratio_min, "tap_max_data": ratio_max})

        self.model.Tap_min =  pyo.Param(self.model.TRANSF, within=pyo.Reals,
                                   initialize=self.trafo_data.tap_min_data[self.model.TRANSF])
        self.model.Tap_max =  pyo.Param(self.model.TRANSF, within=pyo.Reals,
                                   initialize=self.trafo_data.tap_max_data[self.model.TRANSF])

        def trafo_tap_linear_bounds(model, t):
            return model.Tap_min[t], model.Tap[t], model.Tap_max[t]

        self.model.Tap_linear_constr = pyo.Constraint(self.model.TRANSF, rule=trafo_tap_linear_bounds)

        self.unfix_vars('Tap')

    def add_tap_changer_discrete(self):
        """Enable discrete transformer tap ratio optimisation.

        Introduces integer variable Tap_pos and links it to Tap via a
        constraint. Suitable for use with MIP/MINLP solvers (MindtPy, Gurobi).
        """
        tap_neutral = self.net.trafo.tap_neutral
        tap_step = self.net.trafo.tap_step_percent / 100.
        tap_pos_max = self.net.trafo.tap_max
        tap_pos_min = self.net.trafo.tap_min
        tap_side_data = np.where(self.net.trafo.tap_side == 'lv', 1, 0)

        self.trafo_data = self.trafo_data.assign(**{"tap_neutral": tap_neutral, "tap_step": tap_step,
                                                    "tap_pos_max": tap_pos_max, "tap_pos_min": tap_pos_min,
                                                    "tap_side_data": tap_side_data})

        self.model.Tap_pos = pyo.Var(self.model.TRANSF, within=pyo.Integers, initialize=0.)  # transformer tap position
        self.model.Tap_pos_min =  pyo.Param(self.model.TRANSF, within=pyo.Integers,
                                       initialize=self.trafo_data.tap_pos_min[self.model.TRANSF])
        self.model.Tap_pos_max =  pyo.Param(self.model.TRANSF, within=pyo.Integers,
                                       initialize=self.trafo_data.tap_pos_max[self.model.TRANSF])
        self.model.Tap_neutral =  pyo.Param(self.model.TRANSF, within=pyo.Integers,
                                       initialize=self.trafo_data.tap_neutral[
                                           self.model.TRANSF])  # transformer tap neutral position
        self.model.Tap_step =  pyo.Param(self.model.TRANSF, within=pyo.Reals,
                                    initialize=self.trafo_data.tap_step[self.model.TRANSF])  # transformer tap step size
        self.model.Tap_side =  pyo.Param(self.model.TRANSF, initialize=self.trafo_data.tap_side_data[self.model.TRANSF]) # transformer tap side; 0: hv, 1: lv

        def trafo_tap_pos_min_max(model, t):
            return model.Tap_pos_min[t], model.Tap_pos[t], model.Tap_pos_max[t]

        self.model.Tap_pos_constr = pyo.Constraint(self.model.TRANSF, rule=trafo_tap_pos_min_max)

        def trafo_tap_discrete(model, t):
            if model.Tap_side[t]:
                # tap side: lv
                return model.Tap[t] == 1 / (1 + (model.Tap_pos[t] - model.Tap_neutral[t]) * model.Tap_step[t])

            return model.Tap[t] == 1. + (model.Tap_pos[t] - model.Tap_neutral[t]) * model.Tap_step[t]

        self.model.Tap_discrete_constr = pyo.Constraint(self.model.TRANSF, rule=trafo_tap_discrete)

        self.unfix_vars('Tap')
