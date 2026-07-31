"""Multi-period AC OPF combining multi-period AC power flow and OPF
operational limits."""

from pyomo.environ import *
from potpourri.models_multi_period.AC_multi_period import AC_multi_period
from potpourri.models_multi_period.OPF_multi_period import OPF_multi_period
from potpourri.technologies.generator import Generator_multi_period
from potpourri.technologies.demand import Demand_multi_period
from potpourri.technologies.windpower import Windpower_multi_period
from potpourri.technologies.sgens import Sgens_multi_period
from potpourri.technologies.q_control import resolve_grid_code
import numpy as np
import pandas as pd
from loguru import logger


class ACOPF_multi_period(AC_multi_period, OPF_multi_period):
    """Multi-period AC OPF model combining AC power flow and OPF constraints
    over a time horizon."""

    def __init__(self, net, toT, fromT=None, pf=1):
        super().__init__(net, toT, fromT, pf)

    def _calc_opf_parameters(self):
        """Extend OPF parameter calculation with AC-specific limits: voltage
        bounds, Q limits, Q-curve data."""
        super()._calc_opf_parameters()

        max_vm_pu, min_vm_pu = self.get_v_limits()
        self.v_limits = (max_vm_pu, min_vm_pu)

        # create sgen instance
        sgens_object = next(
            (
                obj
                for obj in self.flexibilities
                if isinstance(obj, Sgens_multi_period)
            ),
            None,
        )
        sgens_object.static_generation_reactive_power_limits(
            self.model
        )  # gives the model now instead of the net

        # Populate Q(P)/Q(U) characteristic data for sgens with var_q set.
        # This covers both wind and PV sgens annotated with grid-code variants.
        # The grid code applies model-wide; set it via add_OPF(grid_code=...).
        if "var_q" in self.net.sgen:
            sgens_object.static_generation_q_ctrl_data(
                self.net,
                grid_code=getattr(self, "_grid_code", None),
                qu_deadband=getattr(self, "_qu_deadband", None),
            )

        # Populate inverter S² rating data when sn_mva is present.
        if "sn_mva" in self.net.sgen:
            sgens_object.static_generation_inverter_data(self.net)

        # P(U) active-power curtailment (VDE-AR-N 4105 §8.5).
        if "pu_curtail" in self.net.sgen:
            sgens_object.static_generation_pu_curtail_data(self.net)

        # Fixed cos(φ) equality mode.
        if "fixed_cos_phi" in self.net.sgen:
            sgens_object.static_generation_fixed_cos_phi_data(self.net)

        # cos(φ)(P) profile (quadratic P-Q curve).
        if "cos_phi_p_profile" in self.net.sgen:
            sgens_object.static_generation_cpp_data(self.net)

        # get the object of class 'Windpower' from the 'flexibilities' list
        if "windpot_p_mw" in self.net.bus:
            windpower_object = next(
                (
                    obj
                    for obj in self.flexibilities
                    if isinstance(obj, Windpower_multi_period)
                ),
                None,
            )
            windpower_object.static_generation_wind_var_q(
                self.net, grid_code=getattr(self, "_grid_code", None)
            )

        # create generator instance
        generator_object = next(
            (
                obj
                for obj in self.flexibilities
                if isinstance(obj, Generator_multi_period)
            ),
            None,
        )
        generator_object.generation_reactive_power_limits_acopf()

        # create demand instance
        demand_object = next(
            (
                obj
                for obj in self.flexibilities
                if isinstance(obj, Demand_multi_period)
            ),
            None,
        )
        demand_object.get_demand_reactive_data(self.model)

    def get_v_limits(self):
        """Read per-bus voltage bounds from net.bus, keyed by ppc bus number.

        Returns:
            tuple: (max_vm_pu, min_vm_pu) as :class:`~pandas.Series` indexed
            by **ppc** bus number, covering only the ppc buses that a
            pandapower bus maps onto.

        Note:
            The index is ppc bus numbers, not pandapower bus indices, because
            every consumer looks these up through ``self.bus_lookup``.  Plain
            positional arrays were correct only while the two numbering spaces
            coincided, which fails on grids where pandapower's ppc conversion
            adds auxiliary buses for node-node switches.  Those auxiliary
            buses have no pandapower row and hence no user-supplied limits, so
            they are absent here and their voltage follows from the network
            equations.
        """
        n_bus = len(self.net.bus.index)
        if "max_vm_pu" in self.net.bus:
            vmax = self.net.bus.max_vm_pu.values
        else:
            vmax = np.full(n_bus, 1.1)

        if "min_vm_pu" in self.net.bus:
            vmin = self.net.bus.min_vm_pu.values
        else:
            vmin = np.full(n_bus, 0.9)

        ppc_of_bus = self.pd_bus_to_ppc
        max_vm_pu = pd.Series(np.asarray(vmax, dtype=float), index=ppc_of_bus)
        min_vm_pu = pd.Series(np.asarray(vmin, dtype=float), index=ppc_of_bus)
        # Several pandapower buses can fuse onto one ppc bus; keep the
        # tightest band.
        max_vm_pu = max_vm_pu.groupby(level=0).min()
        min_vm_pu = min_vm_pu.groupby(level=0).max()

        if any(self.net.gen.index):
            self.add_generator_v_limits(max_vm_pu, min_vm_pu)

        return max_vm_pu, min_vm_pu

    def add_generator_v_limits(self, max_vm_pu, min_vm_pu):
        """Apply per-generator voltage limits, overriding bus defaults where
        applicable."""
        # check max_vm_pu / min_vm_pu bus limit violation by gens
        gen_buses = self.bus_lookup[self.net.gen.bus.values]
        if "max_vm_pu" in self.net["gen"].columns:
            v_max_bound = (
                max_vm_pu.loc[gen_buses].to_numpy()
                < self.net["gen"]["max_vm_pu"].values
            )
            if np.any(v_max_bound):
                bound_gens = self.net["gen"].index.values[v_max_bound]
                logger.warning(
                    "gen max_vm_pu > bus max_vm_pu for gens {}. "
                    "Setting bus limit for these gens.",
                    bound_gens,
                )
                # set only vm of gens which do not violate the limits
                max_vm_pu.loc[gen_buses[~v_max_bound]] = self.net["gen"][
                    "max_vm_pu"
                ].values[~v_max_bound]
            else:
                # set vm of all gens
                max_vm_pu.loc[gen_buses] = self.net["gen"]["max_vm_pu"].values

        if "min_vm_pu" in self.net["gen"].columns:
            v_min_bound = (
                self.net["gen"]["min_vm_pu"].values
                < min_vm_pu.loc[gen_buses].to_numpy()
            )
            if np.any(v_min_bound):
                bound_gens = self.net["gen"].index.values[v_min_bound]
                logger.warning(
                    "gen min_vm_pu < bus min_vm_pu for gens {}. "
                    "Setting bus limit for these gens.",
                    bound_gens,
                )
                # set only vm of gens which do not violate the limits
                min_vm_pu.loc[gen_buses[~v_min_bound]] = self.net["gen"][
                    "min_vm_pu"
                ].values[~v_min_bound]
            else:
                # set vm of all gens
                min_vm_pu.loc[gen_buses] = self.net["gen"]["min_vm_pu"].values

        if "controllable" in self.net.gen:
            controllable = self.net["gen"]["controllable"].values
            not_controllable = ~controllable.astype(bool)

            # get voltage setpoints for not controllable generators
            if np.any(not_controllable):
                bus = self.net["gen"]["bus"].values[not_controllable]
                vm_pu = self.net["gen"]["vm_pu"].values[not_controllable]

                not_controllable_buses = self.bus_lookup[bus]
                max_vm_pu[not_controllable_buses] = vm_pu
                min_vm_pu[not_controllable_buses] = vm_pu

        return max_vm_pu, min_vm_pu

    def add_OPF(self, grid_code=None, qu_deadband=None, **kwargs):
        """Extend OPF.add_OPF() with voltage bounds, AC thermal limits, and
        reactive power constraints.

        Args:
            grid_code: Technical connection rule supplying the Q(P)/Q(U)
                capability envelope and the P(U) / cos(phi)(P) thresholds,
                as accepted by
                :func:`~potpourri.technologies.q_control.resolve_grid_code`.
                Applies model-wide and defaults to VDE-AR-N 4120.
            qu_deadband: Replace the Q(U) capability *area* with a Q(U)
                *characteristic* that has a dead band, pinning Q to a curve
                of voltage instead of bounding it.  ``None`` keeps the area;
                ``True`` uses the grid code's own QV plateau; a
                ``(v_low, v_high)`` pair sets the dead band explicitly; a
                :class:`~potpourri.technologies.q_control.QVCurve` is used
                as given.  Applies model-wide.

                The feasible set pinches to Q = 0 inside the dead band and
                is therefore **not convex**, so this builds an integer
                piecewise block per sgen and time step and needs a
                MIP-capable solver (MindtPy, CBC, GLPK, Gurobi).
            **kwargs: Forwarded to the base implementation.
        """
        # Consumed by _calc_opf_parameters, which super().add_OPF() reaches.
        # Resolved here rather than downstream so that ``self._grid_code``
        # means the same thing as it does on the single-period model: the
        # GridCode itself, not whatever selector was passed in.
        self._grid_code = resolve_grid_code(grid_code)
        self._qu_deadband = qu_deadband

        super().add_OPF(**kwargs)

        self.model.name = "ACOPF"

        # run get_all_opf from flexibility instances
        for flex in self.flexibilities:
            flex.get_all_acopf(self.model)

        # voltage limits DONE: make non time dependent
        # Over Bpd only: auxiliary ppc buses have no pandapower row and so
        # no user-supplied limits; their voltage follows from the equations.
        self.model.Vmax = Param(
            self.model.Bpd,
            within=NonNegativeReals,
            initialize=self.v_limits[0][self.model.Bpd],
            mutable=True,
        )  # max voltage (p.u.)
        self.model.Vmin = Param(
            self.model.Bpd,
            within=NonNegativeReals,
            initialize=self.v_limits[1][self.model.Bpd],
            mutable=True,
        )  # min voltage (p.u.)

        # --- line power limits ---
        def line_lim_from_def(model, l, t):
            return (
                model.pLfrom[l, t] ** 2 + model.qLfrom[l, t] ** 2
                <= model.SLmax[l] ** 2 * model.v[model.A[l, 1], t] ** 2
            )

        def line_lim_to_def(model, l, t):
            return (
                model.pLto[l, t] ** 2 + model.qLto[l, t] ** 2
                <= model.SLmax[l] ** 2 * model.v[model.A[l, 2], t] ** 2
            )

        self.model.line_lim_from = Constraint(
            self.model.L, self.model.T, rule=line_lim_from_def
        )
        self.model.line_lim_to = Constraint(
            self.model.L, self.model.T, rule=line_lim_to_def
        )

        # --- power flow limits on transformer lines--- DONE non time dependent
        def transf_lim1_def(model, l, t):
            return (
                model.pThv[l, t] ** 2 + model.qThv[l, t] ** 2
                <= model.SLmaxT[l] ** 2 * model.v[model.AT[l, 1], t] ** 2
            )

        def transf_lim2_def(model, l, t):
            return (
                model.pTlv[l, t] ** 2 + model.qTlv[l, t] ** 2
                <= model.SLmaxT[l] ** 2 * model.v[model.AT[l, 2], t] ** 2
            )

        self.model.transf_lim1 = Constraint(
            self.model.TRANSF, self.model.T, rule=transf_lim1_def
        )
        self.model.transf_lim2 = Constraint(
            self.model.TRANSF, self.model.T, rule=transf_lim2_def
        )

        # voltage bounds are time-dependent
        def v_bounds(model, b, t):
            return model.Vmin[b], model.v[b, t], model.Vmax[b]

        self.model.v_constraint = Constraint(
            self.model.Bpd, self.model.T, rule=v_bounds
        )

    def add_voltage_deviation_objective(self):
        """Set objective to minimise sum of squared bus voltage deviations
        from 1 p.u. over all time steps."""
        self.model.vm = Param(
            self.model.B, initialize=self.bus_data["v_m"][self.model.B]
        )

        def voltage_deviation_objective(model, t):
            return sum(
                (model.v[b, t] - 1.0) ** 2
                for b in model.B - model.b0
                for t in model.T
            ) + sum(
                (model.v[b, t] - model.v_b0[b]) ** 2
                for b in model.b0
                for t in model.T
            )

        self.model.obj_v_deviation = Objective(
            rule=voltage_deviation_objective, sense=minimize
        )

    def add_minimize_power_objective(self):
        """Set objective to minimise total demand served over all loads and
        time steps."""

        def power_minimization_objective(model):
            return sum(model.pD[d, t] for d in model.D for t in model.T)

        self.model.Objective = Objective(
            rule=power_minimization_objective, sense=minimize
        )

    def add_generation_objective(self):
        """Set objective to minimise sum of squared generator real power
        injections."""

        def minimize_generation(model):
            return sum(model.pG[(g, t)] ** 2 for g in model.G for t in model.T)

        self.model.obj = Objective(rule=minimize_generation, sense=minimize)

    def add_weighted_generation_objective(self):
        """Set objective to minimise a weighted sum of external grid and sgen
        power."""

        def weighted_generation_objective(model):
            c1 = 4
            c3 = 1
            return c1 * sum(
                model.pG[(g, t)] for g in model.G for t in model.T
            ) + c3 * sum(model.psG[(g, t)] for g in model.sG for t in model.T)

        self.model.obj = Objective(
            rule=weighted_generation_objective, sense=minimize
        )
