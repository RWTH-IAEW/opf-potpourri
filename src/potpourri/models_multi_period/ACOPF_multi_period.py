# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Multi-period AC OPF: AC power flow plus operating limits."""

import pyomo.environ as pyo
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
    """Multi-period AC OPF over a time horizon.

    Multi-period AC OPF model combining AC power flow and OPF constraints
    over a time horizon.
    """

    def __init__(self, net, toT, fromT=None, pf=1):
        super().__init__(net, toT, fromT, pf)

    def _calc_opf_parameters(self, **kwargs):
        """Extend OPF parameter calculation with AC-specific limits.

        Voltage bounds, Q limits, Q-curve data.

        Args:
            **kwargs: Forwarded to
                :meth:`OPF_multi_period._calc_opf_parameters`, which rejects
                names it does not recognise. Options consumed by
                :meth:`add_OPF` itself (``thermal_limit``, ``free_slack_vm``,
                ``angle_limits``) never reach here.
        """
        super()._calc_opf_parameters(**kwargs)

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
        """Apply per-generator voltage limits over the bus defaults."""
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

    def add_OPF(
        self,
        thermal_limit: str = "current",
        free_slack_vm: bool = True,
        angle_limits: bool = False,
        grid_code=None,
        qu_deadband=None,
        **kwargs,
    ):
        """Add voltage bounds, thermal limits and reactive constraints.

        Extend OPF.add_OPF() with voltage bounds, AC thermal limits, and
        reactive power constraints.

        Args:
            thermal_limit: ``"current"`` enforces ``|S|² ≤ SLmax² · v²``
                (current-limit form, physically meaningful for distribution
                conductors). ``"mva"`` enforces ``|S|² ≤ SLmax²``
                (constant-MVA limit, matches MATPOWER / PGLib-OPF). Defaults
                to ``"current"``, as on the single-period model.
            free_slack_vm: When ``True`` (default), the slack-bus voltage
                magnitude floats within ``[Vmin, Vmax]`` at every time step;
                the reference angle stays fixed. ``AC_multi_period`` pins the
                slack magnitude to its base-case value while building the
                power flow, so leaving this ``False`` reproduces that legacy
                AC-PF behaviour — and gives materially different voltages
                from a single-period AC OPF of the same snapshot, which
                defaults to a free slack.
            angle_limits: When ``True``, enforce branch
                phase-angle-difference constraints
                ``angmin ≤ δ_from − δ_to ≤ angmax`` at every time step, read
                from ``net.line.angmin_degree`` / ``net.line.angmax_degree``
                and the transformer equivalent. Defaults to ``False``.
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

        if thermal_limit not in ("current", "mva"):
            raise ValueError(
                f"thermal_limit must be 'current' or 'mva', got "
                f"{thermal_limit!r}"
            )
        self.thermal_limit_mode = thermal_limit

        super().add_OPF(**kwargs)

        self.model.name = "ACOPF"

        # run get_all_opf from flexibility instances
        for flex in self.flexibilities:
            flex.get_all_acopf(self.model)

        # voltage limits DONE: make non time dependent
        # Over Bpd only: auxiliary ppc buses have no pandapower row and so
        # no user-supplied limits; their voltage follows from the equations.
        self.model.Vmax = pyo.Param(
            self.model.Bpd,
            within=pyo.NonNegativeReals,
            initialize=self.v_limits[0][self.model.Bpd],
            mutable=True,
        )  # max voltage (p.u.)
        self.model.Vmin = pyo.Param(
            self.model.Bpd,
            within=pyo.NonNegativeReals,
            initialize=self.v_limits[1][self.model.Bpd],
            mutable=True,
        )  # min voltage (p.u.)

        # --- line and transformer apparent-power limits ---
        # Same two modes as the single-period model, so a snapshot and a
        # horizon of the same network enforce the same limit.
        if thermal_limit == "current":
            # |S|² ≤ SLmax² · v² (i.e. |I| ≤ I_max). Varies with voltage;
            # physically meaningful for a thermal current rating.
            @self.model.Constraint(self.model.L, self.model.T)
            def line_lim_from(model, l, t):
                r"""Current-based thermal limit at the from end of line `l`.

                $p^2 + q^2 \le S_{max}^2 v^2$; see `ACOPF_base` for why the two
                limit forms exist.

                Args:
                    model: The Pyomo model being extended.
                    l: Branch index.
                    t: Time index from `model.T`.

                Returns:
                    A Pyomo expression.
                """
                return (
                    model.pLfrom[l, t] ** 2 + model.qLfrom[l, t] ** 2
                    <= model.SLmax[l] ** 2 * model.v[model.A[l, 1], t] ** 2
                )

            @self.model.Constraint(self.model.L, self.model.T)
            def line_lim_to(model, l, t):
                """Current-based thermal limit at the to end of line `l`.

                Args:
                    model: The Pyomo model being extended.
                    l: Branch index.
                    t: Time index from `model.T`.

                Returns:
                    A Pyomo expression.
                """
                return (
                    model.pLto[l, t] ** 2 + model.qLto[l, t] ** 2
                    <= model.SLmax[l] ** 2 * model.v[model.A[l, 2], t] ** 2
                )

            @self.model.Constraint(self.model.TRANSF, self.model.T)
            def transf_lim1(model, l, t):
                """Current-based thermal limit at the HV side of trafo `l`.

                Args:
                    model: The Pyomo model being extended.
                    l: Branch index.
                    t: Time index from `model.T`.

                Returns:
                    A Pyomo expression.
                """
                return (
                    model.pThv[l, t] ** 2 + model.qThv[l, t] ** 2
                    <= model.SLmaxT[l] ** 2 * model.v[model.AT[l, 1], t] ** 2
                )

            @self.model.Constraint(self.model.TRANSF, self.model.T)
            def transf_lim2(model, l, t):
                """Current-based thermal limit at the LV side of trafo `l`.

                Args:
                    model: The Pyomo model being extended.
                    l: Branch index.
                    t: Time index from `model.T`.

                Returns:
                    A Pyomo expression.
                """
                return (
                    model.pTlv[l, t] ** 2 + model.qTlv[l, t] ** 2
                    <= model.SLmaxT[l] ** 2 * model.v[model.AT[l, 2], t] ** 2
                )
        else:
            # |S|² ≤ SLmax² (constant-MVA limit, matches MATPOWER /
            # PowerModels' constraint_thermal_limit_* and PGLib-OPF rate_a).
            @self.model.Constraint(self.model.L, self.model.T)
            def line_lim_from(model, l, t):
                r"""Constant-MVA thermal limit at the from end of line `l`.

                $p^2 + q^2 \le S_{max}^2$, the MATPOWER/PGLib convention.

                Args:
                    model: The Pyomo model being extended.
                    l: Branch index.
                    t: Time index from `model.T`.

                Returns:
                    A Pyomo expression.
                """
                return (
                    model.pLfrom[l, t] ** 2 + model.qLfrom[l, t] ** 2
                    <= model.SLmax[l] ** 2
                )

            @self.model.Constraint(self.model.L, self.model.T)
            def line_lim_to(model, l, t):
                """Constant-MVA thermal limit at the to end of line `l`.

                Args:
                    model: The Pyomo model being extended.
                    l: Branch index.
                    t: Time index from `model.T`.

                Returns:
                    A Pyomo expression.
                """
                return (
                    model.pLto[l, t] ** 2 + model.qLto[l, t] ** 2
                    <= model.SLmax[l] ** 2
                )

            @self.model.Constraint(self.model.TRANSF, self.model.T)
            def transf_lim1(model, l, t):
                """Constant-MVA thermal limit at the HV side of trafo `l`.

                Args:
                    model: The Pyomo model being extended.
                    l: Branch index.
                    t: Time index from `model.T`.

                Returns:
                    A Pyomo expression.
                """
                return (
                    model.pThv[l, t] ** 2 + model.qThv[l, t] ** 2
                    <= model.SLmaxT[l] ** 2
                )

            @self.model.Constraint(self.model.TRANSF, self.model.T)
            def transf_lim2(model, l, t):
                """Constant-MVA thermal limit at the LV side of trafo `l`.

                Args:
                    model: The Pyomo model being extended.
                    l: Branch index.
                    t: Time index from `model.T`.

                Returns:
                    A Pyomo expression.
                """
                return (
                    model.pTlv[l, t] ** 2 + model.qTlv[l, t] ** 2
                    <= model.SLmaxT[l] ** 2
                )

        # --- slack voltage magnitude ---
        # AC_multi_period pins v[b0, t] to the base-case magnitude while
        # building the power flow. For a true AC OPF it should float within
        # [Vmin, Vmax] with only the reference angle pinned, which is what the
        # single-period model does by default.
        if free_slack_vm:
            for b0 in self.model.b0:
                for t in self.model.T:
                    self.model.v[b0, t].unfix()

        # voltage bounds are time-dependent
        @self.model.Constraint(self.model.Bpd, self.model.T)
        def v_constraint(model, b, t):
            """Bound the voltage magnitude at bus `b`, time `t`.

            Indexed over `model.Bpd`: auxiliary ppc buses carry no user voltage
            limits.

            Args:
                model: The Pyomo model being extended.
                b: Bus index from `model.B` (a ppc bus number).
                t: Time index from `model.T`.

            Returns:
                A Pyomo expression.
            """
            return model.Vmin[b], model.v[b, t], model.Vmax[b]

        # --- optional branch angle-difference limits ---
        if angle_limits:
            self._add_branch_angle_limits()

    def add_voltage_deviation_objective(self):
        """Minimise the summed squared voltage deviation over time.

        Set objective to minimise sum of squared bus voltage deviations
        from 1 p.u. over all time steps.
        """
        self.model.vm = pyo.Param(
            self.model.B, initialize=self.bus_data["v_m"][self.model.B]
        )

        @self.model.Objective(sense=pyo.minimize)
        def obj_v_deviation(model, t):
            """Summed squared voltage deviation over buses and time.

            Args:
                model: The Pyomo model being extended.
                t: Time index from `model.T`.

            Returns:
                A Pyomo expression.
            """
            return sum(
                (model.v[b, t] - 1.0) ** 2
                for b in model.B - model.b0
                for t in model.T
            ) + sum(
                (model.v[b, t] - model.v_b0[b]) ** 2
                for b in model.b0
                for t in model.T
            )

    def add_minimize_power_objective(self):
        """Minimise the total demand served over all loads and steps.

        Set objective to minimise total demand served over all loads and
        time steps.
        """

        @self.model.Objective(sense=pyo.minimize)
        def Objective(model):
            """Total demand served over all loads and time steps.

            Args:
                model: The Pyomo model being extended.

            Returns:
                A Pyomo expression.
            """
            return sum(model.pD[d, t] for d in model.D for t in model.T)

    def add_generation_objective(self):
        """Minimise the summed squared generator active power.

        Set objective to minimise sum of squared generator real power
        injections.
        """

        @self.model.Objective(sense=pyo.minimize)
        def obj(model):
            """Summed squared generator active power over time.

            Args:
                model: The Pyomo model being extended.

            Returns:
                A Pyomo expression.
            """
            return sum(model.pG[(g, t)] ** 2 for g in model.G for t in model.T)

    def add_weighted_generation_objective(self):
        """Minimise a weighted sum of external-grid and sgen power."""

        @self.model.Objective(sense=pyo.minimize)
        def obj(model):
            """Weighted sum of external-grid and static-generator power.

            The weights trade importing from the grid against running local
            generation; see the method that builds this objective for what they
            mean.

            Args:
                model: The Pyomo model being extended.

            Returns:
                A Pyomo expression.
            """
            c1 = 4
            c3 = 1
            return c1 * sum(
                model.pG[(g, t)] for g in model.G for t in model.T
            ) + c3 * sum(model.psG[(g, t)] for g in model.sG for t in model.T)
