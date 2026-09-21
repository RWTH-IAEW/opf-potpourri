# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Battery mix-in.

Attaches battery storage sets, parameters, variables, and constraints to a
multi-period model.
"""

import math

import pyomo.environ as pyo
from potpourri.technologies.flexibility import Flexibility_multi_period
from potpourri.technologies.q_control import (
    DEFAULT_P_RANGE_PU,
    bus_voltage_range,
    check_var_q,
    resolve_grid_code,
)


class Battery_multi_period(Flexibility_multi_period):
    """Multi-period battery storage device module.

    Randomly places batteries at a fraction of non-slack buses and attaches
    the corresponding Pyomo Sets, Parameters, Variables, and Constraints to
    an existing multi-period model.

    Args:
        net: pandapower network with simbench profiles (already passed to
            ``ACOPF_multi_period``).
        T: Number of time steps (must match the model's time horizon).
        scenario: Predefined penetration scenario 0–3.  Ignored when
            *penetration* is supplied explicitly.

            =========  =================
            scenario   % of non-slack buses
            =========  =================
            0          1.0 %
            1          7.9 %
            2          9.9 %
            3          10.6 %
            =========  =================

        penetration: Percentage of non-slack buses that receive a battery
            (0–100).  Overrides *scenario* when given.
        power_pu: Symmetric charge/discharge power limit in per-unit on the
            system base ``net.sn_mva``.
        soc_max: Maximum state of charge (p.u. of capacity, 0–1).
        soc_min: Minimum state of charge (p.u. of capacity, 0–1).
        capacity_pu_h: Battery energy capacity in per-unit power · hours
            (i.e. ``capacity_MWh / net.sn_mva``).
        efficiency: One-way charge/discharge efficiency (0–1).  Applied as
            ``η`` on charging and ``1/η`` on discharging, so a full cycle
            loses energy and the round-trip efficiency is ``η²``.
        initial_soc_fraction: Initial SOC expressed as a fraction of
            *soc_max* (0–1).  Default 0.5 = 50 % of maximum capacity.
        terminal_soc: State of charge required at the last time step.
            ``"cyclic"`` (default) closes the energy balance over the
            horizon by requiring the final SOC to equal the initial SOC;
            a float pins it to that absolute value; ``None`` leaves it
            free.  Without a terminal condition the optimiser can drain
            the battery and export the stored energy at no cost, which
            biases any objective that values energy or grid exchange.
        s_inv_pu: Apparent-power rating of the battery converter, per-unit on
            ``net.sn_mva``.  Bounds P and Q jointly through the S² circle
            ``BAT_P² + BAT_Q² ≤ s_inv_pu²``.  Defaults to *power_pu*, i.e. a
            converter sized exactly to the active-power limit, which leaves no
            reactive headroom at full charge or discharge.  Oversize it
            (``s_inv_pu = 1.1 * power_pu``) to allow reactive support at rated
            active power.
        q_control: Grid-code reactive-power control mode, as for sgens:

            * ``None`` — no capability area, Q bounded only by the S² circle
              and any ``cos_phi_min`` (default)
            * ``"qp"`` — Q(P) capability area
            * ``"qu"`` — Q(U) capability area (needs an AC model)
            * ``"both"`` — both areas

            The Q(P) area is keyed on the **discharging** leg, since that is
            the mode in which a storage unit acts as a generating unit under
            VDE-AR-N 4105/4110/4120.
        var_q: Grid-code operating variant (0–2), selecting the Q/P envelope
            column, as for sgens.
        cos_phi_min: Power-factor floor.  Adds
            ``|BAT_Q| ≤ tan(arccos(cos_phi_min)) · (BAT_Pchg + BAT_Pdis)``.
            The leg sum is a convex stand-in for ``|BAT_P|``; the two coincide
            whenever only one leg is active, which any loss-making efficiency
            already makes optimal.
        grid_code: Technical connection rule supplying the capability areas,
            as accepted by
            :func:`~potpourri.technologies.q_control.resolve_grid_code`.
        seed: Seed for the placement draw.  Defaults to
            :data:`~potpourri.technologies.flexibility.DEFAULT_PLACEMENT_SEED`,
            so repeated runs equip the same buses.  Vary it to sample
            placements.
        rng: An existing :class:`numpy.random.Generator`, taking precedence
            over *seed*.  Thread one through a Monte Carlo sweep for
            independent scenarios from a run that replays exactly.

    Reactive power is available on AC-style models only.  ``BAT_Q`` and its
    constraints are skipped on the DC formulation, which carries no reactive
    balance.

    Example::

        battery = Battery_multi_period(
            net, T=96, scenario=1
        )
        # or with explicit parameters:
        battery = Battery_multi_period(
            net, T=96,
            penetration=15.0,   # 15 % of buses
            power_pu=0.01,
            capacity_pu_h=0.025,
            efficiency=0.95,
        )
        battery.get_all(model)
    """

    # Default penetration levels per scenario (percentage of non-slack buses)
    SCENARIO_PENETRATION: dict[int, float] = {0: 1.0, 1: 7.9, 2: 9.9, 3: 10.6}

    def __init__(
        self,
        net,
        T=None,
        scenario=None,
        *,
        penetration: float | None = None,
        power_pu: float = 0.006,
        soc_max: float = 1.0,
        soc_min: float = 0.2,
        capacity_pu_h: float = 0.015,
        efficiency: float = 0.9,
        initial_soc_fraction: float = 0.5,
        terminal_soc: float | str | None = "cyclic",
        s_inv_pu: float | None = None,
        q_control: str | None = None,
        var_q: int = 0,
        cos_phi_min: float | None = None,
        grid_code=None,
        seed=None,
        rng=None,
    ):
        super().__init__(net, T, scenario, seed=seed, rng=rng)

        if not 0.0 < efficiency <= 1.0:
            raise ValueError(
                f"efficiency must be in (0, 1]; got {efficiency}. It is "
                "applied as a one-way efficiency, so the round-trip "
                "efficiency is efficiency**2."
            )
        if terminal_soc is not None and terminal_soc != "cyclic":
            terminal_soc = float(terminal_soc)
            if not soc_min <= terminal_soc <= soc_max:
                raise ValueError(
                    f"terminal_soc={terminal_soc} lies outside "
                    f"[soc_min, soc_max] = [{soc_min}, {soc_max}], so the "
                    "model would be infeasible."
                )

        if penetration is not None:
            self.bat_percentage = float(penetration)
        elif scenario is not None:
            self.bat_percentage = self.SCENARIO_PENETRATION[scenario]
        else:
            raise ValueError(
                "Provide either scenario (0–3) or an explicit penetration "
                "percentage via the penetration= argument."
            )

        # Select buses for battery placement, from this device's own generator
        # so the draw is reproducible (see Flexibility_multi_period.rng).
        self.random_indexes = self.draw_placement(self.bat_percentage)

        self.bat_power = power_pu
        self.bat_soc_max = soc_max
        self.bat_soc_min = soc_min
        self.bat_cap = capacity_pu_h
        self.bat_efficiency = efficiency
        self.bat_initial_soc_fraction = initial_soc_fraction
        self.bat_terminal_soc = terminal_soc

        # --- converter apparent-power rating and reactive capability ---
        self.bat_s_inv = float(power_pu if s_inv_pu is None else s_inv_pu)
        if self.bat_s_inv < power_pu:
            raise ValueError(
                f"s_inv_pu={self.bat_s_inv} is below power_pu={power_pu}: the "
                "converter could not deliver the battery's own active-power "
                "limit, so the S**2 circle would cap charging and discharging "
                "instead of the power rating."
            )
        if cos_phi_min is not None and not 0.0 < cos_phi_min <= 1.0:
            raise ValueError(
                f"cos_phi_min must be in (0, 1]; got {cos_phi_min}."
            )
        self.bat_cos_phi_min = cos_phi_min
        self.bat_tan_phi = (
            None if cos_phi_min is None else math.tan(math.acos(cos_phi_min))
        )

        if q_control is not None and q_control not in ("qp", "qu", "both"):
            raise ValueError(
                f"q_control must be None, 'qp', 'qu' or 'both'; got "
                f"{q_control!r}."
            )
        self.bat_q_control = q_control
        self.bat_var_q = int(var_q)
        if q_control is not None:
            self.grid_code = resolve_grid_code(grid_code)
            check_var_q(
                [self.bat_var_q],
                self.grid_code,
                context="Battery_multi_period(var_q=...)",
            )

    def get_all(self, model):
        """Attach the battery to a model and couple it to the balance.

        Attach battery sets, parameters, variables, and constraints to the
        model, and couple its power into the nodal balance.
        """
        self.get_sets(model)
        self.get_parameters(model)
        self.get_variables(model)
        self.get_all_constraints(model)
        self.couple_to_power_balance(model)

    def couple_to_power_balance(self, model):
        """Register the battery's net power with the nodal balance.

        ``BAT_P = BAT_Pchg - BAT_Pdis`` is already in the load sign convention
        — positive while charging, i.e. a load on the grid — so it enters the
        balance with a ``+`` and needs no sign flip. The balance is rebuilt by
        ``add_OPF``, so this works whether the battery is attached before or
        after the OPF constraints.

        ``BAT_Q`` is a generator-convention injection like ``qsG`` — positive
        is capacitive — so it enters the load-convention balance with a ``-``.
        It exists only on model kinds that have a reactive balance.
        """
        self.register_kcl_real(
            model, self.bus_term(model, model.BAT_bus, "BAT_P", sign=1.0)
        )
        if self._has_reactive(model):
            self.register_kcl_reactive(
                model, self.bus_term(model, model.BAT_bus, "BAT_Q", sign=-1.0)
            )

    def get_sets(self, model):
        """Define BAT and BAT_bus sets from randomly placed batteries."""
        super().get_sets(model)
        model.BAT = pyo.Set(initialize=list(range(len(self.random_indexes))))
        model.BAT_bus = pyo.Set(
            initialize=list(enumerate(self.random_indexes))
        )
        return True

    def get_parameters(self, model):
        """Attach the battery power, SOC, capacity and efficiency data."""
        model.BAT_Pmax = pyo.Param(
            model.BAT, within=pyo.Reals, initialize=self.bat_power
        )
        model.BAT_Pmin = pyo.Param(
            model.BAT, within=pyo.Reals, initialize=-self.bat_power
        )
        model.BAT_SOCmax = pyo.Param(
            model.BAT, within=pyo.Reals, initialize=self.bat_soc_max
        )
        model.BAT_SOCmin = pyo.Param(
            model.BAT, within=pyo.Reals, initialize=self.bat_soc_min
        )
        model.BAT_Cap = pyo.Param(
            model.BAT, within=pyo.Reals, initialize=self.bat_cap
        )
        model.BAT_Eff = pyo.Param(
            model.BAT, within=pyo.Reals, initialize=self.bat_efficiency
        )
        model.BAT_SOC_init = pyo.Param(
            model.BAT,
            within=pyo.Reals,
            initialize={
                b: self.bat_initial_soc_fraction * self.bat_soc_max
                for b in range(len(self.random_indexes))
            },
        )
        model.BAT_Sinv = pyo.Param(
            model.BAT, within=pyo.NonNegativeReals, initialize=self.bat_s_inv
        )
        return True

    def _has_reactive(self, model):
        """Whether this model kind carries a reactive power balance.

        The DC formulation has a single real-power balance, so a reactive
        battery variable there would be unconstrained and meaningless.
        """
        return hasattr(model, "KCL_reactive")

    def get_variables(self, model):
        """Create the charge/discharge power and state-of-charge variables.

        Charging and discharging are separate non-negative variables so the
        one-way efficiency can be applied in the right direction on each
        (``η`` charging, ``1/η`` discharging). A single signed power variable
        cannot express that: the same factor would scale both directions and
        a full cycle would return the SOC exactly to its starting value, i.e.
        no round-trip loss for any value of ``efficiency``.

        ``BAT_P`` remains available as the net injection expression
        ``BAT_Pchg - BAT_Pdis`` (positive = charging, i.e. a load on the
        grid), matching ``pSTOR`` in the single-period storage block.

        All three start on the **idle** trajectory: no power either way and the
        SOC held at its initial value. That point satisfies every battery
        constraint, including the terminal condition, so the solver begins from
        a feasible battery and can only improve on the no-battery solution.
        Left at Pyomo's default the SOC would start at 0, below ``soc_min`` and
        away from ``BAT_SOC_init``, i.e. infeasible before the solve begins —
        which cost the optimiser the no-battery solution as a fallback and let
        it settle on a *worse* objective than the same model without a battery.
        """
        idle_soc = self.bat_initial_soc_fraction * self.bat_soc_max
        model.BAT_Pchg = pyo.Var(
            model.BAT, model.T, within=pyo.NonNegativeReals, initialize=0.0
        )
        model.BAT_Pdis = pyo.Var(
            model.BAT, model.T, within=pyo.NonNegativeReals, initialize=0.0
        )
        model.BAT_SOC = pyo.Var(
            model.BAT, model.T, within=pyo.Reals, initialize=idle_soc
        )

        def bat_injection_rule(model, b, t):
            """Net active power of battery `b` at time `t`.

            $P_{chg} - P_{dis}$, in the **load** convention: positive while
            charging. An `Expression`, not a variable, so it carries no
            constraint of its own.

            Args:
                model: The multi-period model being extended.
                b: Battery index from `model.BAT`.
                t: Time index from `model.T`.

            Returns:
                A Pyomo expression in p.u.
            """
            return model.BAT_Pchg[b, t] - model.BAT_Pdis[b, t]

        model.BAT_P = pyo.Expression(
            model.BAT, model.T, rule=bat_injection_rule
        )

        # Reactive power, generator convention (positive = capacitive
        # injection), matching qsG. AC-style models only.
        #
        # The box [-S_inv, +S_inv] is implied by the S**2 circle and so adds no
        # modelling restriction, but stating it explicitly keeps IPOPT's
        # iterates inside the capability region. Without it the solver can
        # wander far outside on the nonconvex AC problem and report a locally
        # infeasible point on models that are feasible.
        if self._has_reactive(model):
            s_inv = self.bat_s_inv
            model.BAT_Q = pyo.Var(
                model.BAT,
                model.T,
                within=pyo.Reals,
                bounds=(-s_inv, s_inv),
                initialize=0.0,
            )
        return True

    def get_all_constraints(self, model):
        """Add power-bound, SOC-bound, SOC-update and terminal constraints."""

        def bat_chg_limit_rule(model, b, t):
            r"""Cap the charging power of battery `b` at time `t`.

            Args:
                model: The multi-period model being extended.
                b: Battery index from `model.BAT`.
                t: Time index from `model.T`.

            Returns:
                A Pyomo inequality, $P_{chg} \le P_{max}$ (p.u.).
            """
            return model.BAT_Pchg[b, t] <= model.BAT_Pmax[b]

        def bat_dis_limit_rule(model, b, t):
            r"""Cap the discharging power of battery `b` at time `t`.

            Args:
                model: The multi-period model being extended.
                b: Battery index from `model.BAT`.
                t: Time index from `model.T`.

            Returns:
                A Pyomo inequality, $P_{dis} \le P_{max}$ (p.u.).
            """
            # BAT_Pmin is the (negative) discharge limit on the signed
            # injection, so its magnitude bounds the discharge leg.
            return model.BAT_Pdis[b, t] <= -model.BAT_Pmin[b]

        model.bat_chg_limit_con = pyo.Constraint(
            model.BAT, model.T, rule=bat_chg_limit_rule
        )
        model.bat_dis_limit_con = pyo.Constraint(
            model.BAT, model.T, rule=bat_dis_limit_rule
        )

        def bat_power_rule(model, b, t):
            r"""Discourage charging and discharging `b` at the same time.

            $P_{chg} + P_{dis} \le P_{max}$: the convex relaxation of the
            charge-or-discharge disjunction. It does not forbid doing both --
            that needs a binary and a MINLP -- so a solver can still split the
            rating if the objective happens to reward the resulting losses.

            Args:
                model: The multi-period model being extended.
                b: Battery index from `model.BAT`.
                t: Time index from `model.T`.

            Returns:
                A Pyomo inequality expression.
            """
            # Converter throughput limit. Also the convex relaxation of "not
            # both legs at once": without it the optimiser can charge and
            # discharge simultaneously to burn energy through the efficiency
            # losses, which is not a physical dispatch.
            return (
                model.BAT_Pchg[b, t] + model.BAT_Pdis[b, t]
                <= model.BAT_Pmax[b]
            )

        model.bat_power_con = pyo.Constraint(
            model.BAT, model.T, rule=bat_power_rule
        )

        def bat_soc_rule(model, b, t):
            """Pin the initial state of charge, bound it afterwards.

            At the first step the SOC is fixed to `BAT_SOC_init`, which is what
            anchors the whole trajectory; at every later step it is only
            required to stay inside $[SOC_{min}, SOC_{max}]$.

            Args:
                model: The multi-period model being extended.
                b: Battery index from `model.BAT`.
                t: Time index from `model.T`.

            Returns:
                An equality at the first time step, otherwise the ranged
                3-tuple `(SOCmin, SOC, SOCmax)`. SOC is a fraction in $[0, 1]$.
            """
            if t == model.T.at(1):
                return model.BAT_SOC[b, t] == model.BAT_SOC_init[b]
            return (
                model.BAT_SOCmin[b],
                model.BAT_SOC[b, t],
                model.BAT_SOCmax[b],
            )

        model.bat_soc_con = pyo.Constraint(
            model.BAT, model.T, rule=bat_soc_rule
        )

        def bat_soc_update_rule(model, b, t):
            r"""Carry the state of charge from one step to the next.

            $$SOC_t = SOC_{t-1} + \frac{\Delta t\,
            (\eta P_{chg,t} - P_{dis,t} / \eta)}{E_{max}}$$

            This is the constraint that couples the time steps, and the reason
            a battery cannot be modelled one snapshot at a time. Note $\eta$ is
            the **one-way** efficiency: a full cycle returns $\eta^2$.

            Units: SOC is a fraction, powers are p.u., $\Delta t$ is in hours,
            so the quotient is dimensionless.

            Args:
                model: The multi-period model being extended.
                b: Battery index from `model.BAT`.
                t: Time index from `model.T`.

            Returns:
                A Pyomo equality expression, or `Constraint.Skip` at the first
                step, which has no predecessor and is pinned by `bat_soc_rule`
                instead.
            """
            if t == model.T.at(1):
                return pyo.Constraint.Skip
            # η on the way in, 1/η on the way out: a charge/discharge cycle
            # of equal grid-side energy loses (1 - η²) of it.
            return (
                model.BAT_SOC[b, t]
                == model.BAT_SOC[b, t - 1]
                + model.deltaT
                * (
                    model.BAT_Eff[b] * model.BAT_Pchg[b, t]
                    - model.BAT_Pdis[b, t] / model.BAT_Eff[b]
                )
                / model.BAT_Cap[b]
            )

        model.bat_soc_update_con = pyo.Constraint(
            model.BAT, model.T, rule=bat_soc_update_rule
        )

        if self.bat_terminal_soc is not None:

            def bat_terminal_soc_rule(model, b):
                """Impose the terminal state of charge of battery `b`.

                Added only when `bat_terminal_soc` was given. `"cyclic"` closes
                the horizon by requiring the final SOC to equal the initial
                one, so the schedule can be repeated; a float pins it to that
                value.

                Without this the optimiser will happily empty the battery by
                the last step, since stored energy has no value at the
                horizon's end.

                Args:
                    model: The multi-period model being extended.
                    b: Battery index from `model.BAT`.

                Returns:
                    A Pyomo equality expression on the last time step.
                """
                last = model.T.last()
                if self.bat_terminal_soc == "cyclic":
                    return model.BAT_SOC[b, last] == model.BAT_SOC_init[b]
                return model.BAT_SOC[b, last] == self.bat_terminal_soc

            model.bat_terminal_soc_con = pyo.Constraint(
                model.BAT, rule=bat_terminal_soc_rule
            )

        if self._has_reactive(model):
            self._add_reactive_constraints(model)
        return True

    def _add_reactive_constraints(self, model):
        """Bound the battery's reactive power by its converter capability.

        Adds, in order of increasing specificity:

        * ``bat_inverter_s2`` — the S² circle ``BAT_P² + BAT_Q² ≤ BAT_Sinv²``.
          Always present. ``BAT_P`` is the affine ``BAT_Pchg − BAT_Pdis``, so
          the constraint stays convex.
        * ``bat_cos_phi_upper`` / ``bat_cos_phi_lower`` — the power-factor
          floor, when ``cos_phi_min`` was given.
        * ``bat_QP_pos`` / ``bat_QP_neg`` — the grid-code Q(P) capability area
          keyed on the discharging leg, when ``q_control`` includes ``"qp"``.
        * ``bat_QU_min`` / ``bat_QU_max`` — the grid-code Q(U) capability area,
          when ``q_control`` includes ``"qu"`` and the model has voltages.

        Each capability area is a piecewise-linear envelope, so it becomes one
        inequality per affine piece: the upper bound is the pointwise minimum
        of its pieces and the lower bound the pointwise maximum. A single line
        cannot express the saturation shelf.
        """

        @model.Constraint(model.BAT, model.T)
        def bat_inverter_s2(model, b, t):
            r"""Converter apparent-power limit of battery `b` at time `t`.

            $P^2 + Q^2 \le S_{inv}^2$, so reactive support competes with active
            throughput for the same rating.

            Args:
                model: The multi-period model being extended.
                b: Battery index from `model.BAT`.
                t: Time index from `model.T`.

            Returns:
                A Pyomo inequality expression.
            """
            return (
                model.BAT_P[b, t] ** 2 + model.BAT_Q[b, t] ** 2
                <= model.BAT_Sinv[b] ** 2
            )

        if self.bat_tan_phi is not None:
            model.BAT_tan_phi = pyo.Param(
                model.BAT, initialize=self.bat_tan_phi
            )

            @model.Constraint(model.BAT, model.T)
            def bat_cos_phi_upper(model, b, t):
                """Upper power-factor bound for battery `b` at time `t`.

                Args:
                    model: The multi-period model being extended.
                    b: Battery index from `model.BAT`.
                    t: Time index from `model.T`.

                Returns:
                    A Pyomo inequality expression.
                """
                return model.BAT_Q[b, t] <= model.BAT_tan_phi[b] * (
                    model.BAT_Pchg[b, t] + model.BAT_Pdis[b, t]
                )

            @model.Constraint(model.BAT, model.T)
            def bat_cos_phi_lower(model, b, t):
                """Lower power-factor bound for battery `b` at time `t`.

                Args:
                    model: The multi-period model being extended.
                    b: Battery index from `model.BAT`.
                    t: Time index from `model.T`.

                Returns:
                    A Pyomo inequality expression.
                """
                return model.BAT_Q[b, t] >= -model.BAT_tan_phi[b] * (
                    model.BAT_Pchg[b, t] + model.BAT_Pdis[b, t]
                )

        if self.bat_q_control is None:
            return

        var_q = self.bat_var_q
        pq_area = self.grid_code.pq_area
        qv_area = self.grid_code.qv_area

        if self.bat_q_control in ("qp", "both"):
            pq_hi = pq_area.upper_pieces(var_q, DEFAULT_P_RANGE_PU)
            pq_lo = pq_area.lower_pieces(var_q, DEFAULT_P_RANGE_PU)

            @model.Constraint(model.BAT, model.T, range(len(pq_hi)))
            def bat_QP_pos(model, b, t, k):
                """Upper Q(P) capability piece `k` for battery `b` at time `t`.

                One affine piece of the grid-code envelope; together the pieces
                form its pointwise minimum.

                Args:
                    model: The multi-period model being extended.
                    b: Battery index from `model.BAT`.
                    t: Time index from `model.T`.
                    k: Piece index.

                Returns:
                    A Pyomo inequality expression.
                """
                m, c = pq_hi[k]
                return (
                    model.BAT_Q[b, t]
                    <= m * model.BAT_Pdis[b, t] + c * model.BAT_Sinv[b]
                )

            @model.Constraint(model.BAT, model.T, range(len(pq_lo)))
            def bat_QP_neg(model, b, t, k):
                """Lower Q(P) capability piece `k` for battery `b` at time `t`.

                Args:
                    model: The multi-period model being extended.
                    b: Battery index from `model.BAT`.
                    t: Time index from `model.T`.
                    k: Piece index.

                Returns:
                    A Pyomo inequality expression.
                """
                m, c = pq_lo[k]
                return (
                    model.BAT_Q[b, t]
                    >= m * model.BAT_Pdis[b, t] + c * model.BAT_Sinv[b]
                )

        if self.bat_q_control in ("qu", "both") and hasattr(model, "v"):
            v_span = bus_voltage_range(self.net)
            qv_hi = qv_area.upper_pieces(var_q, v_span)
            qv_lo = qv_area.lower_pieces(var_q, v_span)
            bat_bus = {
                d: int(self.bus_lookup[int(pd_bus)])
                for d, pd_bus in list(model.BAT_bus)
            }

            @model.Constraint(model.BAT, model.T, range(len(qv_lo)))
            def bat_QU_min(model, b, t, k):
                """Lower Q(U) capability piece `k` for battery `b` at time `t`.

                Evaluated at the battery's own bus voltage, so this piece needs
                an AC model; a DC model has no `v` to read.

                Args:
                    model: The multi-period model being extended.
                    b: Battery index from `model.BAT`.
                    t: Time index from `model.T`.
                    k: Piece index.

                Returns:
                    A Pyomo inequality expression.
                """
                m, c = qv_lo[k]
                bus = bat_bus[b]
                return (
                    model.BAT_Q[b, t]
                    >= (m * model.v[bus, t] + c) * model.BAT_Sinv[b]
                )

            @model.Constraint(model.BAT, model.T, range(len(qv_hi)))
            def bat_QU_max(model, b, t, k):
                """Upper Q(U) capability piece `k` for battery `b` at time `t`.

                Args:
                    model: The multi-period model being extended.
                    b: Battery index from `model.BAT`.
                    t: Time index from `model.T`.
                    k: Piece index.

                Returns:
                    A Pyomo inequality expression.
                """
                m, c = qv_hi[k]
                bus = bat_bus[b]
                return (
                    model.BAT_Q[b, t]
                    <= (m * model.v[bus, t] + c) * model.BAT_Sinv[b]
                )

    def get_all_acopf(self, model):
        """No additional ACOPF components needed for batteries."""

    def get_all_ac(self, model):
        """No additional AC components needed for batteries."""

    def get_all_opf(self, model):
        """No additional OPF components needed for batteries."""
