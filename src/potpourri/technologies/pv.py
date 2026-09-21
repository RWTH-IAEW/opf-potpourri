# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""PV mix-in: attaches PV generation sets, parameters, variables, and
power-bound constraints to a multi-period model."""

import pyomo.environ as pyo
from potpourri.technologies.flexibility import Flexibility_multi_period
from potpourri.technologies.q_control import (
    DEFAULT_P_RANGE_PU,
    bus_voltage_range,
    check_var_q,
    compute_q_curves,
    resolve_grid_code,
)


class PV_multi_period(Flexibility_multi_period):
    """Multi-period PV device module with scenario-based penetration levels.

    PV units are placed randomly at a fraction of non-slack buses.  Their
    available generation follows a time-varying upper bound taken from the
    SimBench renewables profile; curtailment is modelled implicitly (the
    optimiser may dispatch below the available potential).

    Penetration percentages are based on:
    *Anforderungen an aktuelle Verteilnetze und deren zukuenftige
    Versorgungsaufgabe*.

    Args:
        net: pandapower network with simbench profiles.
        T: Number of time steps (must match the model's time horizon).
        scenario: Predefined penetration scenario 0–3.  Ignored when
            *penetration* is supplied.

            =========  =================
            scenario   % of non-slack buses
            =========  =================
            0          13.4 %
            1          22.4 %
            2          24.4 %
            3          25.4 %
            =========  =================

        penetration: Percentage of non-slack buses that receive a PV unit
            (0–100).  Overrides *scenario* when given.
        profile_column: Column name from ``net.pv_load_profiles`` (the
            SimBench renewables table) to use as the PV generation profile.
            Defaults to ``"PV5"``.
        pv_pmin: Minimum PV output in per-unit (curtailment lower bound).
            Defaults to ``0.0`` (full curtailment allowed).
        seed: Seed for the placement draw.  Defaults to
            :data:`~potpourri.technologies.flexibility.DEFAULT_PLACEMENT_SEED`,
            so repeated runs equip the same buses.  Vary it to sample
            placements.
        rng: An existing :class:`numpy.random.Generator`, taking precedence
            over *seed*.  Thread one through a Monte Carlo sweep for
            independent scenarios from a run that replays exactly.

    Example::

        pv = PV_multi_period(net, T=96, scenario=0, profile_column="PV3")
        pv.get_all(model)
    """

    SCENARIO_PENETRATION: dict[int, float] = {
        0: 13.4,
        1: 22.4,
        2: 24.4,
        3: 25.4,
    }

    def __init__(
        self,
        net,
        T=None,
        scenario=None,
        *,
        penetration: float | None = None,
        profile_column: str = "PV5",
        pv_pmin: float = 0.0,
        q_control: str | None = None,
        var_q: int = 0,
        p_inst_mw: float | None = None,
        grid_code=None,
        seed=None,
        rng=None,
    ):
        """
        Args:
            q_control: Reactive-power control mode.  One of:

                * ``None``   — no Q-control (default)
                * ``"qp"``   — Q(P) characteristic only
                * ``"qu"``   — Q(U) droop only (requires AC model with ``v``)
                * ``"both"`` — Q(P) and Q(U) combined

            var_q: VDE-AR-N 4105 operating variant (0–2).  Selects the
                Q/P envelope column from the grid-code table.
            p_inst_mw: Installed PV capacity per unit (MW).  Used as Pn in
                the Q-control characteristic.  Defaults to the peak of the
                generation profile.

        Note:
            ``pPV`` and — when ``q_control`` is set — ``qPV`` are wired into
            the nodal power balance by :meth:`couple_to_power_balance`, which
            ``get_all`` calls. Both use the generator sign convention
            (positive is injection), matching ``psG`` / ``qsG``.
        """
        super().__init__(net, T, scenario, seed=seed, rng=rng)
        self.net = net

        if penetration is not None:
            self.pv_percentage = float(penetration)
        elif scenario is not None:
            self.pv_percentage = self.SCENARIO_PENETRATION[scenario]
        else:
            raise ValueError(
                "Provide either scenario (0–3) or an explicit penetration "
                "percentage via the penetration= argument."
            )

        # SimBench stores the renewables profiles as positive generation in MW.
        # Keep that sign — pPV is a generator-convention injection like psG,
        # which is what the Q(P) / Q(U) characteristics below assume when they
        # bound qPV by `m * pPV + b * p_inst`.
        #
        # Select the column before scaling: net.pv_load_profiles is the raw
        # SimBench renewables table, which carries a non-numeric `time` column
        # alongside the generation ones, so arithmetic on the whole frame
        # raises.
        pv_profiles = self.net.pv_load_profiles
        if profile_column not in pv_profiles.columns:
            available = list(pv_profiles.columns)
            raise ValueError(
                f"profile_column '{profile_column}' not found in "
                f"net.pv_load_profiles. Available columns: {available}"
            )
        # Divided by the system base, because every power quantity in the
        # model is per-unit.
        self.pv_load_profile = pv_profiles[profile_column] / self.baseMVA

        self.pv_pmax = self.pv_load_profile
        self.pv_pmin = pv_pmin

        # Drawn from this device's own generator, so the placement is
        # reproducible (see Flexibility_multi_period.rng).
        self.random_indexes = self.draw_placement(self.pv_percentage)

        self.pv_q_control = q_control
        self.pv_var_q = int(var_q)
        if q_control is not None:
            self.grid_code = resolve_grid_code(grid_code)
            check_var_q(
                [self.pv_var_q],
                self.grid_code,
                context="PV_multi_period(var_q=...)",
            )
            self.q_limit_parameter = compute_q_curves(self.grid_code)
            if p_inst_mw is not None:
                # Supplied in MW, so convert.
                self.pv_p_inst = float(p_inst_mw) / self.baseMVA
            else:
                # pv_load_profile is already per-unit — do not divide twice.
                self.pv_p_inst = float(self.pv_load_profile.abs().max())

    def get_all(self, model):
        """Attach PV sets, parameters, variables, constraints, unfix
        variables, and couple the generation into the nodal balance."""
        self.get_sets(model)
        self.get_parameters(model)
        self.get_variables(model)
        self.get_all_constraints(model)
        self.unfix_variables(model)
        self.couple_to_power_balance(model)

    def couple_to_power_balance(self, model):
        """Register the PV generation with the nodal balance.

        ``pPV`` is a generator-convention injection (``0 <= pPV <= PV_Pmax``),
        so it enters the load-convention balance with a ``-``. ``qPV`` follows
        the same convention as ``qsG`` — positive is capacitive injection —
        and is registered on the reactive balance when Q-control is active.

        Note that ``PV_multi_period`` places units at randomly chosen buses,
        independent of ``net.sgen``. It therefore *adds* generation on top of
        whatever sgens the network already carries rather than describing them,
        which is what a penetration-scenario study wants.
        """
        self.register_kcl_real(
            model, self.bus_term(model, model.PV_bus, "pPV", sign=-1.0)
        )
        if self.pv_q_control is not None and hasattr(model, "qPV"):
            self.register_kcl_reactive(
                model, self.bus_term(model, model.PV_bus, "qPV", sign=-1.0)
            )

    def unfix_variables(self, model):
        """Unfix pPV (and qPV when Q-control is active) for all PV units."""
        for t in model.T:
            for pv in model.PV:
                model.pPV[pv, t].unfix()
                if self.pv_q_control is not None and hasattr(model, "qPV"):
                    model.qPV[pv, t].unfix()

    def get_sets(self, model):
        """Define PV and PV_bus sets from randomly placed PV units."""
        super().get_sets(model)
        model.PV = pyo.Set(initialize=list(range(len(self.random_indexes))))
        model.PV_bus = pyo.Set(initialize=list(enumerate(self.random_indexes)))
        return True

    def get_parameters(self, model):
        """Attach time-indexed PV_Pmax and PV_Pmin parameters."""
        self.PV_Pmax_dict, self.PV_Pmax_tuple = self.make_to_dict(
            model.PV, model.T, self.pv_pmax
        )
        self.PV_Pmin_dict, self.PV_Pmin_tuple = self.make_to_dict(
            model.PV, model.T, self.pv_pmin, False
        )

        model.PV_Pmax = pyo.Param(
            self.PV_Pmax_tuple, within=pyo.Reals, initialize=self.PV_Pmax_dict
        )
        model.PV_Pmin = pyo.Param(
            self.PV_Pmin_tuple, within=pyo.Reals, initialize=self.PV_Pmin_dict
        )

    def get_variables(self, model):
        """Create pPV variable and, when Q-control is active, qPV."""
        self.pPV_data_dict, self.pPV_tuple = self.make_to_dict(
            model.PV, model.T, self.pv_load_profile
        )
        model.pPV = pyo.Var(
            self.pPV_tuple, within=pyo.Reals, initialize=self.pPV_data_dict
        )
        if self.pv_q_control is not None:
            self.qPV_tuple = [(pv, t) for pv in model.PV for t in model.T]
            model.qPV = pyo.Var(
                self.qPV_tuple, within=pyo.Reals, initialize=0.0
            )

    def get_all_constraints(self, model):
        """Add real-power bound constraints for all PV units over all time
        steps."""

        @model.Constraint(model.PV, model.T)
        def PV_real_power_bounds(model, pv, t):
            # (lower, body, upper). The arguments used to be the other way
            # round, which read as Pmax being the floor.
            return model.PV_Pmin[pv, t], model.pPV[pv, t], model.PV_Pmax[pv, t]

    def get_all_acopf(self, model):
        """Add Q(P) and/or Q(U) constraints for PV units when q_control is set.

        Constraint names: ``PV_QP_pos``, ``PV_QP_neg`` (Q(P)) and
        ``PV_QU_min``, ``PV_QU_max`` (Q(U)).  Q(U) constraints require
        ``model.v[bus, t]`` (AC model); they are skipped automatically when
        the model has no voltage variable.
        """
        if self.pv_q_control is None:
            return

        # Each grid-code bound is a piecewise-linear envelope, so it becomes
        # one inequality per affine piece: the upper bound is the pointwise
        # minimum of its pieces, the lower bound the pointwise maximum.  A
        # single line cannot express the saturation shelf.
        v = self.pv_var_q
        p_inst = self.pv_p_inst
        pq_area = self.grid_code.pq_area
        qv_area = self.grid_code.qv_area
        v_span = bus_voltage_range(self.net)
        # list() first: dict() on a scalar Pyomo Set yields {None: <the set>},
        # so every lookup below returned None and every Q(U) constraint was
        # silently skipped. The buses also need mapping into ppc numbering,
        # which is what model.v is indexed over.
        pv_bus_lookup = {
            pv: int(self.bus_lookup[int(pd_bus)])
            for pv, pd_bus in list(model.PV_bus)
        }

        if self.pv_q_control in ("qp", "both"):
            pq_hi = pq_area.upper_pieces(v, DEFAULT_P_RANGE_PU)
            pq_lo = pq_area.lower_pieces(v, DEFAULT_P_RANGE_PU)

            @model.Constraint(model.PV, model.T, range(len(pq_hi)))
            def PV_QP_pos(model, pv, t, k):
                m, b = pq_hi[k]
                return model.qPV[pv, t] <= m * model.pPV[pv, t] + b * p_inst

            @model.Constraint(model.PV, model.T, range(len(pq_lo)))
            def PV_QP_neg(model, pv, t, k):
                m, b = pq_lo[k]
                return model.qPV[pv, t] >= m * model.pPV[pv, t] + b * p_inst

        if self.pv_q_control in ("qu", "both") and hasattr(model, "v"):
            qv_hi = qv_area.upper_pieces(v, v_span)
            qv_lo = qv_area.lower_pieces(v, v_span)

            @model.Constraint(model.PV, model.T, range(len(qv_lo)))
            def PV_QU_min(model, pv, t, k):
                b_bus = pv_bus_lookup.get(pv)
                if b_bus is None:
                    return pyo.Constraint.Skip
                m, b = qv_lo[k]
                return model.qPV[pv, t] >= (m * model.v[b_bus, t] + b) * p_inst

            @model.Constraint(model.PV, model.T, range(len(qv_hi)))
            def PV_QU_max(model, pv, t, k):
                b_bus = pv_bus_lookup.get(pv)
                if b_bus is None:
                    return pyo.Constraint.Skip
                m, b = qv_hi[k]
                return model.qPV[pv, t] <= (m * model.v[b_bus, t] + b) * p_inst

    def get_all_ac(self, model):
        """No additional AC components needed for PV."""

    def get_all_opf(self, model):
        """No additional OPF components needed for PV."""
