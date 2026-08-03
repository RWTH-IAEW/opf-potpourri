"""Static generator (sgen) mix-in: attaches sgen profiles and OPF limits to
a multi-period model."""

import warnings

import numpy as np
import pandas as pd
import pyomo.environ as pyo
from potpourri.technologies.flexibility import Flexibility_multi_period
from potpourri.technologies.q_control import (
    DEFAULT_GRID_CODE,
    DEFAULT_P_RANGE_PU,
    attach_deadband_qu,
    bus_voltage_range,
    check_var_q,
    compute_q_curves,
    resolve_grid_code,
    resolve_qu_curve,
    warn_if_curve_leaves_pq_area,
)


class SgenMinPAboveProfileWarning(UserWarning):
    """Raised when ``net.sgen.min_p_mw`` exceeds the available generation.

    ``sPGmax`` follows the profile, so a constant lower bound can sit above
    it — a PV sgen with ``min_p_mw > 0`` is infeasible at night. The solver
    reports only "infeasible", with nothing pointing at the lower bound.
    """


class Sgens_multi_period(Flexibility_multi_period):
    """Multi-period static generator device module; reads sgen profiles from
    net.profiles."""

    def __init__(self, net, T=None, scenario=None):
        super().__init__(net, T, scenario)

        # SimBench provides q_mvar profiles for sgens; if absent,
        # Basemodel_multi_period.calc_reactive_sgen_power() derives them
        # from the active-power profile and the power-factor argument.
        #
        # Both are divided by the system base, because every other quantity
        # in the model is per-unit: the demand profiles
        # (Flexibility_multi_period),
        # the generator limits (Generator_multi_period), the shunt data, the
        # line and transformer ratings, and the sgen capability data derived
        # from net.sgen further down this class. pyo_to_net_multi_period also
        # multiplies psG by baseMVA on the way out. Taking these two in MW made
        # the power balance mix scales on any network with sn_mva != 1.
        sgen_bus = self.bus_lookup[self.net.sgen.bus.values]
        self.static_generation_data = {
            "p": self.net.profiles[("sgen", "p_mw")] / self.baseMVA,
            "q": self.net.profiles[("sgen", "q_mvar")] / self.baseMVA,
            "in_service": self.net.sgen.in_service.values,
            "bus": sgen_bus,
        }
        self.static_generation_data["gen_bus"] = list(
            enumerate(self.static_generation_data["bus"])
        )

    def get_all(self, model):
        """Attach sets, parameters, and variables; fix all sgen variables to
        profile values."""
        self.get_sets(model)
        self.get_parameters(model)
        self.get_variables(model)
        self.fix_variables(model)

    def get_all_opf(self, model):
        """Attach OPF sets, parameters, and real-power bound constraints for
        controllable sgens."""
        self.get_opf_sets(model)
        self.get_opf_parameters(model)
        self.get_all_Constraints_opf(model)

    def get_all_acopf(self, model):
        """Attach AC-specific OPF parameters (reactive-power limits) and the
        corresponding range constraints. Separated from ``get_all_opf`` so
        that the DC OPF path can ignore Q entirely.

        Also adds the following optional constraint blocks when the matching
        data-population method was called beforehand:

        * Q(P) / Q(U) grid-code constraints (``var_q`` column)
        * Inverter S² circle (``sn_mva`` column)
        * P(U) active-power curtailment (``pu_curtail`` column)
        * Fixed cos(φ) equality (``fixed_cos_phi`` column)
        * cos(φ)(P) profile equality (``cos_phi_p_profile`` column)
        """
        self.get_acopf_parameters(model)
        self.get_all_Constraints_acopf(model)
        if getattr(self, "sgen_qc_indices", []):
            self._add_sgen_q_ctrl_mp(model)
        if getattr(self, "sgen_inv_indices", []):
            self._add_sgen_inverter_s2_mp(model)
        if getattr(self, "sgen_pu_indices", []):
            self._add_sgen_pu_curtail_mp(model)
        if getattr(self, "sgen_fcf_indices", []):
            self._add_sgen_fixed_cos_phi_mp(model)
        if getattr(self, "sgen_cpp_indices", []):
            self._add_sgen_cpp_mp(model)

    def get_sets(self, model):
        super().get_sets(model)
        # list, no set, because list is ordered data source, set is not
        self.sgens_in_service_list = np.where(
            self.static_generation_data["in_service"]
        )[0].tolist()
        model.sG = pyo.Set(
            initialize=self.sgens_in_service_list
        )  # static generators
        model.sGbs = pyo.Set(
            within=model.sG * model.B,
            initialize=self.static_generation_data["gen_bus"],
        )  # set of static generator-bus mapping
        return True

    def get_opf_sets(self, model):
        # list, no set, because list is ordered data source, set is not
        self.sgens_controllable_list = np.where(
            self.static_generation_data["controllable"]
        )[0].tolist()

        model.sGc = pyo.Set(
            within=model.sG,
            initialize=[
                g
                for g in self.sgens_controllable_list
                if g in self.sgens_in_service_list
            ],
        )  # static generators that are controllable and in service

    def get_parameters(self, model):
        self.PsG_data_dict, self.PsG_tuple = self.make_to_dict(
            model.sG, model.T, self.static_generation_data["p"]
        )
        self.QsG_data_dict, self.QsG_tuple = self.make_to_dict(
            model.sG, model.T, self.static_generation_data["q"]
        )
        # --- Parameters ---
        model.PsG = pyo.Param(self.PsG_tuple, initialize=self.PsG_data_dict)
        # reactive generation
        model.QsG = pyo.Param(self.QsG_tuple, initialize=self.QsG_data_dict)
        return True

    def get_opf_parameters(self, model):
        """Attach the OPF parameters that are common to AC and DC: the
        controllable-sgen real-power bounds. The reactive-power bounds
        (`QsGmax` / `QsGmin`) are AC-specific and live in
        :meth:`get_acopf_parameters`; calling ``get_opf_parameters`` from the
        DC OPF path therefore no longer requires `QsGmax_tuple` /
        `QsGmin_tuple` (which are populated only by
        :meth:`static_generation_reactive_power_limits`, an AC step).
        """
        # static generation real power limits
        model.sPGmax = pyo.Param(
            self.PsGmax_tuple, initialize=self.PsGmax_data_dict
        )
        model.sPGmin = pyo.Param(
            self.PsGmin_tuple, initialize=self.PsGmin_data_dict
        )

    def get_acopf_parameters(self, model):
        """Attach AC-only OPF parameters (reactive-power bounds)."""
        model.QsGmax = pyo.Param(
            self.QsGmax_tuple,
            within=pyo.Reals,
            initialize=self.QsGmax_data_dict,
            mutable=True,
        )
        model.QsGmin = pyo.Param(
            self.QsGmin_tuple,
            within=pyo.Reals,
            initialize=self.QsGmin_data_dict,
            mutable=True,
        )

    def static_generation_q_ctrl_data(
        self, net, grid_code=None, qu_deadband=None
    ):
        """Compute Q(P)/Q(U) characteristic data for sgens with ``var_q`` set.

        Populates ``self.q_limit_parameter``, ``self.sgen_var_q``,
        ``self.sgen_p_inst``, and ``self.sgen_qc_indices``.  Call this from
        :meth:`ACOPF_multi_period._calc_opf_parameters` before
        :meth:`get_all_acopf`.

        Args:
            net: pandapower network.  Only processed when ``net.sgen`` has a
                ``var_q`` column.
            grid_code: Technical connection rule supplying the capability
                envelope, as accepted by
                :func:`~potpourri.technologies.q_control.resolve_grid_code`.
                Defaults to VDE-AR-N 4120.
            qu_deadband: Optional dead-band Q(U) characteristic replacing the
                Q(U) capability area, as accepted by
                :func:`~potpourri.technologies.q_control.resolve_qu_curve`.
        """
        if "var_q" not in net.sgen:
            self.sgen_qc_indices = []
            return

        code = resolve_grid_code(grid_code)
        self.grid_code = code
        self.qu_curve = resolve_qu_curve(qu_deadband, code)
        self.q_limit_parameter = compute_q_curves(code)
        self.sgen_var_q = (
            net.sgen.var_q.values
        )  # object array; may contain None

        if "p_inst_mw" in net.sgen:
            p_inst = (
                net.sgen.p_inst_mw.fillna(net.sgen.p_mw).values / self.baseMVA
            )
        else:
            p_inst = net.sgen.p_mw.values / self.baseMVA
        self.sgen_p_inst = p_inst

        self.sgen_qc_indices = [
            g
            for g in self.sgens_in_service_list
            if pd.notna(self.sgen_var_q[g])
        ]
        check_var_q(
            [self.sgen_var_q[g] for g in self.sgen_qc_indices],
            code,
            context="net.sgen.var_q",
        )
        self._reactive_bounds_from_grid_code(code)

    def _add_sgen_q_ctrl_mp(self, model):
        """Add time-indexed Q(P) and Q(U) constraints for controllable sgens.

        Requires :meth:`static_generation_q_ctrl_data` to have been called.
        Creates ``model.sGqc`` and four constraint blocks per (sgen, time):
        ``sG_QP_pos``, ``sG_QP_neg``, ``sG_QU_min``, ``sG_QU_max``.

        Q(U) constraints require a voltage variable ``model.v[bus, t]``;
        they are skipped silently when the model has no such variable (e.g.
        in a linearised or DC formulation).
        """
        qc_list = [g for g in self.sgen_qc_indices if g in set(model.sGc)]
        if not qc_list:
            return

        model.sGqc = pyo.Set(within=model.sGc, initialize=qc_list)
        sGbs_lookup = {g: b for (g, b) in model.sGbs}
        var_q = {g: int(self.sgen_var_q[g]) for g in qc_list}
        p_inst = {g: float(self.sgen_p_inst[g]) for g in qc_list}

        # Each grid-code bound is a piecewise-linear envelope, so it becomes
        # one inequality per affine piece: the upper bound is the pointwise
        # minimum of its pieces, the lower bound the pointwise maximum.  A
        # single line cannot express the saturation shelf.
        pq_area = self.grid_code.pq_area
        qv_area = self.grid_code.qv_area
        v_span = bus_voltage_range(self.net)
        model.sG_QP_PIECE = pyo.RangeSet(
            0, pq_area.max_pieces(DEFAULT_P_RANGE_PU) - 1
        )
        model.sG_QU_PIECE = pyo.RangeSet(0, qv_area.max_pieces(v_span) - 1)

        @model.Constraint(model.sGqc, model.T, model.sG_QP_PIECE)
        def sG_QP_pos(model, g, t, k):
            pieces = pq_area.upper_pieces(var_q[g], DEFAULT_P_RANGE_PU)
            if k >= len(pieces):
                return pyo.Constraint.Skip
            m, b = pieces[k]
            return model.qsG[g, t] <= m * model.psG[g, t] + b * p_inst[g]

        @model.Constraint(model.sGqc, model.T, model.sG_QP_PIECE)
        def sG_QP_neg(model, g, t, k):
            pieces = pq_area.lower_pieces(var_q[g], DEFAULT_P_RANGE_PU)
            if k >= len(pieces):
                return pyo.Constraint.Skip
            m, b = pieces[k]
            return model.qsG[g, t] >= m * model.psG[g, t] + b * p_inst[g]

        if not hasattr(model, "v"):
            return

        # A dead band cannot be expressed as a capability area: the feasible
        # set pinches to Q = 0 around nominal voltage, which is not convex.
        # When one is requested, Q is pinned to the characteristic instead.
        qu_curve = getattr(self, "qu_curve", None)
        if qu_curve is not None:
            # The characteristic assigns Q while the Q(P) area bounds it;
            # where they disagree the model is infeasible with nothing in
            # the solver output pointing here.
            warn_if_curve_leaves_pq_area(
                qu_curve,
                pq_area,
                v_range=v_span,
                context=f"{self.grid_code.title} Q(U) dead band",
            )
            keys = [
                (g, t) for g in qc_list for t in model.T if g in sGbs_lookup
            ]
            v_lo, v_hi = qv_area.exact_range()
            attach_deadband_qu(
                model,
                "sG_qu_db",
                keys,
                q_of=lambda k: model.qsG[k[0], k[1]],
                v_of=lambda k: model.v[sGbs_lookup[k[0]], k[1]],
                pn_of=lambda k: p_inst[k[0]],
                variant_of=lambda k: var_q[k[0]],
                curve=qu_curve,
                v_bounds=(v_lo, v_hi),
            )
            return

        @model.Constraint(model.sGqc, model.T, model.sG_QU_PIECE)
        def sG_QU_min(model, g, t, k):
            if g not in sGbs_lookup:
                return pyo.Constraint.Skip
            pieces = qv_area.lower_pieces(var_q[g], v_span)
            if k >= len(pieces):
                return pyo.Constraint.Skip
            b_bus = sGbs_lookup[g]
            m, b = pieces[k]
            return model.qsG[g, t] >= (m * model.v[b_bus, t] + b) * p_inst[g]

        @model.Constraint(model.sGqc, model.T, model.sG_QU_PIECE)
        def sG_QU_max(model, g, t, k):
            if g not in sGbs_lookup:
                return pyo.Constraint.Skip
            pieces = qv_area.upper_pieces(var_q[g], v_span)
            if k >= len(pieces):
                return pyo.Constraint.Skip
            b_bus = sGbs_lookup[g]
            m, b = pieces[k]
            return model.qsG[g, t] <= (m * model.v[b_bus, t] + b) * p_inst[g]

    def static_generation_inverter_data(self, net):
        """Compute inverter apparent-power ratings for the S² constraint.

        Populates ``self.sgen_s_inv`` (per-unit rating array) and
        ``self.sgen_inv_indices`` (in-service sgen indices with finite
        ``sn_mva``).  Call from
        :meth:`ACOPF_multi_period._calc_opf_parameters` before
        :meth:`get_all_acopf`.

        The apparent-power rating is
        ``S_inv = sn_mva * converter_sizing_pu / baseMVA``.
        ``converter_sizing_pu`` defaults to 1.0 when the column is absent.

        Args:
            net: pandapower network with ``net.sgen.sn_mva`` present.
        """
        if "sn_mva" not in net.sgen:
            self.sgen_inv_indices = []
            return

        conv_sz = (
            net.sgen["converter_sizing_pu"].fillna(1.0)
            if "converter_sizing_pu" in net.sgen
            else pd.Series(1.0, index=net.sgen.index)
        )
        self.sgen_s_inv = (net.sgen["sn_mva"] * conv_sz).fillna(
            0.0
        ).values / self.baseMVA
        self.sgen_inv_indices = [
            g
            for g in self.sgens_in_service_list
            if pd.notna(net.sgen.at[g, "sn_mva"])
            and float(net.sgen.at[g, "sn_mva"]) > 0
        ]

        # cos(φ) cone: tan_phi per sgen if net.sgen["cos_phi_min"] present
        if "cos_phi_min" in net.sgen:
            self.sgen_tan_phi = {
                g: float(
                    np.tan(np.arccos(float(net.sgen.at[g, "cos_phi_min"])))
                )
                for g in self.sgen_inv_indices
                if pd.notna(net.sgen.at[g, "cos_phi_min"])
                and 0 < float(net.sgen.at[g, "cos_phi_min"]) <= 1
            }
        else:
            self.sgen_tan_phi = {}

    def _add_sgen_inverter_s2_mp(self, model):
        """Add time-indexed inverter S² circle and optional cos(φ) cone.

        Requires :meth:`static_generation_inverter_data` to have been called.

        Creates:

        * ``model.sGinv`` — controllable sgens with finite ``sn_mva``
        * ``model.S_inv`` — per-unit apparent-power rating
        * ``model.sgen_inverter_s2`` — ``psG² + qsG² ≤ S_inv²``
        * ``model.sGpf``, ``model.tan_phi``, ``model.sgen_cos_phi_upper/lower``
          — cos(φ) cone constraints, present only when ``cos_phi_min`` data
          was found in ``net.sgen["cos_phi_min"]``.
        """
        inv_list = [g for g in self.sgen_inv_indices if g in set(model.sGc)]
        if not inv_list:
            return

        model.sGinv = pyo.Set(within=model.sGc, initialize=inv_list)
        model.S_inv = pyo.Param(
            model.sGinv,
            initialize={g: float(self.sgen_s_inv[g]) for g in inv_list},
        )

        @model.Constraint(model.sGinv, model.T)
        def sgen_inverter_s2(model, g, t):
            return (
                model.psG[g, t] ** 2 + model.qsG[g, t] ** 2
                <= model.S_inv[g] ** 2
            )

        # cos(φ) cone — only when cos_phi_min data is available
        pf_list = [g for g in inv_list if g in self.sgen_tan_phi]
        if pf_list:
            model.sGpf = pyo.Set(within=model.sGinv, initialize=pf_list)
            model.tan_phi = pyo.Param(
                model.sGpf,
                initialize={g: self.sgen_tan_phi[g] for g in pf_list},
            )

            @model.Constraint(model.sGpf, model.T)
            def sgen_cos_phi_upper(model, g, t):
                return model.qsG[g, t] <= model.tan_phi[g] * model.psG[g, t]

            @model.Constraint(model.sGpf, model.T)
            def sgen_cos_phi_lower(model, g, t):
                return model.qsG[g, t] >= -model.tan_phi[g] * model.psG[g, t]

    def get_all_Constraints_opf(self, model):
        # psG Constraint
        @model.Constraint(model.sGc, model.T)
        def static_generation_real_power_bounds(model, g, t):
            model.psG[(g, t)].unfix()
            return (
                model.sPGmin[(g, t)],
                model.psG[(g, t)],
                model.sPGmax[(g, t)],
            )

    def get_variables(self, model):
        # --- Variables ---
        model.psG = pyo.Var(
            self.PsG_tuple, domain=pyo.NonNegativeReals
        )  # real static generator power
        model.qsG = pyo.Var(
            self.QsG_tuple, domain=pyo.Reals
        )  # reactive power of static generators
        return True

    def unfix_variables(self, model):
        # unfix the static generation values
        for g in model.sG:
            for t in model.T:
                model.psG[(g, t)].unfix()
        return True

    def fix_variables(self, model):
        # fix the static generation values
        for g in model.sG:
            for t in model.T:
                model.psG[(g, t)].fix(model.PsG[(g, t)])

        for g in model.sG:
            for t in model.T:
                model.qsG[(g, t)].fix(model.QsG[(g, t)])

    #     return True

    def static_generation_real_power_limits(self, model):
        if "controllable" in self.net.sgen:
            self.static_generation_data["controllable"] = (
                self.net.sgen.controllable.values
            )
        else:
            self.static_generation_data["controllable"] = np.full(
                len(self.net.sgen), False
            )

        self.PsGmax_data_dict, self.PsGmax_tuple = self.make_to_dict(
            model.sG, model.T, self.static_generation_data["p"]
        )
        self.static_generation_data["min_p"] = self._read_min_p()
        self.PsGmin_data_dict, self.PsGmin_tuple = self.make_to_dict(
            model.sG,
            model.T,
            self.static_generation_data["min_p"],
            False,
        )
        self._warn_if_min_p_above_profile(model)

    def _read_min_p(self):
        """Read the sgen real-power lower bound from ``net.sgen.min_p_mw``.

        Mirrors ``OPF.static_generation_real_power_limits``: a missing column
        or a NaN entry falls back to 0, the distribution-grid convention that
        PV and wind can be curtailed to zero but cannot reverse. Before this
        was read, the bound was hard-coded to 0 for every sgen and time step,
        so ``min_p_mw`` was silently ignored and a single-period study could
        not be reproduced over a horizon.

        Divided by ``baseMVA`` like every other power quantity in the model,
        including ``sPGmax``.
        """
        n_sgen = len(self.net.sgen.index)
        if "min_p_mw" not in self.net.sgen:
            return np.zeros(n_sgen)
        return (
            self.net.sgen.min_p_mw.astype(float).fillna(0.0).values
            / self.baseMVA
        )

    def _warn_if_min_p_above_profile(self, model):
        """Flag a lower bound the profile cannot satisfy.

        ``sPGmin`` is constant over the horizon while ``sPGmax`` follows the
        profile, so the two can cross — most obviously for PV at night. That
        makes the model infeasible with nothing in the solver output naming
        the cause, so say it here instead.

        Only controllable, in-service sgens are checked: those are the ones
        ``get_all_Constraints_opf`` bounds. Everything else has ``psG`` fixed
        to its profile value, so its ``sPGmin`` is never enforced and a
        stray ``min_p_mw`` on it is harmless.
        """
        controllable = np.where(self.static_generation_data["controllable"])[
            0
        ].tolist()
        in_service = getattr(self, "sgens_in_service_list", None)
        bounded = {
            g for g in controllable if in_service is None or g in in_service
        }
        if not bounded:
            return

        conflicts = [
            (g, t)
            for (g, t), lo in self.PsGmin_data_dict.items()
            if g in bounded and lo > self.PsGmax_data_dict[(g, t)]
        ]
        if not conflicts:
            return
        g, t = conflicts[0]
        lo = self.PsGmin_data_dict[(g, t)]
        hi = self.PsGmax_data_dict[(g, t)]
        warnings.warn(
            f"net.sgen.min_p_mw exceeds the available generation for "
            f"{len(conflicts)} (sgen, time step) pairs. First one: sgen {g} "
            f"at t={t} has min_p={lo:.6g} p.u. but the profile offers only "
            f"{hi:.6g} p.u. The model is infeasible there, and the solver "
            f"will report nothing more specific than 'infeasible'. Set "
            f"min_p_mw to 0 (or NaN) for profile-driven sgens, or set "
            f"controllable=False to pin them to the profile.",
            SgenMinPAboveProfileWarning,
            stacklevel=3,
        )

    def static_generation_reactive_power_limits(self, model):
        if "controllable" in self.net.sgen:
            self.static_generation_data["controllable"] = (
                self.net.sgen.controllable.values
            )
        else:
            self.static_generation_data["controllable"] = np.full(
                len(self.net.sgen), False
            )

        self.QsGmax_data_dict, self.QsGmax_tuple = self.make_to_dict(
            model.sG, model.T, abs(self.static_generation_data["q"])
        )
        self.QsGmin_data_dict, self.QsGmin_tuple = self.make_to_dict(
            model.sG, model.T, -abs(self.static_generation_data["q"])
        )
        # self.static_generation_wind_var_q( self.net)
        self.static_generation_data["type"] = self.net.sgen.type.values

    def _reactive_bounds_from_grid_code(self, code):
        """Give Q-controlled sgens the capability the grid code grants them.

        :meth:`static_generation_reactive_power_limits` derives ``QsGmax`` /
        ``QsGmin`` from the ``q_mvar`` profile, which SimBench ships as zero
        for PV.  Those bounds then pin ``qsG`` to zero, and every Q-control
        constraint built on top of it — Q(P), Q(U), the inverter circle —
        becomes vacuous: the model looks Q-controlled and dispatches no
        reactive power at all.  The single-period path has always overridden
        these bounds from the capability table
        (``ACOPF.static_generation_reactive_power_limits``); the
        multi-period one did not, and the call that would have done it sat
        commented out next to the profile-derived defaults.

        The override matches the single-period behaviour: for an sgen with
        ``var_q`` set, the grid code decides the reactive bounds, not the
        profile.  Set ``var_q`` to NaN on sgens that should keep their
        profile-derived limits.
        """
        table = code.vqu_q_max
        if not self.sgen_qc_indices:
            return
        if not hasattr(self, "QsGmax_data_dict"):
            # Silence is what made the original bug survive, so refuse
            # rather than no-op: without the profile-derived bounds to
            # override, the grid-code capability would never reach the
            # model and Q-control would be vacuous again.
            raise RuntimeError(
                "static_generation_reactive_power_limits() must run before "
                "static_generation_q_ctrl_data(): the grid-code reactive "
                "bounds override QsGmax / QsGmin, which do not exist yet. "
                "ACOPF_multi_period._calc_opf_parameters() calls them in "
                "that order; call them in that order too."
            )
        for g in self.sgen_qc_indices:
            variant = int(self.sgen_var_q[g])
            pn = float(self.sgen_p_inst[g])
            hi = float(table[0, variant]) * pn
            lo = float(table[1, variant]) * pn
            for key in self.QsGmax_data_dict:
                if (key[0] if isinstance(key, tuple) else key) == g:
                    self.QsGmax_data_dict[key] = hi
                    self.QsGmin_data_dict[key] = lo

    def get_all_Constraints_acopf(self, model):
        # QsG_Constraint
        @model.Constraint(model.sGc, model.T)
        def static_generation_reactive_power_bounds(model, g, t):
            model.qsG[(g, t)].unfix()
            return (
                model.QsGmin[(g, t)],
                model.qsG[(g, t)],
                model.QsGmax[(g, t)],
            )

    # ------------------------------------------------------------------
    # P(U) active-power curtailment  (VDE-AR-N 4105 §8.5)
    # ------------------------------------------------------------------

    def static_generation_pu_curtail_data(self, net):
        """Compute P(U) curtailment data for sgens with ``pu_curtail`` set.

        Populates ``self.sgen_pu_indices``, ``self.sgen_p_inst_pu_curtail``,
        ``self.sgen_v_curtail``, and ``self.sgen_v_max_curtail``.  Call from
        :meth:`ACOPF_multi_period._calc_opf_parameters`.

        Args:
            net: pandapower network.  ``net.sgen`` must have a ``pu_curtail``
                boolean column.  Per-sgen voltage thresholds are read from
                ``v_curtail_pu`` and ``v_max_curtail_pu``, defaulting to the
                thresholds of the grid code selected in
                :meth:`static_generation_q_ctrl_data`.  Installed capacity is
                read from ``p_inst_mw``, falling back to ``p_mw``.
        """
        if "pu_curtail" not in net.sgen:
            self.sgen_pu_indices = []
            return

        code = getattr(self, "grid_code", DEFAULT_GRID_CODE)

        p_inst = (
            net.sgen.p_inst_mw.fillna(net.sgen.p_mw.abs()).values
            if "p_inst_mw" in net.sgen
            else net.sgen.p_mw.abs().values
        ) / self.baseMVA

        v_curtail = (
            net.sgen.v_curtail_pu.fillna(code.vpu_v_curtail).values
            if "v_curtail_pu" in net.sgen
            else np.full(len(net.sgen), code.vpu_v_curtail)
        )
        v_max_curtail = (
            net.sgen.v_max_curtail_pu.fillna(code.vpu_v_max).values
            if "v_max_curtail_pu" in net.sgen
            else np.full(len(net.sgen), code.vpu_v_max)
        )

        self.sgen_p_inst_pu_curtail = p_inst
        self.sgen_v_curtail = v_curtail
        self.sgen_v_max_curtail = v_max_curtail
        self.sgen_pu_indices = [
            g
            for g in self.sgens_in_service_list
            if pd.notna(net.sgen.at[g, "pu_curtail"])
            and bool(net.sgen.at[g, "pu_curtail"])
        ]

    def _add_sgen_pu_curtail_mp(self, model):
        """Add time-indexed P(U) curtailment constraints.

        Requires :meth:`static_generation_pu_curtail_data`.

        Adds ``model.sGpu``, ``model.P_inst_pu``, ``model.V_curtail``,
        ``model.V_max_curtail``, and the bilinear constraint
        ``model.sgen_pu_curtail``:

            psG[g,t] * (V_max - V_curtail) ≤ P_inst[g] * (V_max - v[bus,t])

        Requires an AC model with ``model.v`` voltage variable; silently
        skipped otherwise.
        """
        if not hasattr(model, "v"):
            return
        pu_list = [g for g in self.sgen_pu_indices if g in set(model.sGc)]
        if not pu_list:
            return

        model.sGpu = pyo.Set(within=model.sGc, initialize=pu_list)
        model.P_inst_pu = pyo.Param(
            model.sGpu,
            initialize={
                g: float(self.sgen_p_inst_pu_curtail[g]) for g in pu_list
            },
        )
        model.V_curtail = pyo.Param(
            model.sGpu,
            initialize={g: float(self.sgen_v_curtail[g]) for g in pu_list},
        )
        model.V_max_curtail = pyo.Param(
            model.sGpu,
            initialize={g: float(self.sgen_v_max_curtail[g]) for g in pu_list},
        )
        sGbs_lookup = {g: b for (g, b) in model.sGbs}

        @model.Constraint(model.sGpu, model.T)
        def sgen_pu_curtail(model, g, t):
            if g not in sGbs_lookup:
                return pyo.Constraint.Skip
            b = sGbs_lookup[g]
            dv = model.V_max_curtail[g] - model.V_curtail[g]
            return model.psG[g, t] * dv <= model.P_inst_pu[g] * (
                model.V_max_curtail[g] - model.v[b, t]
            )

    # ------------------------------------------------------------------
    # Fixed cos(φ) mode
    # ------------------------------------------------------------------

    def static_generation_fixed_cos_phi_data(self, net):
        """Compute fixed-cos(φ) data for sgens with ``fixed_cos_phi`` set.

        Populates ``self.sgen_fcf_indices`` and ``self.sgen_fcf_tan_phi``.
        Call from :meth:`ACOPF_multi_period._calc_opf_parameters`.

        Args:
            net: pandapower network.  ``net.sgen`` must have a
                ``fixed_cos_phi`` column containing per-sgen power factors
                (0 < cos_phi ≤ 1); NaN or missing rows are skipped.
        """
        if "fixed_cos_phi" not in net.sgen:
            self.sgen_fcf_indices = []
            return

        self.sgen_fcf_tan_phi = {
            g: float(np.tan(np.arccos(float(net.sgen.at[g, "fixed_cos_phi"]))))
            for g in self.sgens_in_service_list
            if pd.notna(net.sgen.at[g, "fixed_cos_phi"])
            and 0 < float(net.sgen.at[g, "fixed_cos_phi"]) <= 1
        }
        self.sgen_fcf_indices = list(self.sgen_fcf_tan_phi.keys())

    def _add_sgen_fixed_cos_phi_mp(self, model):
        """Add time-indexed fixed-cos(φ) equality constraints.

        Requires :meth:`static_generation_fixed_cos_phi_data`.

        Adds ``model.sGfcf``, ``model.fixed_tan_phi``, and the equality
        ``model.sgen_fixed_cos_phi``:

            qsG[g, t] == fixed_tan_phi[g] * psG[g, t]
        """
        fcf_list = [g for g in self.sgen_fcf_indices if g in set(model.sGc)]
        if not fcf_list:
            return

        model.sGfcf = pyo.Set(within=model.sGc, initialize=fcf_list)
        model.fixed_tan_phi = pyo.Param(
            model.sGfcf,
            initialize={g: self.sgen_fcf_tan_phi[g] for g in fcf_list},
        )

        @model.Constraint(model.sGfcf, model.T)
        def sgen_fixed_cos_phi(model, g, t):
            return model.qsG[g, t] == model.fixed_tan_phi[g] * model.psG[g, t]

    # ------------------------------------------------------------------
    # cos(φ)(P) profile  (VDE-AR-N 4105 piecewise P-Q curve)
    # ------------------------------------------------------------------

    def static_generation_cpp_data(self, net):
        """Compute cos(φ)(P) profile data for sgens with ``cos_phi_p_profile``.

        Populates ``self.sgen_cpp_indices``, ``self.sgen_cpp_tan_phi``,
        ``self.sgen_cpp_pn``, and ``self.sgen_cpp_p_thresh``.  Call from
        :meth:`ACOPF_multi_period._calc_opf_parameters`.

        The cos(φ)(P) characteristic defines Q = 0 for P ≤ P_thresh and
        Q = P · tan_phi · (P − P_thresh) / (Pn − P_thresh) for P > P_thresh.
        This is implemented as a quadratic equality in the NLP.

        Args:
            net: pandapower network.  ``net.sgen`` must have a truthy
                ``cos_phi_p_profile`` column.  ``cos_phi_min`` sets the
                power factor at full output; ``p_inst_mw`` gives Pn;
                ``cpp_p_threshold_pu`` (optional) sets P_thresh / Pn,
                defaulting to the threshold of the grid code selected in
                :meth:`static_generation_q_ctrl_data`.
        """
        if "cos_phi_p_profile" not in net.sgen:
            self.sgen_cpp_indices = []
            return

        code = getattr(self, "grid_code", DEFAULT_GRID_CODE)

        p_inst = (
            net.sgen.p_inst_mw.fillna(net.sgen.p_mw.abs()).values
            if "p_inst_mw" in net.sgen
            else net.sgen.p_mw.abs().values
        ) / self.baseMVA

        cpp_thresh_pu = (
            net.sgen.cpp_p_threshold_pu.fillna(code.cpp_p_threshold_pu).values
            if "cpp_p_threshold_pu" in net.sgen
            else np.full(len(net.sgen), code.cpp_p_threshold_pu)
        )

        self.sgen_cpp_indices = []
        self.sgen_cpp_tan_phi = {}
        self.sgen_cpp_pn = {}
        self.sgen_cpp_p_thresh = {}

        for g in self.sgens_in_service_list:
            if not (
                pd.notna(net.sgen.at[g, "cos_phi_p_profile"])
                and bool(net.sgen.at[g, "cos_phi_p_profile"])
            ):
                continue
            cos_phi_val = (
                net.sgen.at[g, "cos_phi_min"]
                if "cos_phi_min" in net.sgen
                and pd.notna(net.sgen.at[g, "cos_phi_min"])
                else None
            )
            if cos_phi_val is None or not (0 < float(cos_phi_val) <= 1):
                continue
            pn = float(p_inst[g])
            pt = float(cpp_thresh_pu[g]) * pn
            if pn <= pt:
                continue
            self.sgen_cpp_indices.append(g)
            self.sgen_cpp_tan_phi[g] = float(
                np.tan(np.arccos(float(cos_phi_val)))
            )
            self.sgen_cpp_pn[g] = pn
            self.sgen_cpp_p_thresh[g] = pt

    def _add_sgen_cpp_mp(self, model):
        """Add time-indexed cos(φ)(P) profile equality constraints.

        Requires :meth:`static_generation_cpp_data`.

        Adds ``model.sGcpp``, ``model.cpp_tan_phi``, ``model.cpp_Pn``,
        ``model.cpp_P_thresh``, and the quadratic equality
        ``model.sgen_cpp``:

            qsG[g,t] * (Pn - P_thresh) == cpp_tan_phi * psG[g,t]
                                            * (psG[g,t] - P_thresh)

        This is a smooth quadratic equality tractable by IPOPT.  Q → 0
        as P → 0 or P → P_thresh; Q → Pn·tan_phi at P = Pn.
        """
        cpp_list = [g for g in self.sgen_cpp_indices if g in set(model.sGc)]
        if not cpp_list:
            return

        model.sGcpp = pyo.Set(within=model.sGc, initialize=cpp_list)
        model.cpp_tan_phi = pyo.Param(
            model.sGcpp,
            initialize={g: self.sgen_cpp_tan_phi[g] for g in cpp_list},
        )
        model.cpp_Pn = pyo.Param(
            model.sGcpp,
            initialize={g: self.sgen_cpp_pn[g] for g in cpp_list},
        )
        model.cpp_P_thresh = pyo.Param(
            model.sGcpp,
            initialize={g: self.sgen_cpp_p_thresh[g] for g in cpp_list},
        )

        @model.Constraint(model.sGcpp, model.T)
        def sgen_cpp(model, g, t):
            dPn = model.cpp_Pn[g] - model.cpp_P_thresh[g]
            return model.qsG[g, t] * dPn == model.cpp_tan_phi[g] * model.psG[
                g, t
            ] * (model.psG[g, t] - model.cpp_P_thresh[g])
