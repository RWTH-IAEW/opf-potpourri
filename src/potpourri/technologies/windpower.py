# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Wind power mix-in.

Attaches wind generator sets, parameters, and Q-control constraints to a
multi-period model.
"""

import numpy as np
import pandas as pd
import pyomo.environ as pyo
from loguru import logger
from potpourri.technologies.q_control import (
    DEFAULT_GRID_CODE,
    DEFAULT_P_RANGE_PU,
    DEFAULT_WIND_SGEN_TYPES,
    bus_voltage_range,
    check_var_q,
    compute_q_curves,
    resolve_grid_code,
)
from potpourri.technologies.sgens import Sgens_multi_period


# ---------------------------------------------------------------------------
# Grid-code Q-curve parameters come from the registry in
# potpourri.technologies.q_control, not from a private copy.  This module used
# to carry its own capability table and its own reimplementation of the
# Q-curve maths, which meant `grid_code` never reached the wind or
# hosting-capacity paths and the numbers could drift from the registry
# silently.
#
# The simplified hosting-capacity check uses the widest envelope the selected
# code offers: the largest capacitive entry and the most negative inductive
# one.  For the default VDE-AR-N 4120 that is +0.484322 / -0.410775.
# ---------------------------------------------------------------------------


def _hc_q_bounds(code):
    """Return (max, min) Q/P for the simplified HC check of ``code``."""
    q = code.vqu_q_max
    return float(q[0].max()), float(q[1].min())


# Defaults for the HC keyword arguments below.  Derived from the default
# grid code so they track the registry: +0.484322 / -0.410775 for VDE-AR-N
# 4120, the widest capacitive and widest inductive variant it offers.
_DEFAULT_HC_Q_MAX, _DEFAULT_HC_Q_MIN = _hc_q_bounds(DEFAULT_GRID_CODE)


class Windpower_multi_period(Sgens_multi_period):
    """Multi-period wind, with Q-control and hosting capacity.

    Multi-period wind generator device module, extending sgen with
    Q-control and hosting-capacity (HC) support.

    The Q-control constraints implement the capability area of the selected
    grid code (see :mod:`potpourri.technologies.q_control`), defaulting to
    VDE-AR-N 4120.

    Args:
        net: pandapower network.  If ``net.bus`` contains a ``windpot_p_mw``
            column, it is used as the active-power upper bound for HC
            generators.
        T: Number of time steps.
        scenario: Unused for wind (no penetration scenario); reserved for
            interface compatibility.
        sw_max_mva: Default apparent-power upper limit per HC wind generator
            (MVA).  Can be overridden per network via
            ``_calc_wind_opf_parameters``.
        sw_min_mva: Minimum apparent power for an *active* HC wind generator
            (MVA).  A generator with ``y=1`` must carry at least this much
            apparent power.
        qp_max: Maximum Q/P ratio (capacitive, positive) for the simplified
            HC grid-code Q-P constraint.  Default: 0.484322, the widest
            capacitive variant of VDE-AR-N 4120.
        qp_min: Minimum Q/P ratio (inductive, negative) for the simplified
            HC grid-code Q-P constraint.  Default: -0.410775, the widest
            inductive variant of the same code.
    """

    def __init__(
        self,
        net,
        T=None,
        scenario=None,
        *,
        sw_max_mva: float = 10_000.0,
        sw_min_mva: float = 0.0,
        qp_max: float = _DEFAULT_HC_Q_MAX,
        qp_min: float = _DEFAULT_HC_Q_MIN,
    ):
        super().__init__(net, T, scenario)
        self._sw_max_mva = sw_max_mva
        self._sw_min_mva = sw_min_mva
        self.qp_max = qp_max
        self.qp_min = qp_min

        if "windpot_p_mw" in net.bus:
            # pWmax bounds psG, which is per-unit, so convert from MW.
            self.static_generation_data["windpot"] = (
                net.bus.windpot_p_mw[net.sgen.bus.values].values / self.baseMVA
            )
            self.static_generation_data["type"] = net.sgen.type.values

    def get_all(self, model):
        """No-op. Wind generators are initialised via get_all_opf."""

    def get_all_opf(self, model):
        """Attach OPF sets and parameters for controllable wind generators."""
        self.get_opf_sets(model)
        self.get_opf_parameters(model)

    def get_opf_sets(self, model):
        """Define the WIND_HC, WIND and WINDc sets.

        `WIND_HC` holds the hosting-capacity candidates, `WIND` adds the
        wind units already in the network, and `WINDc` narrows that to the
        ones that are controllable and declare a Q-capability variant.

        The flags are read from `net.sgen` rather than from
        `static_generation_data`, which in the multi-period sgen class is a
        dict of time-indexed arrays and has no `.index` to filter on.

        Args:
            model: The Pyomo model being extended, with `model.sG` and
                `model.sGc` already declared.

        Returns:
            True, once the sets are attached to `model`.
        """
        sgen = self.net.sgen
        # `Windpower_multi_period.get_all` is a no-op and `get_all_opf` calls
        # this directly, so the base class's `get_sets` -- and with it
        # `sgens_in_service_list` -- has not necessarily run.
        in_service = set(sgen.index[sgen.in_service.astype(bool)])

        wind_hc = sgen.get("wind_hc")
        if wind_hc is None:
            hc_candidates = []
        else:
            hc_candidates = [
                g
                for g in sgen.index[wind_hc.fillna(False).astype(bool)]
                if g in in_service
            ]
        model.WIND_HC = pyo.Set(within=model.sG, initialize=hc_candidates)

        # Match every SimBench wind spelling, not just the HV one; see
        # DEFAULT_WIND_SGEN_TYPES. Shared with the single-period path so the
        # two cannot drift apart.
        wind_types = getattr(self, "wind_sgen_types", DEFAULT_WIND_SGEN_TYPES)
        existing_wind = (
            [
                g
                for g in sgen.index[sgen["type"].isin(wind_types)]
                if g in in_service
            ]
            if "type" in sgen
            else []
        )
        model.WIND = model.WIND_HC | pyo.Set(
            within=model.sG, initialize=existing_wind
        )

        var_q = sgen.get("var_q")
        with_variant = (
            []
            if var_q is None
            else [
                g
                for g in sgen.index
                if var_q[g] is not None and not pd.isna(var_q[g])
            ]
        )
        model.WINDc = model.WIND & model.sGc & pyo.Set(initialize=with_variant)
        return True

    def _calc_wind_opf_parameters(
        self,
        model,
        sw_max_mva: float | None = None,
        sw_min_mva: float | None = None,
    ):
        """Derive the apparent-power limits and Q(U) slopes.

        Compute SWmax/SWmin apparent-power limits and Q(U) slope
        parameters for HC wind generators.

        Args:
            model: Pyomo model (requires ``model.WIND_HC`` and ``model.T``).
            sw_max_mva: Override the instance-level ``sw_max_mva`` for this
                call.  Defaults to the value passed at construction.
            sw_min_mva: Override the instance-level ``sw_min_mva`` for this
                call.  Defaults to the value passed at construction.
        """
        if sw_max_mva is None:
            sw_max_mva = self._sw_max_mva
        if sw_min_mva is None:
            sw_min_mva = self._sw_min_mva

        if "windpot_p_mw" in self.net.bus:
            # pWmax bounds psG, which is per-unit, so convert from MW.
            self.windpot = dict(
                zip(
                    self.net.sgen.index,
                    self.net.bus.windpot_p_mw[self.net.sgen.bus.values].values
                    / self.baseMVA,
                )
            )

        wind_hc_set = np.arange(len(self.net.sgen))[
            self.net.sgen.wind_hc.fillna(False).astype(bool)
            & self.net.sgen.in_service
        ]
        # Sizing bounds, so per candidate and not spread over model.T: a
        # candidate is installed once, and only its dispatch varies.
        self.SWmax_data = pd.Series(sw_max_mva / self.baseMVA, wind_hc_set)
        self.SWmin_data = pd.Series(sw_min_mva / self.baseMVA, wind_hc_set)

        # Q(U) slopes from the selected grid code's characteristic.
        # Slope from low-voltage to high-voltage point; intercepts at V3 and V1
        code = resolve_grid_code(getattr(self, "grid_code", None))
        x = code.vqu_v_points
        y = code.vqu_q_max
        hc_max, _hc_min = _hc_q_bounds(code)
        # Narrowest capacitive / widest inductive variant is the last column.
        last = y.shape[1] - 1
        self.m_qu_max = (hc_max + abs(y[1, 0])) / (x[0, 0] - x[1, 0])
        self.qu_max = -self.m_qu_max * x[1, 0] + hc_max
        self.m_qu_min = (abs(y[0, last]) + abs(y[1, last])) / (
            x[0, 0] - x[1, 0]
        )
        self.qu_min = -self.m_qu_min * x[0, 0] + y[0, last]
        return True

    def get_hc_acopf_parameters(self, model, net):
        """Attach the per-candidate sizing bounds, and pWmax if available.

        These are **sizing** bounds, so they carry no time index: a
        candidate is installed once and the same plant is there at every
        step. Only the dispatch varies over the horizon.

        Args:
            model: The Pyomo model being extended, with `model.WIND_HC`.
            net: The pandapower network; `net.bus.windpot_p_mw` caps a
                candidate's active power when the column is present.

        Returns:
            True, once the parameters are attached to `model`.
        """
        model.SWmax = pyo.Param(
            model.WIND_HC,
            initialize=self.SWmax_data.to_dict(),
            mutable=True,
        )
        model.SWmin = pyo.Param(
            model.WIND_HC,
            initialize=self.SWmin_data.to_dict(),
            mutable=True,
        )

        if "windpot_p_mw" in self.net.bus:
            model.pWmax = pyo.Param(
                model.WIND_HC,
                initialize={w: float(self.windpot[w]) for w in model.WIND_HC},
                mutable=True,
            )
        return True

    def get_hc_acopf_variables(self, model):
        """Attach the sizing variables: selection `y` and rating `SW2`.

        Both are per candidate rather than per time step. `SW2` is the
        **squared** installed apparent power, which is what makes the
        per-step limit `p² + q² ≤ SW2` a convex quadratic rather than a
        bilinear one; the installed rating itself is `sqrt(SW2)`.

        Args:
            model: The Pyomo model being extended, with `model.WIND_HC`.

        Returns:
            True, once the variables are attached to `model`.
        """
        model.y = pyo.Var(model.WIND_HC, domain=pyo.Binary, initialize=1.0)
        model.SW2 = pyo.Var(
            model.WIND_HC,
            domain=pyo.NonNegativeReals,
            bounds=lambda m, w: (0.0, float(pyo.value(m.SWmax[w])) ** 2),
            initialize=lambda m, w: float(pyo.value(m.SWmax[w])) ** 2,
        )
        return True

    def get_opf_parameters(self, model):
        """Attach the Q-variant and installed-capacity parameters.

        Both are constant over the horizon -- the grid-code variant a unit
        declares and the capacity it has installed do not change between
        time steps -- but they are declared over `model.T` as well, because
        the capability rules are evaluated per `(w, t)`.

        Args:
            model: The Pyomo model being extended, with `model.WINDc`.

        Returns:
            True, once the parameters are attached to `model`.
        """
        model.var_q = pyo.Param(
            model.WINDc,
            model.T,
            initialize={
                (g, t): int(self.net.sgen.var_q[g])
                for g in model.WINDc
                for t in model.T
            },
        )
        model.PsG_inst = pyo.Param(
            model.WINDc,
            model.T,
            initialize={
                (g, t): float(self._installed_p(g))
                for g in model.WINDc
                for t in model.T
            },
        )
        return True

    def _installed_p(self, g):
        """Installed active power of sgen `g`, in per unit.

        `p_inst_mw` is the nameplate rating where the network carries it;
        otherwise the present `p_mw` is the best available stand-in.

        Args:
            g: Static-generator index.

        Returns:
            The installed active power, divided by `baseMVA`.
        """
        sgen = self.net.sgen
        if "p_inst_mw" in sgen and not pd.isna(sgen.p_inst_mw[g]):
            return float(sgen.p_inst_mw[g]) / self.baseMVA
        return float(sgen.p_mw[g]) / self.baseMVA

    def static_generation_wind_var_q(self, net, grid_code=None):
        """Populate ``static_generation_data`` Q limits from the grid code.

        The capability is parameterised by operating variants (the ``var_q``
        column in ``net.sgen``); each selects one column of the selected grid
        code's capability table.  How many variants exist depends on the
        code — VDE-AR-N 4105 defines two, 4110 one and 4120 three — so
        ``var_q`` is validated against the resolved code.
        """
        code = resolve_grid_code(grid_code)
        self.grid_code = code
        self.q_limit_parameter = compute_q_curves(code)
        # Q/Pn capability table: row 0 capacitive, row 1 inductive;
        # columns are the var_q variants.
        y = code.vqu_q_max

        if "var_q" in self.net.sgen:
            self.static_generation_data["var_q"] = self.net.sgen.var_q.values
            sgens_var_q = self.static_generation_data.index[
                self.static_generation_data.var_q.notna()
            ]
            check_var_q(
                self.static_generation_data.var_q[sgens_var_q],
                code,
                context="net.sgen.var_q",
            )

            try:
                p_inst = self.net.sgen.p_inst_mw.values / self.baseMVA
            except AttributeError:
                logger.warning(
                    "No p_inst_mw attribute found in net.sgen. "
                    "Using p_mw as p_inst for wind generator power limits."
                )
                p_inst = self.static_generation_data["p"]

            self.static_generation_data["p_inst"] = p_inst

            self.static_generation_data["max_q"][sgens_var_q] = [
                y[0, int(self.static_generation_data.var_q[g])]
                * self.static_generation_data["p_inst"][g]
                for g in sgens_var_q
            ]
            self.static_generation_data["min_q"][sgens_var_q] = [
                y[1, int(self.static_generation_data.var_q[g])]
                * self.static_generation_data["p_inst"][g]
                for g in sgens_var_q
            ]

            self.static_generation_data["max_p"][sgens_var_q] = p_inst[
                sgens_var_q
            ]
            self.static_generation_data["min_p"][sgens_var_q] = (
                p_inst[sgens_var_q] * code.qp_p_high
            )

        else:
            self.static_generation_data["var_q"] = None
            self.static_generation_data["p_inst"] = None

        if "wind_hc" in self.net.sgen:
            self.static_generation_data["wind_hc"] = (
                self.net.sgen.wind_hc.values
            )
        else:
            self.static_generation_data["wind_hc"] = False

    def get_objective(self, model):
        """Add a wind-maximisation objective that subtracts network losses.

        Args:
            model: The Pyomo model being extended.

        Returns:
            None. The objective is attached to `model` as `obj`.
        """

        @model.Objective(sense=pyo.maximize)
        def obj(model):
            """Wind energy over the horizon, minus the losses carrying it.

            Maximised, so the objective rewards hosting capacity and charges
            for the losses it causes. Both terms are summed over `model.T`:
            what a multi-period study is asking is how much wind the network
            can absorb across the whole horizon, not at one instant, and a
            single step would let a candidate look free at the hour that
            happens to suit it.

            Args:
                model: The Pyomo model being extended.

            Returns:
                A Pyomo expression.
            """
            infeed = sum(
                model.psG[w, t] for w in model.WIND_HC for t in model.T
            )
            line_losses = sum(
                model.pLfrom[line, t] + model.pLto[line, t]
                for line in model.L
                for t in model.T
            )
            trafo_losses = sum(
                model.pThv[tr, t] + model.pTlv[tr, t]
                for tr in model.TRANSF
                for t in model.T
            )
            return infeed - line_losses - trafo_losses

    def get_constraints(self, model, net):
        """Add the Q(P) and Q(U) capability constraints.

        Add Q(P) and Q(U) constraints for controllable wind and HC
        generators.
        """
        # Each grid-code bound is a piecewise-linear envelope, so it becomes
        # one inequality per affine piece: the upper bound is the pointwise
        # minimum of its pieces, the lower bound the pointwise maximum.  A
        # single line cannot express the saturation shelf.
        pq_area = self.grid_code.pq_area
        qv_area = self.grid_code.qv_area
        v_span = bus_voltage_range(self.net)
        sGbs_lookup = {g: b for (g, b) in model.sGbs}
        model.W_QP_PIECE = pyo.RangeSet(
            0, pq_area.max_pieces(DEFAULT_P_RANGE_PU) - 1
        )
        model.W_QU_PIECE = pyo.RangeSet(0, qv_area.max_pieces(v_span) - 1)

        @model.Constraint(model.WINDc, model.T, model.W_QP_PIECE)
        def QW_pos(model, w, t, k):
            """Upper Q(P) capability piece for wind unit `w` at step `t`.

            Args:
                model: The Pyomo model being extended.
                w: Wind-generator index from `model.WINDc`.
                t: Time step from `model.T`.
                k: Piece index of the piecewise envelope.

            Returns:
                A Pyomo inequality expression, or `Constraint.Skip` when the
                envelope has fewer than `k + 1` pieces.
            """
            pieces = pq_area.upper_pieces(
                int(pyo.value(model.var_q[w, t])), DEFAULT_P_RANGE_PU
            )
            if k >= len(pieces):
                return pyo.Constraint.Skip
            m, b = pieces[k]
            return (
                model.qsG[w, t]
                <= m * model.psG[w, t] + b * model.PsG_inst[w, t]
            )

        @model.Constraint(model.WINDc, model.T, model.W_QP_PIECE)
        def QW_neg(model, w, t, k):
            """Lower Q(P) capability piece for wind unit `w` at step `t`.

            Args:
                model: The Pyomo model being extended.
                w: Wind-generator index from `model.WINDc`.
                t: Time step from `model.T`.
                k: Piece index of the piecewise envelope.

            Returns:
                A Pyomo inequality expression, or `Constraint.Skip` when the
                envelope has fewer than `k + 1` pieces.
            """
            pieces = pq_area.lower_pieces(
                int(pyo.value(model.var_q[w, t])), DEFAULT_P_RANGE_PU
            )
            if k >= len(pieces):
                return pyo.Constraint.Skip
            m, b = pieces[k]
            return (
                model.qsG[w, t]
                >= m * model.psG[w, t] + b * model.PsG_inst[w, t]
            )

        @model.Constraint(model.WINDc, model.T, model.W_QU_PIECE)
        def QV_min(model, w, t, k):
            """Lower Q(U) capability piece for wind unit `w` at step `t`.

            Args:
                model: The Pyomo model being extended.
                w: Wind-generator index from `model.WINDc`.
                t: Time step from `model.T`.
                k: Piece index of the piecewise envelope.

            Returns:
                A Pyomo inequality expression, or `Constraint.Skip` when the
                unit has no bus entry or the envelope has fewer than `k + 1`
                pieces.
            """
            if w not in sGbs_lookup:
                return pyo.Constraint.Skip
            pieces = qv_area.lower_pieces(
                int(pyo.value(model.var_q[w, t])), v_span
            )
            if k >= len(pieces):
                return pyo.Constraint.Skip
            b_bus = sGbs_lookup[w]
            m, b = pieces[k]
            return (
                model.qsG[w, t]
                >= (m * model.v[b_bus, t] + b) * model.PsG_inst[w, t]
            )

        @model.Constraint(model.WINDc, model.T, model.W_QU_PIECE)
        def QV_max(model, w, t, k):
            """Upper Q(U) capability piece for wind unit `w` at step `t`.

            Args:
                model: The Pyomo model being extended.
                w: Wind-generator index from `model.WINDc`.
                t: Time step from `model.T`.
                k: Piece index of the piecewise envelope.

            Returns:
                A Pyomo inequality expression, or `Constraint.Skip` when the
                unit has no bus entry or the envelope has fewer than `k + 1`
                pieces.
            """
            if w not in sGbs_lookup:
                return pyo.Constraint.Skip
            pieces = qv_area.upper_pieces(
                int(pyo.value(model.var_q[w, t])), v_span
            )
            if k >= len(pieces):
                return pyo.Constraint.Skip
            b_bus = sGbs_lookup[w]
            m, b = pieces[k]
            return (
                model.qsG[w, t]
                <= (m * model.v[b_bus, t] + b) * model.PsG_inst[w, t]
            )

        # --- sizing: decided once per candidate, not per time step ---
        @model.Constraint(model.WIND_HC)
        def hc_size_upper(model, w):
            r"""Cap the installed rating of candidate `w`, and gate it on `y`.

            $S^2_w \le S_{max,w}^2 y_w$. A zero selection variable forces the
            rating to zero, which is what makes the problem a MINLP; the
            binary is the only nonconvexity the hosting-capacity layer adds.

            Args:
                model: The Pyomo model being extended.
                w: Candidate index from `model.WIND_HC`.

            Returns:
                A Pyomo inequality expression.
            """
            return model.SW2[w] <= model.SWmax[w] ** 2 * model.y[w]

        @model.Constraint(model.WIND_HC)
        def hc_size_lower(model, w):
            r"""A selected candidate has to be at least `SWmin` in size.

            $S^2_w \ge S_{min,w}^2 y_w$. This is the multi-period home of
            what used to be a per-step minimum *dispatch*, which was wrong
            for wind: a plant that has to produce at every step cannot exist
            in a network whose wind is zero at night. A minimum makes sense
            for the size, not for the output.

            Args:
                model: The Pyomo model being extended.
                w: Candidate index from `model.WIND_HC`.

            Returns:
                A Pyomo inequality expression.
            """
            return model.SW2[w] >= model.SWmin[w] ** 2 * model.y[w]

        # --- dispatch: bound at every step by the size chosen above ---
        @model.Constraint(model.WIND_HC, model.T)
        def SW_max(model, w, t):
            r"""Keep candidate `w` inside its installed rating at step `t`.

            $p_{w,t}^2 + q_{w,t}^2 \le S^2_w$. Convex, because `SW2` carries
            the *squared* rating: writing the same thing against a rating
            variable would put a variable square on the right-hand side.

            Args:
                model: The Pyomo model being extended.
                w: Candidate index from `model.WIND_HC`.
                t: Time step from `model.T`.

            Returns:
                A Pyomo inequality expression.
            """
            return model.psG[w, t] ** 2 + model.qsG[w, t] ** 2 <= model.SW2[w]

        # Simplified HC Q-P bounds: the widest band the grid code offers
        @model.Constraint(model.WIND_HC, model.T)
        def QW_min(model, w, t):
            """Lower Q(P) bound for candidate `w` at step `t`.

            Args:
                model: The Pyomo model being extended.
                w: Candidate index from `model.WIND_HC`.
                t: Time step from `model.T`.

            Returns:
                A Pyomo inequality expression.
            """
            return model.qsG[w, t] >= self.qp_min * model.psG[w, t]

        @model.Constraint(model.WIND_HC, model.T)
        def QW_max(model, w, t):
            """Upper Q(P) bound for candidate `w` at step `t`.

            Args:
                model: The Pyomo model being extended.
                w: Candidate index from `model.WIND_HC`.
                t: Time step from `model.T`.

            Returns:
                A Pyomo inequality expression.
            """
            return model.qsG[w, t] <= self.qp_max * model.psG[w, t]

        @model.Constraint(model.WIND_HC, model.T)
        def QU_min_hc(model, w, t):
            """Lower Q(U) bound for candidate `w` at step `t`.

            Bilinear in the bus voltage and the active power, hence
            nonconvex.

            Args:
                model: The Pyomo model being extended.
                w: Candidate index from `model.WIND_HC`.
                t: Time step from `model.T`.

            Returns:
                A Pyomo inequality expression, or `Constraint.Skip` when the
                candidate has no bus entry in `model.sGbs`.
            """
            if w not in sGbs_lookup:
                return pyo.Constraint.Skip
            b = sGbs_lookup[w]
            return (
                model.qsG[w, t]
                >= (self.m_qu_min * model.v[b, t] + self.qu_min)
                * model.psG[w, t]
            )

        @model.Constraint(model.WIND_HC, model.T)
        def QU_max_hc(model, w, t):
            """Upper Q(U) bound for candidate `w` at step `t`.

            Args:
                model: The Pyomo model being extended.
                w: Candidate index from `model.WIND_HC`.
                t: Time step from `model.T`.

            Returns:
                A Pyomo inequality expression, or `Constraint.Skip` when the
                candidate has no bus entry in `model.sGbs`.
            """
            if w not in sGbs_lookup:
                return pyo.Constraint.Skip
            b = sGbs_lookup[w]
            return (
                model.qsG[w, t]
                <= (self.m_qu_max * model.v[b, t] + self.qu_max)
                * model.psG[w, t]
            )

        if "windpot_p_mw" in net.bus:

            @model.Constraint(model.WIND_HC, model.T)
            def PW_max(model, w, t):
                """Cap candidate `w` at the bus's wind potential at step `t`.

                `windpot_p_mw` records how much wind the site could host, so
                it bounds the output at every step rather than the energy
                over the horizon.

                Args:
                    model: The Pyomo model being extended.
                    w: Candidate index from `model.WIND_HC`.
                    t: Time step from `model.T`.

                Returns:
                    A Pyomo inequality expression.
                """
                return model.psG[w, t] <= model.pWmax[w]

    def unfix_variables(self, model):
        """Free the dispatch of every HC candidate, at every time step.

        The candidates carry a zero profile so the multi-period base class
        can resolve one for them; unfixing is what turns them from that
        placeholder infeed into decision variables.

        Args:
            model: The Pyomo model being extended.

        Returns:
            None. `psG` and `qsG` are unfixed in place.
        """
        for w in model.WIND_HC:
            for t in model.T:
                model.psG[w, t].unfix()
                model.qsG[w, t].unfix()

    def get_all_acopf(self, model):
        """Nothing extra for wind at the AC OPF stage.

        No additional ACOPF components needed for wind (called via
        get_all_opf).
        """
