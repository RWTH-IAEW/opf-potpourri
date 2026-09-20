"""DC Optimal Power Flow model combining linearised DC power flow and OPF
limits."""

import numpy as np
import pyomo.environ as pyo
from potpourri.models.DC import DC
from potpourri.models.OPF import OPF


class DCOPF(DC, OPF):
    """DC OPF model: linearised power flow with generator and thermal limit
    constraints."""

    def __init__(self, net, dc_susceptance: str = "matpower"):
        super().__init__(net, dc_susceptance=dc_susceptance)
        self.create_model()

    def create_model(self):
        """Build the DCOPF Pyomo model (sets model name to 'DCOPF')."""
        super().create_model()
        self.model.name = "DCOPF"

    def add_OPF(self, angle_limits: bool = False, **kwargs):
        """Attach DC-OPF sets, parameters, and constraints to ``self.model``.

        Calls :meth:`OPF.add_OPF` for generator/demand limits and line
        ratings, then adds DC thermal-limit constraints on the from-side
        flow (the DC model is lossless, so ``pLto = −pLfrom`` and a single
        bound suffices).

        Args:
            angle_limits: When ``True``, enforce branch phase-angle-difference
                constraints ``angmin ≤ δ_from − δ_to ≤ angmax`` using
                ``net.line.angmin_degree`` / ``net.line.angmax_degree`` (and
                the transformer equivalent if present). Defaults to ``False``
                to preserve previous behaviour.
            **kwargs: Forwarded to :meth:`OPF.add_OPF`.
        """
        super().add_OPF(**kwargs)

        # --- line power limits (check sending end; DC is approximately
        # lossless) ---
        def line_lim_upper(model, l):
            return model.pLfrom[l] <= model.SLmax[l]

        def line_lim_lower(model, l):
            return model.pLfrom[l] >= -model.SLmax[l]

        self.model.line_lim_from = pyo.Constraint(
            self.model.L, rule=line_lim_upper
        )
        self.model.line_lim_to = pyo.Constraint(
            self.model.L, rule=line_lim_lower
        )

        # --- transformer power limits ---
        def transf_lim_upper(model, l):
            return model.pThv[l] <= model.SLmaxT[l]

        def transf_lim_lower(model, l):
            return model.pThv[l] >= -model.SLmaxT[l]

        self.model.transf_lim1 = pyo.Constraint(
            self.model.TRANSF, rule=transf_lim_upper
        )
        self.model.transf_lim2 = pyo.Constraint(
            self.model.TRANSF, rule=transf_lim_lower
        )

        if angle_limits:
            self._add_dc_branch_angle_limits()

    def _add_dc_branch_angle_limits(self):
        """Attach branch phase-angle-difference constraints (DC variant).

        Reads angle bounds from ``{angmin,angmax}_degree`` on ``net.line``,
        ``net.trafo`` and ``net.impedance`` (impedance rows are the synthetic
        line indices behind ``net.line``) and applies them to
        ``delta[from] − delta[to]`` (radians).
        """

        def _bounds(table, idx_set, hv_col, lv_col):
            if (
                "angmin_degree" not in table.columns
                or "angmax_degree" not in table.columns
            ):
                return {}
            valid = set(table.index)
            imp = self.net.get("impedance")
            n_line = len(self.net.line.index)
            has_imp_bounds = (
                table is self.net.line
                and imp is not None
                and not imp.empty
                and "angmin_degree" in imp.columns
                and "angmax_degree" in imp.columns
            )
            out = {}
            for ix in idx_set:
                if ix in valid:
                    src, row, f_col, t_col = table, ix, hv_col, lv_col
                elif has_imp_bounds and 0 <= ix - n_line < len(imp):
                    # impedance rows sit in model.L behind net.line under
                    # synthetic indices n_line + i (see Basemodel); their
                    # MATPOWER angle bound arrives on net.impedance
                    src = imp
                    row = imp.index[ix - n_line]
                    f_col, t_col = "from_bus", "to_bus"
                else:
                    continue
                amin = float(src.at[row, "angmin_degree"])
                amax = float(src.at[row, "angmax_degree"])
                if (
                    not np.isfinite(amin)
                    or not np.isfinite(amax)
                    or abs(amin) >= 359.0
                    or abs(amax) >= 359.0
                ):
                    continue
                out[ix] = (
                    self.bus_lookup[int(src.at[row, f_col])],
                    self.bus_lookup[int(src.at[row, t_col])],
                    np.deg2rad(amin),
                    np.deg2rad(amax),
                )
            return out

        line_bounds = _bounds(
            self.net.line, list(self.model.L), "from_bus", "to_bus"
        )
        trafo_bounds = _bounds(
            self.net.trafo, list(self.model.TRANSF), "hv_bus", "lv_bus"
        )

        if line_bounds:
            line_idx = list(line_bounds.keys())
            self.model.LineAngleSet = pyo.Set(initialize=line_idx)

            def _line_angle_rule(model, l):
                f, t, amin, amax = line_bounds[l]
                return amin, model.delta[f] - model.delta[t], amax

            self.model.line_angle_diff = pyo.Constraint(
                self.model.LineAngleSet, rule=_line_angle_rule
            )

        if trafo_bounds:
            tr_idx = list(trafo_bounds.keys())
            self.model.TrafoAngleSet = pyo.Set(initialize=tr_idx)

            def _trafo_angle_rule(model, l):
                f, t, amin, amax = trafo_bounds[l]
                return amin, model.delta[f] - model.delta[t], amax

            self.model.trafo_angle_diff = pyo.Constraint(
                self.model.TrafoAngleSet, rule=_trafo_angle_rule
            )
