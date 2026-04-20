"""Full AC OPF model combining AC power flow and OPF operational limits."""

import pandas as pd
import pyomo.environ as pyo
from src.potpourri.models.AC import AC
from src.potpourri.models.OPF import OPF
import numpy as np
from loguru import logger


class ACOPF(AC, OPF):
    """Full AC Optimal Power Flow model.

    Combines AC power flow physics with operational limit constraints via
    multiple inheritance. Adds voltage bounds, reactive power bounds, apparent
    power thermal limits on lines and transformers, and optional wind Q-curve
    constraints for grid-code compliance.

    Args:
        net: A pandapower network with voltage limits and generator data.
    """

    def __init__(self, net):
        super().__init__(net)

    def _calc_opf_parameters(self):
        """Compute all AC-OPF limit data from the network.

        Extends OPF._calc_opf_parameters() with bus voltage limits, reactive
        power limits for static generators and external grids, and reactive
        demand bounds.
        """
        super()._calc_opf_parameters()

        max_vm_pu, min_vm_pu = self.get_v_limits()
        self.v_limits = (max_vm_pu, min_vm_pu)

        self.static_generation_reactive_power_limits()
        self.generation_reactive_power_limits()
        self.get_demand_reactive_data()

    def static_generation_reactive_power_limits(self):
        """Read reactive power limits for static generators from net.sgen.

        Populates static_generation_data['max_q'] and ['min_q'] (per-unit).
        Also calls static_generation_wind_var_q() for wind-specific Q limits.
        """
        if "controllable" in self.net.sgen:
            self.static_generation_data["controllable"] = (
                self.net.sgen.controllable.values
            )
        else:
            self.static_generation_data["controllable"] = False

        lim_q = abs(
            self.net.sgen.q_mvar
        )  # MVAr — divided by baseMVA once below
        if "max_q_mvar" in self.net.sgen:
            self.static_generation_data["max_q"] = (
                self.net.sgen.max_q_mvar.fillna(lim_q).values / self.baseMVA
            )
        else:
            self.static_generation_data["max_q"] = lim_q.values / self.baseMVA

        if "min_q_mvar" in self.net.sgen:
            self.static_generation_data["min_q"] = (
                self.net.sgen.min_q_mvar.fillna(-lim_q).values / self.baseMVA
            )
        else:
            self.static_generation_data["min_q"] = -lim_q.values / self.baseMVA

        self.static_generation_wind_var_q()
        self.static_generation_data["type"] = self.net.sgen.type.values

    def generation_reactive_power_limits(self):
        """Read reactive power limits for external grid generators.

        Populates generation_data['max_q'] and ['min_q'] (per-unit).
        """
        max_q = np.full(len(self.generation_data), 1e9) / self.baseMVA
        min_q = np.full(len(self.generation_data), -1e9) / self.baseMVA

        for element, (f, t) in self.net._gen_order.items():
            if "max_q_mvar" in self.net[element]:
                max_q[f:t] = (
                    self.net[element].max_q_mvar.fillna(1e9).values
                    / self.baseMVA
                )
            if "min_q_mvar" in self.net[element]:
                min_q[f:t] = (
                    self.net[element].min_q_mvar.fillna(-1e9).values
                    / self.baseMVA
                )

        self.generation_data["max_q"] = max_q
        self.generation_data["min_q"] = min_q

    def get_v_limits(self):
        """Read bus voltage bounds from net.bus.

        Returns:
            tuple: (max_vm_pu, min_vm_pu) as numpy arrays indexed by internal
            bus id. Defaults to 1.1 / 0.9 if columns are absent.
            Generator-level limits override bus limits where stricter.
        """
        if "max_vm_pu" in self.net.bus:
            max_vm_pu = self.net.bus.max_vm_pu.values
        else:
            max_vm_pu = np.full(len(self.net.bus.index), 1.1)

        if "min_vm_pu" in self.net.bus:
            min_vm_pu = self.net.bus.min_vm_pu.values
        else:
            min_vm_pu = np.full(len(self.net.bus.index), 0.9)

        if any(self.net.gen.index):
            self.add_generator_v_limits(max_vm_pu, min_vm_pu)

        return max_vm_pu, min_vm_pu

    def add_generator_v_limits(self, max_vm_pu, min_vm_pu):
        """Apply per-generator voltage limits, overriding bus defaults.

        Modifies max_vm_pu and min_vm_pu in-place. Generator limits that exceed
        the bus limit are ignored with a warning; otherwise they override the
        bus-level value.

        Args:
            max_vm_pu: Array of per-bus upper voltage limits (per-unit).
            min_vm_pu: Array of per-bus lower voltage limits (per-unit).
        """
        # check max_vm_pu / min_vm_pu bus limit violation by gens
        gen_buses = self.bus_lookup[self.net.gen.bus.values]
        if "max_vm_pu" in self.net["gen"].columns:
            v_max_bound = (
                max_vm_pu[gen_buses] < self.net["gen"]["max_vm_pu"].values
            )
            if np.any(v_max_bound):
                bound_gens = self.net["gen"].index.values[v_max_bound]
                logger.warning(
                    "gen max_vm_pu > bus max_vm_pu for gens {}. "
                    "Setting bus limit for these gens.",
                    bound_gens,
                )
                # pyo.Set only vm of gens which do not violate the limits
                max_vm_pu[gen_buses[~v_max_bound]] = self.net["gen"][
                    "max_vm_pu"
                ].values[~v_max_bound]
            else:
                # pyo.Set vm of all gens
                max_vm_pu[gen_buses] = self.net["gen"]["max_vm_pu"].values

        if "min_vm_pu" in self.net["gen"].columns:
            v_min_bound = (
                self.net["gen"]["min_vm_pu"].values < min_vm_pu[gen_buses]
            )
            if np.any(v_min_bound):
                bound_gens = self.net["gen"].index.values[v_min_bound]
                logger.warning(
                    "gen min_vm_pu < bus min_vm_pu for gens {}. "
                    "Setting bus limit for these gens.",
                    bound_gens,
                )
                # pyo.Set only vm of gens which do not violate the limits
                min_vm_pu[gen_buses[~v_min_bound]] = self.net["gen"][
                    "min_vm_pu"
                ].values[~v_min_bound]
            else:
                # pyo.Set vm of all gens
                min_vm_pu[gen_buses] = self.net["gen"]["min_vm_pu"].values

        if "controllable" in self.net.gen:
            controllable = self.net["gen"]["controllable"].values
            not_controllable = ~controllable.astype(bool)

            # get voltage pyo.Setpoints for not controllable generators
            if np.any(not_controllable):
                bus = self.net["gen"]["bus"].values[not_controllable]
                vm_pu = self.net["gen"]["vm_pu"].values[not_controllable]

                not_controllable_buses = self.bus_lookup[bus]
                max_vm_pu[not_controllable_buses] = vm_pu
                min_vm_pu[not_controllable_buses] = vm_pu

        return max_vm_pu, min_vm_pu

    def get_demand_reactive_data(self):
        """Read reactive power bounds for loads from net.load.

        Populates self.QDmax_data and self.QDmin_data (per-unit). Falls back
        to active power magnitude if no reactive power values are set.
        """
        # reactive power demand
        # use active power for reactive power limits, if no reactive power
        # given for any sgen
        if self.net.load.q_mvar.sum() == 0:
            lim_q = abs(self.net.load.p_mw)
        else:
            lim_q = abs(self.net.load.q_mvar)

        # add rows with reactive generation limits if not existing
        if "max_q_mvar" not in self.net.load:
            self.net.load["max_q_mvar"] = lim_q

        if "min_q_mvar" not in self.net.load:
            self.net.load["min_q_mvar"] = -lim_q

        # demand limits for loads
        self.QDmax_data = (
            self.net.load.max_q_mvar.fillna(self.net.load.q_mvar)
            / self.baseMVA
        )
        self.QDmin_data = self.net.load.min_q_mvar.fillna(0) / self.baseMVA

    def static_generation_wind_var_q(self):
        """Compute Q-P and Q-U characteristic limits for wind generators.

        Reads sgen.var_q (variant 0–2) and assigns grid-code-compliant reactive
        power bounds. Populates self.q_limit_parameter with slope/intercept
        parameters for the Q-P and Q-U curves used in add_OPF() constraints.
        """
        x = np.array([[96, 103], [120, 127]]) / 110
        y = np.array([[0.48, 0.41, 0.33], [-0.23, -0.33, -0.41]])
        m = (y[1] - y[0]) / (x[0, 1] - x[0, 0])
        b = np.array([y[0] - m * x[i, 0] for i in range(len(x))]).T
        m_qp_max = (0.1 - y[0]) / (0.1 - 0.2)
        m_qp_min = (-0.1 - y[1]) / (0.1 - 0.2)
        b_qp_max = 0.1 - m_qp_max * 0.1
        b_qp_min = -0.1 - m_qp_min * 0.1

        self.q_limit_parameter = pd.DataFrame(
            {
                "m_qv": m,
                "b_qv_min": b[:, 0],
                "b_qv_max": b[:, 1],
                "m_qp_max": m_qp_max,
                "m_qp_min": m_qp_min,
                "b_qp_max": b_qp_max,
                "b_qp_min": b_qp_min,
            }
        )

        if "var_q" in self.net.sgen:
            self.static_generation_data["var_q"] = self.net.sgen.var_q.values
            sgens_var_q = self.static_generation_data.index[
                self.static_generation_data.var_q.notna()
            ]

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
                p_inst[sgens_var_q] * 0.1
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

    def add_OPF(self, **kwargs):
        """Attach AC-OPF sets, parameters, and constraints to self.model.

        Extends OPF.add_OPF() with bus voltage bounds (Vmin, Vmax), apparent
        power limits on lines and transformers, reactive power bounds for
        static generators, external grids, and controllable loads, and wind
        Q-P/Q-U
        curve constraints for sgens with var_q set.

        Args:
            **kwargs: Forwarded to _calc_opf_parameters.
        """
        super().add_OPF(**kwargs)

        self.model.name = "ACOPF"

        # --- pyo.Sets ---
        # generators for hc calculation
        self.model.WIND_HC = pyo.Set(
            within=self.model.sG,
            initialize=self.static_generation_data.index[
                self.static_generation_data["wind_hc"]
                & self.static_generation_data.in_service
            ],
        )
        # all wind generators
        self.model.WIND = self.model.WIND_HC | pyo.Set(
            within=self.model.sG,
            initialize=self.static_generation_data.index[
                (self.static_generation_data["type"] == "Wind")
                & self.static_generation_data.in_service
            ],
        )
        # controllable wind generators, not for hc calculation
        self.model.WINDc = (
            self.model.WIND
            & self.model.sGc
            & pyo.Set(
                initialize=self.static_generation_data.index[
                    self.static_generation_data["var_q"].values != None  # noqa: E711
                ]
            )
        )

        self.model.var_q = pyo.Param(
            self.model.WINDc,
            initialize=self.static_generation_data["var_q"][self.model.WINDc],
        )
        self.model.PsG_inst = pyo.Param(
            self.model.WINDc,
            initialize=self.static_generation_data["p_inst"][self.model.WINDc],
        )

        # voltage limits
        self.model.Vmax = pyo.Param(
            self.model.B,
            within=pyo.NonNegativeReals,
            initialize=self.v_limits[0][self.model.B],
        )  # max voltage (p.u.)
        self.model.Vmin = pyo.Param(
            self.model.B,
            within=pyo.NonNegativeReals,
            initialize=self.v_limits[1][self.model.B],
        )  # min voltage (p.u.)

        # generation reactive power limits
        self.model.QGmax = pyo.Param(
            self.model.G,
            initialize=self.generation_data["max_q"][self.model.G],
        )
        self.model.QGmin = pyo.Param(
            self.model.G,
            initialize=self.generation_data["min_q"][self.model.G],
        )

        # static generation reactive power limits
        self.model.QsGmax = pyo.Param(
            self.model.sGc,
            within=pyo.Reals,
            initialize=self.static_generation_data["max_q"][self.model.sGc],
            mutable=True,
        )
        self.model.QsGmin = pyo.Param(
            self.model.sGc,
            within=pyo.Reals,
            initialize=self.static_generation_data["min_q"][self.model.sGc],
            mutable=True,
        )

        # reactive demand
        self.model.QDmax = pyo.Param(
            self.model.D, initialize=self.QDmax_data[self.model.D]
        )
        self.model.QDmin = pyo.Param(
            self.model.D, initialize=self.QDmin_data[self.model.D]
        )

        # --- line power limits ---
        def line_lim_from_def(model, l):
            return (
                model.pLfrom[l] ** 2 + model.qLfrom[l] ** 2
                <= model.SLmax[l] ** 2 * model.v[model.A[l, 1]] ** 2
            )

        def line_lim_to_def(model, l):
            return (
                model.pLto[l] ** 2 + model.qLto[l] ** 2
                <= model.SLmax[l] ** 2 * model.v[model.A[l, 2]] ** 2
            )

        self.model.line_lim_from = pyo.Constraint(
            self.model.L, rule=line_lim_from_def
        )
        self.model.line_lim_to = pyo.Constraint(
            self.model.L, rule=line_lim_to_def
        )

        # --- power flow limits on transformer lines---
        def transf_lim1_def(model, l):
            return (
                model.pThv[l] ** 2 + model.qThv[l] ** 2
                <= model.SLmaxT[l] ** 2 * model.v[model.AT[l, 1]] ** 2
            )

        def transf_lim2_def(model, l):
            return (
                model.pTlv[l] ** 2 + model.qTlv[l] ** 2
                <= model.SLmaxT[l] ** 2 * model.v[model.AT[l, 2]] ** 2
            )

        self.model.transf_lim1 = pyo.Constraint(
            self.model.TRANSF, rule=transf_lim1_def
        )
        self.model.transf_lim2 = pyo.Constraint(
            self.model.TRANSF, rule=transf_lim2_def
        )

        # --- static generation reactive power limits ---
        def static_generation_reactive_power_bounds(model, g):
            model.qsG[g].unfix()
            return model.QsGmin[g], model.qsG[g], model.QsGmax[g]

        self.model.QsG_pyo = pyo.Constraint(
            self.model.sGc, rule=static_generation_reactive_power_bounds
        )

        # --- reactive generator power limits ---
        def reactive_power_bounds(model, g):
            model.qG[g].unfix()
            return model.QGmin[g], model.qG[g], model.QGmax[g]

        self.model.QG_pyo = pyo.Constraint(
            self.model.G, rule=reactive_power_bounds
        )

        # --- reactive demand limits ---
        def reactive_demand_bounds(model, d):
            model.qD[d].unfix()
            return model.QDmin[d], model.qD[d], model.QDmax[d]

        self.model.QD_pyo = pyo.Constraint(
            self.model.Dc, rule=reactive_demand_bounds
        )

        # --- voltage pyo.Constraints ---
        self.model.v_bPV_setpoint.deactivate()

        def v_bounds(model, b):
            return model.Vmin[b], model.v[b], model.Vmax[b]

        self.model.v_pyo = pyo.Constraint(self.model.B, rule=v_bounds)

        # --- wind generation q requirements variant 3---
        def QW_pos(model, w):
            return (
                model.qsG[w]
                <= self.q_limit_parameter.b_qp_max[model.var_q[w]]
                * model.PsG_inst[w]
                + self.q_limit_parameter.m_qp_max[model.var_q[w]]
                * model.psG[w]
            )

        def QW_neg(model, w):
            return (
                model.qsG[w]
                >= self.q_limit_parameter.b_qp_min[model.var_q[w]]
                * model.PsG_inst[w]
                + self.q_limit_parameter.m_qp_min[model.var_q[w]]
                * model.psG[w]
            )

        self.model.QW_pos_pyo = pyo.Constraint(self.model.WINDc, rule=QW_pos)
        self.model.QW_neg_pyo = pyo.Constraint(self.model.WINDc, rule=QW_neg)

        sGbs_lookup = {g: b for (g, b) in self.model.sGbs}

        def QV_min(model, w):
            if w not in sGbs_lookup:
                return pyo.Constraint.Skip
            b = sGbs_lookup[w]
            return (
                model.qsG[w]
                >= (
                    self.q_limit_parameter.m_qv[model.var_q[w]] * model.v[b]
                    + self.q_limit_parameter.b_qv_min[model.var_q[w]]
                )
                * model.PsG_inst[w]
            )

        self.model.QU_min_pyo = pyo.Constraint(self.model.WINDc, rule=QV_min)

        def QV_max(model, w):
            if w not in sGbs_lookup:
                return pyo.Constraint.Skip
            b = sGbs_lookup[w]
            return (
                model.qsG[w]
                <= (
                    self.q_limit_parameter.m_qv[model.var_q[w]] * model.v[b]
                    + self.q_limit_parameter.b_qv_max[model.var_q[w]]
                )
                * model.PsG_inst[w]
            )

        self.model.QU_max_pyo = pyo.Constraint(self.model.WINDc, rule=QV_max)

    def add_voltage_deviation_objective(self):
        """Set objective to minimise sum of squared voltage deviations from
        1 p.u.

        Minimises Σ (v[b] - 1)² over non-slack buses and
        Σ (v[b] - v_b0[b])² over slack buses.
        """
        self.model.vm = pyo.Param(
            self.model.B, initialize=self.bus_data["v_m"][self.model.B]
        )

        def voltage_deviation_objective(model):
            return sum(
                (model.v[b] - 1.0) ** 2 for b in model.B - model.b0
            ) + sum((model.v[b] - model.v_b0[b]) ** 2 for b in model.b0)

        self.model.obj_v_deviation = pyo.Objective(
            rule=voltage_deviation_objective, sense=pyo.minimize
        )

    def add_reactive_power_flow_objective(self):
        """Set objective to minimise total squared reactive generation.

        Minimises Σ qsG[g]² over all static generators.
        """

        def reactive_objective(model):
            # Minimize the reactive power
            return sum(model.qsG[g] ** 2 for g in model.sG)

        self.model.obj_reactive = pyo.Objective(
            rule=reactive_objective, sense=pyo.minimize
        )
