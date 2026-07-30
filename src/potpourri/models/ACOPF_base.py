"""Full AC OPF model combining AC power flow and OPF operational limits."""

import warnings

import numpy as np
import pandas as pd
import pyomo.environ as pyo
from loguru import logger

from potpourri.models.AC import AC
from potpourri.models.OPF import OPF
from potpourri.technologies.q_control import (
    DEFAULT_WIND_SGEN_TYPES,
    SgenTypeOverlapWarning,
    compute_q_curves,
    resolve_grid_code,
)

# sgen ``type`` values treated as PV for ``pv_q_control``.  SimBench spells PV
# differently per voltage level — its RES dataset uses "PV" in LV (and HV2),
# "PV_MV" in MV and lowercase "pv" in EHV — so matching only "PV" silently
# reaches nothing on any SimBench MV grid.  Matching is exact and
# case-sensitive.  Override per call with ``add_OPF(sgen_types=...)`` to
# include further categories: SimBench also ships "Wind"/"Wind_MV"/"wind
# onshore"/"wind offshore", "Biomass_MV", "Hydro_MV", and the aggregated
# "lv_RES" (low-voltage renewables lumped into one MV element).
DEFAULT_PV_SGEN_TYPES = ("PV", "PV_MV", "pv")


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
                self.net.sgen.max_q_mvar.astype(float)
                .fillna(lim_q.astype(float))
                .values
                / self.baseMVA
            )
        else:
            self.static_generation_data["max_q"] = lim_q.values / self.baseMVA

        if "min_q_mvar" in self.net.sgen:
            self.static_generation_data["min_q"] = (
                self.net.sgen.min_q_mvar.astype(float)
                .fillna(-lim_q.astype(float))
                .values
                / self.baseMVA
            )
        else:
            self.static_generation_data["min_q"] = -lim_q.values / self.baseMVA

        self.static_generation_wind_var_q()
        self.static_generation_data["type"] = self.net.sgen.type.values

    def generation_reactive_power_limits(self):
        """Read reactive power limits for external grid generators.

        Populates ``generation_data['max_q']`` and ``['min_q']`` (per-unit).
        Filters each source table to in-service rows so the slice in
        ``net._gen_order`` matches the broadcast target (see also
        ``generation_real_power_limits``).
        """
        max_q = np.full(len(self.generation_data), 1e9) / self.baseMVA
        min_q = np.full(len(self.generation_data), -1e9) / self.baseMVA

        for element, (f, t) in self.net._gen_order.items():
            if element not in self.net or self.net[element].empty:
                continue
            table = self.net[element]
            if "in_service" in table.columns:
                table = table.loc[table.in_service.astype(bool)]
            if len(table) != t - f:
                continue
            if "max_q_mvar" in table:
                max_q[f:t] = (
                    table.max_q_mvar.astype(float).fillna(1e9).values
                    / self.baseMVA
                )
            if "min_q_mvar" in table:
                min_q[f:t] = (
                    table.min_q_mvar.astype(float).fillna(-1e9).values
                    / self.baseMVA
                )

        self.generation_data["max_q"] = max_q
        self.generation_data["min_q"] = min_q

    def get_v_limits(self):
        """Read bus voltage bounds from net.bus, keyed by ppc bus number.

        Returns:
            tuple: (max_vm_pu, min_vm_pu) as :class:`~pandas.Series` indexed
            by **ppc** bus number, covering only the ppc buses that a
            pandapower bus maps onto (``self.ppc_buses_with_pd``).  Defaults
            to 1.1 / 0.9 if the columns are absent.  Generator-level limits
            override bus limits where stricter.

        Note:
            The index is ppc bus numbers, not pandapower bus indices, because
            every consumer looks these up through ``self.bus_lookup``.  Plain
            positional arrays were correct only while the two numbering spaces
            happened to coincide, which is not the case on grids where
            pandapower's ppc conversion adds auxiliary buses.  Those auxiliary
            buses carry no pandapower row and therefore no user-supplied
            limits, so they are deliberately absent here and their voltage is
            left to follow from the network equations.
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
        # tightest band so the merged node cannot be operated outside the
        # limits of any bus that formed it.
        max_vm_pu = max_vm_pu.groupby(level=0).min()
        min_vm_pu = min_vm_pu.groupby(level=0).max()

        if any(self.net.gen.index):
            self.add_generator_v_limits(max_vm_pu, min_vm_pu)

        return max_vm_pu, min_vm_pu

    def add_generator_v_limits(self, max_vm_pu, min_vm_pu):
        """Apply per-generator voltage limits, overriding bus defaults.

        Modifies max_vm_pu and min_vm_pu in-place. Generator limits that exceed
        the bus limit are ignored with a warning; otherwise they override the
        bus-level value.

        Args:
            max_vm_pu: Per-bus upper voltage limits (per-unit), a
                :class:`~pandas.Series` indexed by ppc bus number.
            min_vm_pu: Per-bus lower voltage limits (per-unit), same indexing.
        """
        # check max_vm_pu / min_vm_pu bus limit violation by gens.
        # gen_buses are ppc bus numbers, matching the Series index.
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
                # pyo.Set only vm of gens which do not violate the limits
                max_vm_pu.loc[gen_buses[~v_max_bound]] = self.net["gen"][
                    "max_vm_pu"
                ].values[~v_max_bound]
            else:
                # pyo.Set vm of all gens
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
                # pyo.Set only vm of gens which do not violate the limits
                min_vm_pu.loc[gen_buses[~v_min_bound]] = self.net["gen"][
                    "min_vm_pu"
                ].values[~v_min_bound]
            else:
                # pyo.Set vm of all gens
                min_vm_pu.loc[gen_buses] = self.net["gen"]["min_vm_pu"].values

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
            self.net.load.max_q_mvar.astype(float).fillna(
                self.net.load.q_mvar.astype(float)
            )
            / self.baseMVA
        )
        self.QDmin_data = (
            self.net.load.min_q_mvar.astype(float).fillna(0.0) / self.baseMVA
        )

    def static_generation_wind_var_q(self):
        """Compute Q-P and Q-U characteristic limits for wind generators.

        Reads sgen.var_q (variant 0–2) and assigns grid-code-compliant reactive
        power bounds. Populates self.q_limit_parameter with slope/intercept
        parameters for the Q-P and Q-U curves used in add_OPF() constraints.

        The curves come from the grid code selected via
        :meth:`add_OPF`'s ``grid_code`` argument (default VDE-AR-N 4105).
        """
        code = resolve_grid_code(getattr(self, "_grid_code", None))
        self.q_limit_parameter = compute_q_curves(code)

        # Q/Pn capability table: row 0 capacitive, row 1 inductive;
        # columns are the var_q variants.
        q_max_table = code.vqu_q_max

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
                q_max_table[0, int(self.static_generation_data.var_q[g])]
                * self.static_generation_data["p_inst"][g]
                for g in sgens_var_q
            ]
            self.static_generation_data["min_q"][sgens_var_q] = [
                q_max_table[1, int(self.static_generation_data.var_q[g])]
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

    def add_OPF(
        self,
        thermal_limit: str = "current",
        free_slack_vm: bool = True,
        fix_hv_buses: bool = False,
        hv_bus_kv: float = 110.0,
        angle_limits: bool = False,
        pv_q_control: "str | bool | None" = None,
        inverter_s2: bool = False,
        cos_phi_min: "float | None" = None,
        pu_curtail: bool = False,
        fixed_cos_phi: "float | None" = None,
        cos_phi_p_profile: bool = False,
        grid_code=None,
        sgen_types=None,
        wind_sgen_types=None,
        **kwargs,
    ):
        """Attach AC-OPF sets, parameters, and constraints to ``self.model``.

        Extends :meth:`OPF.add_OPF` with bus voltage bounds (Vmin, Vmax),
        apparent-power thermal limits on lines and transformers, reactive
        power bounds for static generators, external grids and controllable
        loads, the wind Q-P / Q-U capability constraints for sgens with
        ``var_q`` set, and optionally the same Q(P)/Q(U) grid-code constraints
        for PV-type sgens via ``pv_q_control``.

        Args:
            pv_q_control: VDE-AR-N 4105 Q-control mode for controllable
                PV-type sgens (``type == "PV"`` and ``var_q`` set in
                ``net.sgen``).  Accepted values:

                * ``None`` or ``False`` — no Q-control (default)
                * ``"qp"`` — Q(P) characteristic only
                * ``"qu"`` — Q(U) droop only (requires AC model with ``v``)
                * ``"both"`` or ``True`` — Q(P) and Q(U) combined
            inverter_s2: When ``True``, add the apparent-power circle
                constraint ``psG[g]² + qsG[g]² ≤ S_inv[g]²`` for every
                controllable sgen with a finite ``net.sgen.sn_mva``.  The
                rating ``S_inv = sn_mva * converter_sizing_pu / baseMVA``
                (``converter_sizing_pu`` defaults to 1.0 when absent).
                Defaults to ``False`` (opt-in).
            cos_phi_min: Minimum power factor for the cos(φ) cone constraint
                ``|qsG[g]| ≤ psG[g] · tan(arccos(cos_phi_min))``, applied to
                every sgen in ``sGinv``.  Combined with ``psG ≥ 0`` this
                restricts PV operation to a "pizza-slice" rather than a full
                circle.  Per-sgen overrides can be set via
                ``net.sgen["cos_phi_min"]`` (takes precedence over this
                scalar).  Defaults to ``None`` (unconstrained). Typical
                VDE-AR-N 4105 value: ``0.90``.  Only active when
                ``inverter_s2=True`` and ``sn_mva`` is present.
            thermal_limit: ``"current"`` enforces ``|S|² ≤ SLmax² · v²``
                (current-limit form, physically meaningful for distribution
                conductors). ``"mva"`` enforces ``|S|² ≤ SLmax²``
                (constant-MVA limit, matches MATPOWER / PGLib-OPF). Defaults
                to ``"current"`` for backward compatibility.
            free_slack_vm: When ``True`` (default), the slack-bus voltage
                magnitude floats within ``[Vmin, Vmax]``. The reference angle
                stays fixed. Set ``False`` to reproduce the legacy AC-PF
                behaviour where the slack ``vm`` is pinned to the base-case
                value (e.g. for redispatch studies around a fixed slack).
            fix_hv_buses: When ``True``, pin the voltage magnitude of every
                bus with ``vn_kv == hv_bus_kv`` to its base-case voltage.
                Disabled by default. The historical 110 kV pinning in
                German distribution studies can be re-enabled by setting
                ``fix_hv_buses=True``.
            hv_bus_kv: Voltage level (kV) used by ``fix_hv_buses``.
            pu_curtail: When ``True``, add the P(U) active-power curtailment
                constraint for controllable PV-type sgens (VDE-AR-N 4105
                §8.5).  Above a voltage threshold the allowed active-power
                output is reduced linearly to zero:
                ``psG[g] · ΔV ≤ P_inst[g] · (V_max − v[bus[g]])``.
                Per-sgen thresholds are read from ``net.sgen.v_curtail_pu``
                (default 1.06 p.u.) and ``net.sgen.v_max_curtail_pu``
                (default 1.10 p.u.).  Requires ``net.sgen.p_inst_mw``
                (falls back to ``net.sgen.p_mw``).  Only active in AC models.
            fixed_cos_phi: Fixed power-factor equality
                ``qsG[g] == psG[g] · tan(arccos(cos_phi))`` applied to every
                controllable sgen.  Supply a scalar to apply one value to all
                sgens, or set ``net.sgen["fixed_cos_phi"]`` per-row (takes
                precedence).  Defaults to ``None`` (disabled).
            cos_phi_p_profile: When ``True``, add the VDE-AR-N 4105
                cos(φ)(P) profile as a quadratic equality
                ``qsG · (Pn − Pt) == tan_phi · psG · (psG − Pt)``.
                Reads ``net.sgen.cos_phi_min`` (power factor at full output),
                ``net.sgen.p_inst_mw`` (installed capacity Pn), and
                optionally ``net.sgen.cpp_p_threshold_pu`` (default 0.2).
                Requires IPOPT or another NLP solver.
            angle_limits: When ``True``, enforce branch phase-angle-difference
                constraints ``angmin ≤ δ_from − δ_to ≤ angmax`` using
                ``net.line.angmin_degree`` / ``net.line.angmax_degree`` (and
                the transformer equivalent if present). Defaults disabled to
                preserve previous behaviour.
            sgen_types: sgen ``type`` values treated as PV by
                ``pv_q_control``.  Defaults to
                :data:`DEFAULT_PV_SGEN_TYPES` (``("PV", "PV_MV", "pv")``).
                Matching is exact, so a network whose sgens use other
                category names needs them listed here — for example
                ``sgen_types=("PV", "PV_MV", "pv", "lv_RES")`` to include
                the aggregated LV-renewable units.  Wind categories do not
                belong here: they are handled by ``wind_sgen_types`` below,
                and listing them in both places makes an sgen match both
                paths, which emits
                :class:`~potpourri.technologies.q_control.SgenTypeOverlapWarning`
                and leaves it to the wind path.  Note the multi-period model
                applies no type filter at all (it keys purely off
                ``var_q``).
            wind_sgen_types: sgen ``type`` values treated as wind by the wind
                Q-control path (``model.WIND`` / ``model.WINDc``).  Defaults
                to
                :data:`~potpourri.technologies.q_control.DEFAULT_WIND_SGEN_TYPES`
                (``("Wind", "Wind_MV", "wind onshore", "wind offshore")``),
                covering every SimBench spelling.  Matching is exact.
                Hosting-capacity units selected via ``net.sgen.wind_hc`` are
                included regardless of type.
            grid_code: Technical connection rule supplying the Q(P)/Q(U)
                capability envelope and the P(U) / cos(phi)(P) thresholds.
                Accepts ``None`` (VDE-AR-N 4105, the default), a short name
                such as ``"4105"`` or ``"4110"``, or a
                :class:`~potpourri.technologies.q_control.GridCode`.
                Selecting a grid code whose parameters are still
                placeholders emits a
                :class:`~potpourri.technologies.q_control.ProvisionalGridCodeWarning`.
            **kwargs: Forwarded to :meth:`_calc_opf_parameters`.
        """
        # Resolve before super(), which reaches static_generation_wind_var_q
        # via _calc_opf_parameters and needs the selected code.
        code = resolve_grid_code(grid_code)
        self._grid_code = code

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
        # All wind generators.  SimBench spells wind differently per
        # voltage level ("Wind" in HV, "Wind_MV" in MV, "wind onshore" /
        # "wind offshore" in EHV), so match against a list rather than the
        # single literal that reached nothing outside HV.
        _wind_types = (
            DEFAULT_WIND_SGEN_TYPES
            if wind_sgen_types is None
            else tuple(wind_sgen_types)
        )
        self.model.WIND = self.model.WIND_HC | pyo.Set(
            within=self.model.sG,
            initialize=self.static_generation_data.index[
                self.static_generation_data["type"].isin(_wind_types)
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

        # Voltage limits, over Bpd only: auxiliary ppc buses have no
        # pandapower row and therefore no user-supplied limits.  Their voltage
        # follows from the power-flow equations instead.
        self.model.Vmax = pyo.Param(
            self.model.Bpd,
            within=pyo.NonNegativeReals,
            initialize=self.v_limits[0][self.model.Bpd],
        )  # max voltage (p.u.)
        self.model.Vmin = pyo.Param(
            self.model.Bpd,
            within=pyo.NonNegativeReals,
            initialize=self.v_limits[1][self.model.Bpd],
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

        # --- line and transformer apparent-power limits ---
        if thermal_limit not in ("current", "mva"):
            raise ValueError(
                f"thermal_limit must be 'current' or 'mva', got "
                f"{thermal_limit!r}"
            )

        if thermal_limit == "current":
            # |S|^2 ≤ SLmax^2 · v^2 (i.e. |I| ≤ I_max). Physically meaningful
            # for thermal current rating; varies with voltage.
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
        else:
            # |S|^2 ≤ SLmax^2 (constant-MVA limit, matches MATPOWER /
            # PowerModels' constraint_thermal_limit_* and PGLib-OPF rate_a).
            def line_lim_from_def(model, l):
                return (
                    model.pLfrom[l] ** 2 + model.qLfrom[l] ** 2
                    <= model.SLmax[l] ** 2
                )

            def line_lim_to_def(model, l):
                return (
                    model.pLto[l] ** 2 + model.qLto[l] ** 2
                    <= model.SLmax[l] ** 2
                )

            def transf_lim1_def(model, l):
                return (
                    model.pThv[l] ** 2 + model.qThv[l] ** 2
                    <= model.SLmaxT[l] ** 2
                )

            def transf_lim2_def(model, l):
                return (
                    model.pTlv[l] ** 2 + model.qTlv[l] ** 2
                    <= model.SLmaxT[l] ** 2
                )

        self.model.line_lim_from = pyo.Constraint(
            self.model.L, rule=line_lim_from_def
        )
        self.model.line_lim_to = pyo.Constraint(
            self.model.L, rule=line_lim_to_def
        )
        self.model.transf_lim1 = pyo.Constraint(
            self.model.TRANSF, rule=transf_lim1_def
        )
        self.model.transf_lim2 = pyo.Constraint(
            self.model.TRANSF, rule=transf_lim2_def
        )
        self.thermal_limit_mode = thermal_limit

        # --- static generation reactive power limits ---
        def static_generation_reactive_power_bounds(model, g):
            model.qsG[g].unfix()
            return model.QsGmin[g], model.qsG[g], model.QsGmax[g]

        self.model.QsG_pyo = pyo.Constraint(
            self.model.sGc, rule=static_generation_reactive_power_bounds
        )
        non_ctrl = [g for g in self.model.sG if g not in set(self.model.sGc)]
        self.model.sGnc = pyo.Set(within=self.model.sG, initialize=non_ctrl)

        for g in self.model.sGnc:
            self.model.qsG[g].fix(self.model.QsG[g])

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

        # The base AC model fixes the slack v to its load-flow value. For a
        # true AC OPF, the slack voltage magnitude should float within
        # [Vmin, Vmax]; the reference angle remains pinned at delta_b0.
        if free_slack_vm:
            for b0 in self.model.b0:
                self.model.v[b0].unfix()

        def v_bounds(model, b):
            return model.Vmin[b], model.v[b], model.Vmax[b]

        self.model.v_pyo = pyo.Constraint(self.model.Bpd, rule=v_bounds)

        # Optional opt-in: pin voltage magnitude at every bus whose nominal
        # voltage matches `hv_bus_kv` to the base-case load-flow value. This
        # was the historical default at 110 kV in German distribution
        # studies; off by default for compatibility with generic OPF
        # benchmarks (e.g. PGLib-OPF).
        if fix_hv_buses:
            fixed_buses = list(
                self.net.bus.index[self.net.bus.vn_kv == hv_bus_kv]
            )
        else:
            fixed_buses = []
        self.model.Bfix = pyo.Set(initialize=fixed_buses)

        def fixed_v_rule(model, b):
            return model.v[b] == float(self.bus_data.loc[b, "v_m"])

        self.model.v_fixed = pyo.Constraint(self.model.Bfix, rule=fixed_v_rule)

        # --- optional branch angle-difference limits ---
        # PowerModels.jl convention: angmin ≤ δ_from − δ_to ≤ angmax.
        # MATPOWER stores these as ANGMIN/ANGMAX columns in mpc.branch (deg);
        # we read them from net.line.angmin_degree / angmax_degree and the
        # transformer equivalent if those columns exist. Branches missing
        # angle data fall back to ±π (effectively unconstrained).
        if angle_limits:
            self._add_branch_angle_limits()

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

        # --- optional Q(P) / Q(U) for PV-type sgens ---
        # Normalise legacy bool to string mode; False/None → skip entirely.
        _pv_mode = "both" if pv_q_control is True else pv_q_control
        if _pv_mode and self.static_generation_data["var_q"] is not None:
            _types = (
                DEFAULT_PV_SGEN_TYPES
                if sgen_types is None
                else tuple(sgen_types)
            )
            pv_qctrl_mask = (
                self.static_generation_data["type"].isin(_types)
                & self.static_generation_data.in_service
                & self.static_generation_data["var_q"].notna()
            )
            pv_qctrl_init = list(
                self.static_generation_data.index[pv_qctrl_mask]
            )
            pv_in_sGc = set(self.model.sGc)
            pv_selected = [g for g in pv_qctrl_init if g in pv_in_sGc]

            # The PV and wind paths both constrain qsG with the same
            # grid-code characteristic, so an sgen claimed by both would get
            # two redundant constraint sets.  This can only happen when
            # sgen_types is widened to include a wind category.  Leave those
            # to the wind path, which owns them, and say so rather than
            # silently building either duplicate or no constraints.
            wind_claimed = set(self.model.WINDc)
            overlap = sorted(set(pv_selected) & wind_claimed)
            if overlap:
                warnings.warn(
                    f"sgens {overlap} match both sgen_types and "
                    f"wind_sgen_types, so they are already Q-controlled by "
                    f"the wind path (model.WINDc). Excluding them from PVc "
                    f"to avoid duplicate constraints on the same qsG; drop "
                    f"the wind categories from sgen_types to silence this.",
                    SgenTypeOverlapWarning,
                    stacklevel=2,
                )
                pv_selected = [g for g in pv_selected if g not in wind_claimed]

            self.model.PVc = pyo.Set(
                within=self.model.sGc, initialize=pv_selected
            )
            if list(self.model.PVc):
                self.model.PV_var_q = pyo.Param(
                    self.model.PVc,
                    initialize=self.static_generation_data["var_q"][
                        self.model.PVc
                    ],
                )
                self.model.PV_p_inst = pyo.Param(
                    self.model.PVc,
                    initialize=self.static_generation_data["p_inst"][
                        self.model.PVc
                    ],
                )
                qc = self.q_limit_parameter

                if _pv_mode in ("qp", "both"):

                    def PV_QP_pos(model, g):
                        v = model.PV_var_q[g]
                        return (
                            model.qsG[g]
                            <= qc.b_qp_max[v] * model.PV_p_inst[g]
                            + qc.m_qp_max[v] * model.psG[g]
                        )

                    def PV_QP_neg(model, g):
                        v = model.PV_var_q[g]
                        return (
                            model.qsG[g]
                            >= qc.b_qp_min[v] * model.PV_p_inst[g]
                            + qc.m_qp_min[v] * model.psG[g]
                        )

                    self.model.PV_QP_pos = pyo.Constraint(
                        self.model.PVc, rule=PV_QP_pos
                    )
                    self.model.PV_QP_neg = pyo.Constraint(
                        self.model.PVc, rule=PV_QP_neg
                    )

                if _pv_mode in ("qu", "both"):
                    sGbs_lookup_pv = {g: b for (g, b) in self.model.sGbs}

                    def PV_QU_min(model, g):
                        if g not in sGbs_lookup_pv:
                            return pyo.Constraint.Skip
                        b = sGbs_lookup_pv[g]
                        v = model.PV_var_q[g]
                        return (
                            model.qsG[g]
                            >= (qc.m_qv[v] * model.v[b] + qc.b_qv_min[v])
                            * model.PV_p_inst[g]
                        )

                    def PV_QU_max(model, g):
                        if g not in sGbs_lookup_pv:
                            return pyo.Constraint.Skip
                        b = sGbs_lookup_pv[g]
                        v = model.PV_var_q[g]
                        return (
                            model.qsG[g]
                            <= (qc.m_qv[v] * model.v[b] + qc.b_qv_max[v])
                            * model.PV_p_inst[g]
                        )

                    self.model.PV_QU_min = pyo.Constraint(
                        self.model.PVc, rule=PV_QU_min
                    )
                    self.model.PV_QU_max = pyo.Constraint(
                        self.model.PVc, rule=PV_QU_max
                    )

        # --- inverter apparent-power circle ---
        if inverter_s2 and "sn_mva" in self.net.sgen:
            conv_sz = (
                self.net.sgen["converter_sizing_pu"].fillna(1.0)
                if "converter_sizing_pu" in self.net.sgen
                else pd.Series(1.0, index=self.net.sgen.index)
            )
            s_inv_pu = self.net.sgen["sn_mva"] * conv_sz / self.baseMVA
            inv_idx = [
                g
                for g in self.model.sGc
                if pd.notna(self.net.sgen.at[g, "sn_mva"])
                and float(self.net.sgen.at[g, "sn_mva"]) > 0
            ]
            if inv_idx:
                self.model.sGinv = pyo.Set(
                    within=self.model.sGc, initialize=inv_idx
                )
                self.model.S_inv = pyo.Param(
                    self.model.sGinv,
                    initialize={g: float(s_inv_pu.at[g]) for g in inv_idx},
                )

                def sgen_inverter_s2_rule(model, g):
                    return (
                        model.psG[g] ** 2 + model.qsG[g] ** 2
                        <= model.S_inv[g] ** 2
                    )

                self.model.sgen_inverter_s2 = pyo.Constraint(
                    self.model.sGinv, rule=sgen_inverter_s2_rule
                )

                # cos(φ) cone: |qsG[g]| ≤ psG[g] · tan(arccos(cos_phi_min))
                # Per-sgen "cos_phi_min" column takes precedence over scalar.
                pf_data = {}
                if "cos_phi_min" in self.net.sgen:
                    for g in inv_idx:
                        val = self.net.sgen.at[g, "cos_phi_min"]
                        if pd.notna(val) and 0 < float(val) <= 1:
                            pf_data[g] = float(np.tan(np.arccos(float(val))))
                elif cos_phi_min is not None:
                    tan_val = float(np.tan(np.arccos(cos_phi_min)))
                    pf_data = {g: tan_val for g in inv_idx}

                if pf_data:
                    self.model.sGpf = pyo.Set(
                        within=self.model.sGinv,
                        initialize=list(pf_data.keys()),
                    )
                    self.model.tan_phi = pyo.Param(
                        self.model.sGpf, initialize=pf_data
                    )

                    def sgen_cos_phi_upper(model, g):
                        return model.qsG[g] <= model.tan_phi[g] * model.psG[g]

                    def sgen_cos_phi_lower(model, g):
                        return model.qsG[g] >= -model.tan_phi[g] * model.psG[g]

                    self.model.sgen_cos_phi_upper = pyo.Constraint(
                        self.model.sGpf, rule=sgen_cos_phi_upper
                    )
                    self.model.sgen_cos_phi_lower = pyo.Constraint(
                        self.model.sGpf, rule=sgen_cos_phi_lower
                    )

        # --- P(U) active-power curtailment (VDE-AR-N 4105 §8.5) ---
        if pu_curtail and "p_inst_mw" in self.net.sgen:
            pv_mask = (self.static_generation_data["type"] == "PV") & (
                self.static_generation_data.in_service
            )
            pu_idx = [
                g
                for g in self.model.sGc
                if g in pv_mask.index and pv_mask.loc[g]
            ]
            if pu_idx:
                p_inst_pu = (
                    self.net.sgen["p_inst_mw"].fillna(
                        self.net.sgen["p_mw"].abs()
                    )
                    / self.baseMVA
                )

                def _v_curtail(g):
                    if "v_curtail_pu" in self.net.sgen:
                        v = self.net.sgen.at[g, "v_curtail_pu"]
                        if pd.notna(v):
                            return float(v)
                    return code.vpu_v_curtail

                def _v_max_curtail(g):
                    if "v_max_curtail_pu" in self.net.sgen:
                        v = self.net.sgen.at[g, "v_max_curtail_pu"]
                        if pd.notna(v):
                            return float(v)
                    return code.vpu_v_max

                self.model.sGpu = pyo.Set(
                    within=self.model.sGc, initialize=pu_idx
                )
                self.model.P_inst_pu = pyo.Param(
                    self.model.sGpu,
                    initialize={g: float(p_inst_pu.at[g]) for g in pu_idx},
                )
                self.model.V_curtail = pyo.Param(
                    self.model.sGpu,
                    initialize={g: _v_curtail(g) for g in pu_idx},
                )
                self.model.V_max_curtail = pyo.Param(
                    self.model.sGpu,
                    initialize={g: _v_max_curtail(g) for g in pu_idx},
                )
                sGbs_lookup_pu = {g: b for (g, b) in self.model.sGbs}

                def sgen_pu_curtail_rule(model, g):
                    if g not in sGbs_lookup_pu:
                        return pyo.Constraint.Skip
                    b = sGbs_lookup_pu[g]
                    dv = model.V_max_curtail[g] - model.V_curtail[g]
                    return model.psG[g] * dv <= model.P_inst_pu[g] * (
                        model.V_max_curtail[g] - model.v[b]
                    )

                self.model.sgen_pu_curtail = pyo.Constraint(
                    self.model.sGpu, rule=sgen_pu_curtail_rule
                )

        # --- fixed cos(φ) equality ---
        # Per-sgen column takes precedence over scalar kwarg.
        _fcf_data: dict = {}
        if "fixed_cos_phi" in self.net.sgen:
            for g in self.model.sGc:
                val = self.net.sgen.at[g, "fixed_cos_phi"]
                if pd.notna(val) and 0 < float(val) <= 1:
                    _fcf_data[g] = float(np.tan(np.arccos(float(val))))
        elif fixed_cos_phi is not None:
            tan_val = float(np.tan(np.arccos(fixed_cos_phi)))
            _fcf_data = {g: tan_val for g in self.model.sGc}

        if _fcf_data:
            self.model.sGfcf = pyo.Set(
                within=self.model.sGc, initialize=list(_fcf_data.keys())
            )
            self.model.fixed_tan_phi = pyo.Param(
                self.model.sGfcf, initialize=_fcf_data
            )

            def sgen_fixed_cos_phi_rule(model, g):
                return model.qsG[g] == model.fixed_tan_phi[g] * model.psG[g]

            self.model.sgen_fixed_cos_phi = pyo.Constraint(
                self.model.sGfcf, rule=sgen_fixed_cos_phi_rule
            )

        # --- cos(φ)(P) profile ---
        if cos_phi_p_profile and "cos_phi_min" in self.net.sgen:
            p_inst_arr = (
                self.net.sgen["p_inst_mw"]
                .fillna(self.net.sgen["p_mw"].abs())
                .values
                if "p_inst_mw" in self.net.sgen
                else self.net.sgen["p_mw"].abs().values
            ) / self.baseMVA

            cpp_thresh_pu_arr = (
                self.net.sgen["cpp_p_threshold_pu"]
                .fillna(code.cpp_p_threshold_pu)
                .values
                if "cpp_p_threshold_pu" in self.net.sgen
                else None
            )

            cpp_data: dict = {}
            for g in self.model.sGc:
                cos_v = self.net.sgen.at[g, "cos_phi_min"]
                if not (pd.notna(cos_v) and 0 < float(cos_v) <= 1):
                    continue
                pn = float(p_inst_arr[g])
                thresh_pu = (
                    float(cpp_thresh_pu_arr[g])
                    if cpp_thresh_pu_arr is not None
                    else code.cpp_p_threshold_pu
                )
                pt = thresh_pu * pn
                if pn <= pt:
                    continue
                cpp_data[g] = {
                    "tan_phi": float(np.tan(np.arccos(float(cos_v)))),
                    "pn": pn,
                    "pt": pt,
                }

            if cpp_data:
                self.model.sGcpp = pyo.Set(
                    within=self.model.sGc, initialize=list(cpp_data.keys())
                )
                self.model.cpp_tan_phi = pyo.Param(
                    self.model.sGcpp,
                    initialize={g: cpp_data[g]["tan_phi"] for g in cpp_data},
                )
                self.model.cpp_Pn = pyo.Param(
                    self.model.sGcpp,
                    initialize={g: cpp_data[g]["pn"] for g in cpp_data},
                )
                self.model.cpp_P_thresh = pyo.Param(
                    self.model.sGcpp,
                    initialize={g: cpp_data[g]["pt"] for g in cpp_data},
                )

                def sgen_cpp_rule(model, g):
                    dPn = model.cpp_Pn[g] - model.cpp_P_thresh[g]
                    return model.qsG[g] * dPn == model.cpp_tan_phi[g] * (
                        model.psG[g] * (model.psG[g] - model.cpp_P_thresh[g])
                    )

                self.model.sgen_cpp = pyo.Constraint(
                    self.model.sGcpp, rule=sgen_cpp_rule
                )

    def _add_branch_angle_limits(self):
        """Attach branch phase-angle-difference constraints to ``self.model``.

        Reads per-line / per-transformer angle bounds from
        ``net.line.angmin_degree`` / ``net.line.angmax_degree`` (and the
        transformer equivalent if present), converts to radians, and adds
        ``angmin_rad ≤ delta[from] − delta[to] ≤ angmax_rad`` on every branch
        that has finite bounds.
        """

        def _bounds(table, idx_set, hv_col, lv_col):
            angmin_col = "angmin_degree"
            angmax_col = "angmax_degree"
            if (
                angmin_col not in table.columns
                or angmax_col not in table.columns
            ):
                return {}
            valid = set(table.index)
            out = {}
            for ix in idx_set:
                if ix not in valid:
                    # synthetic impedance indices live in model.L beyond
                    # net.line — they have no MATPOWER angle bound
                    continue
                amin = float(table.at[ix, angmin_col])
                amax = float(table.at[ix, angmax_col])
                if (
                    not np.isfinite(amin)
                    or not np.isfinite(amax)
                    or abs(amin) >= 359.0
                    or abs(amax) >= 359.0
                ):
                    continue
                out[ix] = (
                    self.bus_lookup[int(table.at[ix, hv_col])],
                    self.bus_lookup[int(table.at[ix, lv_col])],
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

    def add_active_change_objective(self):
        for g in self.model.gG:
            self.model.pG[g].fix(self.model.PG[g])  # back to original dispatch
        for g in self.model.eG:
            self.model.pG[g].unfix()

        for g in self.model.eG:
            self.model.pG[g].unfix()

        sgen_type = self.net.sgen["type"] if "type" in self.net.sgen else None

        self.net.sgen["p_avail_mw"] = (
            self.net.sgen.p_mw * self.net.sgen.scaling
        )

        # convenience
        avail_pu = (self.net.sgen["p_avail_mw"] / self.baseMVA).to_dict()

        # define RES set: wind/solar
        def is_res(g):
            if sgen_type is None:
                return False
            t = str(self.net.sgen.at[g, "type"]).lower()
            return ("wind" in t) or ("solar" in t) or ("pv" in t)

        RES = [g for g in self.model.sG if is_res(g)]

        # dispatchable sgenerators = controllable and not RES
        DISPATCH = [
            g
            for g in self.model.sG
            if (
                "controllable" in self.net.sgen.columns
                and bool(self.net.sgen.at[g, "controllable"])
            )
            and g not in RES
        ]

        # non-controllable (fixed) and not RES
        FIXED = [
            g for g in self.model.sG if g not in RES and g not in DISPATCH
        ]

        # res can only go down?
        for g in RES:
            self.model.psG[g].unfix()
            self.model.psG[g].setlb(0.0)
            self.model.psG[g].setub(avail_pu[g])

        for g in DISPATCH:
            self.model.psG[g].unfix()
            if "min_p_mw" in self.net.sgen.columns:
                self.model.psG[g].setlb(
                    float(self.net.sgen.at[g, "min_p_mw"])
                    * self.net.sgen.at[g, "scaling"]
                    / self.baseMVA
                )
            if "max_p_mw" in self.net.sgen.columns:
                self.model.psG[g].setub(
                    float(self.net.sgen.at[g, "max_p_mw"])
                    * self.net.sgen.at[g, "scaling"]
                    / self.baseMVA
                )

        for g in FIXED:
            self.model.psG[g].fix()

        # loads fixed for redispatch
        for d in self.model.D:
            self.model.pD[d].fix()

        self.model.inj_mismatch_pos = pyo.Var(domain=pyo.NonNegativeReals)
        self.model.inj_mismatch_neg = pyo.Var(domain=pyo.NonNegativeReals)

        def active_change_objective(model):
            disp_term = sum(
                (model.psG[g] - model.PsG[g]) ** 2 for g in DISPATCH
            )
            eps_qg = 1e-6
            # qg_pen = eps_qg * sum(model.qG[g] ** 2 for g in model.G)
            qsg_pen = eps_qg * sum(model.qsG[g] ** 2 for g in model.sG)
            eps_qstor = 1e-6
            qstor_pen = eps_qstor * sum(
                model.qSTOR[s] ** 2 for s in model.STOR
            )
            penalty = 1e4
            soft_term = penalty * (
                model.inj_mismatch_pos + model.inj_mismatch_neg
            )
            eps_stor = 1e-4
            stor_term = eps_stor * sum(
                model.STOR_Pchg[s] + model.STOR_Pdis[s] for s in model.STOR
            )

            return disp_term + soft_term + qstor_pen + stor_term + qsg_pen

        self.model.obj_loading = pyo.Objective(
            rule=active_change_objective, sense=pyo.minimize
        )

        def total_injection_mismatch_rule(model):
            total_ref = sum(model.PsG[g] for g in model.sG) - sum(
                model.STOR_P0[s] for s in model.STOR
            )
            total_now = sum(model.psG[g] for g in model.sG) - sum(
                model.pSTOR[s] for s in model.STOR
            )
            # total_now - total_ref = pos - neg
            return (
                total_now - total_ref
            ) == model.inj_mismatch_pos - model.inj_mismatch_neg

        self.model.total_injection_soft = pyo.Constraint(
            rule=total_injection_mismatch_rule
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
