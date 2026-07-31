"""Single-period Basemodel: maps a pandapower network to a Pyomo ConcreteModel
and solves it."""

import copy
from math import pi

import numpy as np
import pandas as pd
import pandapower as pp
import pyomo.environ as pyo
from loguru import logger

from potpourri.models.pyo_to_net import pyo_sol_to_net_res


class Basemodel:
    """Pyomo-based optimization model for single-period power system analysis.

    Extracts buses, lines, transformers, loads, and generators from a
    pandapower network into Pyomo sets and parameters. Subclasses add
    power-flow equations
    and OPF constraints on top.

    Attributes:
        net: Deep-copied pandapower network (pp.runpp already executed).
        model: Pyomo ConcreteModel populated by create_model().
        results: Solver result object populated by solve().
    """

    def __init__(self, net):
        if not isinstance(net, pp.pandapowerNet):
            raise ValueError("Input network must be a pandapower network.")
        # Make sure bus-to-bus switches are handled correctly by merging them
        self.net = preprocess_grid(copy.deepcopy(net))
        pp.runpp(self.net, voltage_depend_loads=False)

        # --- pyo.Sets ---
        # Every ppc bus, not the first len(net.bus) rows.  pandapower's ppc
        # conversion can add auxiliary buses (switch handling), so the ppc bus
        # table is longer than net.bus on some grids — 103 rows for 97
        # pandapower buses on 1-MV-rural--0-sw.  Truncating took the wrong
        # subset there: it kept auxiliary buses while dropping real pandapower
        # buses together with the in-service branches attached to them, and
        # left `_pd2ppc_lookups` able to resolve to a bus the model did not
        # contain (KeyError when writing results back).  The multi-period
        # Basemodel has always used the full table; this matches it.
        bus_set = self.net._ppc["bus"][:, [0, 1, 7, 8]]
        bus_set[:, -1] *= pi / 180
        self.bus_data = pd.DataFrame(
            bus_set[:, 1:],
            index=bus_set[:, 0].astype(int),
            columns=["type", "v_m", "v_a_rad"],
        )

        self.bus_lookup = self.net._pd2ppc_lookups["bus"]
        # ppc bus number carrying each pandapower bus, and the subset of ppc
        # buses that a pandapower bus maps onto.  Auxiliary ppc buses are
        # absent from the latter: they are internal nodes with no pandapower
        # row, so no user-supplied per-bus data (voltage limits in
        # particular) exists for them.
        self.pd_bus_to_ppc = self.bus_lookup[self.net.bus.index.values]
        self.ppc_buses_with_pd = pd.Index(
            sorted({int(b) for b in self.pd_bus_to_ppc})
        )
        self.demand_set = self.net.load.index[self.net.load.in_service]
        self.shunt_set = self.net.shunt.index[self.net.shunt.in_service]
        self.storage_set = self.net.storage.index[self.net.storage.in_service]

        self.bus_demand_set = list(
            zip(
                self.bus_lookup[self.net.load.bus[self.demand_set].values],
                self.demand_set,
            )
        )
        self.bus_shunt_set = list(
            zip(
                self.bus_lookup[self.net.shunt.bus[self.shunt_set].values],
                self.shunt_set,
            )
        )
        self.bus_storage_set = list(
            zip(
                self.bus_lookup[self.net.storage.bus[self.storage_set].values],
                self.storage_set,
            )
        )

        # --- pyo.Param Data ---
        self.baseMVA = self.net.sn_mva

        self.PD_data = (
            self.net.load.p_mw * self.net.load.scaling / self.baseMVA
        )

        self.GB_data = self.net.shunt.p_mw * self.net.shunt.step / self.baseMVA

        # --- generation ---
        pg = self.net._ppc["gen"][:, 1] / self.baseMVA
        ref_gens = self.net._ppc["internal"]["ref_gens"]
        in_service_gens = self.net._ppc["gen"][:, 7].astype(bool)
        gen_bus = self.net._ppc["gen"][:, 0].astype(int)
        self.generation_data = pd.DataFrame(
            {
                "pg": pg,
                "ref": False,
                "in_service": in_service_gens,
                "bus": gen_bus,
            }
        )
        self.generation_data.loc[ref_gens, "ref"] = True
        gen_bus_tuples = list(enumerate(self.generation_data["bus"]))
        self.generation_data["gen_bus"] = gen_bus_tuples

        # --- static generation ---
        psg = self.net.sgen.p_mw * self.net.sgen.scaling / self.baseMVA
        sgen_bus = self.bus_lookup[self.net.sgen.bus.values]
        self.static_generation_data = pd.DataFrame(
            {
                "p": psg.values,
                "in_service": self.net.sgen.in_service.values,
                "bus": sgen_bus,
            }
        )
        self.static_generation_data["gen_bus"] = list(
            enumerate(self.static_generation_data["bus"])
        )

        self.p_ref = (net.sgen.p_mw * net.sgen.scaling) / self.baseMVA
        self.total = (
            (net.sgen.p_mw * net.sgen.scaling).sum() + net.storage.p_mw.sum()
        ) / self.baseMVA
        self.q_ref = net.sgen.q_mvar * net.sgen.scaling / self.baseMVA

        # --- line ---
        hv_bus = self.net._ppc["branch"][:, 0].real
        lv_bus = self.net._ppc["branch"][:, 1].real
        # hv_bus = self.net._ppc['branch'][:, 0].real
        # lv_bus = self.net._ppc['branch'][:, 1].real
        trafo_start = len(self.net.line.index)
        trafo_end = trafo_start + len(self.net.trafo.index)

        hv_bus_line = hv_bus[:trafo_start]
        lv_bus_line = lv_bus[:trafo_start]
        self.line_data = pd.DataFrame(
            {"in_service": self.net.line.in_service.values}
        )
        line_ind = self.line_data.index[self.line_data.in_service]
        self.bus_line_dict = dict(
            zip(
                list(zip(line_ind, [1] * len(line_ind)))
                + list(zip(line_ind, [2] * len(line_ind))),
                np.concatenate([hv_bus_line[line_ind], lv_bus_line[line_ind]]),
            )
        )
        # Build bus_line_dict directly using pandapower net bus indices.
        # We treat ``net.impedance`` rows as additional "lines" in the model.
        # pandapower puts these in a separate table when ``from_vn_kv !=
        # to_vn_kv`` (and the branch has no off-nominal tap or phase shift),
        # but in per-unit on the system base they obey exactly the same
        # ``y = 1/(r + jx) + j·b/2`` model as a normal line, so we include
        # them in ``model.L`` to honour the full network topology. (Skipping
        # them caused PGLib case118_ieee bus 115 to be isolated and case89
        # to find suboptimal dispatches via 15 missing parallel paths.)
        n_line = len(self.net.line.index)
        line_in_service = self.net.line["in_service"].astype(bool).values
        imp_table = self.net.get("impedance")
        has_impedance = imp_table is not None and not imp_table.empty
        if has_impedance:
            imp_in_service = (
                imp_table["in_service"].astype(bool).values
                if "in_service" in imp_table.columns
                else np.ones(len(imp_table), dtype=bool)
            )
            ext_index = pd.Index(
                list(self.net.line.index)
                + [n_line + i for i in range(len(imp_table))]
            )
            ext_in_service = np.concatenate([line_in_service, imp_in_service])
        else:
            ext_index = self.net.line.index
            ext_in_service = line_in_service
        self.line_data = pd.DataFrame(
            {"in_service": ext_in_service}, index=ext_index
        )
        self.bus_line_dict = {}
        for line_idx in self.line_data.index[self.line_data["in_service"]]:
            if int(line_idx) < n_line:
                fb = int(self.net.line.at[line_idx, "from_bus"])
                tb = int(self.net.line.at[line_idx, "to_bus"])
            else:
                imp_row = int(line_idx) - n_line
                fb = int(imp_table.iloc[imp_row]["from_bus"])
                tb = int(imp_table.iloc[imp_row]["to_bus"])
            self.bus_line_dict[(int(line_idx), 1)] = fb
            self.bus_line_dict[(int(line_idx), 2)] = tb
        self._n_native_lines = n_line  # used by subclasses to slice _ppc

        # --- transformer ---
        shift = (
            self.net._ppc["branch"][trafo_start:trafo_end, 9].real * pi / 180
        )
        tap = self.net._ppc["branch"][trafo_start:trafo_end, 8].real
        self.trafo_data = pd.DataFrame(
            {
                "in_service": self.net.trafo.in_service.values,
                "shift_rad": shift,
                "tap": tap,
            }
        )
        # Same for the bus_trafo_dict
        self.trafo_data = pd.DataFrame(
            {
                "in_service": self.net.trafo["in_service"].astype(bool),
                "shift_rad": shift,
                "tap": tap,
            },
            index=self.net.trafo.index,
        )
        self.bus_trafo_dict = {}
        for trafo_idx in self.trafo_data.index[self.trafo_data["in_service"]]:
            self.bus_trafo_dict[(int(trafo_idx), 1)] = int(
                self.net.trafo.at[trafo_idx, "hv_bus"]
            )
            self.bus_trafo_dict[(int(trafo_idx), 2)] = int(
                self.net.trafo.at[trafo_idx, "lv_bus"]
            )

    def create_model(self):
        """Create the Pyomo ConcreteModel with sets, parameters, and fixed
        variables."""
        self.model = pyo.ConcreteModel()

        # Always declare STOR/STOR_bus; KCL sums are 0 when empty
        stor_idx = self.storage_set.tolist()
        self.model.STOR = pyo.Set(initialize=stor_idx)
        if stor_idx:
            stor_bus_ppc = self.bus_lookup[
                self.net.storage.bus[self.storage_set].values
            ]
            stor_bus_dict = dict(zip(stor_idx, stor_bus_ppc))
        else:
            stor_bus_dict = {}
        self.model.STOR_bus = pyo.Param(
            self.model.STOR, initialize=stor_bus_dict, within=pyo.Any
        )
        if not self.storage_set.empty:
            self.add_storage()

        # --- pyo.SetS ---
        self.model.B = pyo.Set(initialize=self.bus_data.index)  # buses
        # Buses that a pandapower bus maps onto, i.e. everything in B except
        # the auxiliary nodes pandapower inserts for node-node switches.
        # Per-bus user data (voltage limits) exists only for these, so
        # constraints derived from net.bus are indexed over Bpd rather than B.
        # On grids without auxiliary nodes Bpd == B and nothing changes.
        self.model.Bpd = pyo.Set(
            within=self.model.B, initialize=self.ppc_buses_with_pd
        )
        self.model.b0 = pyo.Set(
            initialize=self.bus_data.index[self.bus_data.type == 3],
            within=self.model.B,
        )  # reference buses
        self.model.bPV = pyo.Set(
            initialize=self.bus_data.index[self.bus_data.type == 2],
            within=self.model.B,
        )  # PV buses
        self.model.sG = pyo.Set(
            initialize=self.static_generation_data.index[
                self.static_generation_data.in_service
            ]
        )  # static generators
        self.model.G = pyo.Set(
            initialize=self.generation_data.index[
                self.generation_data.in_service
            ]
        )  # external grids and generators
        self.model.eG = pyo.Set(
            initialize=self.generation_data.index[self.generation_data.ref],
            within=self.model.G,
        )  # external grids and slack generators
        self.model.gG = pyo.Set(
            initialize=self.generation_data.index[~self.generation_data.ref],
            within=self.model.G,
        )  # generators (not static) not slack
        self.model.D = pyo.Set(initialize=self.demand_set)
        self.model.L = pyo.Set(
            initialize=self.line_data.index[self.line_data.in_service]
        )
        self.model.SHUNT = pyo.Set(initialize=self.shunt_set)
        self.model.LE = pyo.Set(initialize=[1, 2])
        self.model.TRANSF = pyo.Set(
            initialize=self.trafo_data.index[self.trafo_data.in_service]
        )

        # generators, buses, loads linked to each bus b
        self.model.Dbs = pyo.Set(
            within=self.model.B * self.model.D, initialize=self.bus_demand_set
        )  # pyo.Set of demand-bus mapping
        self.model.SHUNTbs = pyo.Set(
            within=self.model.B * self.model.SHUNT,
            initialize=self.bus_shunt_set,
        )  # pyo.Set of shunt-bus mapping
        self.model.Gbs = pyo.Set(
            within=self.model.G * self.model.B,
            initialize=self.generation_data["gen_bus"][self.model.G],
        )
        self.model.sGbs = pyo.Set(
            within=self.model.sG * self.model.B,
            initialize=self.static_generation_data["gen_bus"][self.model.sG],
        )

        # --- pyo.Parameters ---
        # line and trafo matrix
        self.model.A = pyo.Param(
            self.model.L * self.model.LE, initialize=self.bus_line_dict
        )  # bus-line matrix
        self.model.AT = pyo.Param(
            self.model.TRANSF * self.model.LE, initialize=self.bus_trafo_dict
        )  # bus-transformer matrix

        # generation
        self.model.PsG = pyo.Param(
            self.model.sG,
            initialize=self.static_generation_data.p[self.model.sG],
        )
        self.model.PG = pyo.Param(
            self.model.G, initialize=self.generation_data.pg[self.model.G]
        )

        # demand
        self.model.PD = pyo.Param(
            self.model.D, initialize=self.PD_data[self.model.D]
        )

        # shunt
        self.model.GB = pyo.Param(
            self.model.SHUNT,
            within=pyo.Reals,
            initialize=self.GB_data[self.model.SHUNT],
        )  # shunt conductance

        # trafo
        self.model.shift = pyo.Param(
            self.model.TRANSF,
            within=pyo.Reals,
            initialize=self.trafo_data.shift_rad[self.model.TRANSF],
        )  # transformer phase shift in rad

        # external grid voltage angle
        self.model.delta_b0 = pyo.Param(
            self.model.b0,
            within=pyo.Reals,
            initialize=self.bus_data.v_a_rad[self.model.b0],
        )

        # baseMVA of the net
        self.model.baseMVA = pyo.Param(
            within=pyo.NonNegativeReals, initialize=self.baseMVA
        )

        # --- Variables ---
        delta_init = self.bus_data.v_a_rad.fillna(0.0).to_dict()
        self.model.delta = pyo.Var(
            self.model.B,
            domain=pyo.Reals,
            initialize=delta_init,
            bounds=(-pi, pi),
        )  # voltage phase angle at bus b, rad
        self.model.pD = pyo.Var(
            self.model.D, domain=pyo.Reals
        )  # real power demand delivered
        self.model.psG = pyo.Var(
            self.model.sG, domain=pyo.NonNegativeReals
        )  # real static generator power
        self.model.pG = pyo.Var(
            self.model.G, domain=pyo.Reals, initialize=self.model.PG
        )  # real power injection from static generators
        self.model.pLfrom = pyo.Var(
            self.model.L, domain=pyo.Reals
        )  # real power injected at b onto line
        self.model.pLto = pyo.Var(
            self.model.L, domain=pyo.Reals
        )  # real power injected at b' onto line
        self.model.pThv = pyo.Var(
            self.model.TRANSF, domain=pyo.Reals
        )  # real power injected at b onto transformer
        self.model.pTlv = pyo.Var(
            self.model.TRANSF, domain=pyo.Reals
        )  # real power injected at b' onto transformer
        self.model.Tap = pyo.Var(
            self.model.TRANSF,
            domain=pyo.Reals,
            initialize=self.trafo_data.tap[self.model.TRANSF],
        )  # transformer tap ratio

        # transformer tap ratio
        for t in self.model.TRANSF:
            self.model.Tap[t].fix()

        # --- generator power ---
        for g in self.model.sG:
            self.model.psG[g].fix(self.model.PsG[g])
        for g in self.model.gG:
            self.model.pG[g].fix(self.model.PG[g])

        # --- demand ---
        for d in self.model.D:
            self.model.pD[d].fix(self.model.PD[d])

        # --- reference bus constraint ---
        for b in self.model.b0:
            self.model.delta[b].fix(self.model.delta_b0[b])

    def solve(
        self,
        to_net: bool = True,
        print_solver_output: bool = False,
        solver="ipopt",
        load_solutions: bool = True,
        mip_solver="gurobi",
        max_iter=None,
        time_limit=600,
        init_strategy="rNLP",
        neos_opt="ipopt",
        nlp_solver_args=None,
    ):
        """
        Solves the optimization model using the specified solver.

        Args:
            to_net (bool): Whether to map results back to the pandapower
                network.
            print_solver_output (bool): Whether to print solver output.
            solver (str): The solver to use ('ipopt', 'mindtpy', 'neos',
                'gurobi_direct_minlp', etc.). Use 'gurobi_direct_minlp' to
                send the model — including the polar-form sin/cos power flow
                and any integer variables — straight to Gurobi's global
                spatial branch-and-bound. Requires Pyomo >= 6.10 and
                gurobipy >= 12.
            load_solutions (bool): Whether to load solutions into the model
                after solving.
            mip_solver (str): The mixed-integer programming solver for
                'mindtpy'.
            max_iter (int, optional): Maximum iterations for the solver.
                Mapped to Gurobi's 'IterationLimit' for 'gurobi*' solvers.
            time_limit (int): Time limit for the solver in seconds. Honoured
                by 'mindtpy' and by 'gurobi*' (as 'TimeLimit'); other solvers
                ignore it.
            init_strategy (str): Initialization strategy for 'mindtpy'.
            neos_opt (str): Solver to use with NEOS.
            nlp_solver_args (dict, optional): Extra keyword arguments forwarded
                to the NLP sub-solver when using 'mindtpy' (e.g.
                ``{"max_iter": 10000}`` to raise IPOPT's iteration cap).

        Raises:
            ValueError: If solver settings are invalid or the solver fails.
        """
        logger.info(
            "Solving model '{}' with solver '{}'", self.model.name, solver
        )

        if solver == "mindtpy":
            optimizer = pyo.SolverFactory(solver)
            if not max_iter:
                max_iter = 50

            if mip_solver == "gurobi":
                mip_solver = "gurobi_persistent"

            logger.debug(
                "mindtpy: mip_solver={}, nlp_solver=ipopt, "
                "max_iter={}, init_strategy={}",
                mip_solver,
                max_iter,
                init_strategy,
            )
            try:
                self.results = optimizer.solve(
                    self.model,
                    mip_solver=mip_solver,
                    nlp_solver="ipopt",
                    nlp_solver_args=nlp_solver_args or {},
                    tee=print_solver_output,
                    iteration_limit=max_iter,
                    time_limit=time_limit,
                    init_strategy=init_strategy,
                )
            except ValueError as err:
                logger.error("mindtpy solver error: {}", err)

        elif solver == "neos":
            logger.info("Submitting model to NEOS server (opt={})", neos_opt)
            solver_manager = pyo.SolverManagerFactory("neos")
            self.results = solver_manager.solve(
                self.model, opt=neos_opt, tee=True
            )
        else:
            optimizer = pyo.SolverFactory(solver)

            if solver.startswith("gurobi"):
                # Gurobi uses its own option names; 'max_iter' is IPOPT's and
                # raises GurobiError("Unknown parameter"). A time limit matters
                # here because 'gurobi_direct_minlp' runs a global spatial
                # branch-and-bound that is otherwise unbounded on a nonconvex
                # AC OPF.
                if max_iter:
                    optimizer.options["IterationLimit"] = max_iter
                    logger.debug("Gurobi IterationLimit set to {}", max_iter)
                if time_limit:
                    optimizer.options["TimeLimit"] = time_limit
                    logger.debug("Gurobi TimeLimit set to {} s", time_limit)
            elif max_iter:
                optimizer.options["max_iter"] = max_iter
                logger.debug("Solver max_iter set to {}", max_iter)

            self.results = optimizer.solve(
                self.model,
                load_solutions=load_solutions,
                tee=print_solver_output,
            )

        try:
            if pyo.check_optimal_termination(self.results):
                logger.info(
                    "Optimal solution found for model '{}'", self.model.name
                )
                if to_net:
                    pyo_sol_to_net_res(self.net, self.model)
                    logger.debug("Solution mapped to net.res_*")
            else:
                logger.warning(
                    "Solver did not reach optimal termination for model '{}'"
                    " (condition: {})",
                    self.model.name,
                    self.results.solver.termination_condition,
                )
        except AttributeError as err:
            logger.error("Could not check termination condition: {}", err)

        return self.results

    def change_vals(self, key, value):
        """Set all indices of a named Pyomo component to value."""
        component = self.model.component(key)
        if not component:
            logger.warning(
                "Model '{}' has no component '{}'", self.model.name, key
            )
            return
        try:
            for index in component:
                component[index] = value
        except TypeError as err:
            logger.error("change_vals failed for component '{}': {}", key, err)

    def fix_vars(self, key, value=None):
        """Fix all indices of a named Pyomo variable; optionally set to value
        first."""
        component = self.model.component(key)
        if not component:
            logger.warning(
                "Model '{}' has no component '{}'", self.model.name, key
            )
            return
        try:
            for index in component:
                if value is not None:
                    component[index].fix(value)
                else:
                    component[index].fix()
            logger.debug(
                "Fixed variable '{}' in model '{}'", key, self.model.name
            )
        except AttributeError as err:
            logger.error("fix_vars failed for component '{}': {}", key, err)

    def unfix_vars(self, key, value=None):
        """Unfix all indices of a named Pyomo variable; optionally reset to
        value."""
        component = self.model.component(key)
        if not component:
            logger.warning(
                "Model '{}' has no component '{}'", self.model.name, key
            )
            return
        try:
            for index in component:
                component[index].unfix()
                if value is not None:
                    component[index] = value
            logger.debug(
                "Unfixed variable '{}' in model '{}'", key, self.model.name
            )
        except AttributeError as err:
            logger.error("unfix_vars failed for component '{}': {}", key, err)

    def add_storage(self):
        """Add storage parameters, variables, and constraints.

        STOR set and STOR_bus param are already initialised in
        :meth:`create_model`. Assumes ``net.storage`` has: ``sn_mva``,
        ``p_mw``, ``scaling``, ``max_e_mwh``, ``efficiency_percent`` (0–100),
        ``soc_percent`` (0–100).

        Modelling notes (D12 in the formulation audit):

        * **Symmetric round-trip efficiency**: charge gain and discharge loss
          share the same parameter, modelled as ``η·Pchg − Pdis/η``. This
          corresponds to ``η_chg = η_dis = √η_roundtrip``. PowerModels.jl
          supports independent charge/discharge efficiencies; this model
          does not.
        * **Convex relaxation of complementarity**: the true
          no-simultaneous-charge-discharge constraint ``Pchg·Pdis = 0`` is
          relaxed to ``Pchg + Pdis ≤ Pmax``. This is *stricter* than
          independent ``Pchg, Pdis ∈ [0, Pmax]`` bounds (no double-use of
          the inverter), but *looser* than complementarity (the optimiser
          can still split a small amount of power between both legs to skim
          the round-trip loss). For exact behaviour, use a MILP with a
          binary direction variable or set ``η = 1``.
        * **Single-period energy balance**: the SOC update uses the
          time-step ``STOR_dt`` (15 min by default) and the initial SOC
          ``SOC0``; this is a single-shot snapshot, not a multi-period
          balance.
        """
        if "soc_percent" not in self.net.storage.columns:
            self.net.storage["soc_percent"] = 50.0
        else:
            self.net.storage["soc_percent"] = self.net.storage[
                "soc_percent"
            ].fillna(50.0)

        self.storage_data = copy.deepcopy(
            self.net.storage.loc[self.net.storage.in_service]
        )
        stor_idx = self.storage_data.index.tolist()

        # --- Parameters ---
        self.model.STOR_Pmax = pyo.Param(
            self.model.STOR,
            initialize=(
                (self.storage_data.sn_mva * self.storage_data.scaling)
                / self.baseMVA
            ).to_dict(),
        )
        self.model.STOR_P0 = pyo.Param(
            self.model.STOR,
            initialize=(
                (self.storage_data.p_mw * self.storage_data.scaling)
                / self.baseMVA
            ).to_dict(),
        )
        self.model.STOR_Emax = pyo.Param(
            self.model.STOR,
            initialize=(self.storage_data.max_e_mwh / self.baseMVA).to_dict(),
        )
        # efficiency as fraction (pandapower stores 0–100)
        self.model.STOR_eff = pyo.Param(
            self.model.STOR,
            initialize=(
                self.storage_data.efficiency_percent / 100.0
            ).to_dict(),
        )
        # initial SOC as fraction 0–1 (pandapower stores 0–100 in soc_percent)
        self.model.STOR_SOC0 = pyo.Param(
            self.model.STOR,
            initialize=(self.storage_data.soc_percent / 100.0).to_dict(),
        )
        self.model.STOR_dt = pyo.Param(
            within=pyo.NonNegativeReals, initialize=0.25
        )
        self.model.STOR_SOCmin = pyo.Param(
            self.model.STOR, initialize={s: 0.0 for s in stor_idx}
        )
        self.model.STOR_SOCmax = pyo.Param(
            self.model.STOR, initialize={s: 1.0 for s in stor_idx}
        )

        # --- Variables ---
        self.model.STOR_Pchg = pyo.Var(
            self.model.STOR, domain=pyo.NonNegativeReals
        )
        self.model.STOR_Pdis = pyo.Var(
            self.model.STOR, domain=pyo.NonNegativeReals
        )
        self.model.STOR_SOC = pyo.Var(
            self.model.STOR, domain=pyo.NonNegativeReals
        )
        self.model.qSTOR = pyo.Var(self.model.STOR, domain=pyo.Reals)

        # pSTOR = Pchg - Pdis; negative when discharging → KCL sign: -pSTOR > 0
        def stor_injection_rule(model, s):
            return model.STOR_Pchg[s] - model.STOR_Pdis[s]

        self.model.pSTOR = pyo.Expression(
            self.model.STOR, rule=stor_injection_rule
        )

        # --- Constraints ---
        def stor_chg_limit_rule(model, s):
            return model.STOR_Pchg[s] <= model.STOR_Pmax[s]

        def stor_dis_limit_rule(model, s):
            return model.STOR_Pdis[s] <= model.STOR_Pmax[s]

        self.model.stor_chg_limit = pyo.Constraint(
            self.model.STOR, rule=stor_chg_limit_rule
        )
        self.model.stor_dis_limit = pyo.Constraint(
            self.model.STOR, rule=stor_dis_limit_rule
        )

        # SOC update: single-period energy balance; eff is a fraction (0–1)
        def stor_soc_update_rule(model, s):
            return model.STOR_SOC[s] == (
                model.STOR_SOC0[s]
                + model.STOR_dt
                * (
                    model.STOR_eff[s] * model.STOR_Pchg[s]
                    - model.STOR_Pdis[s] / model.STOR_eff[s]
                )
                / model.STOR_Emax[s]
            )

        self.model.stor_soc_update = pyo.Constraint(
            self.model.STOR, rule=stor_soc_update_rule
        )

        def stor_soc_bounds_rule(model, s):
            return (
                model.STOR_SOCmin[s],
                model.STOR_SOC[s],
                model.STOR_SOCmax[s],
            )

        self.model.stor_soc_bounds = pyo.Constraint(
            self.model.STOR, rule=stor_soc_bounds_rule
        )

        # Convex relaxation of no-simultaneous-charge-discharge
        def stor_no_simul_rule(model, s):
            return (
                model.STOR_Pchg[s] + model.STOR_Pdis[s] <= model.STOR_Pmax[s]
            )

        self.model.stor_no_simul = pyo.Constraint(
            self.model.STOR, rule=stor_no_simul_rule
        )

        # Inverter apparent power limit
        def stor_inverter_cap_rule(model, s):
            return (
                model.pSTOR[s] ** 2 + model.qSTOR[s] ** 2
                <= model.STOR_Pmax[s] ** 2
            )

        self.model.stor_inverter_cap = pyo.Constraint(
            self.model.STOR, rule=stor_inverter_cap_rule
        )


_ELEMENT_BUS_COLUMNS = {
    "line": ("from_bus", "to_bus"),
    "trafo": ("hv_bus", "lv_bus"),
    "trafo3w": ("hv_bus", "mv_bus", "lv_bus"),
    "impedance": ("from_bus", "to_bus"),
    "dcline": ("from_bus", "to_bus"),
    "load": ("bus",),
    "sgen": ("bus",),
    "gen": ("bus",),
    "ext_grid": ("bus",),
    "shunt": ("bus",),
    "storage": ("bus",),
    "ward": ("bus",),
    "xward": ("bus",),
    "motor": ("bus",),
    "asymmetric_load": ("bus",),
    "asymmetric_sgen": ("bus",),
}


def preprocess_grid(grid):
    """Merge zero-impedance bus-bus switches; drop fully orphan buses.

    Updates every element table that carries a bus reference (lines, trafos,
    trafo3w, impedance, dcline, load, sgen, gen, ext_grid, shunt, storage,
    ward/xward, motor, asymmetric loads/sgens) when remapping bus ids, so
    shunt-only / storage-only / trafo-only buses stay consistent after the
    merge. The "referenced buses" computation also includes those tables, so
    a bus that is only connected via, e.g., a shunt or a transformer terminal
    is no longer silently dropped before ``pp.create_continuous_bus_index``.
    """
    grid = copy.deepcopy(grid)

    def _present(name):
        return name in grid and hasattr(grid[name], "loc") and len(grid[name])

    # Iterate over closed bus-bus switches with zero impedance and merge the
    # connected buses across every bus-referencing table.
    for sw_idx, sw in grid.switch.iterrows():
        if (
            sw["closed"]
            and sw["et"] == "b"
            and float(sw.get("z_ohm", 0.0)) == 0.0
        ):
            keep_bus = int(sw["bus"])
            remove_bus = int(sw["element"])
            for name, cols in _ELEMENT_BUS_COLUMNS.items():
                if not _present(name):
                    continue
                table = grid[name]
                for col in cols:
                    if col in table.columns:
                        mask = table[col] == remove_bus
                        if mask.any():
                            table.loc[mask, col] = keep_bus
            grid.switch.drop(sw_idx, inplace=True)

    # Remove self-loop branches/transformers created by the merge.
    if _present("line"):
        self_loop = grid.line.index[
            grid.line["from_bus"] == grid.line["to_bus"]
        ]
        grid.line.drop(self_loop, inplace=True)
    if _present("trafo"):
        self_loop_t = grid.trafo.index[
            grid.trafo["hv_bus"] == grid.trafo["lv_bus"]
        ]
        grid.trafo.drop(self_loop_t, inplace=True)

    for key in list(grid.keys()):
        if key.startswith("res_") and hasattr(grid[key], "drop"):
            grid[key].drop(grid[key].index, inplace=True)
    if "bus_geodata" in grid and len(grid["bus_geodata"]):
        grid["bus_geodata"] = grid["bus_geodata"].loc[
            grid["bus_geodata"].index.intersection(grid.bus.index)
        ]
    # Drop remaining bus-bus switches before reindexing — open bus-bus switches
    # can reference buses that no longer exist in the topology after merging
    # closed ones, causing a KeyError in create_continuous_bus_index.
    if _present("switch"):
        bb = grid.switch.index[grid.switch["et"] == "b"]
        grid.switch.drop(bb, inplace=True)
    pp.create_continuous_bus_index(grid, start=0)

    return grid
