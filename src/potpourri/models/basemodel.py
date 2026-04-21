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

        self.net = copy.deepcopy(net)
        pp.runpp(self.net, voltage_depend_loads=False)

        # --- pyo.Sets ---
        bus_set = self.net._ppc["bus"][:, [0, 1, 7, 8]]
        bus_set[:, -1] *= pi / 180
        self.bus_data = pd.DataFrame(
            bus_set[:, 1:],
            index=bus_set[:, 0].astype(int),
            columns=["type", "v_m", "v_a_rad"],
        )

        self.bus_lookup = self.net._pd2ppc_lookups["bus"]
        self.demand_set = self.net.load.index[self.net.load.in_service]
        self.shunt_set = self.net.shunt.index[self.net.shunt.in_service]

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

        # --- line ---
        hv_bus = self.net._ppc["branch"][:, 0].real
        lv_bus = self.net._ppc["branch"][:, 1].real
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

        hv_bus_trafo = hv_bus[trafo_start:trafo_end]
        lv_bus_trafo = lv_bus[trafo_start:trafo_end]
        trafo_ind = self.trafo_data.index[self.trafo_data.in_service]
        self.bus_trafo_dict = dict(
            zip(
                list(zip(trafo_ind, [1] * len(trafo_ind)))
                + list(zip(trafo_ind, [2] * len(trafo_ind))),
                np.concatenate(
                    [hv_bus_trafo[trafo_ind], lv_bus_trafo[trafo_ind]]
                ),
            )
        )

    def create_model(self):
        """Create the Pyomo ConcreteModel with sets, parameters, and fixed
        variables."""
        self.model = pyo.ConcreteModel()

        # --- pyo.SetS ---
        self.model.B = pyo.Set(initialize=self.bus_data.index)  # buses
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
            solver (str): The solver to use ('ipopt', 'mindtpy', 'neos', etc.).
            load_solutions (bool): Whether to load solutions into the model
                after solving.
            mip_solver (str): The mixed-integer programming solver for
                'mindtpy'.
            max_iter (int, optional): Maximum iterations for the solver.
            time_limit (int): Time limit for the solver in seconds.
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

            if max_iter:
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
