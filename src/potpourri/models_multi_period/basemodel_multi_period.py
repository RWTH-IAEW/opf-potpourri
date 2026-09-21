# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Multi-period Basemodel: extends Basemodel with a time dimension and
simbench profile integration."""

import pandas as pd
from pyomo.environ import *
from math import pi
import copy
import numpy as np
import pandapower as pp
import simbench as sb
import time as ctime
from loguru import logger
from potpourri.models_multi_period.init_pyo_from_pp_res_multi_period import (
    init_pyo_from_pp_res_multi_period,
)
from potpourri.models_multi_period.pyo_to_net_multi_period import (
    pyo_sol_to_net_res,
)
from potpourri.technologies.flexibility import Flexibility_multi_period
from potpourri.technologies.generator import Generator_multi_period
from potpourri.technologies.shunts import Shunts_multi_period
from potpourri.technologies.sgens import Sgens_multi_period
from potpourri.technologies.demand import Demand_multi_period
from potpourri.technologies.windpower import Windpower_multi_period


class Basemodel_multi_period:
    """Multi-period base model that adds a time set and simbench profile
    support to a pandapower network.

    Attributes:
        net: Deep-copied pandapower network with profiles attached.
        model: Pyomo ConcreteModel, built by create_model().
        flexibilities: List of technology mix-in objects
            (Demand, Sgens, Generator, …).
        fromT: First time-step index (inclusive).
        toT: Last time-step index (exclusive).
        T: Number of time steps.
        deltaT: Time-step length in hours (default 0.25 = 15 min).
    """

    def __init__(self, net, toT, fromT=None, pf=1):
        """Initialise the multi-period base model.

        Args:
            net: pandapower network with simbench profiles.
            toT: Model runs up to (exclusive) this time-step index.
            fromT: Model starts from this time-step index (default 0).
            pf: Power factor for reactive power calculation (default 1).
        """
        self.net = copy.deepcopy(net)
        pp.runpp(self.net, voltage_depend_loads=False)

        # --- initialize flexibility list ---
        self.flexibilities = []

        # --- Sets ---
        bus_set = self.net._ppc["bus"][:, [0, 1, 7, 8]]
        bus_set[:, -1] *= pi / 180
        self.bus_data = pd.DataFrame(
            bus_set[:, 1:],
            index=bus_set[:, 0].astype(int),
            columns=["type", "v_m", "v_a_rad"],
        )
        self.bus_lookup = self.net._pd2ppc_lookups["bus"]
        # ppc bus number carrying each pandapower bus, and the subset of ppc
        # buses that a pandapower bus maps onto.  Auxiliary ppc buses added by
        # pandapower's switch handling are absent from the latter: they are
        # internal nodes with no pandapower row, so no user-supplied per-bus
        # data (voltage limits in particular) exists for them.
        self.pd_bus_to_ppc = self.bus_lookup[self.net.bus.index.values]
        self.ppc_buses_with_pd = pd.Index(
            sorted({int(b) for b in self.pd_bus_to_ppc})
        )

        # --- Param Data ---
        self.baseMVA = self.net.sn_mva

        # --- line (+ pandapower impedance branches) ---
        # Mirror of the single-period Basemodel fix: net.impedance carries
        # branches whose from/to vn_kv differ but with no off-nominal tap; in
        # _ppc['branch'] they live in the slice [trafo_end : trafo_end+n_imp]
        # and are otherwise indistinguishable from lines in per-unit. We
        # include them in model.L using synthetic indices >= len(net.line).
        hv_bus = self.net._ppc["branch"][:, 0].real
        lv_bus = self.net._ppc["branch"][:, 1].real
        trafo_start = len(self.net.line.index)
        trafo_end = trafo_start + len(self.net.trafo.index)

        n_line = trafo_start
        imp_table = self.net.get("impedance")
        n_imp = (
            len(imp_table)
            if imp_table is not None and not imp_table.empty
            else 0
        )

        # synthetic indices for impedance rows: n_line, n_line+1, ...
        line_indices = list(self.net.line.index) + [
            n_line + i for i in range(n_imp)
        ]
        line_in_service = np.concatenate(
            [
                self.net.line.in_service.astype(bool).values,
                (
                    imp_table["in_service"].astype(bool).values
                    if n_imp and "in_service" in imp_table.columns
                    else np.ones(n_imp, dtype=bool)
                ),
            ]
        )
        self.line_data = pd.DataFrame(
            {"in_service": line_in_service}, index=line_indices
        )
        # Bus references: native lines from _ppc[:n_line], impedance from
        # _ppc[trafo_end : trafo_end + n_imp].
        if n_imp:
            line_idx_ppc = np.r_[
                np.arange(0, n_line),
                np.arange(trafo_end, trafo_end + n_imp),
            ]
        else:
            line_idx_ppc = np.arange(0, n_line)
        hv_bus_line = hv_bus[line_idx_ppc]
        lv_bus_line = lv_bus[line_idx_ppc]
        line_ind = self.line_data.index[self.line_data.in_service]
        # Map line_ind -> position in line_idx_ppc array
        pos = {idx: p for p, idx in enumerate(self.line_data.index)}
        self.bus_line_dict = {}
        for li in line_ind:
            p = pos[li]
            self.bus_line_dict[(int(li), 1)] = int(hv_bus_line[p])
            self.bus_line_dict[(int(li), 2)] = int(lv_bus_line[p])
        self._n_native_lines = n_line

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

        self.toT = toT
        if fromT is not None:
            self.fromT = fromT
        else:
            self.fromT = 0
        self.T = toT - fromT if fromT else toT
        self.pf = pf  # powerfactor for reactive power calculation
        if toT is not None and not hasattr(self.net, "profiles"):
            raise ValueError(
                "The net object does not have profiles. Please provide a net"
                " object with profiles."
            )
        self.profiles = sb.get_absolute_values(
            self.net, profiles_instead_of_study_cases=True
        )
        # Reactive sgen profiles: present in SimBench data when the network
        # includes a q_mvar column; otherwise derived below from pf.

        # Snapshot the renewables profiles before clearing net.profiles,
        # since PV_multi_period reads them from net.pv_load_profiles later.
        self.net.pv_load_profiles = self.net.profiles["renewables"].iloc[
            self.fromT : self.toT
        ]

        # Clear and refill net.profiles with the [fromT, toT) slice so that
        # all technology modules see the correct time window consistently.
        self.net.profiles.clear()
        # cut off profiles for T, depending on fromT and toT
        if self.fromT is None:
            for profile in self.profiles.keys():
                self.net.profiles[profile] = self.profiles[profile].iloc[
                    : self.toT
                ]  # cut off profiles to T
        elif self.fromT < 0:
            raise ValueError("fromT must be positive integer")
        elif self.toT < 0:
            raise ValueError("toT must be positive integer")
        elif self.toT > len(self.profiles[list(self.profiles.keys())[0]]):
            raise ValueError(
                "toT must be smaller than the length of the profiles"
            )
        elif self.toT < self.fromT:
            raise ValueError("toT must be greater than fromT")
        elif self.fromT is not None:
            for profile in self.profiles.keys():
                self.net.profiles[profile] = self.profiles[profile].iloc[
                    self.fromT : self.toT
                ]

        # print("Model runs from time step " + str(self.fromT+1)
        #       + " to time step " + str(self.toT))

        # calculation of reactive power for static generators with power factor
        if ("sgen", "q_mvar") not in self.net.profiles.keys():
            self.calc_reactive_sgen_power(self.pf)

        # --- create flexibility object and append to list ---

        self.flexibilities.append(Demand_multi_period(self.net))
        self.flexibilities.append(Shunts_multi_period(self.net))
        self.flexibilities.append(Sgens_multi_period(self.net))
        self.flexibilities.append(Generator_multi_period(self.net))

        if "windpot_p_mw" in self.net.bus:
            self.flexibilities.append(Windpower_multi_period(self.net))

    def calc_reactive_sgen_power(self, pf=1):
        """Calculate reactive power profiles for static generators from active
        power and power factor."""
        self.net.profiles[("sgen", "q_mvar")] = self.net.profiles[
            ("sgen", "p_mw")
        ] * np.tan(np.arccos(pf))

    def create_model(self):
        """Create the multi-period Pyomo ConcreteModel with time sets,
        bus/line data, and base variables."""
        logger.info("Creating model at {}", ctime.ctime())
        self.model = ConcreteModel()

        # time dependency to model

        self.model.T = Set(
            initialize=range(self.fromT, self.toT), ordered=True
        )  # time periods
        self.deltaT = 1 / 4  # time step length in hours = 15 minutes
        self.model.deltaT = Param(
            initialize=self.deltaT, within=PositiveReals
        )  # time step length

        # --- iterate through flexibility list ---
        for flex in self.flexibilities:
            flex.get_all(self.model)

        # --- SETS ---
        self.model.b0 = Set(
            initialize=self.bus_data.index[self.bus_data.type == 3],
            within=self.model.B,
        )  # reference buses
        self.model.bPV = Set(
            initialize=self.bus_data.index[self.bus_data.type == 2],
            within=self.model.B,
        )  # PV buses

        self.model.L = Set(
            initialize=self.line_data.index[self.line_data.in_service]
        )
        self.model.LE = Set(initialize=[1, 2])
        self.model.TRANSF = Set(
            initialize=self.trafo_data.index[self.trafo_data.in_service]
        )

        # --- parameters ---
        # line and trafo matrix
        self.model.A = Param(
            self.model.L * self.model.LE, initialize=self.bus_line_dict
        )  # bus-line matrix
        self.model.AT = Param(
            self.model.TRANSF * self.model.LE, initialize=self.bus_trafo_dict
        )  # bus-transformer matrix

        # trafo
        self.model.shift = Param(
            self.model.TRANSF,
            within=Reals,
            initialize=self.trafo_data.shift_rad[self.model.TRANSF],
        )  # transformer phase shift in rad DONE:remove Time dependency

        # external grid voltage angle DONE: remove time depdendency
        self.model.delta_b0 = Param(
            self.model.b0,
            within=Reals,
            initialize=self.bus_data.v_a_rad[self.model.b0],
        )

        # baseMVA of the net
        self.model.baseMVA = Param(
            within=NonNegativeReals, initialize=self.baseMVA
        )

        # --- variables ---
        # Done stay multiperiod
        self.delta_data_dict, self.delta_tuple = self.make_to_dict(
            self.model.B, self.model.T, 0.0, False
        )  # False or true?
        self.pLfrom_tuple = self.make_to_tuple(self.model.L, self.model.T)
        self.pLto_tuple = self.make_to_tuple(self.model.L, self.model.T)
        self.pThv_tuple = self.make_to_tuple(self.model.TRANSF, self.model.T)
        self.pTlv_tuple = self.make_to_tuple(self.model.TRANSF, self.model.T)
        self.Tap_data_dict, self.Tap_tuple = self.make_to_dict(
            self.model.TRANSF,
            self.model.T,
            self.trafo_data.tap[self.model.TRANSF],
            False,
        )  # False or true on time dpendency?

        self.model.delta = Var(
            self.delta_tuple,
            domain=Reals,
            initialize=self.delta_data_dict,
            bounds=(-pi, pi),
        )  # voltage phase angle at bus b, rad
        self.model.pLfrom = Var(
            self.pLfrom_tuple, domain=Reals
        )  # real power injected at b onto line
        self.model.pLto = Var(
            self.pLto_tuple, domain=Reals
        )  # real power injected at b' onto line
        self.model.pThv = Var(
            self.pThv_tuple, domain=Reals
        )  # real power injected at b onto transformer
        self.model.pTlv = Var(
            self.pTlv_tuple, domain=Reals
        )  # real power injected at b' onto transformer
        self.model.Tap = Var(
            self.Tap_tuple, domain=Reals, initialize=self.Tap_data_dict
        )  # transformer tap ratio

        # transformer tap ratio
        for tr in self.model.TRANSF:
            for t in self.model.T:
                self.model.Tap[tr, t].fix()

        # --- reference bus constraint ---
        for b in self.model.b0:
            self.model.delta[b, t].fix(self.model.delta_b0[b])

    def solve(
        self,
        to_net: bool = True,
        print_solver_output: bool = True,
        solver="ipopt",
        load_solutions: bool = True,
        mip_solver="gurobi",
        max_iter=None,
        time_limit=600,
        init_strategy="rNLP",
        neos_opt="bonmin",
        warm_start=True,
    ):
        """Solve the multi-period OPF model with the specified solver.

        Args:
            to_net: Which time step to write into ``net.res_*``. The result
                tables have no time dimension, so exactly one step can be
                mapped. ``True`` maps the **last** step of the horizon; pass
                an int to choose a step from ``model.T``; ``False`` leaves
                ``net.res_*`` untouched. Use :meth:`map_to_net` to map a
                different step afterwards.
            warm_start: Seed every state variable from a per-step pandapower
                power flow before solving (see :meth:`warm_start_from_pf`).
                On by default: a cold start begins with all angles and branch
                flows at zero, which violates the nodal balance everywhere, and
                IPOPT can fail to recover — reporting a locally infeasible
                point on a model that is demonstrably feasible. Set ``False``
                to keep whatever initial values the variables already carry,
                e.g. to re-solve from a previous solution.
            print_solver_output: Whether to stream solver output.
            solver: Solver name ('ipopt', 'mindtpy', 'neos',
                'gurobi_direct_minlp', etc.). 'gurobi_direct_minlp' sends the
                model — including the polar-form sin/cos power flow and any
                integer variables — straight to Gurobi's global spatial
                branch-and-bound. Requires Pyomo >= 6.10 and gurobipy >= 12.
            load_solutions: Whether Pyomo should load the solution after
                solving.
            mip_solver: MIP sub-solver for mindtpy.
            max_iter: Maximum solver iterations. Mapped to Gurobi's
                'IterationLimit' for 'gurobi*' solvers.
            time_limit: Wall-clock time limit in seconds. Honoured by
                'mindtpy' and by 'gurobi*' (as 'TimeLimit'); other solvers
                ignore it.
            init_strategy: Initialization strategy for mindtpy.
            neos_opt: Solver name to use with NEOS.
        """
        if warm_start:
            self.warm_start_from_pf()

        logger.info("Solving model with solver '{}'", solver)
        optimizer = SolverFactory(solver)

        if solver == "mindtpy":
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
                    tee=print_solver_output,
                    iteration_limit=max_iter,
                    time_limit=time_limit,
                    init_strategy=init_strategy,
                )
            except ValueError as err:
                logger.error("mindtpy solver error: {}", err)

        elif solver == "neos":
            logger.info("Submitting model to NEOS server (opt={})", neos_opt)
            solver_manager = SolverManagerFactory("neos")
            self.results = solver_manager.solve(
                self.model, opt=neos_opt, tee=True
            )

        else:
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

        # Only the termination check is guarded: a result object without a
        # solver status is a solver-interface problem, not a modelling one.
        # The mapping call below must NOT be inside the guard — an
        # AttributeError raised while writing net.res_* would otherwise be
        # logged and swallowed, leaving the base-case power flow in place and
        # reporting success.
        try:
            optimal = check_optimal_termination(self.results)
        except AttributeError as err:
            logger.error("Could not check termination condition: {}", err)
            return self.results

        if optimal:
            logger.info("Optimal solution found")
            if to_net is not False:
                self.map_to_net(None if to_net is True else to_net)
        else:
            logger.warning(
                "Solver did not reach optimal termination (condition: {})",
                self.results.solver.termination_condition,
            )
        return self.results

    # --- nodal power balance ---------------------------------------------
    #
    # Each power-flow subclass supplies _kcl_real_rule (and, for the AC-style
    # ones, _kcl_reactive_rule); the construction and the flexibility hook are
    # shared so a device couples the same way whichever formulation is used.

    #: Balance constraints this formulation builds, in construction order.
    KCL_CONSTRAINTS = ("KCL_real", "KCL_reactive")

    def build_kcl(self):
        """Construct the nodal balance constraints, replacing any existing.

        Called once while the power flow is created and again from ``add_OPF``
        via :meth:`rebuild_kcl`, because flexibility devices attach to the
        model *after* the power-flow equations are first built and their
        injections have to reach the balance.
        """
        for name in self.KCL_CONSTRAINTS:
            rule = getattr(self, f"_{name.lower()}_rule", None)
            if rule is None:
                continue
            if hasattr(self.model, name):
                self.model.del_component(getattr(self.model, name))
            # Pyomo materialises an implicit index set for a multi-dimensional
            # constraint; it has to go too, or re-adding clashes on the name.
            index_name = f"{name}_index"
            if hasattr(self.model, index_name):
                self.model.del_component(getattr(self.model, index_name))
            setattr(
                self.model,
                name,
                Constraint(self.model.B, self.model.T, rule=rule),
            )

    def rebuild_kcl(self):
        """Rebuild the balance so late-attached devices are included.

        Idempotent: with no device registered, rebuilding reproduces the same
        constraints.
        """
        self.build_kcl()

    def KCL_flexibility(self, model, b, t, reactive=False):
        """Return the flexible-asset power at bus ``b`` and time ``t``.

        Sums whatever devices registered through
        ``Flexibility_multi_period.register_kcl_real`` and its reactive
        counterpart. The term sits on the consumption side of the balance, so
        it is positive for consumption and negative for injection — the load
        sign convention, matching ``pD`` / ``qD``.

        Returns 0 when nothing registered, which is every model built from
        ``net`` alone.
        """
        terms = Flexibility_multi_period.kcl_terms(model, reactive=reactive)
        if not terms:
            return 0
        return sum(term(model, b, t) for term in terms)

    def warm_start_from_pf(self, curtailment=1.0):
        """Seed every state variable from a power flow at each time step.

        A cold start puts ``v`` at 1.0 and leaves every angle and branch flow
        at zero, so Kirchhoff's laws are violated at every bus by the full
        nodal injection. IPOPT does not always recover from that on a nonconvex
        AC OPF: on a 12-step midday window of ``1-LV-rural1--0-sw`` it reported
        a locally infeasible point even though curtailing the PV to zero is
        both available and feasible. Seeding a *consistent* operating point
        fixes it. The seed need not be near the optimum — an uncurtailed,
        curtailed, or half-curtailed seed all converge to the same solution —
        it only has to satisfy the power flow. Since 0.5.3 that window also
        converges from the cold start (the static-generation lower bound is a
        constraint rather than a variable domain, which changes IPOPT's path);
        the seed remains the safer start and stays the default.

        Args:
            curtailment: Factor applied to the static-generation profile in the
                seeding power flow. Rarely needs changing; see
                :func:`~potpourri.models_multi_period.init_pyo_from_pp_res_multi_period.init_pyo_from_pp_res_multi_period`.

        Returns:
            Number of time steps successfully seeded. A step whose power flow
            does not converge is logged and left at its default values, so a
            partial seed is possible.
        """
        return init_pyo_from_pp_res_multi_period(
            self.net, self.model, self.bus_lookup, curtailment=curtailment
        )

    def map_to_net(self, t=None):
        """Write the solution for one time step into ``self.net.res_*``.

        The pandapower result tables carry no time dimension, so a horizon
        cannot be written whole — one step has to be chosen. Call this once
        per step of interest, reading ``net.res_*`` in between.

        Args:
            t: Time step from ``model.T``. Defaults to the last step of the
                horizon.

        Returns:
            The time step that was written, so callers can label results
            without re-deriving the default.

        Raises:
            ValueError: If ``t`` is not a step in ``model.T``.
        """
        last = self.model.T.last()
        if t is None:
            t = last
        else:
            t = int(t)
            if t not in self.model.T:
                raise ValueError(
                    f"t={t} is not a time step of this model; model.T covers "
                    f"[{self.model.T.first()}, {last}]."
                )
        pyo_sol_to_net_res(self.net, self.model, t)
        logger.debug("Solution for time step {} mapped to net.res_*", t)
        return t

    def change_vals(self, key, value):
        """Set all indices of a named Pyomo component to value."""
        component = self.model.component(key)
        if not component:
            logger.warning("Model has no component '{}'", key)
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
            logger.warning("Model has no component '{}'", key)
            return
        try:
            for index in component:
                if value is not None:
                    component[index].fix(value)
                else:
                    component[index].fix()
            logger.debug("Fixed variable '{}'", key)
        except AttributeError as err:
            logger.error("fix_vars failed for component '{}': {}", key, err)

    def unfix_vars(self, key, value=None):
        """Unfix all indices of a named Pyomo variable; optionally reset to
        value."""
        component = self.model.component(key)
        if not component:
            logger.warning("Model has no component '{}'", key)
            return
        try:
            for index in component:
                component[index].unfix()
                if value is not None:
                    component[index] = value
            logger.debug("Unfixed variable '{}'", key)
        except AttributeError as err:
            logger.error("unfix_vars failed for component '{}': {}", key, err)

    def make_to_dict(self, model_obj, model_time, data, time_dependent=True):
        """
        **make_to_dict** \n
        Args:
            model_obj: Object Indices
            model_time: Time Indices
            data: Data to be converted to dictionary
            time_dependent: should the data be time dependent True or
                constant for all time steps False

        Returns:
            data_dict: Dictionary with the data
            tuple_list: List of tuples with the object and time index

        """

        # Scalar zero sentinel — must check with isinstance to avoid ambiguous
        # truth-value error when `data` is a pandas Series.
        if isinstance(data, (int, float)) and data == 0:
            data_dict = {(o, t): 0 for o in model_obj for t in model_time}
            tuple_list = list([(o, t) for o in model_obj for t in model_time])
            return data_dict, tuple_list

        if isinstance(data, np.ndarray):
            if time_dependent:
                data_dict = {
                    (o, t): data[o][t] for o in model_obj for t in model_time
                }
                tuple_list = [(o, t) for o in model_obj for t in model_time]
            else:
                data_dict = {
                    (o, t): data[o] for o in model_obj for t in model_time
                }
                tuple_list = [(o, t) for o in model_obj for t in model_time]
            return data_dict, tuple_list

        # when data is float put it in dict, correct?
        if isinstance(data, float):
            data_dict = {(o, t): data for o in model_obj for t in model_time}
            tuple_list = list([(o, t) for o in model_obj for t in model_time])
            return data_dict, tuple_list

        # if data is already a dict just put it in data_dict and make
        # tuple_list with time
        if isinstance(data, dict):
            if time_dependent:
                data_dict = {
                    (o, t): data[o][t] for o in model_obj for t in model_time
                }
                tuple_list = list(
                    [(o, t) for o in model_obj for t in model_time]
                )
                return data_dict, tuple_list
            else:
                data_dict = {
                    (o, t): data[o] for o in model_obj for t in model_time
                }
                tuple_list = list(
                    [(o, t) for o in model_obj for t in model_time]
                )
                return data_dict, tuple_list

        # make data_dict with constant values over time if not time dependent
        if time_dependent:
            data_dict = data.to_dict()
            data_dict = {
                (o, t): data_dict[o][t] for o in model_obj for t in model_time
            }
            tuple_list = list([(o, t) for o in model_obj for t in model_time])
        else:
            data_dict = data.to_dict()
            data_dict = {
                (o, t): data_dict[o] for o in model_obj for t in model_time
            }
            tuple_list = list([(o, t) for o in model_obj for t in model_time])

        return data_dict, tuple_list

    def make_to_tuple(self, model_obj, model_time):
        """
        **make_to_tuple** \n
        Args:
            model_obj: Object Indices
            model_time: Time Indices

        Returns:
            tuple_list: List of tuples with the object and time index

        """
        tuple_list = list([(o, t) for o in model_obj for t in model_time])
        return tuple_list
