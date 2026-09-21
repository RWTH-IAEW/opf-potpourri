# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

r"""Nonlinear AC power-flow equations in polar form, as a model layer.

`AC` is the electrical heart of the single-period stack: it takes
the sets, parameters and active-power variables that
[`Basemodel`][potpourri.models.basemodel.Basemodel] extracted from a
pandapower network and adds everything the AC formulation needs on top --
voltage magnitudes, reactive power, the nodal power balance, and the
branch power-flow equations for lines and transformers.

## Where it sits

```text
Basemodel          sets, parameters, active power, base power flow
  └── AC           + voltages, reactive power, AC power flow   <- here
        └── ACOPF (in ACOPF_base.py)  + limits and an objective
              └── HC_ACOPF, LPACOPF, ...
```

`AC` is a **layer, not an application**. It states the physics and
nothing else: it adds no operating limits and no objective. Most users
should reach for [`ACOPF`][potpourri.models.ACOPF_base.ACOPF],
which subclasses `AC` and adds both. Construct `AC` directly only to
inspect or debug the power-flow equations themselves.

## Not a drop-in power-flow solver

The equation system `AC` builds is square -- as many free variables as
constraints -- only when the network has no static generators. Each
in-service `net.sgen` adds a free `qsG` variable that nothing in this
module binds: the reactive set point sits in the `QsG` parameter, but it
is `ACOPF.add_OPF` that either fixes `qsG` to it or replaces it
with a capability constraint. On a two-bus test network the difference
is exactly one degree of freedom:

| network                | vars | fixed | free | constraints | DOF |
|------------------------|------|-------|------|-------------|-----|
| 2 buses, no static gen |   12 |     4 |    8 |           8 |   0 |
| 2 buses, 1 static gen  |   14 |     5 |    9 |           8 |   1 |

For an ordinary power flow use `pandapower.runpp`; the constructor runs
one anyway (see below).

## Formulation

Polar coordinates, per unit on `net.sn_mva`, angles in radians. Buses
carry a magnitude `v[b]` and an angle `delta[b]`; every branch
contributes a power injection at each of its two ends. The model is a
system of nonlinear equalities in `sin`/`cos` -- **nonconvex**, so a
local solver such as IPOPT returns a local optimum with no global
guarantee.

Symbols map to Pyomo components as follows.

| symbol             | component          | meaning                   |
|--------------------|--------------------|---------------------------|
| $v_b$              | `model.v[b]`       | voltage magnitude (p.u.)  |
| $\theta_b$         | `model.delta[b]`   | voltage angle (rad)       |
| $G_{ii}$, $B_{ii}$ | `Gii` / `Bii`      | branch self terms         |
| $G_{ik}$, $B_{ik}$ | `Gik` / `Bik`      | branch mutual terms       |
| $\tau_t$           | `model.Tap[t]`     | transformer tap ratio     |
| $\varphi_t$        | `model.shift[t]`   | phase shift (rad)         |
| $q^{sgen}_g$       | `model.qsG[g]`     | static-generator reactive |

The bus balance (`KCL_real`, `KCL_reactive`) equates generation minus
storage draw to demand, the branch injections at that bus, and the shunt
term. The branch equations (`KVL_real_from` and siblings) are the
standard $\pi$-model relations

$$
p_{ik} = G_{ii} v_i^2
    + v_i v_k (G_{ik}\cos\theta_{ik} + B_{ik}\sin\theta_{ik})
$$

$$
q_{ik} = -B_{ii} v_i^2
    + v_i v_k (G_{ik}\sin\theta_{ik} - B_{ik}\cos\theta_{ik})
$$

with $\theta_{ik} = \theta_i - \theta_k$. Despite the component names,
these are **not** Kirchhoff's voltage law: they are Ohm's law applied to
the branch $\pi$ model, giving the power injected into one end as a
function of both terminal voltages. The `KVL_*` names are historical and
kept because they are public API.

Full derivations, including the transformer equations and the symbol
reference, live in the
[Mathematical Modelling](../mathematical-modelling.md) guide.

Example:
    Build the equations for a two-bus network and inspect them. No
    optimization solver is needed, but note that the constructor does
    run a pandapower power flow internally.

    >>> import pandapower as pp
    >>> from potpourri.models.AC import AC
    >>> net = pp.create_empty_network()
    >>> b0 = pp.create_bus(net, vn_kv=20.0)
    >>> b1 = pp.create_bus(net, vn_kv=20.0)
    >>> _ = pp.create_ext_grid(net, bus=b0, vm_pu=1.02)
    >>> _ = pp.create_line_from_parameters(
    ...     net, b0, b1, length_km=1.0, r_ohm_per_km=0.3,
    ...     x_ohm_per_km=0.1, c_nf_per_km=200.0, max_i_ka=0.4,
    ... )
    >>> _ = pp.create_load(net, bus=b1, p_mw=1.0, q_mvar=0.3)
    >>> ac = AC(net)
    >>> ac.model.name
    'AC'
    >>> ac.model.v[0].fixed          # slack magnitude is fixed
    True
    >>> len(ac.model.KCL_real)       # one balance per bus
    2

    To solve something, use the OPF layer instead:

    ```python
    from potpourri.models.ACOPF_base import ACOPF

    opf = ACOPF(net)
    opf.add_OPF()
    opf.add_active_power_costs_objective()
    results = opf.solve(solver="ipopt")
    # results land on the model's own copy of the network:
    opf.net.res_bus.vm_pu
    ```

See Also:
    - [`Basemodel`][potpourri.models.basemodel.Basemodel]: network
      extraction, sets, active power, the base power flow.
    - [`ACOPF`][potpourri.models.ACOPF_base.ACOPF]: `AC` plus
      operating limits and objectives; the usual entry point.
    - [`DC`][potpourri.models.DC.DC]: linearised alternative that
      replaces this module's equations.
    - `potpourri.models_multi_period.AC_multi_period`: the same
      formulation with a time index.
"""

import numpy as np
import pyomo.environ as pyo
from potpourri.models.basemodel import (
    Basemodel,
    branch_charging_admittance,
)


class AC(Basemodel):
    r"""Nonlinear AC power-flow equations on top of `Basemodel`.

    Adds the reactive-power and voltage-magnitude half of the model: the
    nodal power balance at every bus and the $\pi$-model power-flow
    equations at both ends of every line and transformer. See the module
    docstring for the formulation, the symbol table and a worked example.

    This class states physics only -- no operating limits, no objective.
    Unless you are inspecting the power-flow equations themselves, use
    [`ACOPF`][potpourri.models.ACOPF_base.ACOPF], which subclasses
    `AC` and adds both.

    Constructing an instance immediately builds the Pyomo model: the
    constructor calls `create_model` as its last step, so there is
    no separate build call and `self.model` is ready on return.

    Attributes:
        net: The model's **own deep copy** of the input network,
            preprocessed by `preprocess_grid` (bus-bus switches fused)
            and carrying the base power-flow result. The caller's network
            is never modified, and `solve` writes results
            here -- read `model.net.res_bus`, not your original object.
        model: The `pyomo.environ.ConcreteModel`. Available after
            construction.
        line_data: Per-line frame, extended here with the four derived
            admittance columns `Bii_data`, `Bik_data`, `Gii_data`
            and `Gik_data` (p.u.).
        trafo_data: Per-transformer frame, extended here with
            `BiiT_data`, `BikT_data`, `GiiT_data` and `GikT_data`.
        BB_data: Shunt susceptance per in-service shunt (p.u.), negated
            from `net.shunt.q_mvar` so that a reactive-power-consuming
            shunt gets a negative susceptance.
        QD_data: Reactive demand per load (p.u.), scaled by
            `net.load.scaling`.

    Args:
        net: A pandapower network that `pp.runpp` can handle. It is deep
            copied, so the caller's object is left untouched.

    Raises:
        ValueError: If `net` is not a `pandapowerNet` (from
            `Basemodel`).
    """

    def __init__(self, net):
        # Basemodel deep copies the network, fuses bus-bus switches and
        # runs the base power flow, which is what fills net._ppc below.
        super().__init__(net)

        self.BB_data = (
            -self.net.shunt.q_mvar * self.net.shunt.step / self.baseMVA
        )

        # --- line and transformer admittances (symmetric π model) ---
        # The branch series admittance is y_s = 1/(r + jx) = g + jb. We split
        # the branch into a from/to shunt of half the total line-charging
        # susceptance b_c and a mutual series term y_s.
        #
        # The charging admittance is y_c = g_c + j·b_c. MATPOWER carries
        # only the susceptance (BR_B, column 4); pandapower additionally
        # stores the shunt conductance in BR_G (column 23), which is where
        # the iron-loss (pfe_kw) part of a transformer's magnetising branch
        # lives. Reading only column 4 dropped those losses (0.46 kW on a
        # 160 kVA SimBench LV transformer, 14 kW on a 25 MVA MV one) and
        # made the model disagree with pandapower's own power flow.
        r = self.net._ppc["branch"][:, 2].real
        x = self.net._ppc["branch"][:, 3].real
        y = branch_charging_admittance(self.net._ppc["branch"])
        gt_ik = r / (r**2 + x**2)  # series conductance g
        bt_ik = -x / (r**2 + x**2)  # series susceptance b (b<0 inductive)
        BiiT = bt_ik + y.imag / 2  # self susceptance: b + b_c/2
        BikT = -bt_ik  # mutual susceptance: -b
        GiiT = gt_ik + y.real / 2  # self conductance: g + g_c/2 (g_c==0)
        GikT = -gt_ik  # mutual conductance: -g
        trafo_start = len(self.net.line)
        trafo_end = trafo_start + len(self.net.trafo)
        imp_table = self.net.get("impedance")
        n_imp = (
            len(imp_table)
            if imp_table is not None and not imp_table.empty
            else 0
        )
        # _ppc['branch'] layout (pandapower convention):
        #   [0 : trafo_start)              → lines
        #   [trafo_start : trafo_end)      → trafos
        #   [trafo_end : trafo_end + n_imp)→ impedance branches
        # We treat impedance rows as additional lines in the model — see
        # the matching block in Basemodel.__init__.
        if n_imp:
            line_idx_ppc = np.r_[
                np.arange(0, trafo_start),
                np.arange(trafo_end, trafo_end + n_imp),
            ]
        else:
            line_idx_ppc = np.arange(0, trafo_start)

        self.line_data = self.line_data.assign(
            **{
                "Bii_data": BiiT[line_idx_ppc],
                "Bik_data": BikT[line_idx_ppc],
                "Gii_data": GiiT[line_idx_ppc],
                "Gik_data": GikT[line_idx_ppc],
            }
        )
        self.trafo_data = self.trafo_data.assign(
            **{
                "BiiT_data": BiiT[trafo_start:trafo_end],
                "BikT_data": BikT[trafo_start:trafo_end],
                "GiiT_data": GiiT[trafo_start:trafo_end],
                "GikT_data": GikT[trafo_start:trafo_end],
            }
        )

        # generator and external grids voltage set points
        self.generation_data["v"] = self.net._ppc["gen"][:, 5]

        qsg = (
            self.net.sgen.q_mvar.fillna(0)
            * self.net.sgen.scaling
            / self.baseMVA
        ).values
        self.static_generation_data["q"] = qsg

        self.QD_data = (
            self.net.load.q_mvar * self.net.load.scaling / self.baseMVA
        )

        self.create_model()

    def create_model(self):
        """Add the AC components to `self.model`, in place.

        Called automatically by `__init__`; call it again only to
        rebuild from scratch. It first delegates to
        `create_model`, which **replaces** `self.model`
        with a fresh `ConcreteModel` holding the sets, the active-power
        variables and the base-case fixings. This method then adds the
        reactive and voltage layer on top.

        Returns `None`. The result is the model object itself, reachable
        as `self.model`; nothing is returned to the caller.

        Components added (all indexed over sets from `Basemodel`):

        ===================  ============================================
        Params               `BB` (shunt susceptance), `Bii`/`Bik`/
                             `Gii`/`Gik` (lines), `BiiT`/`BikT`/
                             `GiiT`/`GikT` (transformers), `QsG`,
                             `QD`, `v_bPV`, `v_b0`
        Vars                 `v` (magnitude), `qsG`, `qD`, `qG`,
                             `qLfrom`/`qLto`, `qThv`/`qTlv`
        Constraints          `KCL_real`, `KCL_reactive` (per bus),
                             `KVL_real_from`/`_to`,
                             `KVL_reactive_from`/`_to` (per line),
                             the four `*Transf` analogues,
                             `v_bPV_setpoint`
        ===================  ============================================

        Which quantities end up fixed, bounded or free matters for anyone
        solving the model:

        * **Fixed** here: `qD` at the load set point, and `v` at every
          reference bus. A fixed Pyomo variable is removed from the
          problem, not merely bounded.
        * **Constrained**: `v` at PV buses, via `v_bPV_setpoint`
          rather than a fix, so that the OPF layer can deactivate the
          constraint and let the voltage float.
        * **Bounded**: `v` on `(0.0, 2.0)` p.u. -- a numerical
          safeguard to keep the solver away from the singularity at
          `v = 0`, not a physical voltage band. Real limits come from
          `add_OPF`.
        * **Free**: `qG` (reactive support at PV/reference buses) and
          `qsG`. `qsG` stays free even though its set point sits in
          the `QsG` parameter -- see the module docstring; this is why a
          bare `AC` model is underdetermined once the network has static
          generators.
        * **Initialised**: `v` starts from the base power flow, falling
          back to 1.0 p.u. where that is missing.
        """
        super().create_model()

        self.model.name = "AC"

        # shunt
        self.model.BB = pyo.Param(
            self.model.SHUNT,
            within=pyo.Reals,
            initialize=self.BB_data[self.model.SHUNT],
        )  # shunt susceptance

        # derived line pyo.Parameters
        self.model.Bii = pyo.Param(
            self.model.L,
            within=pyo.Reals,
            initialize=self.line_data.Bii_data[self.model.L],
        )
        self.model.Bik = pyo.Param(
            self.model.L,
            within=pyo.Reals,
            initialize=self.line_data.Bik_data[self.model.L],
        )
        self.model.Gii = pyo.Param(
            self.model.L,
            within=pyo.Reals,
            initialize=self.line_data.Gii_data[self.model.L],
        )
        self.model.Gik = pyo.Param(
            self.model.L,
            within=pyo.Reals,
            initialize=self.line_data.Gik_data[self.model.L],
        )

        ## derived transformer pyo.Parameters
        self.model.BiiT = pyo.Param(
            self.model.TRANSF,
            within=pyo.Reals,
            initialize=self.trafo_data.BiiT_data[self.model.TRANSF],
        )
        self.model.BikT = pyo.Param(
            self.model.TRANSF,
            within=pyo.Reals,
            initialize=self.trafo_data.BikT_data[self.model.TRANSF],
        )
        self.model.GiiT = pyo.Param(
            self.model.TRANSF,
            within=pyo.Reals,
            initialize=self.trafo_data.GiiT_data[self.model.TRANSF],
        )
        self.model.GikT = pyo.Param(
            self.model.TRANSF,
            within=pyo.Reals,
            initialize=self.trafo_data.GikT_data[self.model.TRANSF],
        )

        # reactive generation
        self.model.QsG = pyo.Param(
            self.model.sG,
            initialize=self.static_generation_data["q"][self.model.sG],
        )

        self.model.v_bPV = pyo.Param(
            self.model.bPV,
            within=pyo.NonNegativeReals,
            initialize=self.bus_data.v_m[self.model.bPV],
        )

        # reactive demand
        self.model.QD = pyo.Param(
            self.model.D, initialize=self.QD_data[self.model.D]
        )

        # external grid voltage
        self.model.v_b0 = pyo.Param(
            self.model.b0,
            within=pyo.NonNegativeReals,
            initialize=self.bus_data.v_m[self.model.b0],
        )

        # --- control pyo.Variables ---
        self.model.qsG = pyo.Var(
            self.model.sG, domain=pyo.Reals
        )  # reactive power of static generators

        self.model.qD = pyo.Var(
            self.model.D, domain=pyo.Reals
        )  # reactive power absorbed by demand

        self.model.qLfrom = pyo.Var(
            self.model.L, domain=pyo.Reals
        )  # reactive power injected at b onto line
        self.model.qLto = pyo.Var(
            self.model.L, domain=pyo.Reals
        )  # reactive power injected at b' onto line
        self.model.qThv = pyo.Var(
            self.model.TRANSF, domain=pyo.Reals
        )  # reactive power injected at b onto transformer
        self.model.qTlv = pyo.Var(
            self.model.TRANSF, domain=pyo.Reals
        )  # reactive power injected at b' onto transformer

        v_init = self.bus_data.v_m.fillna(1.0).to_dict()
        self.model.v = pyo.Var(
            self.model.B,
            domain=pyo.NonNegativeReals,
            bounds=(0.0, 2.0),
            initialize=v_init,
        )  # voltage magnitude at bus b (p.u.)

        self.model.qG = pyo.Var(self.model.G, domain=pyo.Reals)

        # --- nodal power balance at each bus b ---
        def KCL_real_def(model, b):
            """Active-power balance at bus `b`.

            `generation - storage = demand + branch outflows + shunt`,
            where the storage and shunt terms follow pandapower's load
            convention (positive = consumption), hence the signs.

            The shunt term `GB * v**2` is the active loss in a shunt's
            conductance and belongs on the consumption side.

            Args:
                model: The Pyomo model being built.
                b: Bus index (an element of `model.B`, i.e. a **ppc**
                    bus number, which may exceed `len(net.bus)` when
                    pandapower inserted auxiliary buses).

            Returns:
                A Pyomo equality expression, or `Constraint.Skip` for a
                bus where every term is constant. That happens on an
                isolated bus with no branches and only fixed injections:
                both sides then evaluate to plain numbers and Python
                collapses `==` to a `bool`, which Pyomo cannot accept
                as a constraint.
            """
            kcl = sum(
                model.psG[g] for g in model.sG if (g, b) in model.sGbs
            ) + sum(model.pG[g] for g in model.G if (g, b) in model.Gbs) - sum(
                model.pSTOR[s] for s in model.STOR if model.STOR_bus[s] == b
            ) == sum(
                model.pD[d] for d in model.D if (b, d) in model.Dbs
            ) + sum(
                model.pLfrom[l] for l in model.L if model.A[l, 1] == b
            ) + sum(
                model.pLto[l] for l in model.L if model.A[l, 2] == b
            ) + sum(
                model.pThv[l] for l in model.TRANSF if model.AT[l, 1] == b
            ) + sum(
                model.pTlv[l] for l in model.TRANSF if model.AT[l, 2] == b
            ) + sum(
                model.GB[s] * model.v[b] ** 2
                for s in model.SHUNT
                if (b, s) in model.SHUNTbs and model.GB[s] != 0
            )
            if isinstance(kcl, (bool, np.bool_)):
                return pyo.Constraint.Skip
            return kcl

        def KCL_reactive_def(model, b):
            """Reactive-power balance at bus `b`.

            Mirrors `KCL_real_def` with `q` in place of `p`. The
            shunt term is `- BB * v**2`: `BB` was negated when it was
            built from `net.shunt.q_mvar`, so a reactive-consuming
            (inductive) shunt has `BB < 0` and the term lands on the
            consumption side with a positive value.

            Args:
                model: The Pyomo model being built.
                b: Bus index from `model.B` (a ppc bus number).

            Returns:
                A Pyomo equality expression, or `Constraint.Skip` when
                the balance degenerates to a constant -- see
                `KCL_real_def`.
            """
            kcl = sum(
                model.qsG[g] for g in model.sG if (g, b) in model.sGbs
            ) + sum(model.qG[g] for g in model.G if (g, b) in model.Gbs) - sum(
                # storage reactive power follows pandapower's load convention
                # (positive = consumption), like pSTOR in KCL_real
                model.qSTOR[s]
                for s in model.STOR
                if model.STOR_bus[s] == b
            ) == sum(
                model.qD[d] for d in model.D if (b, d) in model.Dbs
            ) + sum(
                model.qLfrom[l] for l in model.L if model.A[l, 1] == b
            ) + sum(
                model.qLto[l] for l in model.L if model.A[l, 2] == b
            ) + sum(
                model.qThv[l] for l in model.TRANSF if model.AT[l, 1] == b
            ) + sum(
                model.qTlv[l] for l in model.TRANSF if model.AT[l, 2] == b
            ) - sum(
                model.BB[s] * model.v[b] ** 2
                for s in model.SHUNT
                if (b, s) in model.SHUNTbs and model.BB[s] != 0
            )
            if isinstance(kcl, (bool, np.bool_)):
                return pyo.Constraint.Skip
            return kcl

        self.model.KCL_real = pyo.Constraint(self.model.B, rule=KCL_real_def)
        self.model.KCL_reactive = pyo.Constraint(
            self.model.B, rule=KCL_reactive_def
        )

        # --- branch power flow on each line (both ends) ---
        # Despite the historical `KVL_*` names -- kept because they are
        # public component names -- these are not Kirchhoff's voltage
        # law. Each one is Ohm's law on the branch pi model, giving the
        # power pushed into one end as a function of both terminal
        # voltages. Lines are symmetric, so the same Gii/Bii serve both
        # ends and only the sign of the angle difference flips.
        def KVL_real_fromend(model, l):
            r"""Active power entering line `l` at its *from* bus.

            $p_{ik} = G_{ii} v_i^2 + v_i v_k
            (G_{ik}\cos\theta_{ik} + B_{ik}\sin\theta_{ik})$, with
            $i$ the from bus, $k$ the to bus and
            $\theta_{ik} = \theta_i - \theta_k$.

            Args:
                model: The Pyomo model being built.
                l: Line index from `model.L`. `model.A[l, 1]` and
                    `model.A[l, 2]` give its from and to bus.

            Returns:
                A Pyomo equality expression defining `pLfrom[l]`.
            """
            return model.pLfrom[l] == model.Gii[l] * (
                model.v[model.A[l, 1]] ** 2
            ) + model.v[model.A[l, 1]] * model.v[model.A[l, 2]] * (
                model.Bik[l]
                * pyo.sin(
                    model.delta[model.A[l, 1]] - model.delta[model.A[l, 2]]
                )
                + model.Gik[l]
                * pyo.cos(
                    model.delta[model.A[l, 1]] - model.delta[model.A[l, 2]]
                )
            )

        def KVL_real_toend(model, l):
            r"""Active power entering line `l` at its *to* bus.

            The same relation as `KVL_real_fromend` with the two
            terminals swapped, so the angle difference is
            $\theta_k - \theta_i$. A line is symmetric, so the
            self term reuses the same `Gii`.

            Args:
                model: The Pyomo model being built.
                l: Line index from `model.L`.

            Returns:
                A Pyomo equality expression defining `pLto[l]`.
            """
            return model.pLto[l] == model.Gii[l] * (
                model.v[model.A[l, 2]] ** 2
            ) + model.v[model.A[l, 1]] * model.v[model.A[l, 2]] * (
                model.Bik[l]
                * pyo.sin(
                    model.delta[model.A[l, 2]] - model.delta[model.A[l, 1]]
                )
                + model.Gik[l]
                * pyo.cos(
                    model.delta[model.A[l, 2]] - model.delta[model.A[l, 1]]
                )
            )

        def KVL_reactive_fromend(model, l):
            r"""Reactive power entering line `l` at its *from* bus.

            $q_{ik} = -B_{ii} v_i^2 + v_i v_k
            (G_{ik}\sin\theta_{ik} - B_{ik}\cos\theta_{ik})$. `Bii`
            already carries half the line-charging susceptance, so the
            capacitive generation of the $\pi$ model sits in the
            $-B_{ii} v_i^2$ term.

            Args:
                model: The Pyomo model being built.
                l: Line index from `model.L`.

            Returns:
                A Pyomo equality expression defining `qLfrom[l]`.
            """
            return model.qLfrom[l] == -model.Bii[l] * (
                model.v[model.A[l, 1]] ** 2
            ) + model.v[model.A[l, 1]] * model.v[model.A[l, 2]] * (
                model.Gik[l]
                * pyo.sin(
                    model.delta[model.A[l, 1]] - model.delta[model.A[l, 2]]
                )
                - model.Bik[l]
                * pyo.cos(
                    model.delta[model.A[l, 1]] - model.delta[model.A[l, 2]]
                )
            )

        def KVL_reactive_toend(model, l):
            """Reactive power entering line `l` at its *to* bus.

            `KVL_reactive_fromend` with the terminals swapped.

            Args:
                model: The Pyomo model being built.
                l: Line index from `model.L`.

            Returns:
                A Pyomo equality expression defining `qLto[l]`.
            """
            return model.qLto[l] == -model.Bii[l] * (
                model.v[model.A[l, 2]] ** 2
            ) + model.v[model.A[l, 1]] * model.v[model.A[l, 2]] * (
                model.Gik[l]
                * pyo.sin(
                    model.delta[model.A[l, 2]] - model.delta[model.A[l, 1]]
                )
                - model.Bik[l]
                * pyo.cos(
                    model.delta[model.A[l, 2]] - model.delta[model.A[l, 1]]
                )
            )

        self.model.KVL_real_from = pyo.Constraint(
            self.model.L, rule=KVL_real_fromend
        )
        self.model.KVL_real_to = pyo.Constraint(
            self.model.L, rule=KVL_real_toend
        )
        self.model.KVL_reactive_from = pyo.Constraint(
            self.model.L, rule=KVL_reactive_fromend
        )
        self.model.KVL_reactive_to = pyo.Constraint(
            self.model.L, rule=KVL_reactive_toend
        )

        # --- branch power flow on each transformer (both ends) ---
        # The transformer is the line pi model behind an ideal
        # transformer of complex ratio tau = Tap * exp(j * shift), placed
        # on the *from* (HV) side, as MATPOWER does it. That placement is
        # what makes the two ends asymmetric: the from-end self term
        # picks up 1/Tap**2 and the from-end angle gets - shift, while
        # the to-end self term keeps its bare GiiT/BiiT. Both ends share
        # the mutual factor v_hv * v_lv / Tap.
        #
        # Each rule branches on `if model.shift[l]:`, which is a plain
        # truth test on the Param's value. The two branches are the same
        # equation; the shift == 0 branch just omits the `- 0` terms so
        # the expression tree stays small on the many transformers that
        # are not phase shifters.
        def KVL_real_fromendTransf(model, l):
            r"""Active power entering transformer `l` at its HV bus.

            $p_{hv} = \frac{G_{ii}}{\tau^2} v_{hv}^2 +
            \frac{v_{hv} v_{lv}}{\tau}
            (G_{ik}\cos(\theta_{hv} - \theta_{lv} - \varphi)
            + B_{ik}\sin(\theta_{hv} - \theta_{lv} - \varphi))$.

            Args:
                model: The Pyomo model being built.
                l: Transformer index from `model.TRANSF`.
                    `model.AT[l, 1]` is the HV (from) bus and
                    `model.AT[l, 2]` the LV (to) bus.

            Returns:
                A Pyomo equality expression defining `pThv[l]`.
            """
            if model.shift[l]:
                return model.pThv[l] == model.GiiT[l] / model.Tap[l] ** 2 * (
                    model.v[model.AT[l, 1]] ** 2
                ) + model.v[model.AT[l, 1]] * model.v[
                    model.AT[l, 2]
                ] / model.Tap[l] * (
                    model.GikT[l]
                    * pyo.cos(
                        model.delta[model.AT[l, 1]]
                        - model.delta[model.AT[l, 2]]
                        - model.shift[l]
                    )
                    + model.BikT[l]
                    * pyo.sin(
                        model.delta[model.AT[l, 1]]
                        - model.delta[model.AT[l, 2]]
                        - model.shift[l]
                    )
                )

            return model.pThv[l] == model.GiiT[l] / model.Tap[l] ** 2 * (
                model.v[model.AT[l, 1]] ** 2
            ) + model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] / model.Tap[
                l
            ] * (
                model.GikT[l]
                * pyo.cos(
                    model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]]
                )
                + model.BikT[l]
                * pyo.sin(
                    model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]]
                )
            )

        def KVL_real_toendTransf(model, l):
            r"""Active power entering transformer `l` at its LV bus.

            $p_{lv} = G_{ii} v_{lv}^2 +
            \frac{v_{hv} v_{lv}}{\tau}
            (G_{ik}\cos(\theta_{lv} - \theta_{hv} + \varphi)
            + B_{ik}\sin(\theta_{lv} - \theta_{hv} + \varphi))$.

            The self term has **no** $1/\tau^2$: the ideal
            transformer sits on the HV side, so only that end is
            referred through the ratio.

            Args:
                model: The Pyomo model being built.
                l: Transformer index from `model.TRANSF`.

            Returns:
                A Pyomo equality expression defining `pTlv[l]`.
            """
            if model.shift[l]:
                return model.pTlv[l] == model.GiiT[l] * (
                    model.v[model.AT[l, 2]] ** 2
                ) + model.v[model.AT[l, 1]] * model.v[
                    model.AT[l, 2]
                ] / model.Tap[l] * (
                    model.BikT[l]
                    * pyo.sin(
                        model.delta[model.AT[l, 2]]
                        - model.delta[model.AT[l, 1]]
                        + model.shift[l]
                    )
                    + model.GikT[l]
                    * pyo.cos(
                        model.delta[model.AT[l, 2]]
                        - model.delta[model.AT[l, 1]]
                        + model.shift[l]
                    )
                )

            return model.pTlv[l] == model.GiiT[l] * (
                model.v[model.AT[l, 2]] ** 2
            ) + model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] / model.Tap[
                l
            ] * (
                model.BikT[l]
                * pyo.sin(
                    model.delta[model.AT[l, 2]] - model.delta[model.AT[l, 1]]
                )
                + model.GikT[l]
                * pyo.cos(
                    model.delta[model.AT[l, 2]] - model.delta[model.AT[l, 1]]
                )
            )

        def KVL_reactive_fromendTransf(model, l):
            r"""Reactive power entering transformer `l` at its HV bus.

            The reactive counterpart of `KVL_real_fromendTransf`:
            $q_{hv} = -\frac{B_{ii}}{\tau^2} v_{hv}^2 +
            \frac{v_{hv} v_{lv}}{\tau}
            (G_{ik}\sin(\theta_{hv} - \theta_{lv} - \varphi)
            - B_{ik}\cos(\theta_{hv} - \theta_{lv} - \varphi))$.

            Args:
                model: The Pyomo model being built.
                l: Transformer index from `model.TRANSF`.

            Returns:
                A Pyomo equality expression defining `qThv[l]`.
            """
            if model.shift[l]:
                return model.qThv[l] == -model.BiiT[l] / model.Tap[l] ** 2 * (
                    model.v[model.AT[l, 1]] ** 2
                ) + model.v[model.AT[l, 1]] * model.v[
                    model.AT[l, 2]
                ] / model.Tap[l] * (
                    -model.BikT[l]
                    * pyo.cos(
                        model.delta[model.AT[l, 1]]
                        - model.delta[model.AT[l, 2]]
                        - model.shift[l]
                    )
                    + model.GikT[l]
                    * pyo.sin(
                        model.delta[model.AT[l, 1]]
                        - model.delta[model.AT[l, 2]]
                        - model.shift[l]
                    )
                )

            return model.qThv[l] == -model.BiiT[l] / model.Tap[l] ** 2 * (
                model.v[model.AT[l, 1]] ** 2
            ) + model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] / model.Tap[
                l
            ] * (
                -model.BikT[l]
                * pyo.cos(
                    model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]]
                )
                + model.GikT[l]
                * pyo.sin(
                    model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]]
                )
            )

        def KVL_reactive_toendTransf(model, l):
            """Reactive power entering transformer `l` at its LV bus.

            The reactive counterpart of `KVL_real_toendTransf`,
            again with a bare $-B_{ii} v_{lv}^2$ self term because
            the ideal transformer is on the HV side.

            Args:
                model: The Pyomo model being built.
                l: Transformer index from `model.TRANSF`.

            Returns:
                A Pyomo equality expression defining `qTlv[l]`.
            """
            if model.shift[l]:
                return model.qTlv[l] == -model.BiiT[l] * (
                    model.v[model.AT[l, 2]] ** 2
                ) + model.v[model.AT[l, 1]] * model.v[
                    model.AT[l, 2]
                ] / model.Tap[l] * (
                    -model.BikT[l]
                    * pyo.cos(
                        model.delta[model.AT[l, 2]]
                        - model.delta[model.AT[l, 1]]
                        + model.shift[l]
                    )
                    + model.GikT[l]
                    * pyo.sin(
                        model.delta[model.AT[l, 2]]
                        - model.delta[model.AT[l, 1]]
                        + model.shift[l]
                    )
                )

            return model.qTlv[l] == -model.BiiT[l] * (
                model.v[model.AT[l, 2]] ** 2
            ) + model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] / model.Tap[
                l
            ] * (
                -model.BikT[l]
                * pyo.cos(
                    model.delta[model.AT[l, 2]] - model.delta[model.AT[l, 1]]
                )
                + model.GikT[l]
                * pyo.sin(
                    model.delta[model.AT[l, 2]] - model.delta[model.AT[l, 1]]
                )
            )

        self.model.KVL_real_fromTransf = pyo.Constraint(
            self.model.TRANSF, rule=KVL_real_fromendTransf
        )
        self.model.KVL_real_toTransf = pyo.Constraint(
            self.model.TRANSF, rule=KVL_real_toendTransf
        )
        self.model.KVL_reactive_fromTransf = pyo.Constraint(
            self.model.TRANSF, rule=KVL_reactive_fromendTransf
        )
        self.model.KVL_reactive_toTransf = pyo.Constraint(
            self.model.TRANSF, rule=KVL_reactive_toendTransf
        )

        # --- reactive demand limits ---
        for d in self.model.D:
            self.model.qD[d].fix(self.model.QD[d])

        # --- generator voltage operating point ---
        # A constraint rather than a fix, so that add_OPF() can
        # deactivate it and let the voltage float within its limits.
        def v_bPV_setpoint_rule(model, b):
            """Hold a PV bus at its generator voltage set point.

            Args:
                model: The Pyomo model being built.
                b: Bus index from `model.bPV` (ppc bus type 2).

            Returns:
                A Pyomo equality expression pinning `v[b]` to
                `v_bPV[b]`.
            """
            return model.v[b] == model.v_bPV[b]

        self.model.v_bPV_setpoint = pyo.Constraint(
            self.model.bPV, rule=v_bPV_setpoint_rule
        )

        # --- reference bus voltage pyo.Constraint ---
        for b0 in self.model.b0:
            self.model.v[b0].fix(self.model.v_b0[b0])
