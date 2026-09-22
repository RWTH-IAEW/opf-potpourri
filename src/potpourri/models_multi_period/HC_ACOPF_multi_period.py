# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Multi-period hosting-capacity AC OPF for wind integration."""

import copy
import math

import pyomo.environ as pyo
from potpourri.models_multi_period.ACOPF_multi_period import (
    ACOPF_multi_period,
)
from potpourri.technologies.windpower import Windpower_multi_period
import pandapower as pp

# TODO make multiperiod


class HC_ACOPF_multi_period(ACOPF_multi_period):
    """Multi-period hosting-capacity AC OPF for wind.

    A candidate wind generator is placed at every non-slack bus and the
    model decides how much wind the network can absorb over the horizon,
    with the constraints delegated to
    [`Windpower_multi_period`][potpourri.technologies.windpower.Windpower_multi_period].

    **What is decided once, and what varies.** Building a plant is a single
    decision, so the selection `y[w]` and the squared installed rating
    `SW2[w]` carry no time index. Dispatch does: `psG[w, t]` and
    `qsG[w, t]` move over the horizon, bounded at every step by the rating
    that was chosen, and the grid-code capability limits bind per
    `(w, t)`. The objective sums wind infeed over `model.T` and charges for
    the losses it causes.

    **What the answer means.** The objective rewards *energy*, not
    installed capacity, and nothing prices `SW2`, so a solution's `SW2[w]`
    is an upper envelope rather than a tight rating. The hosting capacity
    a candidate actually realises is the largest apparent power it reaches
    across the horizon, which [`hosting_capacity_mva`][
    potpourri.models_multi_period.HC_ACOPF_multi_period.HC_ACOPF_multi_period.hosting_capacity_mva]
    returns. A formulation that maximises installed capacity directly would
    need a wind-availability profile per candidate, which these candidates
    do not have; see the discussion on issue #22.

    The binary `y` makes this a MINLP. Solve it with a MINLP solver, or
    relax `y` to `UnitInterval` for a bound.
    """

    #: Profile the hosting-capacity candidates are given so that the
    #: multi-period base class can resolve one for every sgen. It is all
    #: zeros and never binds: `Windpower_multi_period.unfix_variables`
    #: frees the candidates' dispatch before the model is solved.
    CANDIDATE_PROFILE = "hc_candidate"

    def __init__(self, net, toT, fromT=None, pf=1):
        if "wind_hc" not in net.sgen:
            net = copy.deepcopy(net)
            buses_excl_extGrids = net.bus.loc[
                ~net.bus.index.isin(net.ext_grid.bus)
            ].index

            pp.create_sgens(net, buses_excl_extGrids, p_mw=0, wind_hc=True)
        net.sgen.wind_hc = net.sgen.wind_hc.fillna(False)
        self._give_candidates_a_profile(net)

        super().__init__(net, toT, fromT, pf)

    @classmethod
    def _give_candidates_a_profile(cls, net):
        """Point every profile-less sgen at a zero profile.

        A hosting-capacity candidate is a decision variable, not an infeed
        read from a time series, but `Basemodel_multi_period` resolves a
        SimBench profile for every sgen before any model exists and raises
        on the ones this class just created. Giving them a column of zeros
        satisfies that pass without inventing generation: the candidates
        are unfixed before the solve, so the zeros never bind.

        Args:
            net: The pandapower network, modified in place. Networks with
                no `profile` column, or none missing, are left alone.

        Returns:
            None.
        """
        if "profile" not in net.sgen or not net.sgen["profile"].isna().any():
            return
        profiles = getattr(net, "profiles", None)
        if not isinstance(profiles, dict) or "renewables" not in profiles:
            return
        renewables = profiles["renewables"]
        if cls.CANDIDATE_PROFILE not in renewables:
            renewables[cls.CANDIDATE_PROFILE] = 0.0
        missing = net.sgen["profile"].isna()
        net.sgen.loc[missing, "profile"] = cls.CANDIDATE_PROFILE

        # noinspection PyProtectedMember

    def _calc_opf_parameters(self, SWmax=10000, SWmin=0, **kwargs):
        """Add the wind hosting-capacity apparent-power bounds.

        Extend AC-OPF parameters with wind HC apparent power bounds via
        Windpower_multi_period.

        Args:
            SWmax: Upper bound on the wind hosting-capacity apparent power.
            SWmin: Lower bound on the wind hosting-capacity apparent power.
            **kwargs: Forwarded up the chain, which rejects names no model
                consumes.
        """
        super()._calc_opf_parameters(**kwargs)

        self._windpower()._calc_wind_opf_parameters(
            self.model, sw_max_mva=SWmax, sw_min_mva=SWmin
        )

    def _windpower(self):
        """The attached wind device, or a clear error saying it is missing.

        Returns:
            The `Windpower_multi_period` in `self.flexibilities`.

        Raises:
            RuntimeError: If none is attached. That used to happen whenever
                `net.bus` had no `windpot_p_mw` column, and the model then
                quietly came out as a plain ACOPF with no hosting-capacity
                variables at all.
        """
        device = next(
            (
                obj
                for obj in self.flexibilities
                if isinstance(obj, Windpower_multi_period)
            ),
            None,
        )
        if device is None:
            raise RuntimeError(
                "no Windpower_multi_period is attached, so this model has no "
                "hosting-capacity layer. Basemodel_multi_period attaches one "
                "when net.sgen carries a wind_hc column or net.bus carries "
                "windpot_p_mw; HC_ACOPF_multi_period.__init__ creates the "
                "former, so reaching this means the network was modified "
                "after construction."
            )
        return device

    def add_OPF(self, **kwargs):
        """Add the hosting-capacity wind constraints and objective.

        Extend ACOPF.add_OPF() with HC wind constraints and objective via
        Windpower_multi_period.
        """
        super().add_OPF(**kwargs)

        self.model.name = "HC_ACOPF"

        # Unconditional: the hosting-capacity layer *is* this model. It used
        # to be built only when net.bus carried the optional windpot_p_mw
        # column, so without it add_OPF() returned a plain ACOPF that looked
        # like a hosting-capacity model and had none of its variables.
        windpower_object = self._windpower()
        windpower_object.get_hc_acopf_parameters(self.model, self.net)
        windpower_object.get_hc_acopf_variables(self.model)
        windpower_object.get_objective(self.model)
        windpower_object.get_constraints(self.model, self.net)
        windpower_object.unfix_variables(self.model)

    def hosting_capacity_mva(self):
        """Apparent power each candidate reaches, in MVA, after a solve.

        The meaningful hosting-capacity figure for this model. `SW2[w]` is
        only an upper envelope — nothing in the objective prices it down —
        so the capacity a candidate realises is the largest apparent power
        it carries at any step of the horizon, scaled by the selection
        variable so unselected sites report zero.

        Returns:
            A dict mapping the candidate's sgen index to its apparent power
            in MVA. Empty if the model has no hosting-capacity layer.

        Raises:
            RuntimeError: If called before the model has been solved, where
                the variables still hold their initial values.
        """
        model = self.model
        if not hasattr(model, "WIND_HC"):
            return {}
        if getattr(self, "results", None) is None:
            raise RuntimeError(
                "solve() first: before a solve the variables still hold "
                "their initial values, and this would report those."
            )
        capacity = {}
        for w in model.WIND_HC:
            selected = pyo.value(model.y[w])
            peak = max(
                math.hypot(
                    pyo.value(model.psG[w, t]), pyo.value(model.qsG[w, t])
                )
                for t in model.T
            )
            capacity[w] = peak * selected * self.baseMVA
        return capacity

    def add_loss_obj(self):
        """Replace the objective with a weighted wind-versus-loss one.

        Replace default objective with weighted wind-vs-loss objective
        using mutable eps parameter.
        """
        self.model.eps = pyo.Param(
            domain=pyo.Reals, initialize=1.0, mutable=True
        )

        # `obj`, not `obj_hc`: the latter is a name nothing in the
        # package ever created, so this raised AttributeError every time.
        self.model.obj.deactivate()

        @self.model.Objective(sense=pyo.maximize)
        def OBJ_with_loss(model):
            """Weighted trade-off between wind energy and network losses.

            Both terms are summed over the horizon, matching the objective
            this replaces: a hosting-capacity answer for a multi-period
            study is about the energy the network can absorb across the
            window, not the power at one instant.

            Args:
                model: The Pyomo model being built.

            Returns:
                A Pyomo expression, weighted by the mutable `eps`.
            """
            infeed = sum(
                model.psG[w, t] for w in model.WIND_HC for t in model.T
            )
            losses = sum(
                model.pLfrom[line, t] + model.pLto[line, t]
                for line in model.L
                for t in model.T
            )
            return model.eps * infeed + (1 - model.eps) * (-losses)
