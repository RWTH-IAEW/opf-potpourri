---
title: 'potpourri: Multi-period Optimal Power Flow for Distribution Grids'
tags:
  - Python
  - power systems
  - optimal power flow
  - distribution grid
  - multi-period optimisation
  - battery storage
  - hosting capacity
authors:
  - name: Steffen Kortmann
    orcid: 0000-XXXX-XXXX-XXXX   # NEEDS REVIEW: add ORCID iD
    corresponding: true
    affiliation: 1
  - name: Andreas Bong
    orcid: 0000-XXXX-XXXX-XXXX   # NEEDS REVIEW: add ORCID iD
    affiliation: 1
  - name: Simon Braun
    orcid: 0000-XXXX-XXXX-XXXX   # NEEDS REVIEW: add ORCID iD
    affiliation: 1
  - name: Alexander Och
    orcid: 0000-XXXX-XXXX-XXXX   # NEEDS REVIEW: add ORCID iD
    affiliation: 1
  - name: Farah Nasr
    orcid: 0000-XXXX-XXXX-XXXX   # NEEDS REVIEW: add ORCID iD
    affiliation: 1
  - name: Philip Kvesic
    orcid: 0000-XXXX-XXXX-XXXX   # NEEDS REVIEW: add ORCID iD
    affiliation: 1
  - name: Nina Stumberger
    orcid: 0000-XXXX-XXXX-XXXX   # NEEDS REVIEW: add ORCID iD
    affiliation: 1
affiliations:
  - name: >-
      Institute for High Voltage Equipment and Grids, Digitalization and
      Energy Economics (IAEW), RWTH Aachen University, Aachen, Germany
    index: 1
date: 2026-04-24
bibliography: references.bib
---

# Summary

`potpourri` is an open-source Python library for **multi-period Optimal Power
Flow (OPF)** in distribution grids. It provides a modular, Pyomo-based
[@bynum2021pyomo; @hart2011pyomo] optimisation layer that operates directly on
pandapower [@thurner2018pandapower] network objects. The library supports both
exact AC and linearised DC power-flow formulations and allows flexible energy
resources — batteries, electric vehicles, heat pumps, photovoltaic systems,
and wind generators — to be attached to a multi-period model as composable
mix-in objects.

A typical workflow has four stages: (1) load a benchmark distribution network
from SimBench [@meinecke2020simbench]; (2) instantiate an OPF model class,
which automatically maps pandapower bus, line, and generator tables to Pyomo
sets and parameters; (3) call `add_OPF()` to activate operational constraints
and an objective function; and (4) invoke `solve()` with a selected solver
(IPOPT [@wachter2006ipopt], GLPK, CBC, or Gurobi). After solving, results are
written back to `net.res_*` DataFrames through a post-processing step, keeping
the interface consistent with the broader pandapower ecosystem.

The library ships with fifteen self-contained example scripts covering
single-period AC/DC OPF, multi-period planning with battery storage, hosting
capacity analysis with binary siting variables, feasible-operation-region
computation, and solver benchmarking across SimBench networks.

# Statement of Need

<!-- NEEDS REVIEW: The following describes the motivation. Please verify that
     the framing accurately reflects the intended contribution and that the
     cited related tools/papers are correct and complete. -->

The increasing penetration of distributed energy resources (DER) in low- and
medium-voltage distribution grids requires planning and operation tools that
can jointly optimise power dispatch and network constraints over multiple time
steps. Single-snapshot AC power-flow solvers — the standard output of tools
such as pandapower — are insufficient for problems where energy storage or
demand flexibility couples decisions across hours or days.

Existing OPF frameworks either target transmission-level networks
(e.g., MATPOWER, PowerModels.jl) or are designed for research prototyping
without a reusable, composable interface for device flexibility. `potpourri`
fills this gap by providing a clean Python library that:

1. Bridges pandapower's network-data model and Pyomo's algebraic modelling
   language, so users need not manually translate bus admittance matrices into
   optimisation variables.
2. Supports multi-period time horizons with configurable step sizes,
   enabling joint optimisation over a full operating day (e.g., 96 × 15 min).
3. Offers a composable device architecture in which flexible resources are
   instantiated as independent objects that attach their own Pyomo constraints
   and variables to an existing model — rather than requiring a monolithic
   class hierarchy.
4. Is validated against pandapower's Newton-Raphson power-flow solver on
   standard SimBench benchmark networks.

The library is primarily intended for power-systems researchers and
distribution-grid engineers who are familiar with pandapower and wish to add
optimisation capabilities without writing Pyomo models from scratch.

# Package Overview

## Class hierarchy

`potpourri` is structured around two parallel class hierarchies:

**Single-period models** (`src/potpourri/models/`):

- `Basemodel` — creates a Pyomo `ConcreteModel`, maps pandapower DataFrames to
  Pyomo sets and parameters, and provides a `solve()` method.
- `AC` / `DC` — extend `Basemodel` with AC (full complex power flow) or DC
  (linearised) power-flow equations.
- `OPF` — adds generator and load limits, line-loading limits, and objective
  functions.
- `ACOPF = AC + OPF`, `DCOPF = DC + OPF` (multiple inheritance).
- `HC_ACOPF` — hosting-capacity variant with binary siting variables and
  grid-code reactive-power constraints.

**Multi-period models** (`src/potpourri/models_multi_period/`):

- `Basemodel_multi_period` — adds a time index and integrates SimBench
  load/generation profiles.
- `ACOPF_multi_period` — full multi-period AC OPF.
- Device modules: `Battery_multi_period`, `HeatPump_multi_period`,
  `PV_multi_period`, `Windpower_multi_period`, `Demand_multi_period`,
  `Sgens_multi_period`, `Generator_multi_period`.

## Composable device interface

Flexible devices are composed, not inherited. Each device module is
instantiated with the pandapower network and a time horizon, and attaches its
own Pyomo `Sets`, `Params`, `Vars`, and `Constraints` to the parent model
object:

```python
import simbench as sb
from potpourri.models_multi_period.ACOPF_multi_period import ACOPF_multi_period
from potpourri.technologies.battery import Battery_multi_period

net = sb.get_simbench_net("1-LV-urban6--0-sw")
opf = ACOPF_multi_period(net, toT=96)      # 96 × 15-min steps = 1 day
battery = Battery_multi_period(opf.net, T=96, scenario=1)
battery.get_all(opf.model)                 # attaches battery constraints
opf.add_OPF()
opf.add_voltage_deviation_objective()
opf.solve(solver="ipopt")
```

## Solver support

`potpourri` uses Pyomo's `SolverFactory` and is therefore solver-agnostic.
Supported solvers include IPOPT [@wachter2006ipopt] for nonlinear AC OPF,
GLPK and CBC for linear/mixed-integer DC OPF and hosting-capacity problems,
Gurobi for commercial MIP/NLP, and NEOS for remote access without a local
solver installation.

# Acknowledgements

<!-- NEEDS REVIEW: Add funding sources, grants, or institutional support if
     applicable. Example: "This work was supported by [grant name / funder]." -->

The authors thank the contributors to pandapower [@thurner2018pandapower],
Pyomo [@bynum2021pyomo], and SimBench [@meinecke2020simbench], on which
`potpourri` depends.

# References
