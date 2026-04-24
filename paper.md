---
title: 'potpourri: A Python-Package for Multi-period Optimal Power Flow in Active Distribution Grids'
tags:
  - Python
  - power systems
  - optimal power flow
  - distribution grid
  - multi-period optimisation
  - battery storage
  - hosting capacity
  - feasible operation region
  - solver benchmarking
  - pandapower
  - Pyomo
  - SimBench
authors:
  - name: Steffen Kortmann
    orcid: 0009-0001-2074-773X
    corresponding: true
    affiliation: 1
  - name: Andreas Bong
    orcid: 0009-0008-0634-1056
    affiliation: 1
  - name: Simon Braun
    orcid: 0009-0008-9080-143X
    affiliation: 1
  - name: Alexander Och
    orcid: 0009-0004-7221-6780
    affiliation: 1
  - name: Farah Nasr
    orcid: 0009-0002-9951-3537
    affiliation: 1
  - name: Philip Kvesic
    orcid: 0009-0005-7514-488X
    affiliation: 1
  - name: Nina Stumberger
    orcid: 0009-0000-4436-4138
    affiliation: 1
  - name: Andreas Ulbig 
    orcid: 0000-0001-5834-1842 
    affiliation: 1
affiliations:
  - name: Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University, Aachen, Germany
    index: 1
date: 24 April 2026
bibliography: paper.bib
---

# Summary

`potpourri` is an open-source Python library for **multi-period Optimal Power
Flow (OPF)** in active distribution grids. It provides a modular, Pyomo-based
[@bynum2021pyomo; @hart2011pyomo] optimisation layer that operates directly on
pandapower [@thurner2018pandapower] network objects. The library supports both
exact AC and linearised DC power-flow formulations and allows flexible energy
resources — including batteries, electric vehicles, heat pumps, photovoltaic
systems, and wind generators — to be attached to a multi-period model as
composable mix-in objects.

A typical workflow has four stages: (1) load a benchmark distribution network
from SimBench [@meinecke2020simbench]; (2) instantiate an OPF model class,
which maps pandapower bus, line, generator, and load tables to Pyomo sets and
parameters; (3) call `add_OPF()` to activate operational constraints and an
objective function; and (4) invoke `solve()` with a selected solver such as
IPOPT [@wachter2006ipopt], GLPK, CBC, or Gurobi. After solving, results are
written back to `net.res_*` DataFrames through a post-processing step, keeping
the interface consistent with the broader pandapower ecosystem.

The library ships with fifteen self-contained example scripts covering
single-period AC/DC OPF, multi-period planning with battery storage, hosting
capacity analysis with binary siting variables, feasible-operation-region
computation, and solver benchmarking across SimBench networks.

# Statement of need

The increasing penetration of distributed energy resources (DERs) in low- and
medium-voltage distribution grids requires planning and operation tools that
can jointly optimise power dispatch and network constraints over multiple time
steps. Single-snapshot power-flow calculations, as commonly provided by tools
such as pandapower, are not sufficient for problems in which energy storage,
demand flexibility, or other inter-temporal constraints couple decisions across
hours or days.

Several established OPF frameworks exist, including MATLAB-based tools such as
`MATPOWER` [@zimmerman2011matpower] and Julia-based frameworks such as
`PowerModels.jl` [@coffrin2018powermodels]. However, there remains a need for a
Python-based framework that combines direct compatibility with pandapower
network objects, algebraic optimisation through Pyomo, and explicit support for
multi-period planning problems. Existing Python- and Pyomo-based OPF tools,
such as `oats` [@bukhsh2020oats], provide valuable modelling concepts but do
not offer the same direct integration with pandapower's network-data model or
the same structured interface for attaching multi-period flexible resources.

`potpourri` addresses this gap by providing a Python library that:

1. Bridges pandapower's network-data model and Pyomo's algebraic modelling
   language, so users do not need to manually translate network tables into
   optimisation variables and constraints.
2. Supports multi-period time horizons with configurable step sizes, enabling
   joint optimisation over full operating days, for example 96 time steps with
   15-minute resolution.
3. Provides a composable device architecture in which flexible resources are
   instantiated as independent objects that attach their own Pyomo constraints
   and variables to an existing model, rather than requiring a monolithic class
   hierarchy.
4. Is validated against pandapower's Newton-Raphson power-flow solver on
   standard SimBench benchmark networks.

The library is primarily intended for power-systems researchers and
distribution-grid engineers who are familiar with pandapower and want to add
optimisation capabilities without writing Pyomo models from scratch.

# Software design

## Class hierarchy

`potpourri` is structured around two parallel class hierarchies.

**Single-period models** (`src/potpourri/models/`):

- `Basemodel` creates a Pyomo `ConcreteModel`, maps pandapower DataFrames to
  Pyomo sets and parameters, and provides a `solve()` method.
- `AC` and `DC` extend `Basemodel` with AC, based on full complex power-flow
  equations, or DC, based on linearised power-flow equations.
- `OPF` adds generator and load limits, line-loading limits, and objective
  functions.
- `ACOPF = AC + OPF` and `DCOPF = DC + OPF` combine these model components
  through multiple inheritance.
- `HC_ACOPF` provides a hosting-capacity variant with binary siting variables
  and grid-code reactive-power constraints.

**Multi-period models** (`src/potpourri/models_multi_period/`):

- `Basemodel_multi_period` adds a time index and integrates SimBench load and
  generation profiles.
- `ACOPF_multi_period` implements full multi-period AC OPF.
- Device modules include `Battery_multi_period`, `HeatPump_multi_period`,
  `PV_multi_period`, `Windpower_multi_period`, `Demand_multi_period`,
  `Sgens_multi_period`, and `Generator_multi_period`.

## Composable device interface

Flexible devices are composed rather than inherited. Each device module is
instantiated with the pandapower network and a time horizon, and attaches its
own Pyomo `Sets`, `Params`, `Vars`, and `Constraints` to the parent model
object:

```python
import simbench as sb
from potpourri.models_multi_period.ACOPF_multi_period import ACOPF_multi_period
from potpourri.technologies.battery import Battery_multi_period

net = sb.get_simbench_net("1-LV-urban6--0-sw")
opf = ACOPF_multi_period(net, toT=96)      # 96 x 15-min steps = 1 day
battery = Battery_multi_period(opf.net, T=96, scenario=1)
battery.get_all(opf.model)                 # attaches battery constraints
opf.add_OPF()
opf.add_voltage_deviation_objective()
opf.solve(solver="ipopt")
```

## Solver support

`potpourri` uses Pyomo's `SolverFactory` and is therefore solver-agnostic.
Supported solvers include IPOPT [@wachter2006ipopt] for nonlinear AC OPF,
GLPK and CBC for linear and mixed-integer DC OPF and hosting-capacity problems,
Gurobi for commercial MIP and NLP workflows, and NEOS for remote solver access
without a local solver installation. While NEOS is not a solver, rather a cloud solver service, 
we discourage the continous use of NEOS for production-level optimisation.
However, it can be very useful for benchmarking and development, when the
installation of a local solver is not feasible.

# Research impact statement

The package has contributed to a range of research use cases in active
distribution grids. These include hosting-capacity calculation for renewable
generation and wind-power integration [@braun2024hcopf; @braun2025hcopf],
electric-vehicle integration and market-oriented charging studies
[@nasr2025ev; @bong2025mopf], uncertainty-aware voltage-control studies
[@kvesic2025mpc], and feasible-operation-region computation for distribution
grid flexibility [@och2025for] and TSO-DSO coordination [@kleinhelmkamp2025ofofor;
@kleinhelmkamp2025safeofo]. Across these applications, `potpourri` has served
as a reusable modelling layer for translating pandapower network data into
single- and multi-period Pyomo optimisation problems.

Beyond individual studies, `potpourri` provides shared research infrastructure
for developing, comparing, and extending OPF formulations for distribution-grid
operation and planning. It has been used in several research and teaching activities at the
Institute for High Voltage Equipment and Grids, Digitalization and Energy
Economics (IAEW), RWTH Aachen University. Since its initial development, 
the package has supported multiple Bachelor’s and Master’s theses as a modelling
framework for optimisation-based distribution-grid studies. Its
pandapower-compatible interface and modular Pyomo formulation have made it
suitable for student research projects as well as peer-reviewed publications
that require reproducible OPF-based analyses across benchmark and applied grid
models.

# AI usage disclosure

Generative AI tools, including ChatGPT by OpenAI and Claude Code, were used
during the development of this software, the preparation of repository
materials, and the writing and editing of this manuscript. In particular,
Claude Code supported the preparation of the public release workflow, including
repository cleanup, publication metadata, and documentation-related tasks. All
scientific content, software functionality, and final manuscript decisions were
reviewed by the authors.

# Acknowledgements

The development of `potpourri` was inspired in part by modelling concepts from
the `oats` package [@bukhsh2020oats]. Early development built on ideas and code
structures from `oats`; over time, the project evolved into a standalone
package focused on providing a pandapower-compatible interface to Pyomo-based
OPF formulations.

The modular and inheritance-based structure of `potpourri` was also influenced
by PowerModels.jl [@coffrin2018powermodels].

The authors thank the contributors to pandapower [@thurner2018pandapower],
Pyomo [@bynum2021pyomo], and SimBench [@meinecke2020simbench], on which
`potpourri` depends.

# References