# Multi-Period OPF

The multi-period models extend single-period OPF with a time dimension, allowing optimisation over load and generation profiles (e.g. one day at 15-minute resolution).

## Overview

The time horizon is controlled by `fromT` and `toT` (0-indexed time step indices into the SimBench profiles). Each time step is 15 minutes (`deltaT = 0.25 h`).

Flexible devices — batteries, EVs, heat pumps — are instantiated separately and automatically attach their Pyomo sets, parameters, variables, and constraints to the parent model.

## Step 1 — Load a SimBench network with profiles

Multi-period models require `net.profiles`, which SimBench provides:

```python
import simbench as sb
net = sb.get_simbench_net("1-LV-urban6--0-sw")
```

## Step 2 — Configure network limits

```python
net.bus["max_vm_pu"] = 1.05
net.bus["min_vm_pu"] = 0.95
net.line["max_loading_percent"] = 100.
net.sgen["controllable"] = True
net.sgen["max_p_mw"] = net.sgen["p_mw"]
net.sgen["min_p_mw"] = 0.
```

## Step 3 — Build the multi-period model

```python
from potpourri.models_multi_period.ACOPF_multi_period import ACOPF_multi_period

fromT = 0
toT = 96        # 1 day at 15-min resolution

opf = ACOPF_multi_period(net, toT=toT, fromT=fromT)
```

The constructor automatically creates device objects for loads (`Demand_multi_period`), static generators (`Sgens_multi_period`), and external grid generators (`Generator_multi_period`). Shunt capacitors and wind power objects are included when the network data supports them.

## Step 4 — Add flexible devices (optional)

Instantiate device classes and pass the model to them. Each device attaches itself to `opf.model`:

```python
# Electric vehicles (V1G + V2G)
# num_vehicles EV profiles are loaded from data/emob_profiles.pkl
opf2 = ACOPF_multi_period(net, toT=toT, fromT=fromT, num_vehicles=5)
```

Battery storage is added similarly — see [Flexible Devices](devices.md) for all options.

## Step 5 — Add OPF constraints and objective

```python
opf.add_OPF()
opf.add_voltage_deviation_objective()     # minimise Σ_t Σ_b (v[b,t] - 1)²
```

Other available objectives:

```python
opf.add_minimize_power_objective()        # minimise total load consumption
opf.add_generation_objective()            # minimise Σ pG²
opf.add_weighted_generation_objective()   # weighted: generators + EVs + sgens
opf.add_arbitrage_objective()             # minimise EV charging cost (day-ahead prices)
opf.add_charging_power_obj()              # maximise EV charging power
opf.add_discharging_power_obj()           # maximise EV discharging power
```

## Step 6 — Solve

```python
opf.solve(solver="ipopt", print_solver_output=False, time_limit=3600)
```

Multi-period problems are large NLPs. IPOPT with a `time_limit` is recommended. For MINLP problems (discrete tap changers), use MindtPy:

```python
opf.solve(solver="mindtpy", mip_solver="gurobi", time_limit=3600)
```

## Step 7 — Access results

Results are indexed by time step. Access Pyomo variables directly:

```python
import pyomo.environ as pyo

for t in opf.model.T:
    p_gen = sum(pyo.value(opf.model.psG[g, t]) * opf.model.baseMVA
                for g in opf.model.sG)
    print(f"t={t}: total sgen output = {p_gen:.3f} MW")
```

Or access EV-specific results:

```python
for v in opf.model.veh:
    for t in opf.model.T:
        print(f"v={v}, t={t}: SOC={pyo.value(opf.model.soc[v, t]):.2f}, "
              f"p={pyo.value(opf.model.p_opf[v, t]) * opf.model.baseMVA:.3f} MW")
```

## Model architecture

The multi-period model is composed from modular device objects. Each device implements a standard interface:

```
Device.__init__(net)
  └─ get_all(model)
       ├─ get_sets(model)
       ├─ get_parameters(model)
       ├─ get_variables(model)
       ├─ get_all_constraints(model)
       └─ get_all_acopf(model)   # AC-OPF specific additions
```

This means you can inspect or extend individual device constraints without touching the core model.
