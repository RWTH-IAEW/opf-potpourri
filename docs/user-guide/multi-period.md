# Multi-Period OPF

The multi-period models extend single-period OPF with a time dimension, allowing optimisation over load and generation profiles (e.g. one day at 15-minute resolution).

## Overview

The time horizon is controlled by `fromT` and `toT` (0-indexed time step indices into the SimBench profiles). Each time step is 15 minutes (`deltaT = 0.25 h`).

Flexible devices — batteries, heat pumps — are instantiated separately and automatically attach their Pyomo sets, parameters, variables, and constraints to the parent model.

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

Instantiate device classes and pass the model to them. Each device attaches itself to `opf.model`. See [Flexible Devices](devices.md) for all options.

## Step 5 — Add OPF constraints and objective

```python
opf.add_OPF()
opf.add_voltage_deviation_objective()     # minimise Σ_t Σ_b (v[b,t] - 1)²
```

Other available objectives:

```python
opf.add_minimize_power_objective()        # minimise total load consumption
opf.add_generation_objective()            # minimise Σ pG²
opf.add_weighted_generation_objective()   # weighted: generators + sgens
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

To read the pandapower result tables instead, map one time step at a time —
`net.res_*` has no time dimension, so it holds a single step:

```python
for t in opf.model.T:
    opf.map_to_net(t)
    print(t, opf.net.res_bus.vm_pu.max())
```

`solve(to_net=True)` (the default) maps the **last** step of the horizon.
See [Solving models](solvers.md).

## Single-period vs multi-period differences

The two model kinds do not expose the same surface. Passing a multi-period
model an option it does not implement raises `TypeError` naming the option
and the class, rather than being silently ignored.

| `add_OPF` option | Single-period | Multi-period | Notes |
|---|---|---|---|
| `thermal_limit` | ✅ `"current"` / `"mva"` | ✅ `"current"` / `"mva"` | AC and LPAC only. The DC model is lossless and carries no reactive power, so it has no current-versus-MVA distinction to make and rejects the option. |
| `free_slack_vm` | ✅ default `True` | ✅ default `True` | The base AC power flow pins the slack magnitude at every time step; `add_OPF` frees it within `[Vmin, Vmax]` unless you pass `False`. |
| `angle_limits` | ✅ | ✅ | Enforced per time step, from `net.line.angmin_degree` / `angmax_degree`. |
| `grid_code`, `qu_deadband` | ✅ | ✅ | |
| `pv_q_control`, `inverter_s2`, `cos_phi_min`, `fixed_cos_phi`, `cos_phi_p_profile`, `pu_curtail` | ✅ as arguments | ⚠️ via `net.sgen` columns | Multi-period reads these from the `net.sgen` table (`var_q`, `sn_mva`, `pu_curtail`, `fixed_cos_phi`, `cos_phi_p_profile`) rather than from `add_OPF`. See [Reactive-power control](reactive-power-control.md). |
| `fix_hv_buses`, `hv_bus_kv` | ✅ | ❌ | No multi-period equivalent. |
| `sgen_types`, `wind_sgen_types` | ✅ | ❌ | No multi-period equivalent. |

Other behavioural differences to know about:

| Behaviour | Single-period | Multi-period |
|---|---|---|
| `solve(to_net=True)` | Writes `net.res_*` | Writes `net.res_*` for one time step (last by default); `map_to_net(t)` for others |
| `net.sgen.min_p_mw` | Honoured | Honoured (constant over the horizon; warns when it exceeds the profile) |
| Storage | `Basemodel.add_storage()` | `Battery_multi_period` device module, with reactive power bounded by the converter's S² circle |
| Storage reactive power | Box bounds via `qSTOR` | `BAT_Q` with S² circle, optional cos φ floor and grid-code Q(P)/Q(U) areas |

Everything in the model is **per-unit on `net.sn_mva`** — profiles, device
ratings, line and transformer limits alike. Device constructor arguments named
`*_pu` (`power_pu`, `capacity_pu_h`, `s_inv_pu`, `power_max_pu`) are already in
that system; convert from MW by dividing by `net.sn_mva`.

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
