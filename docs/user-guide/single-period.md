# Single-Period OPF

This guide walks through a complete AC Optimal Power Flow run using a SimBench low-voltage network.

## Choosing a model

| Class | Formulation | Use case |
|---|---|---|
| `ACOPF` | Full nonlinear AC | General-purpose OPF |
| `DCOPF` | Linearised DC | Fast screening, no reactive power |
| `HC_ACOPF` | AC + binary wind placement | Hosting capacity studies |

## Step 1 — Load a network

Any pandapower network works. SimBench provides realistic distribution grids with pre-configured load and generation profiles:

```python
import simbench as sb
import pandapower as pp

net = sb.get_simbench_net("1-LV-rural1--0-sw")
```

You can also build a network manually with `pandapower.create_*` helpers.

## Step 2 — Configure limits

Operational limits are read directly from the pandapower network's element tables. Set them before constructing the model:

```python
# Voltage limits on buses
net.bus["max_vm_pu"] = 1.05
net.bus["min_vm_pu"] = 0.95

# Generator power limits
net.ext_grid["max_p_mw"] = 1000.
net.ext_grid["min_p_mw"] = -1000.
net.ext_grid["max_q_mvar"] = 1000.
net.ext_grid["min_q_mvar"] = -1000.

# Mark static generators as controllable
net.sgen["controllable"] = True
net.sgen["max_p_mw"] = net.sgen["p_mw"]
net.sgen["min_p_mw"] = 0.

# Line thermal limits (percent of rated current)
net.line["max_loading_percent"] = 100.
```

## Step 3 — Build the model

Pass the network to the model constructor. The constructor runs `pp.runpp()` internally and maps the network data to Pyomo sets and parameters:

```python
from potpourri.models.ACOPF_base import ACOPF

opf = ACOPF(net)
```

## Step 4 — Add OPF constraints and objective

```python
opf.add_OPF()                              # power and thermal limits
opf.add_voltage_deviation_objective()      # minimise (v - 1)²
```

Alternative objectives:

```python
opf.add_reactive_power_flow_objective()    # minimise Σ qsG²
```

## Step 5 — Solve

```python
opf.solve(solver="ipopt", print_solver_output=False)
```

Supported solvers:

| Solver | Type | Notes |
|---|---|---|
| `"ipopt"` | NLP | Default; best for AC OPF |
| `"mindtpy"` | MINLP | For discrete tap changers; requires `mip_solver` kwarg |
| `"neos"` | Remote | Uses NEOS server; set `neos_opt` kwarg |

## Step 6 — Access results

Results are written to the pandapower result tables on `opf.net`:

```python
print(opf.net.res_bus[["vm_pu", "va_degree"]])
print(opf.net.res_line[["p_from_mw", "loading_percent"]])
print(opf.net.res_sgen[["p_mw", "q_mvar"]])
```

## Tap changer optimisation

To include transformer tap ratios as continuous decision variables:

```python
opf.add_tap_changer_linear()
opf.solve(solver="ipopt")
```

For discrete tap positions (requires a MIP solver):

```python
opf.add_tap_changer_discrete()
opf.solve(solver="mindtpy", mip_solver="gurobi", print_solver_output=False)
```

## Hosting capacity analysis

`HC_ACOPF` maximises wind generation subject to Q-curve grid code constraints. If `sgen.wind_hc` is not set, candidate wind generators are automatically placed at every non-external-grid bus:

```python
from potpourri.models.HC_ACOPF import HC_ACOPF

hc = HC_ACOPF(net)
hc.add_OPF()
hc.solve(solver="ipopt")

# Binary y[w] indicates which wind locations are active
for w in hc.model.WIND_HC:
    print(f"Bus {hc.net.sgen.bus[w]}: y={hc.model.y[w].value:.0f}, "
          f"P={hc.model.psG[w].value * hc.model.baseMVA:.3f} MW")
```

To run a weighted wind-generation vs. loss objective:

```python
hc.add_loss_obj()
hc.model.eps.set_value(0.7)   # 70% weight on wind, 30% on loss reduction
hc.solve(solver="ipopt")
```
