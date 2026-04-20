# Flexible Devices

Flexible devices are modular extensions to the multi-period model. Each device class inherits from `Flexibility_multi_period` and attaches its own Pyomo components to an existing `ACOPF_multi_period` (or `DCOPF_multi_period`) model.

## Battery storage

`Battery_multi_period` models grid-connected battery systems. Batteries are placed randomly across network buses according to the `scenario` parameter.

**Key parameters (hardcoded defaults):**

| Parameter | Value | Description |
|---|---|---|
| `bat_power` | 0.006 p.u. | Charge/discharge power limit |
| `bat_cap` | 0.015 p.u. | Energy capacity (MWh) |
| `bat_soc_max` | 1.0 | Maximum state of charge |
| `bat_soc_min` | 0.2 | Minimum state of charge |
| `bat_eff` | 0.9 | Round-trip charging efficiency |

**Variables:**

- `BAT_P[b, t]` — battery power (positive = charging, negative = discharging)
- `BAT_SOC[b, t]` — state of charge (0–1)

**SOC dynamics constraint:**

```
BAT_SOC[b, t] = BAT_SOC[b, t-1] + ΔT · BAT_P[b, t] · η / BAT_Cap[b]
```

**Initial condition:** `BAT_SOC[b, t=0] = 0.5 · BAT_SOCmax[b]`

**Usage:**

```python
from potpourri.technologies.battery import Battery_multi_period

opf = ACOPF_multi_period(net, toT=96)
battery = Battery_multi_period(opf.net)
battery.get_all(opf.model)
opf.add_OPF()
opf.solve(solver="ipopt")
```

---

## Electric vehicles

`EV_multi_period` models a fleet of EVs with per-vehicle driving profiles, charging point assignments, and optional vehicle-to-grid (V2G) capability.

EV profiles are loaded from `data/emob_profiles.pkl`. Day-ahead electricity prices for the arbitrage objective come from `data/da_prices_hourly_2022.xlsx`.

**Vehicle types:**

| Type | Set | Description |
|---|---|---|
| V1G (unidirectional) | `veh_v1g` | Charging only |
| V2G (bidirectional) | `veh_v2g` | Charging and discharging |
| Non-controllable | `veh_nc` | Follows pre-computed charging profile |

**Key variables:**

- `soc[v, t]` — vehicle state of charge (0–1)
- `p_opf[v, t]` — optimised net charging power (per-unit)
- `p_charging[v, t]` — charging power component (≥ 0)
- `p_discharging[v, t]` — discharging power (≤ 0, V2G only)
- `buffer_soc[v, t]` — SOC slack for departure constraint relaxation

**SOC constraints:**

- SOC must stay above `min_soc_parking` while parked
- SOC must reach `min_soc_departing` before each departure
- Net energy over the horizon must be equal (cyclic condition)

**Usage:**

```python
from potpourri.models_multi_period.ACOPF_multi_period import ACOPF_multi_period

# Pass num_vehicles to the constructor
opf = ACOPF_multi_period(net, toT=96, num_vehicles=10)
opf.add_OPF()
opf.add_arbitrage_objective()    # minimise charging cost
opf.solve(solver="ipopt")
```

Enable market constraints (constant aggregated hourly power):

```python
opf.model.ev_object.get_market_constraints(opf.model)
```

---

## Heat pump

`Heatpump_multi_period` models thermal loads with a thermal storage buffer. It uses a simple building thermal model with heat capacity and heat loss parameters.

**Key parameters:**

- `hp_power` — heat pump electrical power (p.u.)
- `COP` — coefficient of performance
- Heat capacity and thermal mass derived from house geometry

**Variables:**

- `hp_p[h, t]` — heat pump electrical power consumption
- `temp[h, t]` — indoor temperature (°C)

**Usage:**

```python
from potpourri.technologies.heat_pump import Heatpump_multi_period

opf = ACOPF_multi_period(net, toT=96)
hp = Heatpump_multi_period(opf.net)
hp.get_all(opf.model)
opf.add_OPF()
opf.solve(solver="ipopt")
```

---

## PV generation

`PV_multi_period` adds controllable PV generators. Generation profiles are taken from `net.profiles`.

**Usage:**

```python
from potpourri.technologies.pv import PV_multi_period

opf = ACOPF_multi_period(net, toT=96)
pv = PV_multi_period(opf.net)
pv.get_all(opf.model)
```

---

## Wind power

`Windpower_multi_period` extends `Sgens_multi_period` with hosting-capacity-style Q-curve constraints for wind generators, taken from grid code requirements. Reactive power control variants (0–3) map to different Q-P and Q-U slopes.

**Usage:**

```python
from potpourri.technologies.windpower import Windpower_multi_period

# Windpower is typically instantiated automatically by ACOPF_multi_period
# when net.bus has a windpot_p_mw column. It can also be added manually:
opf = ACOPF_multi_period(net, toT=96)
wind = Windpower_multi_period(opf.net)
wind.get_all(opf.model)
```

---

## Implementing a custom device

Any class that implements the following interface can be used as a flexible device:

```python
from potpourri.technologies.flexibility import Flexibility_multi_period

class MyDevice(Flexibility_multi_period):
    def __init__(self, net, T=None, scenario=None):
        super().__init__(net, T, scenario)
        # load device-specific data from net

    def get_sets(self, model):
        model.MY_SET = pyo.Set(initialize=...)

    def get_parameters(self, model):
        model.MY_PARAM = pyo.Param(model.MY_SET, initialize=...)

    def get_variables(self, model):
        model.my_var = pyo.Var(model.MY_SET, model.T, domain=pyo.NonNegativeReals)

    def get_all_constraints(self, model):
        def my_constraint_rule(model, i, t):
            return model.my_var[i, t] <= model.MY_PARAM[i]
        model.my_constr = pyo.Constraint(model.MY_SET, model.T, rule=my_constraint_rule)

    def get_all_acopf(self, model):
        pass  # add AC-OPF specific components here

    def get_all(self, model):
        self.get_sets(model)
        self.get_parameters(model)
        self.get_variables(model)
        self.get_all_constraints(model)

    # also implement get_all_opf, get_all_ac if needed
```
