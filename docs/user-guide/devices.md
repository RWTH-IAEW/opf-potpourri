# Flexible Devices

Flexible devices are modular extensions to the multi-period model.  Each class
inherits from `Flexibility_multi_period` and attaches its own Pyomo Sets,
Parameters, Variables, and Constraints to an existing `ACOPF_multi_period`
instance by calling `device.get_all(model)`.

## Placement is reproducible

`Battery`, `PV` and `Heatpump` place their units at a random fraction of the
non-slack buses. The draw comes from a **generator owned by the device**, seeded
with `DEFAULT_PLACEMENT_SEED` (42) unless told otherwise, so a study re-runs to
the same answer:

```python
Battery_multi_period(net, T=96, scenario=1)               # same buses every run
Battery_multi_period(net, T=96, scenario=1, seed=7)       # a different placement
```

Never the process-global `np.random`, so a device's buses do not depend on
whatever else in the program drew a random number first.

For a Monte Carlo sweep, thread one generator through it. Each scenario then
gets an independent placement while the sweep as a whole replays exactly:

```python
import numpy as np

rng = np.random.default_rng(2024)
for scenario in range(50):
    battery = Battery_multi_period(net, T=96, scenario=1, rng=rng)
    ...
```

A single `seed` cannot express that — every scenario would land on the same
buses. `rng` takes precedence over `seed` when both are given.

The drawn buses are on `device.random_indexes`; record them with any result, so
a number carries the placement it came from.

## How a device reaches the power balance

A device attaches to an already-built model, so its injections have to be
registered rather than hard-coded into the power-flow equations:

```
device.get_all(model)
  └─ couple_to_power_balance(model)
       ├─ register_kcl_real(model, term)      # term(model, b, t) -> expression
       └─ register_kcl_reactive(model, term)  # AC-style models only

opf.add_OPF()
  └─ rebuild_kcl()   # re-creates KCL_real / KCL_reactive with the terms
```

Registered terms use the **load sign convention**: positive is consumption,
negative is injection, matching `pD` / `qD`. `Flexibility_multi_period.bus_term`
builds a term from a device's `*_bus` set and maps the placement buses into ppc
numbering, which is what the balance is indexed over.

`rebuild_kcl()` is idempotent — with no device registered it reproduces the
same constraints — so `add_OPF()` calls it unconditionally.

To couple a custom device, register a term in its `get_all`:

```python
def couple_to_power_balance(self, model):
    self.register_kcl_real(
        model, self.bus_term(model, model.MY_bus, "my_p", sign=1.0)
    )
```

Use `sign=-1.0` if the variable is a generator-convention injection.

---

## Battery storage

`Battery_multi_period` models grid-connected battery systems placed randomly
across non-slack network buses.

### Constructor

```python
Battery_multi_period(
    net,
    T=None,
    scenario=None,
    *,
    penetration=None,   # overrides scenario
    power_pu=0.006,
    soc_max=1.0,
    soc_min=0.2,
    capacity_pu_h=0.015,
    efficiency=0.9,
    initial_soc_fraction=0.5,
    terminal_soc="cyclic",
    s_inv_pu=None,          # defaults to power_pu
    q_control=None,
    var_q=0,
    cos_phi_min=None,
    grid_code=None,
)
```

| Parameter | Default | Description |
|---|---|---|
| `scenario` | — | Predefined penetration level (0–3, see table below) |
| `penetration` | `None` | Explicit % of non-slack buses; overrides `scenario` |
| `power_pu` | `0.006` | Symmetric charge/discharge limit (p.u. on `net.sn_mva`) |
| `soc_max` | `1.0` | Maximum state of charge (fraction of capacity) |
| `soc_min` | `0.2` | Minimum state of charge (fraction of capacity) |
| `capacity_pu_h` | `0.015` | Energy capacity (p.u. power × hours) |
| `efficiency` | `0.9` | One-way efficiency: `η` charging, `1/η` discharging, so the round-trip efficiency is `η²` |
| `initial_soc_fraction` | `0.5` | Initial SOC as fraction of `soc_max` |
| `terminal_soc` | `"cyclic"` | SOC at the last step: `"cyclic"` = back to the initial SOC, a float pins an absolute value, `None` leaves it free |
| `s_inv_pu` | `power_pu` | Converter apparent-power rating (p.u.). Bounds P and Q jointly through the S² circle |
| `q_control` | `None` | Grid-code capability area: `None`, `"qp"`, `"qu"`, or `"both"` |
| `var_q` | `0` | Grid-code operating variant (0–2) |
| `cos_phi_min` | `None` | Power-factor floor |
| `grid_code` | `None` | Connection rule supplying the capability areas (default VDE-AR-N 4120) |

### Reactive power

The battery converter supplies reactive power as well as active. `BAT_Q` uses
the generator sign convention — positive is capacitive injection, matching
`qsG` — and is bounded by, in increasing specificity:

| Constraint | Present when | Meaning |
|---|---|---|
| `bat_inverter_s2` | always | `BAT_P² + BAT_Q² ≤ s_inv_pu²`, the converter's apparent-power circle |
| `bat_cos_phi_upper` / `_lower` | `cos_phi_min` set | `\|BAT_Q\| ≤ tan(arccos(cos_phi_min)) · (BAT_Pchg + BAT_Pdis)` |
| `bat_QP_pos` / `bat_QP_neg` | `q_control` includes `"qp"` | Grid-code Q(P) capability area |
| `bat_QU_min` / `bat_QU_max` | `q_control` includes `"qu"` | Grid-code Q(U) capability area |

Notes on the formulation:

- `s_inv_pu` defaults to `power_pu`, a converter sized exactly to the active
  power limit — which leaves **no** reactive headroom at full charge or
  discharge. Oversize it (`s_inv_pu = 1.5 * power_pu`) for reactive support at
  rated active power. It may not be set below `power_pu`.
- The S² circle stays convex: `BAT_P` is the affine `BAT_Pchg − BAT_Pdis`.
  `BAT_Q` also carries the box `[−s_inv_pu, +s_inv_pu]`, which the circle
  already implies but which keeps IPOPT's iterates inside the capability
  region.
- The cos φ limit uses `BAT_Pchg + BAT_Pdis` as a convex stand-in for
  `|BAT_P|`. The two coincide whenever only one leg is active, which a
  loss-making `efficiency` already makes optimal.
- The Q(P) area is keyed on the **discharging** leg, since that is the mode in
  which a storage unit acts as a generating unit under VDE-AR-N
  4105/4110/4120.
- Reactive power exists on AC-style models only. On the DC formulation, which
  has no reactive balance, `BAT_Q` and these constraints are skipped.

**Scenario penetration levels:**

| `scenario` | % of non-slack buses |
|---|---|
| 0 | 1.0 % |
| 1 | 7.9 % |
| 2 | 9.9 % |
| 3 | 10.6 % |

### Pyomo components added

**Parameters:** `BAT_Pmax[b]`, `BAT_Pmin[b]`, `BAT_SOCmax[b]`, `BAT_SOCmin[b]`,
`BAT_Cap[b]`, `BAT_Eff[b]`, `BAT_SOC_init[b]`

**Variables:** `BAT_Pchg[b, t]`, `BAT_Pdis[b, t]` (both ≥ 0), `BAT_SOC[b, t]` (0–1),
`BAT_Q[b, t]` (AC-style models only)

**Expressions:** `BAT_P[b, t] = BAT_Pchg[b,t] − BAT_Pdis[b,t]`, the net
injection (charging positive, i.e. a load on the grid)

**Constraints:**

- `bat_chg_limit_con` / `bat_dis_limit_con` — per-leg bounds: `BAT_Pchg ≤ BAT_Pmax`, `BAT_Pdis ≤ −BAT_Pmin`
- `bat_power_con` — converter throughput: `BAT_Pchg[b,t] + BAT_Pdis[b,t] ≤ BAT_Pmax`, which also keeps the two legs from running at once
- `bat_soc_con` — SOC bounds: `BAT_SOCmin ≤ BAT_SOC[b,t] ≤ BAT_SOCmax`; initial condition `BAT_SOC[b,t₀] = BAT_SOC_init[b]`
- `bat_soc_update_con` — energy balance:

$$
e_{b,t} = e_{b,t-1} + \frac{\Delta t}{C_b}
          \left( \eta_b \, p^\text{chg}_{b,t}
                 - \frac{p^\text{dis}_{b,t}}{\eta_b} \right)
$$

- `bat_terminal_soc_con` — terminal condition, unless `terminal_soc=None`

!!! note "Why two power variables"
    A single signed power variable cannot carry a one-way efficiency: the same
    factor would scale both directions, so a charge/discharge cycle would
    return the SOC exactly to its starting value and `efficiency` would model
    no loss at all. Splitting the legs makes the round-trip efficiency `η²`,
    matching the single-period `Basemodel.add_storage()` block.

### Example

```python
from potpourri.models_multi_period.ACOPF_multi_period import ACOPF_multi_period
from potpourri.technologies.battery import Battery_multi_period

net = sb.get_simbench_net("1-LV-urban6--0-sw")
opf = ACOPF_multi_period(net, toT=96, fromT=0)

# scenario-based placement
battery = Battery_multi_period(opf.net, T=96, scenario=1)

# or with explicit parameters
battery = Battery_multi_period(
    opf.net, T=96,
    penetration=15.0,       # 15 % of buses
    power_pu=0.01,
    capacity_pu_h=0.025,
    efficiency=0.95,
    initial_soc_fraction=0.3,
)

battery.get_all(opf.model)
opf.add_OPF()
opf.solve(solver="ipopt")
```

!!! note "Per-unit values"
    All power and capacity parameters use the per-unit system on the network
    base `net.sn_mva` (MVA).  To convert: `power_pu = power_MW / net.sn_mva`.

---

## Heat pump

`Heatpump_multi_period` models electric heat pumps with a first-order building
thermal model.  Each heat pump consumes electrical power and raises the indoor
temperature; heat is lost proportionally to the load profile.

### Constructor

```python
Heatpump_multi_period(
    net,
    T=None,
    scenario=None,
    *,
    penetration=None,
    power_max_pu=0.005,
    cop=4.0,
    temp_max_c=20.0,
    temp_min_c=15.0,
    avg_house_size_m3=500.0,
    wall_thickness_m=0.2,
    qloss_max=0.01,
)
```

| Parameter | Default | Description |
|---|---|---|
| `scenario` | — | Predefined penetration level (0–3, see table below) |
| `penetration` | `None` | Explicit % of non-slack buses; overrides `scenario` |
| `power_max_pu` | `0.005` | Maximum electrical input power (p.u.) |
| `cop` | `4.0` | Coefficient of performance (heat out / electrical in) |
| `temp_max_c` | `20.0` | Upper indoor temperature bound (°C) |
| `temp_min_c` | `15.0` | Lower indoor temperature bound (°C) |
| `avg_house_size_m3` | `500.0` | Average house volume (m³) for thermal capacity calculation |
| `wall_thickness_m` | `0.2` | Wall thickness (m) for thermal mass |
| `qloss_max` | `0.01` | Peak heat-loss value (p.u. power) mapped to peak load |

**Scenario penetration levels:**

| `scenario` | % of non-slack buses |
|---|---|
| 0 | 6.3 % |
| 1 | 26.4 % |
| 2 | 44.6 % |
| 3 | 59.9 % |

### Pyomo components added

**Parameters:** `HP_Pmax[h]`, `HP_Pmin[h]`, `TempMax[h]`, `TempMin[h]`,
`HP_CoP[h]` *(mutable)*, `HP_ThermCap[h]` *(mutable)*, `Qloss[h, t]`

**Variables:** `hp_p[h, t]` (electrical consumption), `temp[h, t]` (indoor temperature)

**Constraints:**

- `hp_power_con` — power bounds: `HP_Pmin ≤ hp_p[h,t] ≤ HP_Pmax`
- `hp_temp_con` — temperature bounds: `TempMin ≤ temp[h,t] ≤ TempMax`
- `hp_temp_update_con` — thermal dynamics:

$$
\theta_{h,t} = \theta_{h,t-1} + \Delta t \cdot \frac{\text{CoP}_h \cdot p_{h,t}^\text{hp} - Q_{h,t}^\text{loss}}{C_h^\text{therm}}
$$

The thermal capacity $C_h^\text{therm}$ (MWh/K) is derived from the house
geometry using standard air and concrete physical constants.  `HP_CoP` and
`HP_ThermCap` are mutable Pyomo parameters so they can be updated after model
construction without a full rebuild.

### Example

```python
from potpourri.technologies.heat_pump import Heatpump_multi_period

opf = ACOPF_multi_period(net, toT=96, fromT=0)
hp = Heatpump_multi_period(
    opf.net, T=96, scenario=1,
    cop=3.5,
    temp_max_c=22.0,
    temp_min_c=18.0,
)
hp.get_all(opf.model)

# Update COP after construction (e.g. for seasonal sensitivity)
for h in opf.model.HP:
    opf.model.HP_CoP[h].set_value(2.8)  # winter COP
```

---

## PV generation

`PV_multi_period` adds controllable PV units placed randomly at non-slack buses.
Generation is upper-bounded by a time-varying profile from the SimBench
renewables table; curtailment is allowed down to `pv_pmin`.

### Constructor

```python
PV_multi_period(
    net,
    T=None,
    scenario=None,
    *,
    penetration=None,
    profile_column="PV5",
    pv_pmin=0.0,
)
```

| Parameter | Default | Description |
|---|---|---|
| `scenario` | — | Predefined penetration level (0–3, see table below) |
| `penetration` | `None` | Explicit % of non-slack buses; overrides `scenario` |
| `profile_column` | `"PV5"` | Column from `net.pv_load_profiles` (SimBench renewables) |
| `pv_pmin` | `0.0` | Minimum PV output (p.u.); `0.0` = full curtailment allowed |

**Scenario penetration levels:**

| `scenario` | % of non-slack buses |
|---|---|
| 0 | 13.4 % |
| 1 | 22.4 % |
| 2 | 24.4 % |
| 3 | 25.4 % |

### Pyomo components added

**Parameters:** `PV_Pmax[pv, t]` (available irradiance), `PV_Pmin[pv, t]`

**Variables:** `pPV[pv, t]` (dispatched PV power, ≥ 0)

**Constraint:** `PV_real_power_bounds` — `PV_Pmax[pv,t] ≤ pPV[pv,t] ≤ PV_Pmin[pv,t]`

### Example

```python
from potpourri.technologies.pv import PV_multi_period

opf = ACOPF_multi_period(net, toT=96, fromT=0)
pv = PV_multi_period(
    opf.net, T=96, scenario=0,
    profile_column="PV3",   # choose a different irradiance column
)
pv.get_all(opf.model)
```

---

## Wind power

`Windpower_multi_period` extends `Sgens_multi_period` with reactive-power
Q-curve constraints for wind generators based on the selected grid code.
It also supports hosting-capacity (HC) binary placement via
`y[w] ∈ {0, 1}`.

### Constructor

```python
Windpower_multi_period(
    net,
    T=None,
    scenario=None,
    *,
    sw_max_mva=10_000.0,
    sw_min_mva=0.0,
    qp_max=0.48,
    qp_min=-0.41,
)
```

| Parameter | Default | Description |
|---|---|---|
| `sw_max_mva` | `10 000` | Maximum apparent power per HC wind generator (MVA) |
| `sw_min_mva` | `0.0` | Minimum apparent power for an active HC generator (MVA) |
| `qp_max` | `0.484322` | Maximum Q/P ratio (capacitive) — widest variant of the default VDE-AR-N 4120 |
| `qp_min` | `−0.410775` | Minimum Q/P ratio (inductive) — widest variant of the same code |

The Q-curve characteristic is computed from module-level constants
`_VQU_V_POINTS` and `_VQU_Q_MAX` that encode the full three-variant
grid-code table.

### Usage

`Windpower_multi_period` is typically instantiated automatically by
`ACOPF_multi_period` when `net.bus` contains a `windpot_p_mw` column.
It can also be added manually for hosting-capacity studies:

```python
from potpourri.technologies.windpower import Windpower_multi_period

opf = ACOPF_multi_period(net, toT=96)
wind = Windpower_multi_period(
    opf.net,
    sw_max_mva=5.0,    # 5 MVA per candidate site
    qp_max=0.41,       # grid-code variant 1
    qp_min=-0.33,
)
wind.get_all_opf(opf.model)
```

---

## Implementing a custom device

Any class that inherits from `Flexibility_multi_period` and implements the
standard interface can be used as a flexible device:

```python
import pyomo.environ as pyo
from potpourri.technologies.flexibility import Flexibility_multi_period

class MyDevice(Flexibility_multi_period):
    def __init__(self, net, T=None, scenario=None, *, my_param=1.0):
        super().__init__(net, T, scenario)
        self.my_param = my_param
        # load device-specific data from net here

    def get_sets(self, model):
        model.MY_SET = pyo.Set(initialize=...)

    def get_parameters(self, model):
        model.MY_PARAM = pyo.Param(model.MY_SET, initialize=self.my_param)

    def get_variables(self, model):
        model.my_var = pyo.Var(
            model.MY_SET, model.T, domain=pyo.NonNegativeReals
        )

    def get_all_constraints(self, model):
        @model.Constraint(model.MY_SET, model.T)
        def my_constr(model, i, t):
            return model.my_var[i, t] <= model.MY_PARAM[i]

    def get_all(self, model):
        self.get_sets(model)
        self.get_parameters(model)
        self.get_variables(model)
        self.get_all_constraints(model)

    def get_all_acopf(self, model):
        pass  # add AC-OPF-specific components here if needed

    def get_all_opf(self, model):
        pass  # add OPF-specific components here if needed
```

Attach it to a model the same way as any built-in device:

```python
opf = ACOPF_multi_period(net, toT=96)
device = MyDevice(opf.net, T=96, my_param=0.5)
device.get_all(opf.model)
opf.add_OPF()
opf.solve(solver="ipopt")
```
