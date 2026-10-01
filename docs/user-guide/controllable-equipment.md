# Controllable DSO Network Equipment

A distribution system operator does not only dispatch generators: it moves
the taps of its HV/MV transformers and switches capacitor banks and
reactors. This page shows how to make those decisions part of the OPF, so
the solver chooses them *together* with DER curtailment, storage and
reactive power.

Everything here is **opt-in**. A model built without the calls below is the
same NLP as before, with every transformer at its `tap_pos` and every shunt
at its `step`. The mathematics is in
[Mathematical Modelling, Section 8](../mathematical-modelling.md#8-controllable-dso-network-equipment);
the survey, design decisions and references are in
`docs/research/dso_controllable_equipment.md`.

---

## On-load tap changers (OLTC)

### 1. Check which transformers qualify

```python
import simbench as sb
from potpourri.models.oltc import oltc_eligibility

net = sb.get_simbench_net("1-MV-rural--0-sw")
print(oltc_eligibility(net)[["eligible", "reason"]])
```

SimBench networks arrive with complete tap data but **no
`tap_changer_type`**, and pandapower 3.x treats such a transformer as having
no tap changer at all: `tap_pos` is ignored by `pp.runpp`. The report says
so, and the fix is one line:

```python
net.trafo["tap_changer_type"] = "Ratio"
```

A transformer is eligible when it is in service, its changer is a
longitudinal `"Ratio"` changer on the `"hv"` or `"lv"` side without
`tap_step_degree`, without `tap_dependency_table` and without a second
changer, and `tap_neutral`, `tap_min < tap_max`, `tap_step_percent` and
`tap_min <= tap_pos <= tap_max` are integers. The pandapower `oltc` column
is a short-circuit flag and is not consulted.

### 2. Enable the control

```python
from potpourri.models.ACOPF_base import ACOPF

net.bus["min_vm_pu"], net.bus["max_vm_pu"] = 0.95, 1.05
opf = ACOPF(net)
opf.add_OPF()
opf.enable_oltc(mode="discrete")          # every eligible transformer
# opf.enable_oltc(transformers=[0], mode="continuous")   # a selection
opf.add_voltage_deviation_objective()
```

| `mode` | tap position variable | problem class | use it for |
|---|---|---|---|
| `"continuous"` | real in `[tap_min, tap_max]` | NLP (IPOPT) | planning, sensitivity, bounds, warm starts — **a relaxation**, not an implementable position |
| `"discrete"` | integer in `[tap_min, tap_max]` | nonconvex MINLP | the physical schedule |

Controlled and fixed transformers mix freely: the ones you do not name keep
their `tap_pos`.

### 3. Solve

```python
# global MINLP (Gurobi 12+; proves optimality on small cases)
opf.solve(solver="gurobi_direct_minlp", time_limit=120)

# or, with IPOPT only: relax -> round -> fix -> re-solve (a heuristic)
opf.solve_oltc_round_and_fix(solver="ipopt")
```

`solve(solver="ipopt")` on a discrete model raises: IPOPT would silently
solve the continuous relaxation and report a fractional position as if it
were a tap. Pass `relax_integrality=True` if that relaxation is what you
want. MindtPy (`solver="mindtpy"`) and NEOS Bonmin/Couenne are the other
MINLP routes; see [Solvers](solvers.md).

### 4. Read the results

```python
opf.net.res_trafo[["tap_pos", "tap_factor", "loading_percent", "vm_lv_pu"]]
opf.tap_schedule()        # positions (Series; DataFrame over time on a horizon)
opf.tap_operations()      # number of position changes, including the first move
```

`solve()` **does not modify `net.trafo.tap_pos`**; the input network stays
as you built it. To carry a decision over, e.g. for a pandapower power flow
or a time-series run:

```python
import pandapower as pp

checked = opf.apply_tap_positions(net)   # writes net.trafo.tap_pos (rounded)
pp.runpp(checked)
```

It warns if a position had to be rounded (continuous mode) or if pandapower
would ignore it (`tap_changer_type` is None).

### 5. Price tap movement

```python
opf.add_voltage_deviation_objective()
opf.penalize_tap_movement(cost=1e-4)   # per position change, objective units
```

The term is added to the active objective, so add the objective first.
Calling it again only changes the price, which makes a cost sweep cheap.

---

## Multi-period tap scheduling

```python
from potpourri.models_multi_period.ACOPF_multi_period import ACOPF_multi_period

mp = ACOPF_multi_period(net, fromT=40, toT=64)     # 24 x 15 min
mp.add_OPF()
mp.enable_oltc(
    mode="discrete",
    max_change_per_step=1,     # at most one position per 15 min
    max_operations=6,          # at most six operations over the horizon
)
mp.add_voltage_deviation_objective()
mp.penalize_tap_movement(cost=1e-4)
mp.solve_oltc_round_and_fix(solver="ipopt")

mp.tap_schedule()              # DataFrame: time step x transformer
mp.tap_operations()            # Series: operations per transformer
mp.map_to_net(t=50)            # res_trafo.tap_pos etc. for one step
mp.apply_tap_positions(net, t=50)
```

Movement is measured from the network's `tap_pos` into the first period
(`initial_tap="net"`); pass a mapping of positions or `None` to change that.
Both limits and the price act on the split
$k_\tau - k_{\tau-1} = u_\tau - d_\tau$, $u, d \ge 0$, so they stay linear.

On a horizon the global MINLP of SimBench-sized networks is slow (Gurobi
finds incumbents fast but needs many minutes to prove optimality even for
one step of a 15-bus feeder); `solve_oltc_round_and_fix` is the practical
default there and reports the relaxed objective as a bound in
`mp.rounding_info`.

---

## Switchable reactive compensation

pandapower's `net.shunt` is a constant admittance: `p_mw` and `q_mvar` are
per step at 1 p.u., so the bank consumes
$(p + jq)\,\mathrm{step}\,(V_n/\mathtt{vn\_kv})^2\, v^2$ — the voltage-squared
dependence of a real capacitor bank, which the model keeps.

```python
from potpourri.models.shunt_control import shunt_eligibility

print(shunt_eligibility(net)[["eligible", "reason"]])
opf.enable_shunt_control(mode="discrete")     # step in {0, ..., max_step}
opf.penalize_shunt_switching(cost=1e-4)
opf.solve_shunt_round_and_fix(solver="ipopt")
opf.net.res_shunt[["step", "q_mvar", "vm_pu"]]
opf.shunt_schedule(); opf.shunt_operations(); opf.apply_shunt_steps(net)
```

`max_change_per_step`, `max_operations` and `initial_step` work as for the
tap changers. A shunt with `step_dependency_table=True` is rejected rather
than approximated. Tap positions and shunt steps enabled on the same model
are rounded in one pass by either `solve_*_round_and_fix` method.

STATCOM/SVC-like continuous reactive support is already an `sgen` with
`p_mw = 0` and Q limits; inverter Volt/VAR capability is covered in
[Reactive-Power Control](reactive-power-control.md).

---

## OPF control versus a pandapower controller

pandapower's `DiscreteTapControl` is a *local* rule: set the tap, run a
power flow, step one position towards a voltage band, repeat. The OPF makes
the position a decision variable chosen against the whole objective and
every network limit at once, together with curtailment and reactive power.
They answer different questions, and the two can disagree legitimately: a
controller holding the LV bus of the substation at 1.0 p.u. knows nothing
about the feeder end that the OPF is keeping below 1.05 p.u. while it
curtails as little PV as possible.

A DSO workflow that uses both:

1. optimise the schedule with `ACOPF_multi_period` + `enable_oltc`;
2. `apply_tap_positions(net, t)` for each step and check it with `pp.runpp`
   (`opf.diagnose()` does this replay for you);
3. read the voltage the optimised position produces at the controller's
   bus — that is the set point a `DiscreteTapControl` or
   `ContinuousTapControl` would need to reproduce the schedule.

`scripts/oltc_voltage_control_demo.py` runs this comparison on the SimBench
LV rural1 network and its regulated distribution transformer.

---

## What is not supported

| configuration | behaviour |
|---|---|
| `tap_changer_type` None (SimBench default) | not eligible; set `"Ratio"` |
| `"Ratio"` with `tap_step_degree != 0`, `"Symmetrical"`, `"Ideal"` (phase shifters) | not eligible; see the research note for the designed extensions |
| `"Tabular"` / `tap_dependency_table=True` (tap-dependent impedance) | not eligible; a finite-state model is designed, not built |
| second tap changer (`tap2_*`), three-winding transformers | not eligible / not modelled |
| `step_dependency_table=True` on a shunt | not eligible |
| DC OPF | `enable_oltc` / `enable_shunt_control` raise `NotImplementedError` |
| network reconfiguration (`net.switch`) | not an OPF variable; design note only |

The legacy `add_tap_changer_linear()` / `add_tap_changer_discrete()` methods
still work unchanged but emit a `DeprecationWarning`: they free every
transformer regardless of its data and build the LV-side ratio differently
from pandapower.
