# Load-Case OPF Analysis

**Script:** `scripts/loadcases.py`

Demonstrates two successive AC-OPF analyses on the SimBench LV rural1 network:
reactive-power minimisation at a single operating point, followed by
voltage-deviation minimisation for representative load cases.

---

## Purpose

Power system operators assess grid behaviour under standardised loading
scenarios called **load cases** (e.g. high load, low wind). This script
shows how to:

- Fix a snapshot from a time-series profile and run an AC-OPF that
  minimises reactive power exchange.
- Sweep through SimBench load cases, scale loads and generation accordingly,
  and solve the OPF to find the reactive dispatch that minimises voltage
  deviations from 1 p.u.

---

## Workflow

### Step 1 — Reactive-power minimisation (single snapshot)

```
load network (1-LV-rural1--0-sw)
select profile index 1190
set PV reactive-power envelope from 0.95 power factor
fix sgen active power (max_p_mw = min_p_mw = p_mw)
pp.runpp() → reference power flow
deep-copy network → constrain ext_grid P to pp result, Q to ±0.05 MVAR
ACOPF.add_OPF() + add_reactive_power_flow_objective()
solve → minimise Σ qsG²
print initial and optimised sgen set-points
```

The external grid active power is **pinned** to the pandapower result,
restricting the ext_grid to a small reactive-power band (±0.05 MVAR).
This forces the AC-OPF to dispatch reactive power from the PV generators
to satisfy the network reactive balance, and the objective drives their
reactive output towards zero.

### Step 2 — Voltage-deviation minimisation per load case

```
for case in ["lW", "hL"]:
    deep-copy profile-loaded network
    scale loads and generation with loadcase factors
    set slack voltage to loadcase value
    ACOPF.add_voltage_deviation_objective()  (no add_OPF())
    solve → minimise Σ_b (v[b] - 1)²
    print ext_grid dispatch and line 0 active power
```

`add_OPF()` is intentionally **omitted** for the case loop. The optimiser
is free to dispatch unlimited reactive power from all generators; it finds
the reactive set-point that best flattens the voltage profile. This is an
exploratory analysis, not an operational planning problem.

---

## Load cases

Standard SimBench load-case identifiers used in this script:

| Key | Meaning | Load | Wind | PV |
|---|---|---|---|---|
| `lW` | Low wind | 10 % | 100 % | 80 % |
| `hL` | High load | 100 % | varies | varies |

Access all available cases for any network:

```python
import simbench as sb
net = sb.get_simbench_net("1-LV-rural1--0-sw")
print(net.loadcases)
```

---

## Usage

```bash
conda activate potpourri_env
export NEOS_EMAIL="your@email.address"
python scripts/loadcases.py
```

To use a local solver, change both `solve()` calls from `solver="neos"` to
`solver="ipopt"`.

!!! note
    Step 1 may report **infeasible** with IPOPT on snapshot index 1190.
    At this high-PV operating point the tight ext_grid reactive-power limit
    (±0.05 MVAR) conflicts with the network's reactive balance requirements.
    The problem is feasible on NEOS where the solver explores a wider region
    from a different initialisation. Step 2 converges reliably with both
    IPOPT and NEOS.

---

## Example output

```
=== Step 1: reactive-power minimisation ===
  Converged: True
  Objective sum(qsG^2): 0.000017
  sgen P [pu] / Q [pu] after OPF:
    sgen 0: psG=0.00766  qsG=-0.00253
    sgen 1: psG=0.00784  qsG=-0.00259
    sgen 2: psG=0.00190  qsG=-0.00063
    sgen 3: psG=0.00230  qsG=-0.00076

=== Step 2: voltage-deviation per load case ===
  Case lW: converged=True  obj_v_deviation=0.000000
    ext_grid 0: P=-0.01192  Q=0.25487 [pu]
    line 0 pLfrom = 0.00012 [pu]
    sgen Q dispatch [pu]:
      sgen 0: qsG=-0.01683
      ...
  Case hL: converged=True  obj_v_deviation=0.000003
    ext_grid 0: P=0.02304  Q=-0.15417 [pu]
    line 0 pLfrom = 0.00116 [pu]
```

---

## Key API calls

| Call | Effect |
|---|---|
| `ACOPF(net)` | Build AC model; runs `pp.runpp()` internally |
| `add_OPF()` | Add voltage bounds, Q limits, line loading limits |
| `add_reactive_power_flow_objective()` | Minimise Σ qsG² → `model.obj_reactive` |
| `add_voltage_deviation_objective()` | Minimise Σ (v−1)² → `model.obj_v_deviation` |
| `solve(solver="ipopt")` | Solve and write results to `net.res_*` |
