# AC Power Flow Validation

**Script:** `scripts/validation.py`

Validates the `potpourri` AC power flow solver against
[pandapower](https://pandapower.readthedocs.io/)'s Newton-Raphson power flow
on six standard SimBench benchmark networks.

---

## Purpose

The `AC` class in `potpourri` formulates AC power flow as a Pyomo
optimisation problem and solves it with an NLP solver (IPOPT or NEOS).
This script provides a quantitative comparison of the results against
pandapower's reference solver.

---

## Workflow

```
for each network in NETS
    ┌─ load simbench network
    ├─ solve with potpourri AC (NEOS / IPOPT)
    ├─ solve with pandapower pp.runpp()
    ├─ compute error metrics on res_bus and res_line
    └─ collect results
print pivot tables and error ranking
```

1. **Load the network** via `simbench.get_simbench_net()`.
2. **Build and solve the `AC` model** — the constructor runs
   `pp.runpp()` internally to initialise voltages and angles; the NLP
   solver then finds the feasible AC power flow solution.
3. **Run the reference power flow** with `pp.runpp()` on the original
   network using voltage-independent loads.
4. **Compare** the following fields:

   | Table | Fields |
   |---|---|
   | `res_bus` | `vm_pu`, `va_degree` |
   | `res_line` | `pl_mw`, `ql_mvar` |

5. **Report** mean absolute error (MAE), maximum absolute error, RMSE,
   and standard deviation per field; identify the bus or line with the
   largest deviation.

---

## Networks tested

| Network ID | Voltage level | Notes |
|---|---|---|
| `1-HV-mixed--0-sw` | HV | mixed urban/rural, switches |
| `1-HV-urban--0-sw` | HV | urban, switches |
| `1-MV-rural--0-sw` | MV | rural |
| `1-MV-urban--0-sw` | MV | urban |
| `1-LV-urban6--0-sw` | LV | urban |
| `1-LV-rural1--0-sw` | LV | rural |

---

## Usage

Set the `NEOS_EMAIL` environment variable to a registered NEOS address,
then run:

```bash
conda activate potpourri_env
export NEOS_EMAIL="your@email.address"
python scripts/validation.py
```

To use a local IPOPT solver instead of NEOS, change line 87 of the script:

```python
# change
results = ac.solve(solver="neos")
# to
results = ac.solve(solver="ipopt")
```

---

## Example output

```
Detailed comparison:
                    net   element  variable  mean_abs   max_abs      rmse   std_abs  worst_index  solver_converged  vm_pu_close
0    1-HV-mixed--0-sw  res_bus     vm_pu  0.023850  0.042419  0.026011  0.010122            5              True        False
1    1-HV-mixed--0-sw  res_bus  va_degree  0.234032  0.739218  0.289415  0.173017           12              True        False
...

Compact pivot table (mean absolute error):
variable              pl_mw   ql_mvar  va_degree     vm_pu
net
1-HV-mixed--0-sw   0.012336  0.051528   0.234032  0.023850
1-LV-rural1--0-sw  0.002670  0.001040   7.052760  0.280162
1-MV-rural--0-sw   0.050133  0.027593   2.634818  0.035633

Networks ranked by total mean absolute error:
net
1-HV-mixed--0-sw     0.321746
1-MV-rural--0-sw     2.748177
1-LV-rural1--0-sw    7.336632
```

---

## Interpreting the results

The `AC` model treats static generator reactive power (`qsG`) as a **free
optimisation variable**, whereas pandapower fixes it to the value specified in
`net.sgen.q_mvar`. This structural difference means the NLP solver finds a
different (but equally valid) AC power flow solution in terms of reactive
power dispatch, leading to measurable voltage angle and magnitude differences
relative to pandapower.

Key observations:

- **HV networks** (larger impedances, stronger voltage coupling) show smaller
  relative errors than LV networks.
- **`vm_pu_close`** (within 1×10⁻³ p.u. of pandapower) is `False` for all
  networks due to the free reactive dispatch — this is expected behaviour.
- For pure power flow validation (fixed reactive dispatch), initialise the
  `qsG` variables from a prior `pp.runpp()` result using
  `init_pyo_from_pp_res`.

---

## Key functions

```python
calculate_error_metrics(reference, candidate)
```
Returns a dict with `mean_abs`, `max_abs`, `rmse`, `std_abs` for a pair of
`pd.Series`.

```python
compare_results(pp_net, ac_net, delta_keys)
```
Iterates over `delta_keys = {result_table: [column, ...]}`, calls
`calculate_error_metrics`, and records the index of the worst-case element.
