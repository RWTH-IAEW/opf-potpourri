# PGLib-OPF Benchmark Validation

**Script:** `scripts/pglib_benchmark.py`

Benchmarks `potpourri` DC-OPF and AC-OPF against the published reference
objective values from the
[IEEE PES Power Grid Library (PGLib-OPF)](https://github.com/power-grid-lib/pglib-opf),
solved by [PowerModels.jl](https://github.com/lanl-ansi/PowerModels.jl) with
IPOPT.

---

## Purpose

PGLib-OPF is the standard benchmark suite for comparing OPF solvers and
formulations. Each case ships as a MATPOWER `.m` file with a published
reference objective value (minimum generation cost in $/h) for DC and AC
OPF under three scenario groups:

| Group | Abbreviation | Description |
|---|---|---|
| Typical Operations | TYP | Near-nominal loading, baseline scenario |
| Congested Operations | API | Elevated loading, active constraints |
| Small Angle Difference | SAD | Tight phase-angle bounds |

The benchmark script runs `potpourri` DC and AC OPF on a configurable subset
of TYP cases and reports the objective gap relative to the PowerModels.jl
reference.

---

## Setup

The PGLib-OPF case files are included as a git submodule. Initialise it once:

```bash
git submodule update --init benchmarks/pglib-opf
```

Or clone with submodules from the start:

```bash
git clone --recurse-submodules https://github.com/RWTH-IAEW/opf-potpourri.git
```

---

## Workflow

```
for each selected PGLib case (.m file)
    ┌─ load via load_pglib_case()
    │    (from_mpc → pandapower net; drop out-of-service generators together
    │     with their cost rows; keep negative-demand sgens fixed; attach
    │     ANGMIN/ANGMAX; put transformer taps on MATPOWER's from bus;
    │     rebalance the initial dispatch)
    ├─ DC-OPF: DCOPF + add_poly_cost_objective + angle_limits=True
    ├─ AC-OPF: ACOPF + add_poly_cost_objective
    │    + thermal_limit='mva'
    │    + free_slack_vm=True
    │    + fix_hv_buses=False
    │    + angle_limits=True
    └─ report objective vs BASELINE.md reference, write scripts/results/
```

The PGLib-compatible flags are passed directly to `ACOPF.add_OPF`; no
monkey-patching of the model is required. Building either model runs one
flat-start pandapower power flow first (`Basemodel` needs the `ppc` tables).
Where Newton-Raphson diverges, as it does on `case162_ieee_dtc`,
`case240_pserc` and `case300_ieee` with the shipped setpoints (PYPOWER
diverges on them too), `Basemodel` falls back to a DC power flow for the
tables and the starting point and logs a warning.

---

## Usage

The script has no command-line options. Its behaviour is set by the
configuration block at the top of the file:

```python
SOLVER = "ipopt"
MAX_BUSES = 300  # skip cases with more buses than this
RUN_DC = True    # include DC-OPF column
RUN_AC = True    # include AC-OPF column
CASES = None     # None → all TYP cases within MAX_BUSES; or ["case5_pjm", "case14_ieee"]
RESULTS_DIR = Path(__file__).parent / "results"
```

Edit the constants and run:

```bash
conda activate potpourri_env
python scripts/pglib_benchmark.py
```

Results are written to `scripts/results/pglib_benchmark.csv` and
`scripts/results/pglib_benchmark.md` (the directory is git-ignored).

!!! note "Which IPOPT runs"
    Pyomo resolves the `ipopt` executable through its own configuration
    directory (`~/.pyomo/bin`, filled by `pyomo download-extensions`) before
    it searches `PATH`. A copy there is used even inside `potpourri_env`. To
    force the environment's pinned IPOPT, point `PYOMO_CONFIG_DIR` at an empty
    directory for the run.

---

## Results

Run on 2026-09-20 with the default configuration (all TYP cases up to 300
buses, DC and AC).

| Component | Version |
|---|---|
| potpourri | 0.5.1 plus the loader and power-flow fixes below (unreleased) |
| pandapower | 3.5.4 |
| Pyomo | 6.10.1 |
| IPOPT | 3.14.20 (conda-forge) |
| Python | 3.12.0 |
| PGLib-OPF | v23.07 |

| Case | Nodes | DC ($/h) | DC ref ($/h) | DC gap (%) | AC ($/h) | AC ref ($/h) | AC gap (%) |
|---|---|---|---|---|---|---|---|
| case3_lmbd | 3 | 5 693.80 | 5 695.90 | −0.04 | 5 812.64 | 5 812.60 | 0.00 |
| case5_pjm | 5 | 17 479.90 | 17 480.00 | 0.00 | 17 551.89 | 17 552.00 | 0.00 |
| case14_ieee | 14 | 2 051.53 | 2 051.50 | 0.00 | 2 178.08 | 2 178.10 | 0.00 |
| case24_ieee_rts | 24 | 61 001.24 | 61 001.00 | 0.00 | 63 352.20 | 63 352.00 | 0.00 |
| case30_as | 30 | 767.60 | 767.60 | 0.00 | 803.13 | 803.13 | 0.00 |
| case30_ieee | 30 | 7 506.48 | 7 472.80 | +0.45 | 8 208.52 | 8 208.50 | 0.00 |
| case39_epri | 39 | 136 816.15 | 136 890.00 | −0.05 | 138 415.56 | 138 420.00 | 0.00 |
| case57_ieee | 57 | 34 772.95 | 34 773.00 | 0.00 | 37 589.34 | 37 589.00 | 0.00 |
| case60_c | 60 | 90 700.00 | 90 700.00 | 0.00 | 92 693.67 | 92 694.00 | 0.00 |
| case73_ieee_rts | 73 | 183 003.72 | 183 000.00 | 0.00 | 189 764.08 | 189 760.00 | 0.00 |
| case89_pegase | 89 | 105 117.82 | 105 040.00 | +0.07 | 107 285.67 | 107 290.00 | 0.00 |
| case118_ieee | 118 | 93 152.38 | 93 101.00 | +0.06 | 97 213.61 | 97 214.00 | 0.00 |
| case162_ieee_dtc | 162 | 101 505.86 | 101 460.00 | +0.05 | 108 075.68 | 108 080.00 | 0.00 |
| case179_goc | 179 | 751 888.13 | 751 880.00 | 0.00 | 754 266.41 | 754 270.00 | 0.00 |
| case197_snem | 197 | 1.47 | 1.47 | 0.00 | 1.50 | 1.50 | −0.02 |
| case200_activ | 200 | 27 479.64 | 27 480.00 | 0.00 | 27 557.57 | 27 558.00 | 0.00 |
| case240_pserc | 240 | 3 270 857.31 | 3 271 400.00 | −0.02 | 3 329 670.06 | 3 329 700.00 | 0.00 |
| case300_ieee | 300 | 517 352.61 | 517 850.00 | −0.10 | 565 220.11 | 565 220.00 | 0.00 |

Gaps are `(potpourri − reference) / reference`; `0.00` means below 0.005 %
in magnitude.

What the table says:

- All 18 cases solve for both DC and AC. Every AC objective is within
  0.02 % of PowerModels.jl (the PGLib references are rounded to five
  significant digits, so this is the reference's own resolution). The DC
  objectives are within 0.1 % except `case30_ieee` at +0.45 %.
- Every individual `solve()` took less than 1.1 s of wall time except
  `case300_ieee` AC at 4.1 s (Pyomo translation included); the complete run
  took about 20 s on one core of an Intel Xeon Platinum 8568Y+.

### What it took to get there

The first rerun on 0.5.1 solved 15 of the 18 cases, with `case200_activ`
14 % and `case240_pserc` 5 % above the reference, and three cases that did
not even build a model. None of it was the OPF formulation; all of it was in
how the MATPOWER file became a pandapower network, or in the power flow
that model construction runs. Each item below is pinned by a unit test in
`tests/unit_tests/test_pglib_loader.py` and
`tests/unit_tests/test_base_powerflow_fallback.py`.

| Symptom | Cause | Fix |
|---|---|---|
| `case162_ieee_dtc`, `case240_pserc`, `case300_ieee` raised `LoadflowNotConverged` before any OPF | With the shipped setpoints the flat-start Newton-Raphson diverges; PYPOWER's NR, fast-decoupled and Gauss-Seidel solvers diverge on the same files | `Basemodel` falls back to a DC power flow for the `ppc` tables and the start and logs a warning |
| `case200_activ` +14 % on DC and AC alike | The loader dropped 11 out-of-service generators but left `net.poly_cost` untouched: cost curves shifted onto the wrong units, constant terms of dead units kept being charged | Cost rows are dropped and renumbered together with the generator rows |
| `case240_pserc` +4.9 % | `from_mpc` turns the two buses with negative demand (4.6 GW) into sgens; the loader's dispatch rebalance overwrote their setpoint with 0 because they have no `max_p_mw` | Only sgens with an active-power limit are rebalanced; negative-demand sgens keep their setpoint and are no longer marked controllable |
| `case162_ieee_dtc`, `case300_ieee` AC-OPF locally infeasible within `[Vmin, Vmax]` | MATPOWER's `TAP` acts on the from bus, `from_mpc` always encodes it on the high-voltage side; for 18 and 16 transformers the from bus is the low-voltage side, so the wrong diagonal of the admittance matrix was divided by `TAP²` (differences of 5–12 p.u.) | `load_pglib_case` sets `tap_side="lv"` on those transformers; the admittance matrix then matches MATPOWER's to 0.01 p.u., and the small residual gaps on `case24`, `case73` and `case89` vanished as well |

The console output for the first cases looks like this:

```
case3_lmbd                (    3 buses)  DC:     5693.80 (ref     5695.90,  -0.04%,  0.25s) ✓  AC:     5812.64 (ref     5812.60,  +0.00%,  0.35s) ✓
case5_pjm                 (    5 buses)  DC:    17479.90 (ref    17480.00,  -0.00%,  0.30s) ✓  AC:    17551.89 (ref    17552.00,  -0.00%,  0.28s) ✓
case14_ieee               (   14 buses)  DC:     2051.53 (ref     2051.50,  +0.00%,  0.24s) ✓  AC:     2178.08 (ref     2178.10,  -0.00%,  0.54s) ✓
```

---

## Key API

### `potpourri.benchmarks.load_pglib_case`

```python
from potpourri.benchmarks import load_pglib_case

net = load_pglib_case("case14_ieee")
```

Converts a PGLib `.m` file into a `pandapowerNet` with:

- Out-of-service generators removed, with their cost rows, so
  `net.poly_cost.element` still addresses the right units
- All generators flagged `controllable=True`; sgens that `from_mpc` created
  from negative demand stay fixed injections
- Polynomial cost coefficients in `net.poly_cost`
- Phase-angle limits (`angmin_degree` / `angmax_degree`) on `net.line` and
  `net.trafo`, read from the MATPOWER `ANGMIN` / `ANGMAX` fields
- Transformer taps on the side MATPOWER puts them (`align_tap_sides=True`)
- Initial dispatch rebalanced so the flat-start power flow that model
  construction runs has a chance to converge

### `potpourri.models.cost_objective.add_poly_cost_objective`

```python
from potpourri.models.cost_objective import add_poly_cost_objective

acopf = ACOPF(net)
acopf.add_OPF(thermal_limit="mva", free_slack_vm=True, angle_limits=True)
add_poly_cost_objective(acopf)
acopf.solve(solver="ipopt")
```

Wires `net.poly_cost` coefficients (`c2·P² + c1·P + c0`) into a Pyomo
objective over all controllable generators (`ext_grid`, `gen`, `sgen`).

### Reference baselines

```python
from potpourri.benchmarks import PGLIB_BASELINE_TYP, PGLIB_BASELINE_API, PGLIB_BASELINE_SAD

ref = PGLIB_BASELINE_TYP["pglib_opf_case14_ieee"]
print(ref["dc"], ref["ac"])   # 2051.5, 2178.1
```

---

## Known limitations

- **`case30_ieee` DC-OPF is +0.45 % above the reference** while its AC-OPF
  matches to five digits, and `case300_ieee` DC is −0.10 %. The remaining DC
  differences are formulation details of the two DC models (PowerModels'
  `DCPPowerModel` versus `potpourri`'s `DCOPF`, for example how branch
  limits and quadratic costs enter), not data conversion; the AC agreement
  shows the networks match.
- Cases above 300 buses are not part of the default run and have not been
  validated.
- Only the **Typical Operations (TYP)** baseline is used by default. API and
  SAD cases can be loaded via `load_pglib_case` but are not included in the
  automated benchmark table.
