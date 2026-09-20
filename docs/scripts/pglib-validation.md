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
    │    (from_mpc → pandapower net, attach ANGMIN/ANGMAX,
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
flat-start pandapower power flow first (`Basemodel` needs the `ppc` tables),
which is why the initial dispatch is rebalanced by the loader.

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
| potpourri | 0.5.1 (commit `d3a384a`) |
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
| case24_ieee_rts | 24 | 61 001.24 | 61 001.00 | 0.00 | 63 384.08 | 63 352.00 | +0.05 |
| case30_as | 30 | 767.60 | 767.60 | 0.00 | 803.13 | 803.13 | 0.00 |
| case30_ieee | 30 | 7 506.48 | 7 472.80 | +0.45 | 8 208.52 | 8 208.50 | 0.00 |
| case39_epri | 39 | 136 816.15 | 136 890.00 | −0.05 | 138 415.56 | 138 420.00 | 0.00 |
| case57_ieee | 57 | 34 772.95 | 34 773.00 | 0.00 | 37 589.34 | 37 589.00 | 0.00 |
| case60_c | 60 | 90 700.00 | 90 700.00 | 0.00 | 92 693.67 | 92 694.00 | 0.00 |
| case73_ieee_rts | 73 | 183 003.72 | 183 000.00 | 0.00 | 189 843.46 | 189 760.00 | +0.04 |
| case89_pegase | 89 | 105 117.82 | 105 040.00 | +0.07 | 107 259.47 | 107 290.00 | −0.03 |
| case118_ieee | 118 | 93 152.38 | 93 101.00 | +0.06 | 97 213.61 | 97 214.00 | 0.00 |
| case162_ieee_dtc | 162 | — | 101 460.00 | — | — | 108 080.00 | — |
| case179_goc | 179 | 751 888.13 | 751 880.00 | 0.00 | 754 236.69 | 754 270.00 | 0.00 |
| case197_snem | 197 | 1.47 | 1.47 | 0.00 | 1.50 | 1.50 | −0.01 |
| case200_activ | 200 | 31 474.92 | 27 480.00 | +14.54 | 31 474.92 | 27 558.00 | +14.21 |
| case240_pserc | 240 | — | 3 271 400.00 | — | — | 3 329 700.00 | — |
| case300_ieee | 300 | — | 517 850.00 | — | — | 565 220.00 | — |

Gaps are `(potpourri − reference) / reference`; `0.00` means below 0.005 %
in magnitude. A dash marks a case where no model could be built (see below).

What the table says:

- 15 of the 18 cases solve for both DC and AC. On 14 of them the AC objective
  is within ±0.05 % of PowerModels.jl; the DC objective is within ±0.07 % on
  13, with `case30_ieee` at +0.45 %.
- Every individual `solve()` took less than one second of wall time (Pyomo
  translation included); the complete run took about 12 s on one core of an
  Intel Xeon Platinum 8568Y+.
- The objectives are insensitive to the IPOPT patch level: a second pass with
  IPOPT 3.14.6 reproduced every value to 1e-9. They also match the run made in
  August 2026 on pandapower 3.4.0 with the 0.5.0 code, so none of the 0.5.1
  model fixes (which concern switches, transformer iron losses and storage
  units) touch these transmission cases.

The console output for the first cases looks like this:

```
case3_lmbd                (    3 buses)  DC:     5693.80 (ref     5695.90,  -0.04%,  0.20s) ✓  AC:     5812.64 (ref     5812.60,  +0.00%,  0.25s) ✓
case5_pjm                 (    5 buses)  DC:    17479.90 (ref    17480.00,  -0.00%,  0.28s) ✓  AC:    17551.89 (ref    17552.00,  -0.00%,  0.25s) ✓
case14_ieee               (   14 buses)  DC:     2051.53 (ref     2051.50,  +0.00%,  0.19s) ✓  AC:     2178.08 (ref     2178.10,  -0.00%,  0.25s) ✓
```

---

## Key API

### `potpourri.benchmarks.load_pglib_case`

```python
from potpourri.benchmarks import load_pglib_case

net = load_pglib_case("case14_ieee")
```

Converts a PGLib `.m` file into a `pandapowerNet` with:

- All generators flagged `controllable=True`
- Polynomial cost coefficients in `net.poly_cost`
- Phase-angle limits (`angmin_degree` / `angmax_degree`) on `net.line` and
  `net.trafo`, read from the MATPOWER `ANGMIN` / `ANGMAX` fields
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

- **Three cases do not get past model construction.** For `case162_ieee_dtc`,
  `case240_pserc` and `case300_ieee` the flat-start Newton–Raphson power flow
  that `Basemodel` runs to obtain the `ppc` tables stops with
  `Power Flow nr did not converge after 10 iterations!`, so no OPF is
  attempted. The rebalanced initial dispatch is not enough for these
  networks; a better starting point (for example a DC power-flow
  initialisation or more iterations) is the obvious next step.
- **`case200_activ` is 14 % above the reference for both DC and AC.** The two
  objectives coincide exactly, so the AC power-flow equations are not the
  cause; the difference sits in the data translation or the generator and
  cost bounds. The case has many `net.impedance` rows after `from_mpc`, and
  PowerModels.jl's input preprocessing (topology simplification, merging of
  low-impedance branches) has no counterpart in `potpourri`. Not yet
  explained.
- **`case30_ieee` DC-OPF is +0.45 % above the reference** while its AC-OPF
  matches to five digits. DC-OPF in PGLib uses linear costs only;
  `add_poly_cost_objective(..., allow_quadratic=True)` keeps quadratic terms,
  which is one candidate for the offset.
- Cases above 300 buses are not part of the default run and have not been
  validated.
- Only the **Typical Operations (TYP)** baseline is used by default. API and
  SAD cases can be loaded via `load_pglib_case` but are not included in the
  automated benchmark table.
