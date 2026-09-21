# PGLib-OPF Benchmark Validation

**Script:** `scripts/pglib_benchmark.py`

Benchmarks `potpourri` DC-OPF and AC-OPF against the published reference
objective values from the
[IEEE PES Power Grid Library (PGLib-OPF)](https://github.com/power-grid-lib/pglib-opf),
solved by [PowerModels.jl](https://github.com/lanl-ansi/PowerModels.jl) with
IPOPT, for all three operating conditions and every case size.

---

## Purpose

PGLib-OPF is the standard benchmark suite for comparing OPF solvers and
formulations. Each case ships as a MATPOWER `.m` file with a published
reference objective value (minimum generation cost in $/h) for DC and AC
OPF under three scenario groups, 66 cases each:

| Group | Abbreviation | Description | What it exercises in potpourri |
|---|---|---|---|
| Typical Operations | TYP | Near-nominal loading | power-flow equations, generator limits |
| Congested Operations | API | Loads scaled up until branches bind | thermal limits |
| Small Angle Difference | SAD | Tight phase-angle-difference bounds | branch angle limits |

The API and SAD groups matter because the constraints they make binding stay
slack in the TYP cases; a TYP-only table validates neither the thermal nor the
angle-limit code paths.

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
for each group (typ, api, sad) and each case (.m file), in parallel workers
    ┌─ load via load_pglib_case()
    │    (from_mpc → pandapower net; out-of-service generators dropped with
    │     their cost rows; negative-demand sgens kept fixed; ANGMIN/ANGMAX
    │     on lines, transformers and impedance elements; transformer taps
    │     on MATPOWER's from bus; reference bus guaranteed; dispatch rebalanced)
    ├─ DC-OPF: DCOPF(dc_convention='powermodels') + add_poly_cost_objective
    │    + angle_limits=True
    ├─ AC-OPF: ACOPF + add_poly_cost_objective
    │    + thermal_limit='mva' + free_slack_vm=True + fix_hv_buses=False
    │    + angle_limits=True
    ├─ a solve that does not reach optimality is repeated once from a
    │    flat, mid-range start (`RETRY_STARTS`)
    └─ report objective vs BASELINE.md reference, write scripts/results/
```

The AC-OPF starts from the network's own dispatch. The congested (API)
files leave that setpoint far from anything feasible, and IPOPT then
converges to a locally infeasible point although the problem has a
solution: `case179_goc__api` reaches the published 1 883 400 $/h from a
flat, mid-range start. A solve that does not reach optimality is therefore
repeated once from that start, and the `ac_start` column of the CSV says
which of the two produced the reported result. Nothing about the
optimisation problem changes, only the point the solver starts from.

The PGLib-compatible flags are passed directly to `ACOPF.add_OPF`; no
monkey-patching of the model is required. Building either model runs one
flat-start pandapower power flow first (`Basemodel` needs the `ppc` tables).
Where Newton-Raphson diverges with the shipped setpoints (PYPOWER diverges on
the same files), `Basemodel` falls back to a DC power flow for the tables and
the starting point and logs a warning; where no power flow can run at all
because a branch has zero reactance (`case1803_snem`), it takes pandapower's
network tables with a flat start.

The DC column follows PowerModels' `DCPPowerModel`, which produced the PGLib DC
references: branch susceptance `−x/(r² + x²)` and a flow `p = −b (θ_from − θ_to)`
in which neither the tap nor the transformer phase shift appears. potpourri's
default is MATPOWER's `−1/x` with the phase shift; the two agree where r ≪ x and
no phase shifters exist, and differ by several percent otherwise (up to −2.8 %
on API cases, +6 % on the RTE cases with their 17 phase shifters).
`DCOPF(net, dc_convention="powermodels")` selects the convention; the model's
default is unchanged.

---

## Usage

The script has no command-line options. Its behaviour is set by the
configuration block at the top of the file:

```python
SOLVER = "ipopt"
GROUPS = ("typ", "api", "sad")  # PGLib operating conditions to run
MAX_BUSES = None  # None → every case; an int skips cases with more buses
RUN_DC = True
RUN_AC = True
CASES = None  # None → all cases of each group; or ["case5_pjm", "case14_ieee"]
N_WORKERS = 8  # parallel worker processes, one solver thread each
TIME_LIMIT_S = 3600  # IPOPT wall-time limit per solve
RETRY_STARTS = ("midrange",)  # retry starts for a failed AC solve; "dc" also available
DC_CONVENTION = "powermodels"  # DC linearisation convention, see DCOPF
```

Edit the constants and run:

```bash
conda activate potpourri_env
python scripts/pglib_benchmark.py
```

Results are written to `scripts/results/pglib_benchmark.csv` and one
`scripts/results/pglib_benchmark_<group>.md` per group (git-ignored). The
largest cases are dispatched first so they overlap with the many small ones.
A full run takes about four and a half hours with 24 workers (the run below raised
`N_WORKERS` to 24 on a 64-core host; the shipped default is 8);
`MAX_BUSES = 300` gives a
two-minute pass. The largest cases need a worker of their own for up to an
hour of model construction before the solver starts.

!!! note "Which IPOPT runs"
    Pyomo resolves the `ipopt` executable through its own configuration
    directory (`~/.pyomo/bin`, filled by `pyomo download-extensions`) before
    it searches `PATH`. A copy there is used even inside `potpourri_env`. To
    force the environment's pinned IPOPT, point `PYOMO_CONFIG_DIR` at an empty
    directory for the run.

---

## Results

Run on 2026-09-21 with the default configuration: all three groups, every case
size, DC and AC.

| Component | Version |
|---|---|
| potpourri | 0.5.3 (commit `c37e259`) |
| pandapower | 3.5.4 |
| Pyomo | 6.10.1 |
| IPOPT | 3.14.20 (conda-forge), one thread per solve, 24 solves in parallel |
| Python | 3.12.0 |
| PGLib-OPF | v23.07 |
| Host | Intel Xeon Platinum 8568Y+, 393 GB |

Every row comes from one code version. An earlier pass over the same 198
cases exposed three further defects (the solver time limit never reached
IPOPT, the DC convention kept pandapower's tap-referred transformer
impedance, and static generators were clipped at zero), which between them
changed results for about a third of the cases, so the table was recomputed
from scratch rather than patched case by case. The findings are in the
table at the end of this section.

### Accounting

`case78484_epigrids` (TYP, API and SAD) is not in the tables below. Its DC-OPF solves, but
building the AC-OPF of 78 484 buses takes over an hour and the solve did not
reach a conclusion within the time budget of this run, in any of the three
groups; the run was stopped rather than reported with a guess. The three
cases are the only ones missing: 195 of the 198 are below.

Five cases do not reach optimality from the network's own setpoint and solve
from the mid-range start: `case179_goc__api`, `case4837_goc__api`,
`case9591_goc__api`, `case10000_goc__api` and `case10480_goc`. A run with the
retry switched off misses all five. The benchmark also tries the DC-OPF
solution as a second start; it rescued nothing here and costs a full extra
model build on the largest cases, so `RETRY_STARTS` ships with the mid-range
start alone.

Every one of the 195 cases in the table produced a DC and an AC result row.
`✓` marks a
solve that agrees with the reference: an optimal solution where PowerModels
reports a value, or an infeasible outcome where PowerModels reports `inf.`.

| Group | Model | Solved | Infeasible, reference `inf.` | Solved, reference `inf.` | Not solved, finite reference | Time limit | Errors | median gap | max gap |
|---|---|---|---|---|---|---|---|---|---|
| TYP | DC | 65 | 0 | 0 | 0 | 0 | 0 | 0.001 % | 0.28 % |
| TYP | AC | 64 | 0 | 0 | 1 | 0 | 0 | 0.001 % | 1.79 % |
| API | DC | 65 | 0 | 0 | 0 | 0 | 0 | 0.001 % | 1.01 % |
| API | AC | 63 | 0 | 0 | 2 | 0 | 0 | 0.002 % | 2.47 % |
| SAD | DC | 20 | 45 | 0 | 0 | 0 | 0 | 0.001 % | 0.00 % |
| SAD | AC | 64 | 0 | 0 | 1 | 0 | 0 | 0.002 % | 3.04 % |

Cases that do not agree with the reference, or differ by more than 0.1 %:

| Group | Model | Case | Finding |
|---|---|---|---|
| TYP | DC | case197_snem | gap −0.28 % |
| TYP | AC | case2742_goc | not solved, reference 275 710.00 |
| TYP | AC | case197_snem | gap −0.11 % |
| TYP | AC | case10192_epigrids | gap +1.79 % |
| API | DC | case1803_snem__api | gap +0.56 % |
| API | DC | case5658_epigrids__api | gap +0.89 % |
| API | DC | case10192_epigrids__api | gap +1.01 % |
| API | AC | case2742_goc__api | not solved, reference 609 960.00 |
| API | AC | case24464_goc__api | not solved, reference 2 684 000.00 |
| API | AC | case2736sp_k__api | gap +0.24 % |
| API | AC | case2746wop_k__api | gap −0.15 % |
| API | AC | case3375wp_k__api | gap +0.37 % |
| API | AC | case4837_goc__api | gap +0.24 % |
| API | AC | case5658_epigrids__api | gap +0.89 % |
| API | AC | case6470_rte__api | gap +0.14 % |
| API | AC | case7336_epigrids__api | gap +0.12 % |
| API | AC | case9591_goc__api | gap −2.47 % |
| API | AC | case10192_epigrids__api | gap +0.33 % |
| API | AC | case10480_goc__api | gap +0.10 % |
| SAD | AC | case2742_goc__sad | not solved, reference 275 710.00 |
| SAD | AC | case2383wp_k__sad | gap +0.33 % |
| SAD | AC | case2736sp_k__sad | gap +0.41 % |
| SAD | AC | case2737sop_k__sad | gap +0.45 % |
| SAD | AC | case4837_goc__sad | gap −0.18 % |
| SAD | AC | case9591_goc__sad | gap −0.40 % |
| SAD | AC | case10192_epigrids__sad | gap +3.04 % |
| SAD | AC | case30000_goc__sad | gap +0.56 % |

Across all groups the AC objective agrees with PowerModels.jl to a median of
0.001 % (maximum 3.04 %) on 191 solved cases;
the PGLib references are rounded to five significant digits, so agreement
below 0.01 % is at the reference's own resolution.
Summed over all models the run used 41.0 core-hours. The longest AC solve
was case2742_goc at 16132 s, which is the sum of its three
attempts from the different starts rather than one solve; the longest model
build was case24464_goc__sad at 1334 s. Pyomo construction scales
roughly linearly with the bus count, at 6 to 40 ms per bus.

### Typical Operating Conditions (TYP)

| Case | Nodes | DC ($/h) | DC ref | DC gap (%) | DC | AC ($/h) | AC ref | AC gap (%) | AC | AC build (s) | AC solve (s) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| case3_lmbd | 3 | 5 695.90 | 5 695.90 | 0.00 | ✓ | 5 812.64 | 5 812.60 | 0.00 | ✓ | <1 | <1 |
| case5_pjm | 5 | 17 479.90 | 17 480.00 | 0.00 | ✓ | 17 551.89 | 17 552.00 | 0.00 | ✓ | <1 | <1 |
| case14_ieee | 14 | 2 051.53 | 2 051.50 | 0.00 | ✓ | 2 178.08 | 2 178.10 | 0.00 | ✓ | <1 | <1 |
| case24_ieee_rts | 24 | 61 001.24 | 61 001.00 | 0.00 | ✓ | 63 352.20 | 63 352.00 | 0.00 | ✓ | <1 | <1 |
| case30_as | 30 | 767.60 | 767.60 | 0.00 | ✓ | 803.13 | 803.13 | 0.00 | ✓ | <1 | <1 |
| case30_ieee | 30 | 7 472.81 | 7 472.80 | 0.00 | ✓ | 8 208.52 | 8 208.50 | 0.00 | ✓ | <1 | <1 |
| case39_epri | 39 | 136 889.69 | 136 890.00 | 0.00 | ✓ | 138 415.56 | 138 420.00 | 0.00 | ✓ | <1 | <1 |
| case57_ieee | 57 | 34 772.95 | 34 773.00 | 0.00 | ✓ | 37 589.34 | 37 589.00 | 0.00 | ✓ | <1 | <1 |
| case60_c | 60 | 90 700.00 | 90 700.00 | 0.00 | ✓ | 92 693.67 | 92 694.00 | 0.00 | ✓ | <1 | <1 |
| case73_ieee_rts | 73 | 183 003.72 | 183 000.00 | 0.00 | ✓ | 189 764.08 | 189 760.00 | 0.00 | ✓ | <1 | <1 |
| case89_pegase | 89 | 105 044.27 | 105 040.00 | 0.00 | ✓ | 107 285.67 | 107 290.00 | 0.00 | ✓ | <1 | <1 |
| case118_ieee | 118 | 93 100.73 | 93 101.00 | 0.00 | ✓ | 97 213.61 | 97 214.00 | 0.00 | ✓ | <1 | <1 |
| case162_ieee_dtc | 162 | 101 462.26 | 101 460.00 | 0.00 | ✓ | 108 075.68 | 108 080.00 | 0.00 | ✓ | <1 | <1 |
| case179_goc | 179 | 751 881.01 | 751 880.00 | 0.00 | ✓ | 754 266.41 | 754 270.00 | 0.00 | ✓ | <1 | <1 |
| case197_snem | 197 | 1.47 | 1.47 | −0.28 | ✓ | 1.50 | 1.50 | −0.11 | ✓ | <1 | <1 |
| case200_activ | 200 | 27 479.64 | 27 480.00 | 0.00 | ✓ | 27 557.57 | 27 558.00 | 0.00 | ✓ | <1 | <1 |
| case240_pserc | 240 | 3 271 437.38 | 3 271 400.00 | 0.00 | ✓ | 3 329 670.06 | 3 329 700.00 | 0.00 | ✓ | <1 | <1 |
| case300_ieee | 300 | 517 851.08 | 517 850.00 | 0.00 | ✓ | 565 220.11 | 565 220.00 | 0.00 | ✓ | <1 | 1 |
| case500_goc | 500 | 440 549.16 | 440 550.00 | 0.00 | ✓ | 454 936.40 | 454 950.00 | 0.00 | ✓ | 1 | <1 |
| case588_sdet | 588 | 310 125.52 | 310 130.00 | 0.00 | ✓ | 313 139.78 | 313 140.00 | 0.00 | ✓ | 1 | <1 |
| case793_goc | 793 | 258 307.88 | 258 310.00 | 0.00 | ✓ | 260 198.57 | 260 200.00 | 0.00 | ✓ | 2 | 1 |
| case1354_pegase | 1354 | 1 218 182.02 | 1 218 200.00 | 0.00 | ✓ | 1 258 843.99 | 1 258 800.00 | 0.00 | ✓ | 4 | 3 |
| case1803_snem | 1803 | 87 703.26 | 87 696.00 | +0.01 | ✓ | 98 312.40 | 98 335.00 | −0.02 | ✓ | 6 | 4 |
| case1888_rte | 1888 | 1 352 871.74 | 1 352 900.00 | 0.00 | ✓ | 1 402 630.44 | 1 402 500.00 | +0.01 | ✓ | 6 | 5 |
| case1951_rte | 1951 | 2 031 627.90 | 2 031 600.00 | 0.00 | ✓ | 2 085 582.80 | 2 085 600.00 | 0.00 | ✓ | 7 | 4 |
| case2000_goc | 2000 | 943 042.20 | 943 040.00 | 0.00 | ✓ | 973 432.47 | 973 430.00 | 0.00 | ✓ | 8 | 5 |
| case2312_goc | 2312 | 440 328.37 | 440 330.00 | 0.00 | ✓ | 441 329.99 | 441 330.00 | 0.00 | ✓ | 14 | 5 |
| case2383wp_k | 2383 | 1 804 091.75 | 1 804 100.00 | 0.00 | ✓ | 1 868 531.53 | 1 868 200.00 | +0.02 | ✓ | 9 | 5 |
| case2736sp_k | 2736 | 1 276 033.65 | 1 276 000.00 | 0.00 | ✓ | 1 307 853.17 | 1 308 000.00 | −0.01 | ✓ | 11 | 3 |
| case2737sop_k | 2737 | 764 008.57 | 764 010.00 | 0.00 | ✓ | 777 737.62 | 777 730.00 | 0.00 | ✓ | 12 | 3 |
| case2742_goc | 2742 | 259 696.43 | 259 700.00 | 0.00 | ✓ | — | 275 710.00 | — | not solved ✗ | 13 | 16132 |
| case2746wop_k | 2746 | 1 178 163.94 | 1 178 200.00 | 0.00 | ✓ | 1 208 284.08 | 1 208 300.00 | 0.00 | ✓ | 12 | 4 |
| case2746wp_k | 2746 | 1 581 425.01 | 1 581 400.00 | 0.00 | ✓ | 1 631 800.12 | 1 631 700.00 | +0.01 | ✓ | 12 | 4 |
| case2848_rte | 2848 | 1 267 731.65 | 1 267 700.00 | 0.00 | ✓ | 1 286 608.00 | 1 286 600.00 | 0.00 | ✓ | 12 | 5 |
| case2853_sdet | 2853 | 2 036 957.85 | 2 037 000.00 | 0.00 | ✓ | 2 052 386.72 | 2 052 400.00 | 0.00 | ✓ | 16 | 24 |
| case2868_rte | 2868 | 1 966 683.71 | 1 966 700.00 | 0.00 | ✓ | 2 009 599.86 | 2 009 600.00 | 0.00 | ✓ | 12 | 6 |
| case2869_pegase | 2869 | 2 386 379.35 | 2 386 400.00 | 0.00 | ✓ | 2 462 790.43 | 2 462 800.00 | 0.00 | ✓ | 19 | 9 |
| case3012wp_k | 3012 | 2 509 001.09 | 2 509 000.00 | 0.00 | ✓ | 2 600 844.49 | 2 600 800.00 | 0.00 | ✓ | 18 | 5 |
| case3022_goc | 3022 | 599 221.11 | 599 220.00 | 0.00 | ✓ | 601 383.81 | 601 380.00 | 0.00 | ✓ | 14 | 8 |
| case3120sp_k | 3120 | 2 087 975.32 | 2 088 000.00 | 0.00 | ✓ | 2 147 971.57 | 2 148 000.00 | 0.00 | ✓ | 14 | 6 |
| case3375wp_k | 3374 | 7 317 010.90 | 7 317 000.00 | 0.00 | ✓ | 7 435 703.75 | 7 438 200.00 | −0.03 | ✓ | 19 | 7 |
| case3970_goc | 3970 | 934 219.29 | 934 220.00 | 0.00 | ✓ | 960 985.26 | 960 990.00 | 0.00 | ✓ | 27 | 16 |
| case4020_goc | 4020 | 795 061.59 | 795 060.00 | 0.00 | ✓ | 822 247.28 | 822 250.00 | 0.00 | ✓ | 28 | 26 |
| case4601_goc | 4601 | 793 813.76 | 793 810.00 | 0.00 | ✓ | 826 241.53 | 826 240.00 | 0.00 | ✓ | 38 | 13 |
| case4619_goc | 4619 | 457 436.33 | 457 440.00 | 0.00 | ✓ | 476 703.32 | 476 700.00 | 0.00 | ✓ | 39 | 15 |
| case4661_sdet | 4661 | 2 216 303.42 | 2 216 300.00 | 0.00 | ✓ | 2 251 344.07 | 2 251 300.00 | 0.00 | ✓ | 35 | 353 |
| case4837_goc | 4837 | 850 396.85 | 850 400.00 | 0.00 | ✓ | 872 261.14 | 872 260.00 | 0.00 | ✓ | 38 | 21 |
| case4917_goc | 4917 | 1 383 655.47 | 1 383 700.00 | 0.00 | ✓ | 1 387 798.35 | 1 387 800.00 | 0.00 | ✓ | 38 | 25 |
| case5658_epigrids | 5658 | 1 195 466.11 | 1 195 500.00 | 0.00 | ✓ | 1 207 362.67 | 1 207 300.00 | +0.01 | ✓ | 54 | 14 |
| case6468_rte | 6468 | 1 982 818.37 | 1 982 800.00 | 0.00 | ✓ | 2 069 795.37 | 2 069 700.00 | 0.00 | ✓ | 60 | 16 |
| case6470_rte | 6470 | 2 136 095.35 | 2 136 100.00 | 0.00 | ✓ | 2 237 615.31 | 2 237 600.00 | 0.00 | ✓ | 65 | 20 |
| case6495_rte | 6495 | 2 561 787.20 | 2 561 800.00 | 0.00 | ✓ | 3 068 512.98 | 3 067 800.00 | +0.02 | ✓ | 62 | 43 |
| case6515_rte | 6515 | 2 559 328.97 | 2 559 300.00 | 0.00 | ✓ | 2 825 714.97 | 2 825 500.00 | +0.01 | ✓ | 69 | 113 |
| case7336_epigrids | 7336 | 1 855 899.64 | 1 855 900.00 | 0.00 | ✓ | 1 882 408.93 | 1 882 400.00 | 0.00 | ✓ | 86 | 23 |
| case8387_pegase | 8387 | 2 502 808.00 | 2 502 800.00 | 0.00 | ✓ | 2 771 427.26 | 2 771 400.00 | 0.00 | ✓ | 140 | 3060 |
| case9241_pegase | 9241 | 6 028 744.84 | 6 028 700.00 | 0.00 | ✓ | 6 243 094.60 | 6 243 100.00 | 0.00 | ✓ | 200 | 1805 |
| case9591_goc | 9591 | 1 030 939.10 | 1 030 900.00 | 0.00 | ✓ | 1 061 901.71 | 1 061 700.00 | +0.02 | ✓ | 153 | 61 |
| case10000_goc | 10000 | 1 346 113.34 | 1 346 100.00 | 0.00 | ✓ | 1 354 031.31 | 1 354 000.00 | 0.00 | ✓ | 143 | 46 |
| case10192_epigrids | 10192 | 1 666 679.87 | 1 665 600.00 | +0.06 | ✓ | 1 717 033.11 | 1 686 900.00 | +1.79 | ✓ | 199 | 90 |
| case10480_goc | 10480 | 2 215 825.63 | 2 215 800.00 | 0.00 | ✓ | 2 316 320.66 | 2 314 600.00 | +0.07 | ✓ | 216 | 386 |
| case13659_pegase | 13659 | 8 769 893.12 | 8 769 900.00 | 0.00 | ✓ | 8 948 048.92 | 8 948 000.00 | 0.00 | ✓ | 346 | 4258 |
| case19402_goc | 19402 | 1 897 793.81 | 1 897 800.00 | 0.00 | ✓ | 1 977 815.40 | 1 977 800.00 | 0.00 | ✓ | 721 | 2446 |
| case20758_epigrids | 20758 | 2 572 324.67 | 2 572 300.00 | 0.00 | ✓ | 2 618 592.44 | 2 618 600.00 | 0.00 | ✓ | 835 | 87 |
| case24464_goc | 24464 | 2 512 810.93 | 2 512 800.00 | 0.00 | ✓ | 2 629 660.98 | 2 629 500.00 | +0.01 | ✓ | 1313 | 10016 |
| case30000_goc | 30000 | 1 092 247.12 | 1 092 100.00 | +0.01 | ✓ | 1 142 331.57 | 1 142 300.00 | 0.00 | ✓ | 1147 | 616 |

### Congested Operating Conditions (API)

| Case | Nodes | DC ($/h) | DC ref | DC gap (%) | DC | AC ($/h) | AC ref | AC gap (%) | AC | AC build (s) | AC solve (s) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| case3_lmbd__api | 3 | 10 444.36 | 10 444.00 | 0.00 | ✓ | 11 242.12 | 11 242.00 | 0.00 | ✓ | <1 | <1 |
| case5_pjm__api | 5 | 78 025.19 | 78 025.00 | 0.00 | ✓ | 78 949.91 | 78 950.00 | 0.00 | ✓ | <1 | <1 |
| case14_ieee__api | 14 | 4 797.60 | 4 797.60 | 0.00 | ✓ | 5 999.36 | 5 999.40 | 0.00 | ✓ | <1 | <1 |
| case24_ieee_rts__api | 24 | 148 845.53 | 148 850.00 | 0.00 | ✓ | 161 222.58 | 161 220.00 | 0.00 | ✓ | <1 | <1 |
| case30_as__api | 30 | 3 092.12 | 3 092.10 | 0.00 | ✓ | 4 996.20 | 4 996.20 | 0.00 | ✓ | <1 | <1 |
| case30_ieee__api | 30 | 16 145.05 | 16 145.00 | 0.00 | ✓ | 18 036.59 | 18 037.00 | 0.00 | ✓ | <1 | <1 |
| case39_epri__api | 39 | 252 754.58 | 252 750.00 | 0.00 | ✓ | 256 769.34 | 256 770.00 | 0.00 | ✓ | <1 | <1 |
| case57_ieee__api | 57 | 34 081.06 | 34 081.00 | 0.00 | ✓ | 36 242.46 | 36 242.00 | 0.00 | ✓ | <1 | <1 |
| case60_c__api | 60 | 176 378.43 | 176 380.00 | 0.00 | ✓ | 185 002.89 | 185 000.00 | 0.00 | ✓ | <1 | <1 |
| case73_ieee_rts__api | 73 | 472 183.13 | 472 180.00 | 0.00 | ✓ | 509 847.89 | 509 850.00 | 0.00 | ✓ | <1 | <1 |
| case89_pegase__api | 89 | 118 630.14 | 118 630.00 | 0.00 | ✓ | 129 568.36 | 129 570.00 | 0.00 | ✓ | <1 | <1 |
| case118_ieee__api | 118 | 231 291.90 | 231 290.00 | 0.00 | ✓ | 249 614.51 | 249 610.00 | 0.00 | ✓ | <1 | <1 |
| case162_ieee_dtc__api | 162 | 111 570.83 | 111 570.00 | 0.00 | ✓ | 120 881.19 | 120 880.00 | 0.00 | ✓ | <1 | <1 |
| case179_goc__api | 179 | 1 813 207.61 | 1 813 200.00 | 0.00 | ✓ | 1 883 405.36 | 1 883 400.00 | 0.00 | ✓ | <1 | 2 |
| case197_snem__api | 197 | 15 689.71 | 15 690.00 | 0.00 | ✓ | 16 366.45 | 16 363.00 | +0.02 | ✓ | <1 | <1 |
| case200_activ__api | 200 | 40 129.76 | 40 130.00 | 0.00 | ✓ | 40 700.07 | 40 700.00 | 0.00 | ✓ | <1 | <1 |
| case240_pserc__api | 240 | 4 624 849.61 | 4 624 800.00 | 0.00 | ✓ | 4 692 231.48 | 4 692 200.00 | 0.00 | ✓ | <1 | 4 |
| case300_ieee__api | 300 | 659 835.46 | 659 840.00 | 0.00 | ✓ | 686 040.79 | 686 040.00 | 0.00 | ✓ | <1 | <1 |
| case500_goc__api | 500 | 646 879.66 | 646 870.00 | 0.00 | ✓ | 688 327.13 | 688 290.00 | +0.01 | ✓ | 2 | 1 |
| case588_sdet__api | 588 | 392 953.70 | 392 950.00 | 0.00 | ✓ | 398 761.70 | 398 760.00 | 0.00 | ✓ | 1 | 2 |
| case793_goc__api | 793 | 372 180.72 | 372 180.00 | 0.00 | ✓ | 379 813.58 | 379 800.00 | 0.00 | ✓ | 2 | 1 |
| case1354_pegase__api | 1354 | 1 558 525.15 | 1 558 500.00 | 0.00 | ✓ | 1 608 226.87 | 1 608 200.00 | 0.00 | ✓ | 5 | 4 |
| case1803_snem__api | 1803 | 62 065.88 | 61 723.00 | +0.56 | ✓ | 80 251.35 | 80 240.00 | +0.01 | ✓ | 6 | 10 |
| case1888_rte__api | 1888 | 1 961 568.54 | 1 961 600.00 | 0.00 | ✓ | 2 019 737.49 | 2 019 700.00 | 0.00 | ✓ | 7 | 40 |
| case1951_rte__api | 1951 | 2 411 503.63 | 2 411 500.00 | 0.00 | ✓ | 2 490 372.71 | 2 490 300.00 | 0.00 | ✓ | 7 | 34 |
| case2000_goc__api | 2000 | 1 409 993.10 | 1 410 000.00 | 0.00 | ✓ | 1 483 886.47 | 1 483 900.00 | 0.00 | ✓ | 11 | 8 |
| case2312_goc__api | 2312 | 615 231.89 | 615 220.00 | 0.00 | ✓ | 662 992.21 | 663 440.00 | −0.07 | ✓ | 10 | 9 |
| case2383wp_k__api | 2383 | 279 125.82 | 279 130.00 | 0.00 | ✓ | 279 125.82 | 279 130.00 | 0.00 | ✓ | 10 | 3 |
| case2736sp_k__api | 2736 | 978 236.08 | 977 820.00 | +0.04 | ✓ | 1 020 221.70 | 1 017 800.00 | +0.24 | ✓ | 12 | 7 |
| case2737sop_k__api | 2737 | 755 516.89 | 755 320.00 | +0.03 | ✓ | 788 745.13 | 788 310.00 | +0.06 | ✓ | 11 | 7 |
| case2742_goc__api | 2742 | 505 399.17 | 505 400.00 | 0.00 | ✓ | — | 609 960.00 | — | not solved ✗ | 14 | 135 |
| case2746wop_k__api | 2746 | 531 413.46 | 531 550.00 | −0.03 | ✓ | 549 633.55 | 550 480.00 | −0.15 | ✓ | 14 | 6 |
| case2746wp_k__api | 2746 | 581 827.85 | 581 830.00 | 0.00 | ✓ | 581 829.91 | 581 830.00 | 0.00 | ✓ | 12 | 10 |
| case2848_rte__api | 2848 | 1 488 984.49 | 1 489 000.00 | 0.00 | ✓ | 1 531 112.80 | 1 531 100.00 | 0.00 | ✓ | 13 | 7 |
| case2853_sdet__api | 2853 | 2 456 065.06 | 2 456 100.00 | 0.00 | ✓ | 2 484 262.86 | 2 484 300.00 | 0.00 | ✓ | 14 | 9 |
| case2868_rte__api | 2868 | 2 277 459.40 | 2 277 500.00 | 0.00 | ✓ | 2 343 858.79 | 2 343 900.00 | 0.00 | ✓ | 12 | 104 |
| case2869_pegase__api | 2869 | 2 966 625.83 | 2 966 600.00 | 0.00 | ✓ | 3 062 988.91 | 3 063 000.00 | 0.00 | ✓ | 18 | 11 |
| case3012wp_k__api | 3012 | 851 751.87 | 851 750.00 | 0.00 | ✓ | 914 591.92 | 914 590.00 | 0.00 | ✓ | 15 | 9 |
| case3022_goc__api | 3022 | 666 194.32 | 666 190.00 | 0.00 | ✓ | 687 774.70 | 687 360.00 | +0.06 | ✓ | 22 | 10 |
| case3120sp_k__api | 3120 | 1 326 456.07 | 1 326 500.00 | 0.00 | ✓ | 1 403 484.63 | 1 403 500.00 | 0.00 | ✓ | 20 | 10 |
| case3375wp_k__api | 3374 | 6 272 119.47 | 6 272 100.00 | 0.00 | ✓ | 6 387 877.89 | 6 364 100.00 | +0.37 | ✓ | 17 | 8 |
| case3970_goc__api | 3970 | 1 227 776.26 | 1 227 800.00 | 0.00 | ✓ | 1 749 410.63 | 1 749 400.00 | 0.00 | ✓ | 31 | 14 |
| case4020_goc__api | 4020 | 1 082 421.41 | 1 082 400.00 | 0.00 | ✓ | 1 281 691.37 | 1 281 700.00 | 0.00 | ✓ | 29 | 53 |
| case4601_goc__api | 4601 | 796 409.64 | 796 410.00 | 0.00 | ✓ | 871 237.69 | 871 240.00 | 0.00 | ✓ | 36 | 15 |
| case4619_goc__api | 4619 | 1 010 784.39 | 1 010 800.00 | 0.00 | ✓ | 1 068 791.38 | 1 068 800.00 | 0.00 | ✓ | 42 | 23 |
| case4661_sdet__api | 4661 | 2 670 622.41 | 2 670 600.00 | 0.00 | ✓ | 2 731 463.16 | 2 731 500.00 | 0.00 | ✓ | 45 | 53 |
| case4837_goc__api | 4837 | 1 209 605.67 | 1 209 600.00 | 0.00 | ✓ | 1 295 011.17 | 1 291 900.00 | +0.24 | ✓ | 40 | 142 |
| case4917_goc__api | 4917 | 1 703 463.42 | 1 703 500.00 | 0.00 | ✓ | 1 717 372.57 | 1 717 400.00 | 0.00 | ✓ | 42 | 14 |
| case5658_epigrids__api | 5658 | 1 317 875.61 | 1 306 300.00 | +0.89 | ✓ | 1 338 548.66 | 1 326 800.00 | +0.89 | ✓ | 56 | 20 |
| case6468_rte__api | 6468 | 2 343 288.88 | 2 343 300.00 | 0.00 | ✓ | 2 452 903.11 | 2 452 700.00 | +0.01 | ✓ | 59 | 20 |
| case6470_rte__api | 6470 | 2 604 136.28 | 2 604 100.00 | 0.00 | ✓ | 2 724 691.73 | 2 720 900.00 | +0.14 | ✓ | 62 | 44 |
| case6495_rte__api | 6495 | 2 924 493.91 | 2 924 500.00 | 0.00 | ✓ | 3 129 972.26 | 3 130 400.00 | −0.01 | ✓ | 70 | 123 |
| case6515_rte__api | 6515 | 2 959 743.12 | 2 959 700.00 | 0.00 | ✓ | 3 127 215.45 | 3 129 200.00 | −0.06 | ✓ | 68 | 210 |
| case7336_epigrids__api | 7336 | 1 980 515.20 | 1 980 100.00 | +0.02 | ✓ | 2 038 032.18 | 2 035 500.00 | +0.12 | ✓ | 97 | 35 |
| case8387_pegase__api | 8387 | 4 975 882.04 | 4 975 900.00 | 0.00 | ✓ | 5 242 827.16 | 5 242 800.00 | 0.00 | ✓ | 128 | 3775 |
| case9241_pegase__api | 9241 | 6 816 325.13 | 6 816 300.00 | 0.00 | ✓ | 7 068 768.32 | 7 067 900.00 | +0.01 | ✓ | 210 | 82 |
| case9591_goc__api | 9591 | 1 470 244.15 | 1 470 200.00 | 0.00 | ✓ | 1 531 521.29 | 1 570 300.00 | −2.47 | ✓ | 172 | 252 |
| case10000_goc__api | 10000 | 2 499 130.85 | 2 499 100.00 | 0.00 | ✓ | 2 678 659.41 | 2 678 700.00 | 0.00 | ✓ | 146 | 857 |
| case10192_epigrids__api | 10192 | 1 875 183.93 | 1 856 500.00 | +1.01 | ✓ | 1 984 160.49 | 1 977 700.00 | +0.33 | ✓ | 175 | 113 |
| case10480_goc__api | 10480 | 2 712 411.58 | 2 712 400.00 | 0.00 | ✓ | 2 866 506.41 | 2 863 500.00 | +0.10 | ✓ | 195 | 88 |
| case13659_pegase__api | 13659 | 9 126 339.96 | 9 126 300.00 | 0.00 | ✓ | 9 385 903.56 | 9 385 800.00 | 0.00 | ✓ | 399 | 4736 |
| case19402_goc__api | 19402 | 2 458 062.23 | 2 458 100.00 | 0.00 | ✓ | 2 583 662.75 | 2 583 700.00 | 0.00 | ✓ | 697 | 8965 |
| case20758_epigrids__api | 20758 | 3 035 073.94 | 3 034 800.00 | +0.01 | ✓ | 3 126 647.95 | 3 126 500.00 | 0.00 | ✓ | 773 | 1874 |
| case24464_goc__api | 24464 | 2 531 140.79 | 2 531 100.00 | 0.00 | ✓ | — | 2 684 000.00 | — | not solved ✗ | 1265 | 12176 |
| case30000_goc__api | 30000 | 1 701 584.82 | 1 700 900.00 | +0.04 | ✓ | 1 778 080.54 | 1 777 900.00 | +0.01 | ✓ | 1136 | 9031 |

### Small Angle Difference Conditions (SAD)

| Case | Nodes | DC ($/h) | DC ref | DC gap (%) | DC | AC ($/h) | AC ref | AC gap (%) | AC | AC build (s) | AC solve (s) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| case3_lmbd__sad | 3 | 5 855.99 | 5 856.00 | 0.00 | ✓ | 5 959.31 | 5 959.30 | 0.00 | ✓ | <1 | <1 |
| case5_pjm__sad | 5 | — | inf. | — | inf. ✓ | 26 108.84 | 26 109.00 | 0.00 | ✓ | <1 | <1 |
| case14_ieee__sad | 14 | — | inf. | — | inf. ✓ | 2 776.79 | 2 776.80 | 0.00 | ✓ | <1 | <1 |
| case24_ieee_rts__sad | 24 | 78 122.48 | 78 122.00 | 0.00 | ✓ | 76 917.96 | 76 918.00 | 0.00 | ✓ | <1 | <1 |
| case30_as__sad | 30 | — | inf. | — | inf. ✓ | 897.35 | 897.35 | 0.00 | ✓ | <1 | <1 |
| case30_ieee__sad | 30 | — | inf. | — | inf. ✓ | 8 208.52 | 8 208.50 | 0.00 | ✓ | <1 | <1 |
| case39_epri__sad | 39 | 150 669.88 | 150 670.00 | 0.00 | ✓ | 148 340.50 | 148 340.00 | 0.00 | ✓ | <1 | <1 |
| case57_ieee__sad | 57 | — | inf. | — | inf. ✓ | 38 663.28 | 38 663.00 | 0.00 | ✓ | <1 | <1 |
| case60_c__sad | 60 | — | inf. | — | inf. ✓ | 113 498.76 | 113 500.00 | 0.00 | ✓ | <1 | <1 |
| case73_ieee_rts__sad | 73 | 232 679.11 | 232 680.00 | 0.00 | ✓ | 227 603.74 | 227 600.00 | 0.00 | ✓ | <1 | <1 |
| case89_pegase__sad | 89 | — | inf. | — | inf. ✓ | 107 285.67 | 107 290.00 | 0.00 | ✓ | <1 | <1 |
| case118_ieee__sad | 118 | — | inf. | — | inf. ✓ | 105 155.04 | 105 160.00 | 0.00 | ✓ | <1 | <1 |
| case162_ieee_dtc__sad | 162 | 106 287.85 | 106 290.00 | 0.00 | ✓ | 108 690.73 | 108 690.00 | 0.00 | ✓ | <1 | <1 |
| case179_goc__sad | 179 | — | inf. | — | inf. ✓ | 762 532.53 | 762 530.00 | 0.00 | ✓ | <1 | <1 |
| case197_snem__sad | 197 | — | inf. | — | inf. ✓ | 1.51 | 1.51 | −0.02 | ✓ | <1 | <1 |
| case200_activ__sad | 200 | — | inf. | — | inf. ✓ | 27 557.57 | 27 558.00 | 0.00 | ✓ | <1 | <1 |
| case240_pserc__sad | 240 | — | inf. | — | inf. ✓ | 3 405 363.60 | 3 405 400.00 | 0.00 | ✓ | <1 | 1 |
| case300_ieee__sad | 300 | 527 294.20 | 527 290.00 | 0.00 | ✓ | 565 704.43 | 565 700.00 | 0.00 | ✓ | <1 | 2 |
| case500_goc__sad | 500 | — | inf. | — | inf. ✓ | 487 421.83 | 487 400.00 | 0.00 | ✓ | 1 | <1 |
| case588_sdet__sad | 588 | — | inf. | — | inf. ✓ | 329 356.03 | 329 360.00 | 0.00 | ✓ | 1 | 1 |
| case793_goc__sad | 793 | — | inf. | — | inf. ✓ | 285 799.98 | 285 800.00 | 0.00 | ✓ | 2 | 1 |
| case1354_pegase__sad | 1354 | — | inf. | — | inf. ✓ | 1 258 848.04 | 1 258 800.00 | 0.00 | ✓ | 7 | 3 |
| case1803_snem__sad | 1803 | — | inf. | — | inf. ✓ | 106 282.37 | 106 340.00 | −0.05 | ✓ | 6 | 4 |
| case1888_rte__sad | 1888 | 1 353 178.46 | 1 353 200.00 | 0.00 | ✓ | 1 414 113.32 | 1 413 900.00 | +0.02 | ✓ | 9 | 5 |
| case1951_rte__sad | 1951 | — | inf. | — | inf. ✓ | 2 092 521.53 | 2 092 400.00 | +0.01 | ✓ | 7 | 4 |
| case2000_goc__sad | 2000 | — | inf. | — | inf. ✓ | 992 879.79 | 992 880.00 | 0.00 | ✓ | 8 | 8 |
| case2312_goc__sad | 2312 | — | inf. | — | inf. ✓ | 462 355.31 | 462 350.00 | 0.00 | ✓ | 9 | 5 |
| case2383wp_k__sad | 2383 | — | inf. | — | inf. ✓ | 1 917 421.76 | 1 911 200.00 | +0.33 | ✓ | 10 | 6 |
| case2736sp_k__sad | 2736 | — | inf. | — | inf. ✓ | 1 332 057.63 | 1 326 600.00 | +0.41 | ✓ | 13 | 5 |
| case2737sop_k__sad | 2737 | — | inf. | — | inf. ✓ | 794 478.12 | 790 950.00 | +0.45 | ✓ | 12 | 4 |
| case2742_goc__sad | 2742 | 259 696.43 | 259 700.00 | 0.00 | ✓ | — | 275 710.00 | — | not solved ✗ | 14 | 2147 |
| case2746wop_k__sad | 2746 | — | inf. | — | inf. ✓ | 1 233 541.35 | 1 233 700.00 | −0.01 | ✓ | 12 | 6 |
| case2746wp_k__sad | 2746 | — | inf. | — | inf. ✓ | 1 667 753.55 | 1 666 900.00 | +0.05 | ✓ | 13 | 5 |
| case2848_rte__sad | 2848 | — | inf. | — | inf. ✓ | 1 288 966.02 | 1 289 000.00 | 0.00 | ✓ | 13 | 9 |
| case2853_sdet__sad | 2853 | — | inf. | — | inf. ✓ | 2 069 157.72 | 2 069 200.00 | 0.00 | ✓ | 16 | 15 |
| case2868_rte__sad | 2868 | — | inf. | — | inf. ✓ | 2 021 266.72 | 2 021 300.00 | 0.00 | ✓ | 16 | 6 |
| case2869_pegase__sad | 2869 | — | inf. | — | inf. ✓ | 2 468 676.30 | 2 468 700.00 | 0.00 | ✓ | 20 | 10 |
| case3012wp_k__sad | 3012 | — | inf. | — | inf. ✓ | 2 619 458.51 | 2 619 500.00 | 0.00 | ✓ | 14 | 6 |
| case3022_goc__sad | 3022 | 599 221.11 | 599 220.00 | 0.00 | ✓ | 601 434.00 | 601 430.00 | 0.00 | ✓ | 16 | 11 |
| case3120sp_k__sad | 3120 | — | inf. | — | inf. ✓ | 2 174 869.90 | 2 174 900.00 | 0.00 | ✓ | 15 | 8 |
| case3375wp_k__sad | 3374 | 7 319 636.49 | 7 319 600.00 | 0.00 | ✓ | 7 435 703.75 | 7 438 200.00 | −0.03 | ✓ | 18 | 7 |
| case3970_goc__sad | 3970 | — | inf. | — | inf. ✓ | 965 550.65 | 965 550.00 | 0.00 | ✓ | 29 | 16 |
| case4020_goc__sad | 4020 | — | inf. | — | inf. ✓ | 889 686.87 | 889 690.00 | 0.00 | ✓ | 31 | 32 |
| case4601_goc__sad | 4601 | 1 195 549.95 | 1 195 500.00 | 0.00 | ✓ | 878 179.90 | 878 180.00 | 0.00 | ✓ | 36 | 14 |
| case4619_goc__sad | 4619 | — | inf. | — | inf. ✓ | 484 325.65 | 484 350.00 | −0.01 | ✓ | 39 | 18 |
| case4661_sdet__sad | 4661 | — | inf. | — | inf. ✓ | 2 260 951.32 | 2 261 000.00 | 0.00 | ✓ | 36 | 36 |
| case4837_goc__sad | 4837 | — | inf. | — | inf. ✓ | 875 520.25 | 877 120.00 | −0.18 | ✓ | 42 | 23 |
| case4917_goc__sad | 4917 | 1 384 324.78 | 1 384 300.00 | 0.00 | ✓ | 1 389 056.32 | 1 389 000.00 | 0.00 | ✓ | 42 | 22 |
| case5658_epigrids__sad | 5658 | — | inf. | — | inf. ✓ | 1 237 030.37 | 1 235 800.00 | +0.10 | ✓ | 66 | 20 |
| case6468_rte__sad | 6468 | 1 982 818.37 | 1 982 800.00 | 0.00 | ✓ | 2 069 795.37 | 2 069 700.00 | 0.00 | ✓ | 64 | 18 |
| case6470_rte__sad | 6470 | 2 139 077.30 | 2 139 100.00 | 0.00 | ✓ | 2 241 665.97 | 2 241 600.00 | 0.00 | ✓ | 69 | 20 |
| case6495_rte__sad | 6495 | 2 561 803.78 | 2 561 800.00 | 0.00 | ✓ | 3 068 512.98 | 3 067 800.00 | +0.02 | ✓ | 68 | 53 |
| case6515_rte__sad | 6515 | 2 559 522.72 | 2 559 500.00 | 0.00 | ✓ | 2 869 895.95 | 2 869 800.00 | 0.00 | ✓ | 65 | 179 |
| case7336_epigrids__sad | 7336 | — | inf. | — | inf. ✓ | 1 888 377.56 | 1 888 400.00 | 0.00 | ✓ | 100 | 28 |
| case8387_pegase__sad | 8387 | 2 597 699.48 | 2 597 700.00 | 0.00 | ✓ | 2 803 912.40 | 2 803 900.00 | 0.00 | ✓ | 137 | 3742 |
| case9241_pegase__sad | 9241 | — | inf. | — | inf. ✓ | 6 318 490.72 | 6 318 500.00 | 0.00 | ✓ | 201 | 2207 |
| case9591_goc__sad | 9591 | — | inf. | — | inf. ✓ | 1 162 687.25 | 1 167 400.00 | −0.40 | ✓ | 161 | 57 |
| case10000_goc__sad | 10000 | — | inf. | — | inf. ✓ | 1 490 209.58 | 1 490 200.00 | 0.00 | ✓ | 147 | 42 |
| case10192_epigrids__sad | 10192 | — | inf. | — | inf. ✓ | 1 772 570.61 | 1 720 200.00 | +3.04 | ✓ | 194 | 84 |
| case10480_goc__sad | 10480 | — | inf. | — | inf. ✓ | 2 316 373.94 | 2 314 700.00 | +0.07 | ✓ | 212 | 345 |
| case13659_pegase__sad | 13659 | — | inf. | — | inf. ✓ | 9 042 198.38 | 9 042 200.00 | 0.00 | ✓ | 380 | 3781 |
| case19402_goc__sad | 19402 | 1 910 223.69 | 1 910 200.00 | 0.00 | ✓ | 1 983 808.73 | 1 983 800.00 | 0.00 | ✓ | 704 | 2513 |
| case20758_epigrids__sad | 20758 | — | inf. | — | inf. ✓ | 2 638 068.12 | 2 638 200.00 | 0.00 | ✓ | 857 | 105 |
| case24464_goc__sad | 24464 | 2 554 863.90 | 2 554 900.00 | 0.00 | ✓ | 2 654 098.46 | 2 654 000.00 | 0.00 | ✓ | 1334 | 10350 |
| case30000_goc__sad | 30000 | 1 297 526.48 | 1 297 500.00 | 0.00 | ✓ | 1 317 270.99 | 1 310 000.00 | +0.56 | ✓ | 1079 | 2498 |

Gaps are `(potpourri − reference) / reference`; `0.00` means below 0.005 % in
magnitude; `inf.` is a reference problem PowerModels.jl found infeasible.

### What the full run exposed

The first pass over all groups solved most cases but could not build twelve
of the 66 networks (the eight RTE cases, `case2746wop_k`, `case500_goc`,
`case1803_snem` and `case78484_epigrids`),
found a 3.7 % *cheaper* solution than the reference on `case60_c__sad`, and
carried DC gaps of up to 9 % on the RTE cases. None of it was the OPF
formulation; each item below is pinned by a unit test in
`tests/unit_tests/test_pglib_loader.py`, `test_base_powerflow_fallback.py` or
`test_dcopf.py`.

| Symptom | Cause | Fix |
|---|---|---|
| `case60_c__sad` AC-OPF 3.7 % below the reference with six angle limits violated; several SAD DC cases solved although PowerModels reports them infeasible | A MATPOWER branch between different voltage levels without a tap becomes a pandapower `impedance` element; the loader attached `ANGMIN`/`ANGMAX` to lines and transformers only, and both models skipped the impedance rows | Limits attached to `net.impedance`, both models enforce them (27 of the 88 branches in `case60_c` were unconstrained before) |
| The eight RTE cases raised `No reference bus is available` before any model was built | Their MATPOWER slack bus carries no generator; `from_mpc` derives `net.ext_grid` from that generator | The loader adds a zero-capacity external grid at the slack bus: an angle reference without dispatchable power |
| `case2746wop_k` and `case500_goc` raised the same error | Their first slack-bus unit is out of service, so pandapower created the external grid out of service (the remaining units at the bus are kept as `sgen` rows and were fine) | The external grid is reduced to an in-service, zero-capacity angle reference and its cost row dropped |
| DC objectives up to 9 % above the reference on the RTE cases, −2.8 % on `case14_ieee__api`, SAD DC cases solving where PowerModels reports `inf.` | potpourri's DC uses MATPOWER's `−1/x` with the transformer phase shift in the angle difference; the PGLib DC references come from PowerModels' `−x/(r² + x²)` with neither tap nor shift | `DCOPF(net, dc_convention="powermodels")`, used by the benchmark; the default is unchanged |
| 230 000 Pyomo warnings about angle starts outside (−π, π) on the 20 000+ bus cases | The DC-fallback start took the raw DC angles, which reach −22 rad on heavily loaded networks | The fallback wraps the start into (−180°, 180°]; branch flows depend on angle differences only, so the starting residuals are unchanged |
| `case588_sdet` 1.6 % (DC) and 1.5 % (AC) above the reference, `case4661_sdet` 0.4 %, with every network limit slack | Four units in those cases have `PMIN = −200 MW`; pandapower keeps them as `sgen` rows and potpourri's `psG` variable was declared over the non-negative reals, so the lower bound in the data was silently replaced by zero | `psG` is declared over the reals; the bound comes from `net.sgen.min_p_mw`, still defaulting to 0 where the column is missing or `NaN` |
| DC objectives off by a few tenths of a percent on the 120 case files whose transformers have a tap on the low-voltage side, `case24_ieee_rts__sad` infeasible where PowerModels reports 78 122 | pandapower refers a transformer's series impedance to the tapped LV voltage, so the `ppc` reactance carries a factor (vn_trafo_lv / vn_lv_kv)²; PowerModels reads the MATPOWER reactance, which no tap enters | The `powermodels` convention divides that factor out; the MATPOWER convention keeps pandapower's value |
| AC solves of the largest cases ran for up to 2.6 h although the script asked for a one-hour limit | `solve(time_limit=...)` was honoured by Gurobi and mindtpy only and silently dropped for IPOPT | IPOPT receives it as `max_wall_time`; the benchmark's `TIME_LIMIT_S` now applies |
| `case78484_epigrids` could not be built: Pyomo refused a `KCL_const` rule ("does not have a proper value. Found bool 'True'") | Six buses have every branch out of service; their balance holds constants only, numpy floats, so the equality is a numpy bool that the `isinstance(kcl, bool)` guard in the KCL rules did not recognise | The DC and AC rules skip numpy bools as well |
| `case1803_snem` could not be built: `FloatingPointError` from every pandapower power flow | Two purely resistive ties (x = 0); pandapower's DC initialisation and its DC power flow divide by `1/x`, so neither the AC power flow nor the DC fallback could run | `Basemodel` falls back a level further and takes pandapower's network tables without a power flow, with a flat start |

---

## Key API

### `potpourri.benchmarks.load_pglib_case`

```python
from potpourri.benchmarks import load_pglib_case

net = load_pglib_case("case14_ieee")          # TYP
net = load_pglib_case("case14_ieee__api")     # API and SAD variants by suffix
```

Converts a PGLib `.m` file into a `pandapowerNet` with:

- Out-of-service generators removed, with their cost rows, so
  `net.poly_cost.element` still addresses the right units
- All generators flagged `controllable=True`; sgens that `from_mpc` created
  from negative demand stay fixed injections
- Polynomial cost coefficients in `net.poly_cost`
- Phase-angle limits (`angmin_degree` / `angmax_degree`) on `net.line`,
  `net.trafo` and `net.impedance`, read from the MATPOWER `ANGMIN` / `ANGMAX`
  fields
- Transformer taps on the side MATPOWER puts them (`align_tap_sides=True`)
- An in-service reference bus even when the MATPOWER slack bus has no usable
  generator
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
PGLIB_BASELINE_SAD["pglib_opf_case14_ieee__sad"]["dc"]   # inf: PowerModels found the DC-OPF infeasible
```

---

## Known limitations

- The DC column is a comparison of two linear approximations under the same
  convention, not a validation of a physical model. The residual DC
  differences listed in the accounting table above are documented as they
  are; the AC agreement shows that the networks themselves match.
- IPOPT returns a local optimum of the non-convex AC-OPF; PowerModels' values
  are local optima too. A potpourri objective below the reference is not an
  error by itself, but every such case in the table was checked for violated
  limits before being accepted (that check is how the impedance-branch angle
  limits were found).
- The run is single-threaded per solve; times are wall time of `solve()`
  including Pyomo's model translation, measured with 24 solves sharing the
  host.
- The largest case, `case78484_epigrids`, is not in the tables. Its AC model
  takes over an hour to build and its solve did not conclude within the time
  budget of this run in any of the three groups. Its DC-OPF does solve.
- `TIME_LIMIT_S` reaches IPOPT as `max_wall_time`, which IPOPT checks between
  iterations. On the 20 000-bus cases and larger a single iteration can take
  many minutes, so a solve can overrun the limit considerably before it stops.
