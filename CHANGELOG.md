# Changelog

All notable changes to `potpourri` are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added

- **The PGLib benchmark covers all three operating conditions and every case
  size.** `scripts/pglib_benchmark.py` runs the Typical (TYP), Congested (API)
  and Small Angle Difference (SAD) groups, all 66 cases each by default
  (`MAX_BUSES = None`), in parallel worker processes with the largest cases
  dispatched first, an IPOPT time limit per solve and build/solve times per
  model; one results table per group. `parse_baseline_md` now keeps the rows
  whose reference reads `inf.` (PowerModels.jl found the problem infeasible;
  45 SAD cases have no DC-OPF solution), returning `inf`, and the script counts
  an infeasible potpourri result on such a row as agreement.

### Fixed

- **Angle limits reach the branches pandapower stores as impedance
  elements.** A MATPOWER branch between different voltage levels without a
  tap becomes a `net.impedance` row. The PGLib loader attached
  `ANGMIN`/`ANGMAX` to lines and transformers only, and both the AC and the DC
  model skipped the synthetic line indices that carry impedance rows, so those
  branches had no angle limit at all. On `case60_c__sad` 27 of 88 branches
  were unconstrained: the AC-OPF settled 3.7 % below the PGLib reference with
  six limits violated, and DC cases PowerModels reports infeasible solved. The
  loader now writes the limits to `net.impedance` and both models enforce
  them; `case60_c__sad` matches the reference and its DC-OPF is infeasible as
  published.

## [0.5.2] — 2026-09-20

A patch release. The single-period model no longer fails to build when the
flat-start power flow diverges, the PGLib loader converts MATPOWER cases
faithfully (all 18 Typical-Operations cases up to 300 buses now match
PowerModels.jl), and the declared Python floor moves to 3.10, which the
pinned dependencies already required. The only API addition is the
`align_tap_sides` keyword of `load_pglib_case`.

8 new regression tests, suite at 360.

### Fixed

- **Model construction no longer fails when the base power flow diverges.**
  `Basemodel` runs a flat-start pandapower power flow only to obtain the
  `ppc` tables and a starting point, yet a `LoadflowNotConverged` there
  aborted the whole model. On heavily loaded transmission cases such as
  PGLib `case162_ieee_dtc`, `case240_pserc` and `case300_ieee` Newton-Raphson
  diverges with the shipped setpoints (PYPOWER diverges on the same files),
  so no OPF could be attempted. `Basemodel` now falls back to a DC power
  flow, which yields identical bus, branch and generator tables with the DC
  angles and flat voltage magnitudes as the start, and logs a warning. The
  multi-period model is unchanged.
- **PGLib cases are loaded the way MATPOWER describes them.** Four defects in
  `load_pglib_case` distorted the benchmark networks: out-of-service
  generators were dropped without renumbering `net.poly_cost`, so cost curves
  landed on the wrong units and the constant terms of dead units kept being
  charged (`case200_activ` was 14 % above the PGLib reference); the sgens
  that `from_mpc` creates from negative demand were rebalanced to zero
  because they carry no `max_p_mw` (`case240_pserc` lost 4.6 GW of injection,
  +4.9 %) and were marked controllable although they are fixed injections;
  and every transformer tap was encoded on the high-voltage side although
  MATPOWER's `TAP` acts on the from bus, which for 18 transformers in
  `case162_ieee_dtc` and 16 in `case300_ieee` is the low-voltage side, so the
  wrong diagonal of the admittance matrix was divided by `TAP²` and the
  AC-OPF was locally infeasible within `[Vmin, Vmax]`. With the loader
  fixed (new keyword `align_tap_sides`, default on) all 18 TYP cases up to
  300 buses solve and every AC objective is within 0.02 % of PowerModels.jl;
  see `docs/scripts/pglib-validation.md`.

### Changed

- **Python 3.10 is now the minimum.** `requires-python` declared `>=3.9,<3.13`,
  but the pinned pyomo 6.10 and simbench 1.6 both require 3.10, so a 3.9
  install resolved to older releases of both that the suite has never been run
  against. Continuous integration exercises 3.10–3.12 only; the declaration,
  the classifiers, `environment.yaml` and the docs now match. The upper bound
  is unchanged.

## [0.5.1] — 2026-09-20

A patch release: three model fixes found while cross-checking the single-period
models against pandapower power flows on SimBench networks, and an IPOPT pin.
No API changes. Results change only on networks that hit one of the bugs below;
each entry names the affected networks.

10 new regression tests, suite at 352.

### Fixed

- **Single-period models wired branches to the wrong buses on networks with
  switches or orphan buses.** `Basemodel` keyed its line and transformer
  endpoint maps by pandapower bus index while `model.B` and the load and
  generator maps use ppc bus numbers. pandapower moves orphan buses to the end
  of its numbering and routes branches with an open switch to auxiliary buses,
  so on `1-MV-rural--0-sw` every line was attached to the wrong bus pair and the
  AC-OPF was infeasible. Endpoints now come from the ppc branch table, matching
  the multi-period model, and open line switches are honoured. `preprocess_grid`
  also left the merged-away bus of a closed bus-bus switch behind as an orphan
  (the cause of the shifted numbering), did not move switches sitting on it,
  could pick an already removed bus as survivor in chained merges and skipped
  switches with `NaN` `z_ohm`; the survivor now inherits the tightest voltage
  band of the pair. Networks without switches or orphan buses are unaffected.
- **Transformer iron losses were dropped from the AC models.** The AC
  power-flow equations read the charging admittance from ppc column `BR_B` only
  and treated it as purely imaginary, while pandapower stores the shunt
  conductance from `pfe_kw` in column `BR_G`. The models therefore disagreed
  with pandapower's own power flow at the slack by the iron loss: 0.46 kW on a
  160 kVA SimBench LV transformer, 14 kW on a 25 MVA MV one. A new
  `branch_charging_admittance()` in `basemodel` reads both columns (purely
  imaginary for MATPOWER-style arrays without `BR_G`) and feeds the single- and
  multi-period AC models.
- **Storage reactive power entered the AC reactive balance with the generator
  sign.** pandapower treats storage `q_mvar` like a load (positive =
  consumption), and the active-power balance already subtracted `pSTOR`
  accordingly; the reactive balance added `qSTOR`. A storage unit with
  `q_mvar = 0.04` therefore injected 40 kvar instead of drawing it. The
  reactive balance now uses the load convention too. Only networks with
  storage units carrying a non-zero `q_mvar` are affected.

### Changed

- **IPOPT pinned to 3.14.20** (released 2026-08-27) in `environment.yaml` and
  the Dockerfile, up from 3.14.19. The README, `CLAUDE.md` and the
  getting-started guide still quoted 3.14.16 and now name the same version.
  3.14.20 makes the `max_filter_resets` limit effective (it was silently
  ignored before), so a run that relied on unlimited filter resets may take a
  different iteration path; set `max_filter_resets` to a large value to
  restore the old behaviour.

## [0.5.0] — 2026-08-03

A minor rather than a patch release: it adds battery reactive-power capability,
warm starting, seeded device placement and three `add_OPF` options, and changes
five defaults. Each behaviour change is reversible with one argument; they are
collected under *Changed*.

74 new regression tests, suite at 351.

### Fixed

- **pandapower 3.5 broke every single-period model constructor.** 3.5 removed
  `create_continuous_bus_index` from its top-level namespace, and
  `preprocess_grid` — which every single-period constructor runs — called it
  there, so a fresh install produced a package that could not build a model at
  all. It also removed `clear_result_tables`, the first call in **both** result
  mappers, so on 3.5 `solve(to_net=True)` silently left the base-case power flow
  in `net.res_*`. Both are now imported from `pandapower.toolbox`, which carries
  them on every supported version, and the `pandapower` requirement gained the
  upper bound it was missing.
- **`solve()` reported success when result mapping failed.** Both `solve()`
  methods wrapped the termination check *and* the mapping in one
  `except AttributeError`, so a failure while writing `net.res_*` was logged
  into a logger this package disables by default and then discarded. That is
  what hid the breakage above. The guard now covers only the termination check.
- **Multi-period `solve(to_net=True)` never wrote `net.res_*`.** It logged
  `"Solution mapped to net.res_*"` and called no mapper. Because the constructor
  runs `pp.runpp`, callers read **base-case** values as if they were the
  optimisation result. `solve()` now maps, and a new `map_to_net(t)` writes any
  chosen step — `net.res_*` has no time dimension, so `to_net=True` writes the
  last step of the horizon and an `int` selects another.
- **The multi-period battery efficiency produced no round-trip loss.** The
  one-way efficiency was applied to a single *signed* power variable, so the same
  factor scaled both directions and a charge/discharge cycle returned the state
  of charge exactly to its start — an implied round trip of 100 % for any value
  of `efficiency`. Charging and discharging are now separate non-negative legs
  with `η·Pchg − Pdis/η`, matching the single-period storage block, so the round
  trip costs `η²`. Nothing constrained the final state of charge either, so the
  optimiser could drain the battery and export the stored energy for free.
- **Multi-period sgen power parameters were in MW while demand was per-unit.**
  Every other quantity in the model — demand, generator limits, shunt data, line
  and transformer ratings, and the capability data derived from `net.sgen` — is
  per-unit, and the result mapper multiplies by the base on the way out. On any
  network with `sn_mva != 1` the nodal balance therefore added quantities on two
  scales, with no error and no warning. A fully pinned model missed pandapower's
  own power flow by 6.3e-2 p.u. at `sn_mva = 10` and was driven to a locally
  infeasible point; it now reproduces it at any base. Every SimBench network has
  `sn_mva == 1`, which is why this stayed invisible. The same slip is fixed for
  the PV profile, the heat pump's heat-loss profile and the wind potential.
- **`Battery`, `Heatpump` and `PV` were not coupled to the power balance.**
  `KCL_flexibility()` returned `0`, so the devices attached their variables and
  internal dynamics but their power never entered the Kirchhoff constraints —
  feasible trajectories with no effect on voltages, line loading or the
  objective. Devices now register a coupling term and `add_OPF` rebuilds the
  balance, so a device attached after the power-flow equations still reaches
  them. The DC formulation, which had no hook at all, gets the same treatment.
- **PV's active power had the opposite sign to its own Q characteristics**, so
  the Q(P) envelope was mirrored; its bound arguments were also ordered so that
  Pyomo read the maximum as the lower bound. And **every PV Q(U) constraint was
  silently skipped**, because a lookup built with `dict()` over a scalar Pyomo
  Set yields no members and the rule took its skip branch.
- **A multi-period AC OPF could report infeasible on a feasible model.** A cold
  start puts every voltage at 1.0, every angle at 0 and leaves the branch flows
  with no value at all, so Kirchhoff's laws are violated at every bus by the full
  nodal injection, and IPOPT does not always recover on a nonconvex problem. On a
  midday window it reported a locally infeasible point even though curtailing the
  PV to zero is available and feasible. `solve()` now seeds a consistent
  operating point from a power flow at each step first; the seed need not be near
  the optimum, only satisfy the power flow. Also much faster — the unit suite
  went from roughly eight minutes to three.
- **Device placement drew from the unseeded global NumPy generator**, so two runs
  of the same script equipped different buses and gave different results. It
  changed conclusions rather than only bus numbers: the same converter
  configuration in the battery example came out infeasible on one run and optimal
  on the next.
- **The AC model could not be built on a network with shunts.** The conductance
  and susceptance parameters were declared with opposite index sets, and each
  consumer used the wrong one for one of them, so construction raised
  `KeyError`. The LPAC formulation declared the susceptance a third way. All of
  them now agree. Invisible because every SimBench network has no shunts, so the
  loop over them never iterated.
- **Filling `net.res_*` mutated temporaries.** Both mappers used chained
  `inplace=True` calls, which operate on the object the column access returns.
  pandas warns today and makes it a silent no-op under copy-on-write, which
  would leave `NaN` where a value was meant to be — and in one case would have
  replaced a curtailed sgen dispatch with its uncurtailed profile. Replaced with
  assignment, the result frames are built with a float dtype so filling no longer
  downcasts, and `pandas` gained an upper bound.
- **`potpourri.__version__` was hard-coded** and had drifted three releases
  behind `pyproject.toml`, sending bug reports to the wrong changelog entries. It
  now comes from the installed distribution metadata and cannot drift.
- **Documented names and capabilities that do not exist.** EV support was
  advertised in the README, the architecture diagram and the docs index; there is
  no EV module. The README named a heat-pump class whose real name differs in
  case, so copying it raised `ImportError`, omitted the shunt device, and showed
  the wind device inheriting from the wrong base. The multi-period `add_OPF`
  surface also diverged from single-period: `thermal_limit` and `angle_limits`
  raised `TypeError` naming a private method, `net.sgen.min_p_mw` was silently
  ignored, and the slack voltage was pinned with no way to free it — so the same
  snapshot gave materially different voltages depending on the period kind.
- **`performance_test_solver.py` could not run as shipped.** Every case is
  submitted to the NEOS remote service with a placeholder address, which NEOS
  rejects; the rejection surfaced as a bare Pyomo `ActionManagerError` saying
  nothing about the cause. It now explains what to set and skips.

### Added

- **Battery converter reactive power.** `BAT_Q` in the generator sign
  convention, bounded by the apparent-power circle `P² + Q² ≤ S_inv²`, an
  optional power-factor floor, and the grid-code Q(P) and Q(U) capability areas.
  The Q(P) area is keyed on the discharging leg, the mode in which a storage unit
  acts as a generating unit under VDE-AR-N 4105/4110/4120. New arguments
  `s_inv_pu`, `q_control`, `var_q`, `cos_phi_min` and `grid_code`; `s_inv_pu`
  defaults to the active-power limit, which deliberately leaves no reactive
  headroom at full charge or discharge.
- **`warm_start_from_pf()`** on the multi-period models, and `warm_start` on
  `solve()`.
- **`seed` and `rng` on every placing device.** `rng` takes precedence, so one
  generator threaded through a Monte Carlo sweep gives independent scenarios from
  a run that replays exactly — which a single seed cannot express.
- **`thermal_limit`, `angle_limits` and `free_slack_vm`** on the multi-period
  models, matching the single-period surface, plus `angle_limits` on the DC one.
- **`map_to_net(t)`** for reading a chosen time step into `net.res_*`.
- **`terminal_soc` on the battery**, and a single-period-versus-multi-period
  differences table in the user guide.

### Changed

Five defaults move. Each is reversible with one argument:

- `solve(warm_start=True)` — every multi-period solve now runs one power flow per
  time step first. Pass `False` to keep existing initial values.
- `Battery_multi_period(terminal_soc="cyclic")` — the final state of charge
  returns to its initial value. Always feasible, since holding both legs at zero
  satisfies it. `None` restores the previous freedom.
- PV's `pPV` is a generator-convention injection rather than a negated load, so
  its sign flips for anyone reading it.
- Device placement is seeded, so default placements differ from previous runs.
  Nobody can have depended on a specific draw, because there was no way to
  obtain one twice.
- Unsupported `add_OPF` options raise `TypeError` naming the option and the model
  class, instead of being silently discarded.

`environment.yaml` moves to pandapower 3.5.4, which brings `scipy`, `tqdm` and
`geojson` with it, so CI exercises the version a fresh install resolves to.

## [0.4.2] — 2026-07-31

### Fixed

- **A dead-band Q(U) characteristic could silently conflict with the Q(P)
  area.**  Both are imposed on the same reactive power: the area *bounds* Q
  from active power, the characteristic *assigns* it from voltage.  Where
  the assigned value falls outside the bound there is no feasible Q at all
  and the solver reports a plain infeasibility.  For VDE-AR-N 4110 with its
  default dead band the curve assigns +0.484 Pn at 0.90 p.u., which the
  Q(P) area permits only at rated output.  Building such a model now emits
  a `QuCurveOutsidePqAreaWarning` naming the voltage and the active power
  required, so the conflict is diagnosable before solving.  Widening the
  dead band does not help — the curve still assigns its full reactive limit
  at its outermost breakpoints.
- **The multi-period reactive-bound override failed silently when called out
  of order.**  `static_generation_q_ctrl_data()` mutates the dicts built by
  `static_generation_reactive_power_limits()`; running it first left them
  absent and the override became a no-op, restoring exactly the defect
  0.4.1 fixed.  It now raises `RuntimeError` naming the required order.

## [0.4.1] — 2026-07-31

### Added

- **Operating-area figures in the reactive-power-control user guide** —
  one figure per constraint group, showing the region each permits and what it
  leaves unbounded: the Q(P) envelope per `var_q` variant, the Q(U) droop
  across the code's own voltage breakpoints, the inverter S² circle and
  cos(φ) cone separately and then intersected into the PV operating region,
  P(U) curtailment, and the two equality modes whose locus is a line rather
  than an area.  Generated by `scripts/generate_qcontrol_figures.py` directly
  from the grid-code registry, so they cannot drift from the implementation;
  rendered with rwthplots (RWTH style, LaTeX, IEEE column width).

- **Global MINLP solving via `gurobi_direct_minlp`.**  Pyomo's
  `gurobi_direct_minlp` interface (Pyomo >= 6.10, gurobipy >= 12) maps
  sin/cos/exp/log/sqrt onto `gurobipy.nlfunc`, so the polar-form AC power
  flow reaches Gurobi's global spatial branch-and-bound unchanged rather
  than going through IPOPT or MindtPy.  The older `gurobi`,
  `gurobi_direct` and `gurobi_persistent` interfaces are limited to
  expressions of degree 2 and reject the model with `DegreeError`; they
  remain usable for DCOPF and LPACOPF.  `max_iter` and `time_limit` now
  route to Gurobi's own parameter names in both `solve()` implementations
  — `time_limit` was previously dropped for every solver except MindtPy,
  leaving the branch-and-bound unbounded, and `max_iter` was forwarded
  under IPOPT's name and raised `GurobiError`.

- **Q(U) characteristic with a dead band** (`qu_deadband` on both
  `add_OPF` methods).  Where a capability *area* bounds Q and leaves the
  optimiser free inside it, a characteristic *assigns* Q from voltage,
  which is what makes a dead band expressible — a voltage span around
  nominal over which Q is held at zero.  Accepts `True` (the grid code's
  own QV plateau), a `(v_low, v_high)` pair, or a `QVCurve` built with
  `GridCode.deadband_curve()`.  The feasible set pinches to a point in the
  dead band and is therefore not convex, so the constraint is built with
  `pyomo.Piecewise` and needs a MIP-capable solver; `gurobi_direct_minlp`
  solves the documented examples in about a second.  MindtPy's outer
  approximation is unsound on the AC power flow's non-convex nonlinear
  equalities and can report `infeasible` for a feasible model.
- **Capability areas aligned with pandapower.**  The registry now carries
  `Envelope` objects — piecewise-linear Q bounds over active power
  (`pq_area`) and voltage (`qv_area`) — that reproduce
  `pandapower.control.controller.DERController`'s `PQArea4105/4110/4120`
  and `QVArea4105/4110/4120` to machine precision.  pandapower models these
  as shapely polygons (LV/MV) or explicit branch logic (HV); both collapse
  to the same `numpy.interp` envelope, which is the form the Pyomo
  constraints need.  `Envelope.q_flexibility()` mirrors pandapower's method
  of the same name.
- **VDE-AR-N 4110 has real values.**  It was a placeholder copy of the
  (mislabelled) 4105 entry; it now carries the medium-voltage rule's own
  parameters, including the 0.05 p.u. active-power threshold that neither
  other code shares.  Nothing in the registry is provisional any more, so
  `ProvisionalGridCodeWarning` no longer fires for a shipped code.
- **VDE-AR-N 4120** (high voltage, 110 kV, three variants) is registered
  under `"4120"`.
- `check_var_q()` validates `net.sgen.var_q` against the selected code's
  variant count, and `bus_voltage_range()` reads the operating band a
  network permits.
- A dead-band figure in the user guide, drawn against the capability area
  it sits inside, and the Q(P) / Q(U) figures regenerated from the
  corrected envelopes — the previous pair drew the unclipped bug.

### Fixed

- **Multi-period Q-control dispatched no reactive power at all**
  (`src/potpourri/technologies/sgens.py`).  `QsGmax` / `QsGmin` were derived
  from the `q_mvar` profile, which SimBench ships as **zero** for PV, so
  `qsG` was pinned to zero and every multi-period Q-control constraint —
  Q(P), Q(U), the inverter S² circle — was satisfied trivially.  The model
  looked Q-controlled and was not: on `1-LV-rural1` the solver could reach
  only 5.1e-27 p.u. of reactive power where the grid code permits 3.8e-2.
  The single-period path had always overridden these bounds from the
  capability table; the multi-period one now does the same, and the two
  produce identical bounds for the same network.  **Any multi-period
  Q-control result from an earlier version should be rerun.**
- **Q(P) and Q(U) bounds ran past the grid code**
  (`src/potpourri/technologies/q_control.py` and the five constraint sites).
  Each bound was a single unclipped straight line, so the saturation shelves
  the standards define were missing.  Q(P) reached **+3.52·Pn** at rated
  output against a limit of +0.484, and the Q(U) band was **3.4× too wide at
  every voltage**, including nominal, where it spanned [−0.940, +1.494]
  instead of [−0.228, +0.484].  Both are now piecewise-linear envelopes,
  imposed as one linear inequality per affine piece.
- **No reactive power was feasible below P = 0.061·Pn.**  Because the
  pieces extrapolated below the standard's first active-power breakpoint,
  the lower bound rose above the upper one and the two crossed, so any sgen
  curtailed below that point made the model infeasible with nothing in the
  solver output pointing at the cause.  Each bound is now taken over the
  operating range, which substitutes its concave majorant / convex minorant
  where the exact area is non-convex.  The substitution can only widen the
  feasible band, never narrow it, and is exact from the reference point
  (0.2·Pn) to rated output.
- **The default grid code was mislabelled.**  The parameter set named
  `VDE_AR_N_4105` carried the **VDE-AR-N 4120** values: voltage breakpoints
  at 96 / 103 / 120 / 127 kV on the 110 kV base and the three 4120 variant
  pairs.  It is renamed `VDE_AR_N_4120` and remains the default, so models
  that relied on those numbers are unaffected beyond the saturation fix.
  `grid_code="4105"` now selects the real low-voltage rule.
- **Fractional `var_q` was silently truncated.**  `var_q=0.9` was accepted
  and mapped to variant 0, selecting the wrong capability column with no
  error anywhere.  Fractional values are now rejected; integral floats keep
  working, since pandas stores the column as `float64` whenever it holds
  NaN.  Out-of-range values name the code and its variant count.
- **The Q(U) pieces demanded reactive power the standard does not.**  Below
  a code's own voltage span the end segment was extrapolated rather than
  relaxed: for VDE-AR-N 4105 at 0.85 p.u. the model required Q ≥ +0.329·Pn
  where the standard requires nothing, forbidding dispatch the code
  permits.  The bus voltage range is now passed through to the Q(U) pieces,
  as it already was for Q(P).
- **`_grid_code` meant different things on the two models.**  The
  single-period model stored the resolved `GridCode`, the multi-period one
  the raw selector, so `mp._grid_code.pq_area` raised `AttributeError`.
  Both now store the resolved object.

## [0.4.0] — 2026-07-30

### Added

- **Selectable grid codes** (`src/potpourri/technologies/q_control.py`): the
  connection rules are now `GridCode` parameter sets in a `GRID_CODES`
  registry instead of module-level constants, selected model-wide via
  `ACOPF.add_OPF(grid_code=…)` or `ACOPF_multi_period.add_OPF(grid_code=…)`.
  Accepts `None` (VDE-AR-N 4105, the default, so existing models are
  unaffected), a short name (`"4105"`, `"4110"`), or a `GridCode`; unknown
  names raise `ValueError`. A `GridCode` carries the Q(U) breakpoints, the
  Q/Pn capability table and its variants, the Q(P) breakpoints, and the P(U)
  and cos(φ)(P) thresholds.
  - **`VDE_AR_N_4110` (medium voltage) is provisional**: it is registered but
    reuses the VDE-AR-N 4105 values as a placeholder, because its normative
    medium-voltage figures have not been entered yet. Results obtained with
    `grid_code="4110"` are **not** 4110-compliant. Selecting it emits a
    `ProvisionalGridCodeWarning` instead of failing.
  - `ACOPF.static_generation_wind_var_q()` previously reimplemented the
    Q-curve maths inline with hard-coded constants, so the single-period path
    would have ignored the selected code; it now routes through
    `compute_q_curves()`. Numerically identical for VDE-AR-N 4105.
  - The module-level constants (`VQU_Q_MAX`, `QP_P_HIGH`, `VPU_V_CURTAIL`, …)
    are retained as aliases of the default grid code, so existing imports
    keep working.
  - `windpower.py` no longer keeps a private copy of the VDE-AR-N 4105 table,
    nor its own reimplementation of the Q-curve maths -- the third
    byte-identical copy in the codebase. Both now come from the registry, so
    `grid_code` reaches the wind and hosting-capacity paths too. The
    simplified HC check derives its Q/P bounds from the selected code (widest
    capacitive and most negative inductive entry), and the public `qp_max` /
    `qp_min` keyword defaults on `Windpower_multi_period` follow the default
    code. Verified bit-identical for VDE-AR-N 4105: all four Q(U)
    hosting-capacity slopes match the previous hard-coded values to 1e-12, and
    the defaults still resolve to +0.48 / -0.41.
- **`ACOPF.add_OPF(sgen_types=…)`** — the sgen `type` values that
  `pv_q_control` treats as PV, defaulting to `DEFAULT_PV_SGEN_TYPES`
  (`("PV", "PV_MV", "pv")`). Previously the filter matched `type == "PV"` exactly,
  which reached **zero** sgens on every SimBench medium-voltage grid (they are
  labelled `PV_MV`, `Wind_MV`, `lv_RES`, …), so `pv_q_control` silently built
  no constraints there and raised no warning. The new default reaches 2–5
  units per MV grid; pass a wider list such as
  `("PV", "PV_MV", "Wind_MV", "lv_RES")` to include the aggregated
  LV-renewable and MV wind units, raising the reach to 87–133. LV grids are
  unaffected. Note that the wind Q-control path is selected separately and
  still matches `type == "Wind"` exactly, and that the multi-period model
  applies no type filter at all.
- **`ACOPF.add_OPF(wind_sgen_types=…)`** — the sgen `type` values treated as
  wind by the wind Q-control path (`model.WIND` / `model.WINDc`), defaulting to
  `DEFAULT_WIND_SGEN_TYPES`
  (`("Wind", "Wind_MV", "wind onshore", "wind offshore")`). The path previously
  matched `type == "Wind"` exactly — the same defect fixed for PV — so it
  reached **zero** sgens on every SimBench MV and EHV grid. Measured reach is
  now 6 / 5 / 0 / 3 on the four MV grids, previously 0 throughout. The constant
  lives in `q_control` so the single- and multi-period wind paths cannot drift
  apart; both now use it.
  - Because both paths impose the same characteristic on the same `qsG`, an
    sgen matching `sgen_types` *and* `wind_sgen_types` would get duplicate
    constraints. That now raises `SgenTypeOverlapWarning`, names the sgens and
    leaves them to the wind path. The defaults are disjoint, so it can only
    arise from a widened `sgen_types`; the documented widening example no
    longer suggests `Wind_MV` for that reason.
- **`tests/unit_tests/test_q_control.py`** — 81 solver-free unit tests for the
  Q-control code, which previously had no automated coverage at all. Covers
  the grid-code registry (capability curves against the VDE-AR-N 4105 table,
  variant ordering, alias resolution, unknown-code `ValueError`, the
  provisional-value warning contract), single-period activation for every
  `add_OPF` argument, and multi-period column-driven activation. Includes
  regression tests for the two failure modes that are otherwise silent: the
  cos(φ) cone quietly not being added without `inverter_s2=True`, and the
  `var_q` reactive-bound path, which a refactor once broke while the whole
  suite still passed.
- **`scripts/grid_code_q_strategies.py`** — worked example of both selection
  mechanisms: one snapshot solved under each registered grid code (surfacing
  the provisional-values warning rather than silencing it), then Q(P)/Q(U),
  fixed cos(φ), cos(φ)(P) and P(U) curtailment assigned to different PV units
  within a single multi-period model.
- **Documentation** for per-sgen strategy assignment in
  `docs/user-guide/reactive-power-control.md`, including which strategies may
  share an sgen: fixed cos(φ) and the cos(φ)(P) profile must not (both are
  equalities on the same reactive power), while P(U) curtailment may be
  combined with a Q rule since it constrains active power. Per-row assignment
  is a multi-period feature; in the single-period model `pu_curtail` and
  `cos_phi_p_profile` are model-wide switches.
- **VDE-AR-N 4105 / BDEW reactive-power control for PV and wind sgens**
  (`src/potpourri/technologies/q_control.py`,
  `src/potpourri/technologies/pv.py`,
  `src/potpourri/technologies/sgens.py`):
  - `compute_q_curves()` returns the slope/intercept coefficients for all
    three Q-control variants (`var_q` 0–2).
  - `ACOPF.add_OPF(pv_q_control=…)` accepts `None` (off), `"qp"` (Q(P)
    only), `"qu"` (Q(U) only), or `"both"`. Adds `PV_QP_pos/neg` and
    `PV_QU_min/max` on set `PVc`.
  - `ACOPF_multi_period` detects `var_q` in `net.sgen` automatically and
    adds the time-indexed `sG_QP_pos/neg` / `sG_QU_min/max` constraints.
- **Inverter S² apparent-power circle**: `ACOPF.add_OPF(inverter_s2=True)`
  (opt-in; defaults to `False`) reads `net.sgen.sn_mva` and
  `net.sgen.converter_sizing_pu` and adds
  `psG[g]² + qsG[g]² ≤ S_inv[g]²`. The multi-period model enables the
  time-indexed form automatically when `sn_mva` is present.
- **cos(φ) cone** completing the PV operating region (P ≥ 0 ∩ S² circle ∩
  cone): `|qsG[g]| ≤ psG[g] · tan(arccos(cos_phi_min))`, from
  `net.sgen["cos_phi_min"]` or the `cos_phi_min` keyword. In the
  single-period model the cone sits inside the S² block, so it also
  requires `inverter_s2=True` and a usable `sn_mva`.
- **P(U) active-power curtailment** (VDE-AR-N 4105 §8.5): above a voltage
  threshold the allowed active output falls linearly to zero.
  `add_OPF(pu_curtail=True)`, thresholds from `net.sgen.v_curtail_pu` /
  `v_max_curtail_pu` (defaults 1.06 / 1.10 p.u.). Bilinear — needs an NLP
  solver.
- **Fixed cos(φ) mode**: `qsG[g] = psG[g] · tan(arccos(cos_phi))` as an
  equality, via `add_OPF(fixed_cos_phi=…)` or `net.sgen["fixed_cos_phi"]`.
- **cos(φ)(P) profile**: quadratic equality
  `qsG · (Pn − Pt) = tan_phi · psG · (psG − Pt)`, via
  `add_OPF(cos_phi_p_profile=True)`.
- **Multi-period OLTC tap optimisation**: `add_tap_changer_linear()` and
  `add_tap_changer_discrete()` now use the time-indexed `Tap[tr, t]`
  variable, and `add_tap_changer_linear()` accepts an optional
  `max_tap_change_per_step` rate limit.
- **`scripts/q_control_opf.py`** — compares the Q-control modes in the
  single-period case, then the PV inverter controller modes, then runs a
  24-step multi-period AC OPF with automatic Q-control detection.
- **`docs/user-guide/reactive-power-control.md`** — user-guide page for all
  of the above, including an activation-reference table giving the exact
  precondition for each constraint group in both the single- and
  multi-period model, since a missing precondition is silent (the
  constraint is simply not added, with no warning).
- **Docs CI**: `.github/workflows/ci.yml` gained a `docs` job running
  `mkdocs build --strict`, so broken internal links and invalid navigation
  now fail the build.

### Fixed

- **Bus numbering: the auxiliary buses that pandapower inserts for node-node
  switches are handled correctly.** Previously this left every medium-voltage
  grid either wrong or unusable.

  `Basemodel.__init__` truncated the ppc bus table to `len(net.bus)` rows and
  used those ppc bus numbers as `model.B`. Where the ppc table is longer —
  `1-MV-rural--0-sw` has 103 ppc buses for 97 pandapower buses — that slice
  kept auxiliary buses while **dropping real pandapower buses along with the
  in-service branches attached to them**, so the optimisation ran on an
  incomplete network. It also let `net._pd2ppc_lookups["bus"]` resolve to a bus
  the model did not contain, raising `KeyError` in `pyo_to_net` on write-back
  (observed on `1-MV-rural--0-sw`, latent on `1-MV-comm--0-sw`). The
  multi-period model already used the full ppc table but failed earlier still,
  with `IndexError` in `add_OPF`, so multi-period AC OPF had never run on any
  SimBench MV grid.

  SimBench documents these nodes: in its CSV format all switches are modelled
  as node-node switches, which inserts `auxiliary`-type nodes between busbars
  and the edge elements attached through them. SimBench's `no_sw` variant reduces its
  own auxiliary nodes, as its documentation describes — 95 rather than 97
  pandapower buses for `1-MV-rural` — but pandapower still derives the
  same number of auxiliary ppc buses either way (6 in both), so choosing
  `no_sw` does not avoid this defect.

  `model.B` now spans the whole ppc bus table, and a new set `model.Bpd` holds
  the buses that a pandapower bus maps onto. Voltage limits and their bound
  constraints are indexed over `Bpd`, because auxiliary nodes have no
  pandapower row and hence no user-supplied limits; their voltage follows from
  the power-flow equations. Where the ppc has no auxiliary buses `Bpd == B` and
  nothing changes, which covers every LV grid.

  `get_v_limits()` now returns `pandas.Series` keyed by ppc bus number instead
  of arrays in pandapower positional order. That also repairs
  `add_generator_v_limits()`, which indexed a positional array with ppc bus
  numbers — correct only while the two numbering spaces coincided. Where
  several pandapower buses fuse onto one ppc bus (the multi-period path, which
  does not run `preprocess_grid`), the tightest band is kept.

  Eight regression tests cover the ppc/pandapower split, including the
  compatibility guarantee that `Bpd == B` without auxiliary buses.

- `Basemodel` no longer deletes buses aggressively during network
  preparation, which raised `KeyError` on HV/MV grids.
- `mkdocs.yml` now enables the `admonition` markdown extension. It was
  missing, so every `!!! note` block in the documentation rendered as
  literal text.
- `mkdocs.yml` disables `mkdocs-bibtex` inline citations. Inline (bare)
  `@key` parsing runs on the raw markdown before code fences are handled, so
  the Pyomo decorator `@model.Constraint(...)` in the device-development
  example was read as a citation key and failed a strict docs build.
  Bracketed `[@key]` citations are unaffected.

### Changed

- Example scripts in `scripts/` use a uniform style: tunable parameters sit
  in a configuration block of module-level constants instead of
  command-line arguments. Some previously infeasible example setups were
  corrected.

## [0.3.1] — 2026-05-23

### Added

- **PGLib-OPF git submodule** at `benchmarks/pglib-opf/` — the
  [IEEE PES Power Grid Library](https://github.com/power-grid-lib/pglib-opf)
  benchmark dataset is now included as a submodule. Clone with
  `git clone --recurse-submodules` to enable benchmark runs. The submodule
  is not part of the PyPI package; `pip install opf-potpourri` is unaffected.
- Documentation updates: README "Benchmarking against PGLib-OPF" section,
  `docs/scripts/examples.md` entries for `pglib_benchmark.py` and
  `timeseries_acopf.py`, developer clone command updated to
  `--recurse-submodules`.

## [0.3.0] — 2026-05-23

### Added

- `scripts/timeseries_acopf.py` — pandapower `run_timeseries` demo that
  drives loads and sgens from SimBench profiles via
  `ConstControl` / `DFData`, runs the potpourri AC OPF at every step via
  a custom `run=run_acopf` callback, and logs results with
  `OutputWriter`.
- **PGLib-OPF benchmark integration**
  - New `potpourri.benchmarks` package with `load_pglib_case`, which parses
    a PGLib `.m` file (via `pandapower.converter.matpower.from_mpc`) into a
    pandapower network ready for OPF, attaches MATPOWER `ANGMIN` / `ANGMAX`
    onto `net.line` / `net.trafo`, and rebalances the initial dispatch so
    `pp.runpp` converges on PGLib snapshots.
  - `PGLIB_BASELINE_TYP` / `PGLIB_BASELINE_API` / `PGLIB_BASELINE_SAD` dicts
    parsed from upstream `BASELINE.md`, plus `parse_baseline_md` helper.
  - `scripts/pglib_benchmark.py` runs DC and AC OPF on a configurable
    PGLib subset and writes `results/pglib_benchmark.{csv,md}` with a
    BASELINE.md-style markdown table.
- **Polynomial generator cost objective**
  - `potpourri.models.cost_objective.add_poly_cost_objective` wires
    `net.poly_cost` polynomial coefficients (`c2·P² + c1·P + c0`) into a
    Pyomo objective over `ext_grid + gen + sgen`. Raises a clear error if
    `net.pwl_cost` is non-empty (PWL costs not yet supported).
- **Configurable AC-OPF flags** on `ACOPF.add_OPF`:
  - `thermal_limit="current"|"mva"` — choose between the current-based
    constraint `|S|² ≤ SLmax²·v²` (default, physically meaningful for
    distribution conductors) and the MATPOWER / PowerModels-compatible
    constant-MVA limit `|S|² ≤ SLmax²`.
  - `free_slack_vm=True` (new default) — the slack-bus voltage magnitude
    floats within `[Vmin, Vmax]` for AC OPF, while the reference angle
    remains fixed. Set `False` to reproduce the legacy AC-PF behaviour
    where the slack `vm` is pinned to its base-case value.
  - `fix_hv_buses=False` (new default) / `hv_bus_kv=110.0` — replace the
    hardcoded 110 kV bus pinning with an explicit opt-in flag.
  - `angle_limits=False` — optional enforcement of branch
    phase-angle-difference constraints `angmin ≤ δ_from − δ_to ≤ angmax`
    from `net.line.angmin_degree` / `net.line.angmax_degree`. Available
    on both AC and DC OPF.
- **`net.impedance` branch support** in single- and multi-period AC / DC
  OPF. pandapower's third branch table — used for branches whose
  from/to buses have different nominal voltage but no off-nominal tap —
  is now treated as additional lines (per-unit physics is identical) and
  enters `model.L` with synthetic indices `≥ len(net.line)`. Without
  this, several PGLib cases (e.g. `case118_ieee`, `case89_pegase`,
  `case200_activ`) had electrically isolated buses or under-modelled
  parallel paths.
- **Regression test suite for the formulation audit fixes** in
  `tests/unit_tests/test_audit_fixes.py` — covers D1–D14 plus the
  impedance-branch fix and the multi-period DC OPF construction path.
- **CHANGELOG.md** (this file).

### Changed

- **Slack-bus voltage in AC OPF** is now free within `[Vmin, Vmax]` by
  default; the reference angle stays fixed. The previous pinning is
  available via `free_slack_vm=False`.
- **110 kV bus voltage pinning** is no longer applied by default. Re-enable
  with `fix_hv_buses=True` (and optionally `hv_bus_kv=...`).
- **Storage docstring** in `Basemodel.add_storage` now states explicitly
  the round-trip-efficiency convention (`η·Pchg − Pdis/η`), the convex
  relaxation `Pchg + Pdis ≤ Pmax` used in place of complementarity, and
  the single-period nature of the SOC update.
- **AC branch admittance comment** in `AC.py` now documents the explicit
  assumption that branch shunt conductance is zero (the standard MATPOWER
  convention).
- **`Basemodel.preprocess_grid`** now references a comprehensive list of
  element tables (lines, trafos, trafo3w, impedance, dcline, load, sgen,
  gen, ext_grid, shunt, storage, ward, xward, motor, asymmetric loads
  and sgens) for both the bus-bus switch merge and the orphan-bus prune.
- **PGLib benchmark script** no longer monkey-patches the model;
  PGLib-compatible behaviour is now driven entirely by the new
  `ACOPF.add_OPF` flags (`thermal_limit="mva"`, `free_slack_vm=True`,
  `fix_hv_buses=False`, `angle_limits=True`).
- **Generator-cost helper** writes the polynomial cost in MW units
  (per-unit pyomo variables multiplied by `baseMVA`), matching the
  MATPOWER `mpc.gencost` polynomial convention.

### Fixed

- **`OPF.generation_real_power_limits`** and
  **`ACOPF_base.generation_reactive_power_limits`** now filter source
  tables to in-service rows before broadcasting into the
  `net._gen_order[element]` slice. Out-of-service generators no longer
  trigger `ValueError: could not broadcast …` (e.g. PGLib
  `case200_activ` with 11 OOS gens).
- **Degenerate generator P-range** (`Pmin == Pmax`) — synchronous
  condensers in PGLib (`case30_ieee`, `case60_c`, `case240_pserc`, …)
  no longer drive IPOPT into `TOO_FEW_DEGREES_OF_FREEDOM`. The variable
  is pinned via tight ε-padded bounds while the redundant range
  constraint is skipped; the variable remains in the NL file so IPOPT
  doesn't lose a degree of freedom to Pyomo's fixed-variable elimination.
- **`Basemodel.preprocess_grid`** no longer silently drops buses that
  are referenced only by `shunt` or `storage` (the prior
  `referenced_buses` union omitted those tables, which crashed
  `pp.create_continuous_bus_index` on, e.g., PGLib `case200_activ`).
- **Bus-bus switch merging** now updates `shunt`, `storage`, `trafo`,
  `trafo3w`, `ward`, `xward`, and the asymmetric load/sgen tables when
  remapping a removed bus, not only `line` / `load` / `sgen` / `gen` /
  `ext_grid`.
- **`pyo_to_net._generation_results_to_net`** now uses key-aware
  iteration (`model.G` ∩ `_gen_order[et]` slice) rather than positional
  `range(f, t)` indexing, making it robust to non-contiguous generator
  indices.
- **`pyo_to_net._line_results_to_net`** now writes back line flows via
  `net.line.index`-keyed lookup with `dict.get(idx, 0.0)`. Out-of-service
  lines and impedance-only synthetic indices in `model.L` no longer
  produce length mismatches.
- **DC line-susceptance lookup** (`DC.BL_data`, `DC_multi_period.BL_data`)
  is now a `pd.Series` keyed by the actual `model.L` indices, so cases
  with out-of-service lines (e.g. PGLib `case2000_goc` with 6 OOS lines)
  no longer fail with `KeyError: Index 'k' is not valid for indexed
  component 'BL'`.
- **`DC_multi_period.create_model`** has been rewritten:
  - The dead `if self.T is None: ... else: ...` split has been removed
    (`self.T` is always an int set by `Basemodel_multi_period.__init__`).
  - Constraint indexing now uses the Pyomo Set `self.model.T` rather
    than the bare int `self.T` (the cause of
    `TypeError: Cannot create a Set from data that does not support
    __contains__. ... received 'int'`).
  - Single-period parameters (`BL[l]`, `BLT[l]`, `shift[l]`, `GB[s]`)
    are no longer accessed with a spurious time index.
  - Time-indexed variables (`delta[bus, t]`, `deltaL[l, t]`,
    `deltaLT[l, t]`, `pLfrom/pLto/pThv/pTlv[l, t]`) are correctly
    addressed throughout.
- **`DC_multi_period.__init__`** now accepts a `pf=1` keyword so
  subclasses (`DCOPF_multi_period`) can forward it through the MRO.
- **`DCOPF_multi_period`** rewritten to follow the single-period DCOPF
  pattern: `__init__` no longer triggers a duplicate `create_model()`,
  and an explicit `add_OPF` method attaches the line/transformer
  thermal limits using the correct variable names (`pLfrom`, `pThv`
  instead of the previously non-existent `pL`, `pLT`).
- **`Sgens_multi_period.get_opf_parameters`** no longer references
  `self.QsGmax_tuple` / `self.QsGmin_tuple`, which are populated only
  by the AC-only step `static_generation_reactive_power_limits`.
  Reactive-power parameter setup has moved into a new
  `get_acopf_parameters` hook called only via `get_all_acopf`. Without
  this, any multi-period DC OPF call raised
  `AttributeError: 'Sgens_multi_period' object has no attribute
  'QsGmax_tuple'`.

### Notes / unresolved

- For several large PGLib transmission cases with many `net.impedance`
  rows (`case200_activ`, `case588_sdet`, `case793_goc`, `case2312_goc`),
  the AC and DC objectives match the PGLib reference within a consistent
  ±8 % to ±32 % band rather than to machine precision. case197_snem
  (55 impedance rows) and case118_ieee (2 impedance rows) match to
  ≤ 0.01 %, indicating the impedance-branch fix is mathematically
  correct; the remaining gap on the larger cases is the most plausibly
  explained by PowerModels.jl's input-preprocessing pipeline
  (e.g. topology simplification, low-impedance branch merging) that
  potpourri does not perform. Side-by-side investigation against
  PowerModels.jl is the recommended next step.
- The pandapower import warning `some BR_B of transformers in
  ppc['branch'] are positive, but transformers are inductive` was
  investigated and confirmed to be **benign**: it reflects pandapower's
  internal `i0_percent` interpretation. potpourri reads `BR_B` directly
  from `_ppc['branch']` and treats it as the generic line-charging
  susceptance, matching MATPOWER and PowerModels conventions. The
  underlying data is roundtripped to float epsilon.

## [0.2.1]

### Changed
- README: split installation into a user section (`pip install opf-potpourri`)
  and a developer section (conda environment + editable install).
- README: added PyPI badge and direct links to the Read the Docs sections.
- README: moved solver installation guidance into a dedicated Solvers section
  and noted that PyPI users must install solvers separately.

## [0.2.0] — initial public release

### Added
- Single-period AC OPF (`ACOPF`, `HC_ACOPF`) and DC OPF (`DCOPF`) models.
- Multi-period AC OPF (`ACOPF_multi_period`) with a configurable time horizon
  and SimBench profile integration.
- Flexible-resource mix-in modules: `Battery_multi_period`,
  `HeatPump_multi_period`, `PV_multi_period`, `Windpower_multi_period`,
  `Demand_multi_period`, `Sgens_multi_period`, `Generator_multi_period`.
- Hosting-capacity OPF (`HC_ACOPF`, `HC_ACOPF_multi_period`) with binary
  placement variables and VDE-AR-N 4105 reactive-power constraints.
- Post-processing utilities (`pyo_to_net`, `pyo_to_net_multi_period`) that
  write Pyomo solutions back to `net.res_*` DataFrames.
- Warm-start helpers (`init_pyo_from_pp_res`,
  `init_pyo_from_pp_res_multi_period`) based on pandapower power-flow results.
- Runnable example scripts covering single-period OPF, multi-period planning,
  hosting-capacity analysis, feasible-operation-region computation, and
  solver benchmarking.
- MkDocs documentation with API reference, mathematical formulations, and
  class architecture diagrams.
- GitHub Actions CI: lint (Ruff), unit tests (Python 3.10–3.12), and
  installation smoke test.
- MIT licence.
