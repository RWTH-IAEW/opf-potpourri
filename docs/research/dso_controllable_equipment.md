# Controllable DSO network equipment in the OPF: audit, research and design

**Status:** design note and modelling survey for the opt-in controllable-equipment
feature set (on-load tap changers first, switched shunts second). Written before
the implementation and updated with what was actually built; the mathematics
here is the one the code implements (`src/potpourri/models/oltc.py`,
`src/potpourri/models/shunt_control.py`) and the one the user guide
(`docs/user-guide/controllable-equipment.md`) and the mathematical-modelling
page (Section 8) document. Section 10 holds the measured validation results.

**Date:** 2026-10-01. **Verified against:** pandapower 3.5.4 (source and the
3.5.5 rendered documentation), Pyomo 6.10.1, simbench 1.6.2.

---

## 1. Problem statement

`potpourri` models network physics (AC, DC) and flexible resources
(batteries, PV, wind, heat pumps, demand, inverter Q-control). It did not
expose the equipment a distribution system operator (DSO) actually operates as
OPF decision variables in a pandapower-aligned, validated way. The question
this note answers for each candidate device is:

1. Is it relevant to distribution-system operation?
2. How does pandapower represent it, exactly (3.5.x semantics)?
3. What does `potpourri` already do with it?
4. Which variables and constraints does co-optimisation need, in the AC, DC and
   multi-period formulations?
5. Which solver class results?
6. Verdict: implement now / design now, defer implementation / out of scope.

The hard constraint on everything below: **default behaviour must not change.**
`ACOPF(net); add_OPF(); solve()` has to produce the same numbers as before, and
no integer variable may appear unless a user asks for one.

---

## 2. Repository audit (architecture note)

### 2.1 Data flow and where the transformer lives

```
pandapowerNet
  -> Basemodel.__init__            deep copy, preprocess_grid (bus-bus switch fusion),
                                   pp.runpp -> net._ppc tables (bus, branch, gen)
  -> Basemodel.create_model        Sets B, L, TRANSF, ...; Params A, AT, shift, ...;
                                   Vars delta, pLfrom/pLto, pThv/pTlv, Tap (FIXED)
  -> AC.__init__ / create_model    per-branch admittance Params Gii/Bii/Gik/Bik (lines)
                                   and GiiT/BiiT/GikT/BikT (transformers), Vars v, q*,
                                   KCL_real/KCL_reactive, KVL_* branch equations
  -> OPF / ACOPF.add_OPF           limits, thermal ratings, Q capability, objectives
  -> solve()                        Pyomo SolverFactory
  -> pyo_to_net.pyo_sol_to_net_res  writes net.res_*
```

Files: `src/potpourri/models/basemodel.py`, `AC.py`, `OPF.py`,
`ACOPF_base.py`, `pyo_to_net.py`; multi-period mirrors in
`src/potpourri/models_multi_period/` (`basemodel_multi_period.py`,
`AC_multi_period.py`, `OPF_multi_period.py`, `pyo_to_net_multi_period.py`).

### 2.2 How transformers are represented today

* **There is no Y-bus in the model.** Every branch has its own four admittance
  parameters and its own four power-flow equalities (`KVL_real_fromTransf`,
  `KVL_real_toTransf`, `KVL_reactive_fromTransf`, `KVL_reactive_toTransf`).
  The admittance parameters are read from `net._ppc["branch"]`
  (`BR_R`, `BR_X`, `BR_B`, pandapower's `BR_G`) and split into self and mutual
  terms, `G_ii = g + g_c/2`, `B_ii = b + b_c/2`, `G_ik = -g`, `B_ik = -b`.
* **The tap ratio is an explicit Pyomo variable**, `model.Tap[t]`, initialised
  from the ppc `TAP` column and **fixed** in `Basemodel.create_model`. The
  transformer equations use it as MATPOWER does: the ideal transformer of
  ratio `Tap * exp(j*shift)` sits on the from (HV) side, so the HV self term
  carries `1/Tap**2`, the mutual terms `1/Tap`, and the LV self term nothing.
* **Consequence for this work:** the tap is *not* frozen into a fixed
  admittance matrix. Unfixing `Tap` is enough to make the HV-side ratio a
  decision variable. What is *not* available is a second ratio on the LV side,
  which pandapower's LV-side tap semantics need (Section 3.1).
* **Phase shift** is a fixed parameter `model.shift[t]` (ppc `SHIFT`, radians)
  that enters the angle difference. Each transformer rule branches on
  `if model.shift[l]:` at construction time to keep the expression small, so a
  variable phase shift would need the rules to reference a variable instead
  (deferred, Section 9.1).
* **Three-winding transformers** are not in the model. pandapower converts
  each `net.trafo3w` row into three ppc branches placed after the two-winding
  transformers; `Basemodel` slices only `[n_line, n_line + n_trafo)` for
  `TRANSF` and the impedance rows for `L`, so trafo3w branches and their
  auxiliary star bus are silently dropped from the equations. Documented as a
  pre-existing limitation; not changed here.
* **Verified numerically** (`scratchpad/verify_trafo.py`, 110/21 kV
  transformer on 110/20 kV buses, 30 kW iron losses, 0.1 % magnetising
  current, Dyn5 shift, taps -9, 0, +4, +9 on both sides): the residuals of all
  four transformer equations at pandapower's own power-flow solution are
  below `2e-15` p.u. The fixed-tap formulation reproduces pandapower exactly on
  both tap sides.
* **One pre-existing gap found:** pandapower's default `trafo_model="t"`
  converts the T equivalent to a pi equivalent (`build_branch._wye_delta`) and
  can produce *asymmetric* shunt admittances (`BR_G_ASYM`, `BR_B_ASYM`) when the
  leakage split `leakage_resistance_ratio_hv` / `leakage_reactance_ratio_hv`
  is not 0.5. `potpourri` reads only `BR_G`/`BR_B` and uses the same half at
  both ends. With the default 0.5 split the asymmetric columns are exactly
  zero and nothing is lost; with a 0.3 split the LV-end residual is
  `1.2e-4` p.u. This is independent of the OLTC work and is recorded as a
  limitation (Section 11).

### 2.3 Legacy tap methods already in the package

`OPF.add_tap_changer_linear()` / `add_tap_changer_discrete()` and their
multi-period twins exist, are mentioned in `docs/architecture.md`,
`docs/user-guide/single-period.md` and `docs/user-guide/reactive-power-control.md`,
and are used by two example scripts. The audit found:

| aspect | `add_tap_changer_linear` | `add_tap_changer_discrete` |
|---|---|---|
| selection | every transformer in `TRANSF`, no eligibility check; NaN tap data produces NaN bounds | same |
| HV-side ratio | bounds from a copy of pandapower's `_calc_tap_from_dataframe` (correct, incl. nominal mismatch) | `Tap = 1 + (k - k_n) s` — **wrong whenever `vn_hv/vn_lv` differs from the bus-voltage ratio** (e.g. a 110/21 kV unit on 110/20 kV buses gives 1.060 instead of 1.0095 at k = 3) |
| LV-side ratio | bounds correct | `Tap = 1/(1 + (k - k_n) s)` — same nominal-ratio defect, and **ignores that pandapower refers the impedance to the tapped LV winding** (admittances scale with `1/n²`) |
| integer start | — | `Tap_pos` initialised at 0 regardless of `tap_pos` |
| multi-period | per-step ratio bounds, optional ratio rate limit | no movement/count constraints |
| tests | none | none |

They are **kept unchanged and deprecated** (a `DeprecationWarning` points to
the new API) so existing callers keep getting exactly what they got; the new
`enable_oltc` is the supported path.

### 2.4 Shunts today

`net.shunt` is a fixed admittance: `GB = p_mw*step/baseMVA`,
`BB = -q_mvar*step/baseMVA`, entering the balance as `GB*v²` and `-BB*v²`
(voltage-squared dependence already correct). Two details differ from
pandapower: the `(V_bus,n / vn_kv,shunt)²` factor pandapower applies when a
shunt's rated voltage differs from its bus voltage is not applied (identical
when `vn_kv` is left at the bus default, which `create_shunt` does), and
`step_dependency_table` / `shunt_characteristic_table` are ignored. The
multi-period `Shunts_multi_period` module only carries the same constants over
time. No shunt is controllable.

### 2.5 DER reactive power today

Already covered, and reused rather than duplicated:

* box limits `QsGmin/QsGmax` from `net.sgen.min_q_mvar/max_q_mvar`, generator
  and external-grid limits from `min_q_mvar/max_q_mvar`;
* apparent-power circle `psG² + qsG² <= S²` (`add_OPF(inverter_s2=True)`,
  `net.sgen.sn_mva`), cos(phi) cone (`cos_phi_min`), fixed cos(phi),
  cos(phi)(P) profile;
* grid-code Q(P) / Q(U) capability areas and Q(U) dead-band characteristics
  (VDE-AR-N 4105/4110/4120) in `technologies/q_control.py`, with parity tests
  against pandapower's `PQVArea` classes;
* battery converter circle `BAT_inverter_s2`, storage circle `stor_inverter_cap`.

Verdict for item 14 of the brief: adequate; no change. The OLTC and shunt
controls simply join these resources in the same balance.

### 2.6 Switches and topology today

`preprocess_grid` fuses closed zero-impedance bus-bus switches, drops
self-loops and the merged-away buses, and renumbers; open line/transformer
switches are left to pandapower, which routes the branch to an auxiliary ppc
bus. Topology is therefore **fixed at model construction**. There is no branch
status variable anywhere.

### 2.7 Integer variables and solver handling today

* `HC_ACOPF.y[w]` (binary placement), Q(U) dead-band piecewise blocks, and the
  legacy `Tap_pos` are the only integer variables. The hosting-capacity script
  uses `gurobi_direct_minlp`; `docs/user-guide/solvers.md` documents
  `gurobi_direct_minlp` and `mindtpy`.
* `solve(solver="ipopt")` on a model with free integer variables **silently
  solves the continuous relaxation**: Pyomo's NL writer flags the variables as
  integer, IPOPT ignores the flag. Nothing warned.
* Available in the development environment: IPOPT 3.14.20, GLPK 5.0,
  Gurobi 13 (`gurobi`, `gurobi_direct`, `gurobi_direct_minlp`), MindtPy.
  Not available: CBC, SCIP, HiGHS, Bonmin, Couenne. CI has IPOPT and GLPK only.

### 2.8 What a DSO cannot do with the package today (gap list)

1. Optimise a tap position that pandapower would apply the same way
   (LV-side units and non-unity nominal ratios are mis-modelled by the legacy
   discrete method; nothing validates against `pp.runpp`).
2. Select *which* transformers are controllable; mix fixed and controllable
   units.
3. Limit tap movement per step, count tap operations over a horizon, or price
   them.
4. Read an optimised tap schedule back as pandapower data.
5. Switch a capacitor bank or reactor.
6. Be told that a chosen solver ignores the integrality of the model.
7. Reconfigure the network, or control a phase-shifting transformer.

Items 1–6 are addressed here; 7 is designed and deferred (Section 9).

### 2.9 What the audit changed beyond the new feature

* `Tap_lv`, the LV-side ratio, exists on every model (fixed at 1), so the
  two-sided transformer model is available without rebuilding equations.
* `solve()` guards integrality (Section 5.5).
* `opf.diagnose()`'s replay copies optimised positions and steps onto the
  check network, so a discrete OLTC solution replays exactly; the
  diagnostics metadata names every new component in pandapower terms.
* The legacy tap methods are deprecated, not removed (Section 2.3).

---

## 3. pandapower mapping (verified against 3.5.4 source and documentation)

### 3.1 Two-winding transformer tap changer

Source: `pandapower/build_branch.py` (`_calc_tap_from_dataframe`,
`_get_trafo_shift`, `_calc_nominal_ratio_from_dataframe`,
`_calc_r_x_from_dataframe`, `_calc_y_from_dataframe`, `_wye_delta`) and the
element documentation (`elements/trafo.html`).

**Tap factor.** For `tap_changer_type == "Ratio"` with `tap_step_degree`
zero or NaN (a *longitudinal* regulator):

$$
n(k) = 1 + (k - k_{\text{neutral}})\,\frac{s}{100},
\qquad s = \texttt{tap\_step\_percent},
$$

applied to the rated voltage of the tapped winding:
`tap_side="hv"`: $V_{n,\text{HV}}^{\text{trafo}} = v_{n,\text{hv}}\,n(k)$;
`tap_side="lv"`: $V_{n,\text{LV}}^{\text{trafo}} = v_{n,\text{lv}}\,n(k)$.

**Off-nominal ratio in the ppc** (`TAP` column):

$$
\tau = \frac{V_{n,\text{HV}}^{\text{trafo}} / V_{n,\text{LV}}^{\text{trafo}}}
            {V^{\text{bus}}_{\text{HV}} / V^{\text{bus}}_{\text{LV}}}
= r_0 \cdot \begin{cases} n(k) & \text{hv-side tap} \\ 1/n(k) & \text{lv-side tap}\end{cases},
\qquad
r_0 = \frac{v_{n,\text{hv}}/v_{n,\text{lv}}}{V^{\text{bus}}_{\text{HV}}/V^{\text{bus}}_{\text{LV}}}.
$$

$r_0$ is the nominal mismatch (1 when the rated voltages equal the bus
voltages, 20/21 for a 110/21 kV unit on 110/20 kV buses).

**Impedance referral.** pandapower refers the short-circuit impedance to the
LV side using the *tapped* LV rated voltage:
`z_pu = vk/100 * (V_n,LV^trafo / V_bus,LV)^2 * sn_mva/sn_trafo`
and divides the magnetising admittance by `(V_n,LV^trafo / vn_lv_kv)^2`.
Hence for an LV-side tap **every** per-unit admittance of the branch scales:

$$
y_s(k) = \frac{y_s^{(1)}}{n(k)^2}, \qquad y_c(k) = \frac{y_c^{(1)}}{n(k)^2},
$$

where $y^{(1)}$ is the value at the neutral position; for an HV-side tap the
admittances do not depend on the tap. The wye-delta conversion of the default
T model is homogeneous in this scaling, so it holds for the pi-equivalent
values the ppc actually carries. Verified numerically: at $k = 9$,
$s = 2\,\%$, LV-side, `r/r0 = x/x0 = 1.3924 = n²` and `g/g0 = b/b0 = 0.71818 = 1/n²`;
HV-side, all ratios 1.000000.

**`tap_changer_type` semantics (3.x):**

| value | effect | OLTC support here |
|---|---|---|
| `None` / NaN | **no tap changer**: `tap_pos` is ignored by the power flow even if `tap_step_percent` is set. SimBench 1.6.2 networks arrive like this (verified: `1-MV-urban--0-sw` has `tap_pos = -1` but a ppc ratio of exactly 1.0; setting `"Ratio"` gives 0.985) | not eligible; the eligibility report tells the user to set `"Ratio"` |
| `"Ratio"`, `tap_step_degree` 0/NaN | longitudinal regulator, formulas above | **supported** |
| `"Ratio"`, `tap_step_degree != 0` | cross regulator: magnitude $\sqrt{(1 + d\cos\theta)^2 + (d\sin\theta)^2}$ and an angle $\arctan(\pm d\sin\theta/(1 + d\cos\theta))$ per step | not eligible (nonlinear ratio *and* tap-dependent shift) |
| `"Symmetrical"` | cross regulator with $\theta = 90°$ | not eligible |
| `"Ideal"` | pure phase shifter, $\theta_{tp} = (k - k_n)\,\texttt{tap\_step\_degree}$ or $2\arcsin(\tfrac12 s/100)(k - k_n)$ | not eligible; design in Section 9.1 |
| `"Tabular"` / `tap_dependency_table=True` | ratio, angle and impedance per step from `net.trafo_characteristic_table` | not eligible; a finite-state (one-hot) design is given in Section 9.3 |

**Second tap changer** (`tap2_*`): applied on top of the first by the same
code path. A transformer with a non-NaN `tap2_pos` is not eligible.

**`oltc` column:** the documentation and `create_transformer_from_parameters`
say "(short circuit only)"; it feeds the IEC 60909 correction factors and
says nothing about operation. It is *not* used for eligibility.

**`tap_min` / `tap_max`:** the documentation states they are *not* considered
by the power flow ("the user is responsible to ensure that
tap_min < tap_pos < tap_max"). They are exactly the bounds an OPF needs.

**Loading:** pandapower's `trafo_loading="current"` uses the *nameplate*
`vn_hv_kv`/`vn_lv_kv`, not the tapped values; `pyo_to_net` does the same, so
loading stays comparable at any tap.

### 3.2 Eligibility rules implemented (`oltc_eligibility(net)`)

A transformer is eligible when **all** of the following hold; otherwise the
report names the first failing reason:

1. `in_service` is True (and the row survived `preprocess_grid`);
2. `tap_side` in {`"hv"`, `"lv"`};
3. `tap_changer_type == "Ratio"`;
4. `tap_step_degree` is NaN or 0;
5. `tap_dependency_table` is False/NaN and `tap2_pos` is NaN (no second
   changer);
6. `tap_step_percent` finite and non-zero, `tap_neutral` finite,
   `tap_min` and `tap_max` finite integers with `tap_min < tap_max`;
7. `tap_pos` finite and within `[tap_min, tap_max]`;
8. the ratio pandapower actually built (`ppc TAP`) equals the one the
   "Ratio" formula predicts for `tap_pos` (self-consistency check; guards
   against data the formula does not describe).

`enable_oltc(transformers=None)` takes every eligible transformer and logs
them; an explicit list is validated and a `ValueError` names each rejected
index with its reason. Nothing is inferred from `oltc`, from controllers in
`net.controller`, or from the presence of `tap_min/tap_max` alone.

### 3.3 Three-winding transformers

`net.trafo3w` taps (`tap_side` in hv/mv/lv, `tap_at_star_point`) are converted
by pandapower into the equivalent two-winding branches. Since `potpourri`
does not model trafo3w at all (Section 2.2), they are out of scope here.

### 3.4 Shunts

`net.shunt` rows carry `p_mw`, `q_mvar` **per step at 1 p.u.**, an integer
`step` (>= 1 in the schema, 0 allowed by the controllers) and `max_step`.
pandapower builds a constant admittance
$y = (p + jq)\,\text{step}\,(V^{\text{bus}}_n/\texttt{vn\_kv})^2 / S_N$ into
the bus (`GS`, `BS`), so consumption is $S = y\,v^2$ — the voltage-squared
dependence a capacitor bank has. `step_dependency_table=True` replaces the
linear-in-step law by a table (`shunt_characteristic_table`: `step`,
`q_mvar`, `p_mw`). `DiscreteShuntController` moves `step` by `increment`
within `[0, max_step]` to hold a bus voltage within `vm_set_pu ± tol`.

### 3.5 Controllers are not OPF controls

pandapower's `DiscreteTapControl` runs the loop *set tap → power flow →
check voltage band → step one position → repeat* until the controlled bus
(default: LV side) is inside `[vm_lower_pu, vm_upper_pu]` or a tap limit is
reached; `ContinuousTapControl` computes a fractional position from the
voltage error and `t_nom`. Both are *local, single-criterion* rules with a
dead band. The OPF formulation here makes the position a decision variable
chosen jointly with DER dispatch, storage and reactive power against the
model's objective and all network limits. The controllers are used as

* semantic reference (the sign convention `tap_side_coeff * tap_sign`
  encodes the same physics as the ratio formulas above),
* a validation oracle (set the optimised position, run `pp.runpp`, compare),
* the comparison case in the demonstration script (Section 10).

---

## 4. Literature and implementation survey

Every entry below was inspected (publisher metadata via Crossref, abstract or
full text via the publisher, author page or arXiv, software sources directly);
where a detail was *not* visible in the inspected text it is not claimed.
Reference numbers refer to Section 12.

### 4.1 Transformer branch models in established tools

* **MATPOWER** [13] (manual v8.1, Section 3.2): every branch is a $\pi$ line in
  series with an ideal transformer of ratio magnitude $\tau$ and phase shift
  $\theta_{\text{shift}}$ *at the from end*,
  $Y_{ff} = (y_s + jb_c/2)/\tau^2$, $Y_{ft} = -y_s/(\tau e^{-j\theta})$,
  $Y_{tf} = -y_s/(\tau e^{j\theta})$, $Y_{tt} = y_s + jb_c/2$. The DC model
  (Section 3.7, eq. 3.24) uses $b_i = 1/(x_s \tau)$ and a shift injection
  $P_{\text{shift}} = \theta_{\text{shift}} b_i [-1, 1]^\top$.
* **PowerModels.jl** [14]: the same from-side complex ratio
  $T_{ij} = \text{tap}\, e^{j\,\text{shift}}$ (`calc_branch_t`), with
  $S_{ij} = (Y_{ij} + Y^c_{ij})^* |V_i|^2/|T_{ij}|^2 - Y_{ij}^* V_i V_j^*/T_{ij}$
  and no ratio on the to-side shunt term; the branch-flow variant writes
  $V_i/T_{ij} = V_j + z_{ij} I^s_{ij}$. The polar implementation
  (`constraint_ohms_yt_from`) is, term for term, the expression `potpourri`
  builds in `KVL_real_fromTransf` (with `tr = tap·cos(shift)`,
  `ti = tap·sin(shift)`).
* **pandapower** [15] adds what neither of the two has: the nameplate
  conversion ($v_k$, $v_{kr}$, $p_{fe}$, $i_0$), the T-model with a leakage
  split, and the *side-dependent* tap semantics of Section 3.1 (impedance
  referred to the tapped LV winding). This is why a model that merely scales
  `Tap` is not enough for LV-side changers.

`potpourri` already followed the MATPOWER/PowerModels convention exactly; the
contribution here is the second ratio that lets the same equations express
pandapower's LV-side convention without rescaling parameters (Section 5.1).

### 4.2 OLTC formulations

| family | representative sources | what they do | relevance here |
|---|---|---|---|
| Nonlinear AC with a continuous ratio | classical NR/OPF practice; Capitanescu et al. [10] treat the OLTC ratio among mostly discrete controls in a centralised voltage-management MINLP; Daratha et al. [16] OLTC + SVC as an MINLP that did not solve in 8 h and needed a two-stage scheme | ratio (or position) is a continuous variable in the AC equations | the `continuous` mode; also the relaxation step of the rounding heuristic |
| Exact linearisation in branch-flow/SOCP models | Wu, Tian & Zhang [1] (binary expansion of the discrete tap + big-M, keeps the SOCP convex → MISOCP); Tian et al. [2] (MISOCP VAR optimisation + reconfiguration with tap positions, SOC relaxation, big-M, piecewise linearisation); Ding et al. [20] (MISOCP reactive-power optimisation coordinating discrete/continuous compensators and tap ratios, robust to wind) | needed because $n^2 v^2$ is a product of a discrete and a continuous quantity in a *linear/conic* model | not needed in the polar AC model, which carries $1/a^2$ natively; recorded as the route for a future LinDistFlow/SOCP layer |
| Convex relaxations with continuous taps | Robbins, Zhu & Domínguez-García [3] (rank-constrained SDP with a virtual secondary bus, tap read off as a voltage ratio and rounded); Bazrafshan, Gatsis & Zhu [A4] (branch-flow SDP for step-voltage regulators, trilinear voltage–tap terms via McCormick); Ayyagari et al. [A3] (LinDist3Flow with a continuous effective SVR ratio as an LP, ~1 % gap) | relax, solve, round | the rounding heuristic `solve_oltc_round_and_fix` is the same idea on the exact AC model |
| Linear models with discrete taps | Borghetti [11] (single-period VVO as a MILP over tap changers, switchable capacitor banks and DG reactive output); Li, Disfani, Haghi & Kleissl [A2] (MILP coordination of OLTC and smart inverters minimising voltage deviation *and the number of tap operations*); Liu, Li & Wu [A1] (mixed-integer chordal SDP with binaries for switching status and discrete taps, continuous vs discrete tap models compared) | tractable MILP/MISDP | confirm that counting tap operations is a standard objective term, and that integer positions with linear coupling are the usual discrete representation |

**Choice made:** general integer position + affine tap factor inside the exact
polar AC equations (Section 5.2). It is the most compact exact representation
for a longitudinal regulator, its continuous relaxation is the published
continuous-ratio model, and it needs no big-M. The MISOCP/MILP literature is
the right template if a convex layer is added later.

### 4.3 Switched capacitors, reactors and continuous compensation

Baran & Wu [6] formulate capacitor placement/sizing as a mixed-integer
program with DistFlow voltage constraints; Borghetti [11], Capitanescu et al.
[10], Ding et al. [20] and Chen, Strothers & Benigni [B2] all carry switched
capacitor states as discrete controls next to the OLTC. The physical model in
every case is a constant admittance ($Q \propto v^2$), which is also
pandapower's shunt model (Section 3.4). Continuous compensation
(STATCOM/SVC as a Q source) is a reactive injection with limits — Daratha et
al. [16] coordinate an SVC with the OLTC that way — and in pandapower it is
represented by an `sgen` with `p_mw = 0` and Q limits (or by `net.svc`, which
`potpourri` does not read); the existing sgen Q limits cover it.

### 4.4 DER inverter reactive capability

Turitsyn et al. [12] compare local and centralised PV-inverter Q control
schemes; IEEE Std 1547-2018 [21] requires reactive capability of 44 % of
nameplate apparent power (injection; absorption 44 % for Category B, 25 %
for Category A), VDE-AR-N 4105/4110/4120 define the German Q(P)/Q(U) areas
the package already implements. Nothing to add.

### 4.5 Network reconfiguration

Baran & Wu [5] (branch exchange with DistFlow), Jabr, Singh & Pal [7]
(mixed-integer conic reconfiguration), Taylor & Hover [8] (MIQP/MIQCP/MISOCP
convex models), Lavorato et al. [9] and Ahmadi & Martí [18] (radiality
constraints; the latter shows "$n-1$ branches plus one parent per node" is
*not* sufficient and proposes spanning-tree constraints via the dual graph),
Liu, Li & Wu [A1] (reconfiguration + regulator + DER in one MISDP). The
consensus: reconfiguration belongs in a convex branch-flow model with explicit
radiality/connectivity constraints, not in a bus-injection polar AC model
with line-status binaries. Deferred (Section 9.2).

### 4.6 Phase-shifting transformers

Verboomen et al. [19] classify PST technologies; in optimisation models the
PST is the angle $\theta_{\text{shift}}$ of the MATPOWER/PowerModels branch,
linear in the DC model (MATPOWER eq. 3.24) and inside $\sin/\cos$ in AC.
Discrete PST steps are an integer variable on that angle. Deferred
(Section 9.1).

### 4.7 Multi-period scheduling and tap wear

Agalgaonkar, Pal & Jabr [B1] minimise the *number* of OLTC/regulator tap
operations under PV; Chen, Strothers & Benigni [B2] constrain "the maximum
allowable number of operations for OLTCs and SCs … in predefined limits",
making the day-ahead problem time-coupled; Wang et al. [17] schedule OLTC
taps and capacitor states by model-predictive control (MINLP); Xu, Dong,
Zhang & Hill [B3] schedule taps and capacitors hourly and inverter Q every
15 min (mixed-integer QP); Li et al. [A2] put the operation count into the
objective. This is the basis of Section 5.3: per-step movement limit,
horizon operation limit, and a priced movement term, all linear in the
position variable.

### 4.8 Research table

| Device | Typical decision variable | Continuous/discrete | Resulting model class (with polar AC) | Distribution relevance | pandapower representation | Recommended support |
|---|---|---|---|---|---|---|
| OLTC, longitudinal (HV/MV, MV/LV regulated units) | tap position $k$ / ratio $n(k)$ | both (discrete physical, continuous relaxation) | NLP / MINLP | high (HV/MV substations; increasingly MV/LV) | `net.trafo` `tap_*`, `tap_changer_type="Ratio"`, `tap_side` | **implement now** |
| OLTC with cross regulation (`Ratio` + angle, `Symmetrical`) | position | discrete | MINLP with tap-dependent shift | low in DSO grids | same fields + `tap_step_degree` | design, defer (finite-state) |
| Tabular tap changer / tap-dependent impedance | position (one-hot) | discrete | MINLP with admittance variables | medium (CIM/PowerFactory imports) | `tap_dependency_table`, `trafo_characteristic_table` | design, defer; detect and reject |
| Second tap changer, three-winding transformers | positions | discrete | MINLP | low / medium (HV) | `tap2_*`, `net.trafo3w` | out of scope (trafo3w not modelled) |
| Switched capacitor bank / reactor | step count $m$ | discrete (continuous relaxation) | MINLP / NLP | medium (MV compensation, reactors in cable grids) | `net.shunt` `step`, `max_step`, `q_mvar` per step | **implement now** |
| Table-driven shunt | step (one-hot) | discrete | MINLP | low | `step_dependency_table`, `shunt_characteristic_table` | detect and reject |
| STATCOM / SVC / continuous Q source | $q$ | continuous | NLP | medium | `sgen` with Q limits (`net.svc` not read) | already covered by sgen Q limits |
| DER inverter Volt/VAR | $q$ within capability | continuous (dead bands → integer) | NLP (MINLP with dead band) | high | `net.sgen` + grid-code columns | already covered |
| Network reconfiguration / tie switches | branch status $z$ | discrete | MISOCP (branch-flow) / MINLP | high | `net.switch` (`et` b/l/t) | design, defer |
| Phase-shifting transformer | shift angle / position | both | MILP (DC), MINLP (AC) | low (transmission) | `tap_changer_type="Ideal"` | design, defer |
| Voltage regulator set point / dead band | controller parameters | — | not an OPF decision | — | `DiscreteTapControl` etc. | document the workflow, do not model |

---

## 5. Mathematical formulation

### 5.1 AC transformer with ratios on both sides

Transformer $t$ joins HV bus $i$ (ppc from bus) and LV bus $j$ (to bus).
From the ppc: series admittance $y_s = g + jb$ with
$g = r/(r^2 + x^2)$, $b = -x/(r^2 + x^2)$, charging admittance
$y_c = g_c + jb_c$, shift $\varphi$ (rad). `potpourri` stores
$G_{ii} = g + g_c/2$, $B_{ii} = b + b_c/2$, $G_{ik} = -g$, $B_{ik} = -b$.

Place an ideal transformer of real ratio $a_{hv}$ on the HV side and one of
ratio $a_{lv}$ on the LV side of the pi section, the phase shift with the HV
one (MATPOWER/pandapower place it on the from side). The branch admittance
matrix becomes

$$
Y_{ff} = \frac{y_s + y_c/2}{a_{hv}^2}, \quad
Y_{tt} = \frac{y_s + y_c/2}{a_{lv}^2}, \quad
Y_{ft} = -\frac{y_s}{a_{hv} a_{lv}} e^{j\varphi}, \quad
Y_{tf} = -\frac{y_s}{a_{hv} a_{lv}} e^{-j\varphi},
$$

and the polar power-flow equations, with $\theta_{ij} = \delta_i - \delta_j$:

$$
\begin{aligned}
p^{hv}_t &= \frac{G_{ii}}{a_{hv}^2} v_i^2
  + \frac{v_i v_j}{a_{hv} a_{lv}}\big(G_{ik}\cos(\theta_{ij} - \varphi) + B_{ik}\sin(\theta_{ij} - \varphi)\big) \\
q^{hv}_t &= -\frac{B_{ii}}{a_{hv}^2} v_i^2
  + \frac{v_i v_j}{a_{hv} a_{lv}}\big(G_{ik}\sin(\theta_{ij} - \varphi) - B_{ik}\cos(\theta_{ij} - \varphi)\big) \\
p^{lv}_t &= \frac{G_{ii}}{a_{lv}^2} v_j^2
  + \frac{v_i v_j}{a_{hv} a_{lv}}\big(G_{ik}\cos(\theta_{ji} + \varphi) + B_{ik}\sin(\theta_{ji} + \varphi)\big) \\
q^{lv}_t &= -\frac{B_{ii}}{a_{lv}^2} v_j^2
  + \frac{v_i v_j}{a_{hv} a_{lv}}\big(G_{ik}\sin(\theta_{ji} + \varphi) - B_{ik}\cos(\theta_{ji} + \varphi)\big)
\end{aligned}
$$

With $a_{lv} = 1$ and $a_{hv} = \tau$ these are exactly the equations the
package has always used (`Tap` = $a_{hv}$). The Pyomo variables are
`Tap[t]` ($a_{hv}$) and the new `Tap_lv[t]` ($a_{lv}$); both are fixed by
default (`Tap` at the ppc ratio, `Tap_lv` at 1.0), so the default model is
numerically identical to the previous one — the only change is a division by
the constant 1.0 that the NL writer folds away.

**Mapping to pandapower, relative to the base case.** Let $k_0$ be
`tap_pos` when the model was built, $n_0 = n(k_0)$, $\tau_0$ the ppc ratio and
$Y_0$ the ppc admittances (Section 3.1). For a controlled transformer with
tap factor $n = n(k)$:

| tap side | $a_{hv}$ | $a_{lv}$ | admittances | why |
|---|---|---|---|---|
| hv | $\tau_0\, n/n_0 = r_0\, n$ | $1$ | $Y_0$ (tap-independent) | pandapower only rescales the HV rated voltage |
| lv | $\tau_0 = r_0/n_0$ | $n/n_0$ | $Y_0 = Y^{(1)}/n_0^2$ | substituting $\tau(k) = r_0/n$, $Y(k) = Y^{(1)}/n^2$ into the MATPOWER form gives $Y_{ff} = Y^{(1)}/r_0^2$ (constant), $Y_{ft} = -y_s^{(1)}/(r_0 n)$, $Y_{tt} = Y^{(1)}/n^2$, which is the two-sided form with $a_{hv} = r_0$, $a_{lv} = n$ when written with $Y^{(1)}$, or $a_{hv} = \tau_0$, $a_{lv} = n/n_0$ when written with the stored $Y_0$ |

At $k = k_0$ both rows reduce to the stored fixed-tap model. The
"lv" row is the statement that **for an LV-side tap the pi section referred to
the HV side is tap-invariant and the ideal transformer sits on the LV side**;
pandapower's "impedance scales with $n^2$" is the same physics seen from the
LV side.

$n_0$ is not taken from `tap_pos` but recovered from the ppc
($n_0 = \tau_0/r_0$ for hv, $r_0/\tau_0$ for lv) and compared with $n(k_0)$;
a mismatch fails eligibility rule 8.

### 5.2 Tap position, tap factor and the fixed / continuous / discrete hierarchy

For $t \in$ `TRANSF_OLTC` (and each period $\tau$ in the multi-period model):

$$
k_{t} \in [k^{\min}_t, k^{\max}_t], \qquad
n_t = 1 + (k_t - k^{\text{neutral}}_t)\, s_t/100
\quad (\texttt{trafo\_tap\_factor\_def}),
$$

$$
\text{hv: } \texttt{Tap}_t = r_{0,t}\, n_t \;(\texttt{trafo\_tap\_ratio\_hv\_def}), \qquad
\text{lv: } \texttt{Tap\_lv}_t = n_t / n_{0,t} \;(\texttt{trafo\_tap\_ratio\_lv\_def}).
$$

| mode | domain of `trafo_tap_position` | model class with the polar AC equations | meaning |
|---|---|---|---|
| fixed (default) | no variable; `Tap`, `Tap_lv` fixed | NLP (as before) | pandapower's `tap_pos` |
| `continuous` | $\mathbb{R}$, bounds $[k^{\min}, k^{\max}]$ | NLP (one extra bilinear coupling $v_i v_j / a$) | **relaxation**: the convex hull of the discrete positions in the position variable; not implementable as such; the exact LP/NLP relaxation of the discrete model, hence a valid bound and a warm start |
| `discrete` | $\mathbb{Z}$, same bounds | MINLP (nonconvex) | physical OLTC |

The discrete model uses **one general integer per transformer and period**.
Alternatives were weighed (Section 8): a one-hot selection is only needed when
the ratio is *not* affine in the position (tabular/symmetrical changers,
deferred); a binary expansion and big-M products are the tools for
*linear* models (LinDistFlow/SOCP, Wu et al. 2017), where $n^2 v^2$ has to be
linearised — the polar AC model already carries $1/a^2$ natively, so
introducing big-M would only loosen it. The LP relaxation of the integer
model is exactly the continuous mode, which makes relax-round-resolve a
well-defined heuristic (`solve_oltc_round_and_fix`).

All derived variables get tight bounds from the position bounds
(`trafo_tap_factor` in $[\min(n(k^{\min}), n(k^{\max})), \max(\cdot)]$,
`Tap` / `Tap_lv` scaled accordingly).

### 5.3 Multi-period tap scheduling

With `T` ordered and $k_{t,\tau^-}$ the position of the previous period
($k_{t,0} := $ `trafo_tap_pos_init`, the `tap_pos` of the network, for the
first period):

$$
k_{t,\tau} - k_{t,\tau^-} = u_{t,\tau} - d_{t,\tau}, \quad
u, d \ge 0 \quad (\texttt{trafo\_tap\_movement\_def}),
$$

$$
u_{t,\tau} + d_{t,\tau} \le \Delta k^{\max}_t \quad (\texttt{trafo\_tap\_change\_limit}), \qquad
\sum_\tau (u_{t,\tau} + d_{t,\tau}) \le N^{\max}_t \quad (\texttt{trafo\_tap\_operations\_limit}),
$$

$$
C^{\text{switch}} = \sum_t c_t \sum_\tau (u_{t,\tau} + d_{t,\tau}) \quad (\texttt{trafo\_tap\_movement\_cost}),
$$

added to the active objective by `penalize_tap_movement(cost)`. $u + d \ge
|\Delta k|$ always; it equals $|\Delta k|$ whenever the cost is positive or the
operation limit binds, and the reported operation count is computed from the
positions themselves, never from $u, d$. Both are bounded by
$k^{\max} - k^{\min}$, which is tight. The single-period model carries the
same movement variables relative to `tap_pos`, so "do not move the tap unless
it pays" is expressible there too. Dwell-time constraints are not built; the
up/down split is the structure they would attach to.

### 5.4 Switchable reactive compensation (implemented in Section 7)

For $s \in$ `SHUNT_CTRL`: step variable $m_s \in [0, m^{\max}_s]$, integer in
`discrete` mode, real in `continuous` mode, and the balance terms

$$
p^{sh}_s = \frac{p_s\,\rho_s}{S_N}\, m_s\, v_b^2, \qquad
q^{sh}_s = \frac{q_s\,\rho_s}{S_N}\, m_s\, v_b^2, \qquad
\rho_s = \big(V^{\text{bus}}_{n,b} / \texttt{vn\_kv}_s\big)^2,
$$

in the load (consumption) sign convention of `net.shunt` (positive `q_mvar`
is an inductive reactor, negative a capacitor), replacing the constant
$GB_s v^2$ / $-BB_s v^2$ terms for the controlled shunts only. The step
movement variables, change limit, operation limit and switching cost mirror
the OLTC ones (`shunt_step_up/down`, `shunt_step_change_limit`,
`shunt_step_operations_limit`, `shunt_switching_cost`). A shunt with
`step_dependency_table=True` is rejected (table-driven values would need a
one-hot finite-state model, Section 9.3).

### 5.5 Solver classes

| formulation | class | solvers in this project |
|---|---|---|
| AC OPF, fixed taps/steps | NLP | IPOPT (default), Gurobi 12+ (`gurobi_direct_minlp`), NEOS |
| AC OPF + continuous tap / continuous shunt | NLP (extra bilinear terms) | same |
| AC OPF + discrete tap and/or discrete shunt | nonconvex MINLP | `gurobi_direct_minlp` (global), MindtPy (local, needs a MIP sub-solver), NEOS Bonmin/Couenne; **IPOPT refused** unless `relax_integrality=True`; `solve_oltc_round_and_fix` gives a feasible rounded solution with IPOPT alone |
| DC OPF | LP | GLPK, CBC, Gurobi — no OLTC (no voltage magnitude) |
| DC OPF + discrete phase shifter | MILP (deferred) | GLPK/CBC/Gurobi |
| LinDistFlow / SOCP + discrete tap | MILP / MISOCP (not in this package) | — |

`solve()` now counts the free integer variables of the model before calling
the solver and raises a `ValueError` naming them when the solver is IPOPT
(or NEOS with a continuous solver), unless `relax_integrality=True` is passed
explicitly. Previously the relaxation happened silently.

---

## 6. API and Pyomo component design

### 6.1 Entry points

```python
opf = ACOPF(net)                      # or ACOPF_multi_period(net, toT=24)
opf.add_OPF()
opf.enable_oltc(                      # opt-in; nothing happens without this call
    transformers=None,                # None = every eligible transformer; or indices
    mode="discrete",                  # or "continuous"
    max_change_per_step=1,            # multi-period: |Δk| ≤ 1 per step   (optional)
    max_operations=6,                 # Σ|Δk| over the horizon            (optional)
)
opf.add_voltage_deviation_objective()
opf.penalize_tap_movement(cost=1e-3)  # adds c·Σ(u+d) to the active objective (optional)
opf.solve(solver="gurobi_direct_minlp")          # MINLP
# or, with IPOPT only:
opf.solve_oltc_round_and_fix(solver="ipopt")     # relax → round → fix → re-solve
```

`oltc_eligibility(net)` returns the per-transformer eligibility table and is
what `enable_oltc(transformers=None)` uses; it is public so a user can see
*why* a transformer was skipped. `enable_shunt_control(shunts=None,
mode=..., ...)` and `penalize_shunt_switching(cost)` follow the same pattern.

Both methods raise `NotImplementedError` on DC models (class
attribute `OLTC_SUPPORTED` / `SHUNT_CONTROL_SUPPORTED` on the AC layers).

### 6.2 Components (single-period; multi-period adds the time index)

| kind | name | index | meaning |
|---|---|---|---|
| Var (existing) | `Tap` | TRANSF | HV-side ratio $a_{hv}$, fixed unless hv-side OLTC |
| Var (new, fixed 1.0) | `Tap_lv` | TRANSF | LV-side ratio $a_{lv}$, fixed unless lv-side OLTC |
| Set | `TRANSF_OLTC` | — | controlled transformers |
| Param | `trafo_tap_pos_min`, `trafo_tap_pos_max`, `trafo_tap_pos_neutral`, `trafo_tap_pos_init`, `trafo_tap_step`, `trafo_tap_factor_base`, `trafo_tap_ratio_nominal` | TRANSF_OLTC | $k^{\min}, k^{\max}, k^{n}, k_0, s/100, n_0, r_0$ |
| Param | `trafo_tap_switching_cost` (mutable) | TRANSF_OLTC | $c_t$, 0 until `penalize_tap_movement` |
| Var | `trafo_tap_position` | TRANSF_OLTC (× T) | $k$; Integers or Reals |
| Var | `trafo_tap_factor` | TRANSF_OLTC (× T) | $n$ |
| Var | `trafo_tap_up`, `trafo_tap_down` | TRANSF_OLTC (× T) | $u, d \ge 0$ |
| Constraint | `trafo_tap_factor_def` | TRANSF_OLTC (× T) | $n = 1 + (k - k^n) s$ |
| Constraint | `trafo_tap_ratio_hv_def` | hv-side subset (× T) | $\texttt{Tap} = r_0 n$ |
| Constraint | `trafo_tap_ratio_lv_def` | lv-side subset (× T) | $\texttt{Tap\_lv} = n/n_0$ |
| Constraint | `trafo_tap_movement_def` | TRANSF_OLTC (× T) | $k - k^- = u - d$ |
| Constraint | `trafo_tap_change_limit` | TRANSF_OLTC (× T) | $u + d \le \Delta k^{\max}$ (if requested) |
| Constraint | `trafo_tap_operations_limit` | TRANSF_OLTC | $\sum (u + d) \le N^{\max}$ (if requested) |
| Expression | `trafo_tap_movement_cost` | — | $\sum c (u + d)$ |

Shunt control: `SHUNT_CTRL`, `shunt_step_max`, `shunt_step_init`,
`shunt_p_step`, `shunt_q_step`, `shunt_switching_cost_coeff`; Vars
`shunt_step`, `shunt_step_up`, `shunt_step_down`; Constraints
`shunt_step_movement_def`, `shunt_step_change_limit`,
`shunt_step_operations_limit`; Expression `shunt_switching_cost`. The balance
constraints `KCL_real` / `KCL_reactive` are rebuilt to carry the variable
terms (the single-period AC layer gains the `build_kcl` / `rebuild_kcl` pair
the multi-period layer already had).

All names are registered in `potpourri.diagnostics.metadata` so
`opf.diagnose()` reports them in pandapower terms.

### 6.3 Result mapping

* `net.res_trafo["tap_pos"]` and `net.res_trafo["tap_factor"]` are written for
  every transformer whenever OLTC control is enabled (uncontrolled units show
  their fixed `tap_pos` and factor 1 + (k − k_n)s, i.e. exactly what pandapower
  applied). Flows, losses and loading come from the same mapper as before.
* `net.trafo.tap_pos` is **not** modified by `solve()`. `apply_tap_positions(net=None,
  t=None)` writes the (rounded) optimised positions into a network's `trafo`
  table on request, and warns when a continuous solution is rounded or when a
  transformer's `tap_changer_type` is None (pandapower would ignore the value).
* Multi-period: `tap_schedule()` returns a DataFrame (index = time steps,
  columns = transformer indices); `tap_operations()` returns the number of
  position changes per transformer, including the move away from the initial
  position; `map_to_net(t)` writes one step.
* Shunts: `net.res_shunt["step"]`, `p_mw`, `q_mvar` at the solved voltage;
  `shunt_schedule()`, `apply_shunt_steps()`.

---

## 7. Switched shunts and coordinated reactive resources

Implemented as described in Sections 5.4 and 6 (`enable_shunt_control`,
`penalize_shunt_switching`, `shunt_schedule`, `shunt_operations`,
`apply_shunt_steps`, `solve_shunt_round_and_fix`; `shunt_eligibility(net)`
for the report). The single-period AC balance gained the `build_kcl` /
`rebuild_kcl` pair the multi-period layer already had, so the step variable
can replace the constant admittance after construction. The coordinated
reactive-power resources (inverters, batteries, generators) already existed
and needed no change — the new controls enter the same `KCL_reactive` and
are traded off by the same objective, which is the Volt/VAR coordination of
Section 4. The DC formulation gets no shunt or Volt/VAR control: it has no
voltage magnitude and no reactive power, so such controls would be
meaningless there.

---

## 8. Rejected alternatives

| alternative | why not |
|---|---|
| Re-run `pp.runpp` with a changed `tap_pos` inside the optimisation / call `DiscreteTapControl` | that is a simulation loop, not a decision variable; the OPF must see the tap's effect through its own equations |
| Keep a single `Tap` variable and rescale the admittance *parameters* by $n_0^2/n^2$ for LV-side taps | admittances would become variables (or mutable parameters mutated after construction); the two-sided ratio expresses the same physics with constants |
| Rebuild the four transformer constraints when `enable_oltc` is called, so the default model has no `Tap_lv` | more code paths (AC, AC_multi_period) producing the same equations; the fixed `Tap_lv = 1` is folded by the NL writer and is verified numerically identical (Section 10) |
| One-hot binaries $z_{t,k}$, $\sum_k z = 1$, $n = \sum_k n_k z_k$ | $(k^{\max} - k^{\min} + 1)$ binaries per transformer and period for no gain while $n$ is affine in $k$; kept as the design for tabular changers |
| Binary expansion $k = k^{\min} + \sum_j 2^j b_j$ with big-M / McCormick products | the exact-linearisation tool for LinDistFlow/SOCP models; the polar AC model has no product to linearise, and big-M would only weaken it |
| Generalized Disjunctive Programming | adds a transformation step and no modelling power beyond the integer position here |
| Optimising a voltage set point and dead band instead of the position | a controller parameterisation; the OPF decides the position directly, and the set point a DSO wants can be read off the optimised schedule (Section 10) |
| Using `net.trafo.oltc` to select transformers | short-circuit flag, unrelated to operation |
| Writing the optimised position into `net.trafo.tap_pos` automatically | destroys the input state; done only on explicit request |

---

## 9. Designed, deferred

### 9.1 Phase-shifting transformers (`tap_changer_type="Ideal"`)

Relevance: transmission (380/220 kV) and some 110 kV interconnections; rare in
DSO operation. pandapower: shift $\varphi(k) = \varphi_0 + \theta_{tp}(k)$,
linear in $k$ for `tap_step_degree`, otherwise $2\arcsin(\tfrac12 s/100)(k -
k_n)$ (also linear in $k$). Model: replace the Param `shift[t]` by a fixed Var
`shift_var[t]` (same pattern as `Tap_lv`) that the transformer rules
reference; `enable_phase_shifter(transformers, mode)` adds
`trafo_shift_position` (integer/continuous) and `trafo_shift_def`. AC: the
angle enters $\sin/\cos(\theta_{ij} - \varphi)$ → NLP/MINLP. DC: $p = -B(\delta_i -
\delta_j - \varphi)$ stays linear → LP/MILP, the one case where a tap control
belongs in `DCOPF`. Not implemented now: no pandapower network in the test
data has a PST, and the transformer rules' construction-time branch on
`shift` would have to be rewritten first.

### 9.2 Network reconfiguration

Relevance: high (MV loss reduction, restoration), but a different model class.
pandapower: `net.switch` with `et` in {`b`, `l`, `t`} and `closed`; the OPF
variable would be the status $z_e \in \{0,1\}$ of each switchable line /
transformer. Requirements that a naive line-status binary does not meet:

* branch equations must vanish when open — in the polar AC model via
  disjunctions / big-M on $p, q$ **and** on the voltage coupling (the
  equality $p = f(v_i, v_j, \theta)$ cannot simply be multiplied by $z$ without
  creating products of a binary with trigonometric terms);
* radiality (spanning tree: $\sum_e z_e = n_{\text{bus}} - n_{\text{sources}}$ plus
  single-commodity or parent-child flow constraints to exclude cycles and
  islands — Lavorato et al. 2012, Ahmadi & Martí 2015) and connectivity of
  every load;
* `preprocess_grid` currently *fuses* closed bus-bus switches and relies on
  pandapower's auxiliary buses for open ones, so switch candidates would have
  to be kept as explicit zero-impedance branches instead.

The convex route (MISOCP branch-flow, Jabr et al. 2012; Taylor & Hover 2012)
is the robust one and would be a new formulation family, not an extension of
`AC`. Deferred.

### 9.3 Tabular / symmetrical tap changers and table-driven shunts

Finite-state formulation: binaries $z_{t,k}$ for $k \in [k^{\min}, k^{\max}]$,
$\sum_k z_{t,k} = 1$, and **every** tap-dependent quantity as a convex
combination of its tabulated values: $a_{hv} = \sum_k a_{hv,k} z_k$,
$\varphi = \sum_k \varphi_k z_k$, $G_{ii} = \sum_k G_{ii,k} z_k$, … The
transformer equations then contain products $z_k v_i^2$ and $z_k v_i v_j$
(bilinear, MINLP). Exact for any `trafo_characteristic_table`, including
tap-dependent impedance; the same construction covers `shunt_characteristic_table`.
Requires admittance *variables* in the branch rules, hence deferred; today
these configurations are detected and rejected with a message rather than
approximated.

### 9.4 DC and linearised AC with OLTC

DC has no voltage magnitude and raises `NotImplementedError` from
`enable_oltc`. A linearised AC formulation is linear only while `Tap` is
constant; a variable tap needs a linearisation around the base tap
(first-order in $n$), which is a research item.

---

## 10. Validation results (`tests/unit_tests/test_oltc.py`, `test_shunt_control.py`)

All numbers below were measured on pandapower 3.5.4 / IPOPT 3.14.20 /
Gurobi 13.0.2 on 2026-10-01.

1. **Pandapower round trip, fixed positions.** Feeder: 40 MVA 110/21 kV
   unit on 110/20 kV buses ($r_0 = 20/21$), Dyn5 (150°), 30 kW iron
   losses, 0.1 % magnetising current, ±9 × 2 % changer; model built at the
   neutral position, control enabled, position fixed at $k$ and compared
   with `pp.runpp` at $k$.

   | tap side | positions | residuals of the 7 transformer/tap equations at pandapower's solution (no solver) | IPOPT solve vs `res_bus.vm_pu` | `res_trafo` P/Q both ends and losses | `loading_percent` |
   |---|---|---|---|---|---|
   | hv | −9, −4, 0, +5, +9 | ≤ 2e-15 p.u. | ≤ 1.1e-9 p.u. | ≤ 4.2e-7 MW/Mvar | ≤ 5e-7 % |
   | lv | −9, −4, 0, +5, +9 | ≤ 2e-15 p.u. | ≤ 6.2e-11 p.u. | ≤ 2.2e-8 MW/Mvar | ≤ 2.6e-8 % |

   The LV-bus voltage spans 0.88–1.29 p.u. (hv side) and 0.86–1.24 p.u.
   (lv side) over the range, so the agreement is not a small-signal
   artefact. Test tolerances are 1e-6 p.u., 1e-5 MW, 1e-4 %.
2. **Default path unchanged.** Six reference cases solved with the
   pre-change code (`simple_four_bus_system`, `1-LV-rural1--0-sw`,
   `1-MV-rural--0-sw` with both taps at −2 and typed, the feeder above on
   each tap side as a square AC power flow, and a 4-step multi-period
   `1-LV-rural1--0-sw`) were re-solved after the change: objective,
   every voltage magnitude, every angle and every transformer flow agree
   **exactly** (difference 0.0), as expected from a division by the
   constant 1.0 that the NL writer folds. Structurally the default model
   has no `TRANSF_OLTC`, both ratios fixed, `Tap_lv == 1` and no free
   integer variable.
3. **Optimisation behaviour** (30 MW load behind the feeder, slack pinned
   at 1.0 p.u., band 0.98–1.06; the fixed-tap model is infeasible):
   continuous mode moves the tap to −1.109 (hv) / +1.199 (lv) and lifts the
   LV bus to exactly 0.98; `solve_oltc_round_and_fix` and Gurobi's global
   MINLP both pick −1 (hv) / +1 (lv) with objectives 2.62e-6 / 2.66e-6 (hv)
   and 7.16e-6 / 7.23e-6 (lv), the rounded re-solve lying slightly below
   the global solution because the slack-side relaxation differs; the
   relaxed objective (≈5e-15) bounds both. `apply_tap_positions` +
   `pp.runpp` reproduces the OPF's LV voltage to 1e-6. With 10 MW PV at the
   end of a 12 km MV line the fixed tap curtails to hold 1.04 p.u. and the
   free tap moves up and admits more than 0.5 MW extra. A priced tap
   (cost 1.0 per operation) stays at its initial position where the
   unpriced one moves; a slack at 0.75 p.u. drives the position to the
   −9 bound (−8.9999985 at IPOPT's interior-point tolerance), never beyond.
4. **Multi-period** (`1-LV-rural1--0-sw`, 3 steps, band 0.97–1.03):
   rounded schedules stay within ±2, change by at most one position per
   step and use at most the allowed operations; a switching cost of 10
   removes every operation; objectives order continuous ≤ discrete ≤
   fixed.
5. **Shunts** (two-bus 20 kV feeder, four-step 0.5 Mvar bank, with and
   without a 21 kV rating on a 20 kV bus): fixed steps 0, 1, 3, 4 reproduce
   `pp.runpp` voltages to 7e-11 p.u. and the bank's P/Q to 3e-11. With the
   slack pinned at 1.0 the load bus sits at 0.950 p.u. and each step lifts
   it by about 0.0043 p.u.; a 0.96 p.u. limit is infeasible with the bank
   off and both the rounding heuristic and Gurobi switch in three or four
   steps. With a free slack a 0.99 p.u. limit is met by all four steps.
6. **Solver guard**: IPOPT on a free integer position raises the
   `ValueError`; `relax_integrality=True` solves the relaxation and returns
   a fractional position.
7. **Global MINLP on a horizon** (`1-LV-rural1--0-sw`, Gurobi 13,
   `gurobi_direct_minlp`, 300 s): one step already stops at the time limit
   with a 10.9 % gap (incumbent found within seconds); four steps with a
   120 s limit produced no incumbent at all. The rounding heuristic solves
   the same models in a few seconds. Hence the documented guidance: global
   MINLP for single-period cases, `solve_oltc_round_and_fix` (or MindtPy)
   for horizons.
8. **Demonstration** — `scripts/oltc_voltage_control_demo.py` on SimBench
   `1-LV-rural1--0-sw` (160 kVA 20/0.4 kV regulated distribution
   transformer, ±2 × 2.5 %), 24 × 15 min of the sunniest day with the
   installed PV doubled and a 1.03 p.u. upper limit: fixed vs. continuous
   vs. discrete (rounded) taps against a curtailment-plus-voltage objective,
   and pandapower's `DiscreteTapControl` for comparison. Results in
   Section 10.1. The MV rural network was the first choice (two 25 MVA
   110/20 kV units, ±9 × 1.5 %) but the multi-period AC model does not
   converge on it — see Section 11; the single-period tests cover
   110/20 kV units with ±9 positions.

### 10.1 Demonstration results

`scripts/oltc_voltage_control_demo.py`, SimBench `1-LV-rural1--0-sw`, day
146 (the sunniest noon of the SimBench year), 10:00–16:00 in 24 × 15 min
steps, installed PV and inverter ratings doubled, band 0.95–1.03 p.u., the
20 kV slack pinned at the SimBench value (1.025 p.u.), objective = curtailed
energy + 0.01·Σ(v−1)² + 0.002 per tap operation, at most one position per
step and four per horizon. Figure: `results/oltc_voltage_control_demo.png`.

| case | objective | curtailed PV [MWh] | tap operations | v_min / v_max [p.u.] | transformer max [%] | line max [%] | IPOPT time |
|---|---|---|---|---|---|---|---|
| A fixed tap (neutral) | 0.28292 | 0.280 | — | 1.0249 / 1.0300 | 51.1 | 30.1 | 4.1 s |
| B continuous tap | 0.00237 | 0.000 | 0.84 (position 0.84) | 1.0038 / 1.0250 | 94.0 | 55.2 | 2.1 s |
| C discrete tap, rounded | 0.00242 | 0.000 | 1 (position +1 all day) | 0.9999 / 1.0250 | 93.9 | 55.4 | 2.7 s (both stages) |
| E pandapower `DiscreteTapControl` (0.99–1.02 at the 0.4 kV busbar) | — | 0 (cannot curtail) | 1 (position +1) | 0.9999 / 1.0250 | 94.2 | 55.4 | 24 power flows |

Reading: at the neutral tap the 1.03 p.u. limit is the binding constraint
and the OPF has to curtail 0.28 MWh of PV over six hours; one tap position
up lowers the whole LV feeder by 2.5 % and the same PV fits without any
curtailment, with the transformer at 94 % and the lines at 55 %. The
continuous relaxation (0.84 positions) bounds the discrete objective from
below by 2 %; rounding to +1 costs nothing in curtailment. The local
controller arrives at the same position in this case because the PV lifts
the 0.4 kV busbar above its 1.02 p.u. band; it would not have reacted had
the busbar stayed inside the band while a feeder end exceeded 1.03 p.u.,
which is the situation the OPF covers and the controller cannot.

Two data-preparation lessons the script encodes: the SimBench tap changer
needs `tap_changer_type = "Ratio"` before pandapower applies it, and the
static generators carry an inverter rating `sn_mva` that the multi-period
model enforces, so scaling installed PV without the rating caps the infeed
and looks like curtailment.

### 10.2 Test-suite status

Baseline before any change: 477 unit tests, 476 passed, 1 failed
(`test_dcopf_solves_with_neos`, which needs the remote NEOS service and
fails in this environment regardless of the code). After the change: 550
unit tests, 549 passed, the same single environment-dependent failure;
every previous test passes unchanged. The two new modules add 55 OLTC
tests (parametrised over tap side and position) and 17 shunt-control
tests, including three doctests; the Gurobi-dependent ones skip without
the solver. `ruff check .`, `ruff format --check .`,
`interrogate src/potpourri` (100 %), `python tools/check_license_headers.py`
and `mkdocs build --strict` pass.

---

## 11. Limitations recorded

* Tap changer types `Symmetrical`, `Ideal`, `Tabular`, `Ratio` with
  `tap_step_degree != 0`, `tap_dependency_table=True`, second tap changers and
  three-winding transformers are detected and rejected, not approximated.
* Transformers whose `tap_changer_type` is None (every SimBench network as
  delivered) are not eligible because pandapower itself ignores their
  `tap_pos`; set `net.trafo["tap_changer_type"] = "Ratio"` first.
* Non-default leakage splits (`leakage_*_ratio_hv != 0.5`) produce
  asymmetric magnetising admittances the package does not read (pre-existing,
  ~1e-4 p.u.).
* `net.shunt.vn_kv` different from the bus voltage: the new controllable-shunt
  path applies pandapower's $(V_n^{\text{bus}}/\texttt{vn\_kv})^2$ factor; the
  pre-existing fixed-shunt path does not (pre-existing).
* Discrete modes are nonconvex MINLPs: `gurobi_direct_minlp` is the only
  global solver wired in; MindtPy is local; the rounding heuristic is a
  heuristic.
* No network reconfiguration, no phase-shifter control, no OLTC in DC.
* **Pre-existing, found while preparing the demonstration:** the
  multi-period AC model reports a locally infeasible point on the SimBench
  MV networks `1-MV-rural--0-sw` and `1-MV-rural--0-no_sw` (two steps,
  ±5 % band, PV curtailable, no tap control, pre-change source tree) even
  where `pp.runpp` at the same steps gives 1.02–1.06 p.u.; the diagnostics
  show the solver stopping with violated branch-flow equations. The
  single-period `ACOPF` on the same networks solves. The multi-period base
  model does not run `preprocess_grid` and was never exercised on an MV
  SimBench network by the test suite; this is a follow-up item independent
  of the controls.
* The multi-period base model indexes its transformer data positionally
  (`trafo_data` has a RangeIndex while `net.trafo` may not), as the
  multi-period result mapper always has; the controls inherit that
  pre-existing assumption that `net.trafo.index` is `0..n-1`, which holds
  for every SimBench network.
* `net.trafo` tables from pandapower < 3.0 (no `tap_changer_type`, a
  `tap_phase_shifter` flag) are read as "Ratio" when the flag is False; the
  ppc consistency check guards the inference.

---

## 12. References

1. W. Wu, Z. Tian, B. Zhang, "An Exact Linearization Method for OLTC of
   Transformer in Branch Flow Model", *IEEE Trans. Power Syst.* 32(3),
   2475–2476, 2017. doi:10.1109/TPWRS.2016.2603438
2. Z. Tian, W. Wu, B. Zhang, A. Bose, "Mixed-integer second-order cone
   programing model for VAR optimisation and network reconfiguration in active
   distribution networks", *IET Gener. Transm. Distrib.* 10(8), 1938–1946,
   2016. doi:10.1049/iet-gtd.2015.1228
3. B. A. Robbins, H. Zhu, A. D. Domínguez-García, "Optimal Tap Setting of
   Voltage Regulation Transformers in Unbalanced Distribution Systems", *IEEE
   Trans. Power Syst.* 31(1), 256–267, 2016. doi:10.1109/TPWRS.2015.2392693
4. M. Farivar, S. H. Low, "Branch Flow Model: Relaxations and
   Convexification—Part I", *IEEE Trans. Power Syst.* 28(3), 2554–2564, 2013.
   doi:10.1109/TPWRS.2013.2255317
5. M. E. Baran, F. F. Wu, "Network reconfiguration in distribution systems for
   loss reduction and load balancing", *IEEE Trans. Power Del.* 4(2),
   1401–1407, 1989. doi:10.1109/61.25627
6. M. E. Baran, F. F. Wu, "Optimal capacitor placement on radial distribution
   systems", *IEEE Trans. Power Del.* 4(1), 725–734, 1989. doi:10.1109/61.19265
7. R. A. Jabr, R. Singh, B. C. Pal, "Minimum Loss Network Reconfiguration Using
   Mixed-Integer Convex Programming", *IEEE Trans. Power Syst.* 27(2),
   1106–1115, 2012. doi:10.1109/TPWRS.2011.2180406
8. J. A. Taylor, F. S. Hover, "Convex Models of Distribution System
   Reconfiguration", *IEEE Trans. Power Syst.* 27(3), 1407–1413, 2012.
   doi:10.1109/TPWRS.2012.2184307
9. M. Lavorato, J. F. Franco, M. J. Rider, R. Romero, "Imposing Radiality
   Constraints in Distribution System Optimization Problems", *IEEE Trans.
   Power Syst.* 27(1), 172–180, 2012. doi:10.1109/TPWRS.2011.2161349
10. F. Capitanescu, I. Bilibin, E. Romero Ramos, "A Comprehensive Centralized
    Approach for Voltage Constraints Management in Active Distribution Grid",
    *IEEE Trans. Power Syst.* 29(2), 933–942, 2014.
    doi:10.1109/TPWRS.2013.2287897
11. A. Borghetti, "Using mixed integer programming for the volt/var
    optimization in distribution feeders", *Electric Power Systems Research*
    98, 39–50, 2013. doi:10.1016/j.epsr.2013.01.003
12. K. Turitsyn, P. Šulc, S. Backhaus, M. Chertkov, "Options for Control of
    Reactive Power by Distributed Photovoltaic Generators", *Proc. IEEE* 99(6),
    1063–1073, 2011. doi:10.1109/JPROC.2011.2116750
13. R. D. Zimmerman, C. E. Murillo-Sánchez, R. J. Thomas, "MATPOWER:
    Steady-State Operations, Planning, and Analysis Tools for Power Systems
    Research and Education", *IEEE Trans. Power Syst.* 26(1), 12–19, 2011.
    doi:10.1109/TPWRS.2010.2051168. Branch model: *MATPOWER User's Manual*
    v8.1, Section 3.2 and 3.7, https://matpower.org/docs/MATPOWER-manual.pdf
14. C. Coffrin, R. Bent, K. Sundar, Y. Ng, M. Lubin, "PowerModels.jl: An
    Open-Source Framework for Exploring Power Flow Formulations", *2018 Power
    Systems Computation Conference (PSCC)*, 2018.
    doi:10.23919/PSCC.2018.8442948. Mathematical model:
    https://lanl-ansi.github.io/PowerModels.jl/stable/math-model/ and
    `src/form/acp.jl`, `src/core/data.jl`.
15. L. Thurner, A. Scheidler, F. Schäfer, J.-H. Menke, J. Dollichon, F. Meier,
    S. Meinecke, M. Braun, "pandapower—An Open-Source Python Tool for
    Convenient Modeling, Analysis, and Optimization of Electric Power Systems",
    *IEEE Trans. Power Syst.* 33(6), 6510–6521, 2018.
    doi:10.1109/TPWRS.2018.2829021. Transformer element documentation:
    https://pandapower.readthedocs.io/en/latest/elements/trafo.html;
    source `pandapower/build_branch.py` (3.5.4).
16. N. Daratha, B. Das, J. Sharma, "Coordination Between OLTC and SVC for
    Voltage Regulation in Unbalanced Distribution System Distributed
    Generation", *IEEE Trans. Power Syst.* 29(1), 289–299, 2014.
    doi:10.1109/TPWRS.2013.2280022
17. Z. Wang, J. Wang, B. Chen, M. M. Begovic, Y. He, "MPC-Based Voltage/Var
    Optimization for Distribution Circuits With Distributed Generators and
    Exponential Load Models", *IEEE Trans. Smart Grid* 5(5), 2412–2420, 2014.
    doi:10.1109/TSG.2014.2329842
18. H. Ahmadi, J. R. Martí, "Mathematical representation of radiality
    constraint in distribution system reconfiguration problem", *Int. J.
    Electr. Power Energy Syst.* 64, 293–299, 2015.
    doi:10.1016/j.ijepes.2014.06.076
19. J. Verboomen, D. Van Hertem, P. H. Schavemaker, W. L. Kling, R. Belmans,
    "Phase shifting transformers: principles and applications", *2005
    International Conference on Future Power Systems*, 2005.
    doi:10.1109/FPS.2005.204302
20. T. Ding, S. Liu, W. Yuan, Z. Bie, B. Zeng, "A Two-Stage Robust Reactive
    Power Optimization Considering Uncertain Wind Power Integration in Active
    Distribution Networks", *IEEE Trans. Sustain. Energy* 7(1), 301–311, 2016.
    doi:10.1109/TSTE.2015.2494587
21. IEEE Std 1547-2018, *IEEE Standard for Interconnection and
    Interoperability of Distributed Energy Resources with Associated Electric
    Power Systems Interfaces*, 2018. doi:10.1109/IEEESTD.2018.8332112
- A1. Y. Liu, J. Li, L. Wu, "Coordinated Optimal Network Reconfiguration and
  Voltage Regulator/DER Control for Unbalanced Distribution Systems", *IEEE
  Trans. Smart Grid* 10(3), 2912–2922, 2019. doi:10.1109/TSG.2018.2815010
- A2. C. Li, V. R. Disfani, H. V. Haghi, J. Kleissl, "Coordination of OLTC and
  smart inverters for optimal voltage regulation of unbalanced distribution
  networks", *Electric Power Systems Research* 187, 106498, 2020.
  doi:10.1016/j.epsr.2020.106498
- A3. K. S. Ayyagari, S. A. Abraham, Y. Yao, S. Ghosh, F. Flores-Espino,
  A. Nagarajan, N. Gatsis, "Assessing the Optimality of LinDist3Flow for
  Optimal Tap Selection of Step Voltage Regulators in Unbalanced Distribution
  Networks", *2022 IEEE 61st Conference on Decision and Control (CDC)*,
  3116–3122, 2022. doi:10.1109/CDC51059.2022.9992433
- A4. M. Bazrafshan, N. Gatsis, H. Zhu, "Optimal Power Flow With Step-Voltage
  Regulators in Multi-Phase Distribution Networks", *IEEE Trans. Power Syst.*
  34(6), 4228–4239, 2019. doi:10.1109/TPWRS.2019.2915795
- B1. Y. P. Agalgaonkar, B. C. Pal, R. A. Jabr, "Distribution Voltage Control
  Considering the Impact of PV Generation on Tap Changers and Autonomous
  Regulators", *IEEE Trans. Power Syst.* 29(1), 182–192, 2014.
  doi:10.1109/TPWRS.2013.2279721
- B2. Y. Chen, M. Strothers, A. Benigni, "All-day coordinated optimal
  scheduling in distribution grids with PV penetration", *Electric Power
  Systems Research* 164, 112–122, 2018. doi:10.1016/j.epsr.2018.07.028
- B3. Y. Xu, Z. Y. Dong, R. Zhang, D. J. Hill, "Multi-Timescale Coordinated
  Voltage/Var Control of High Renewable-Penetrated Distribution Systems",
  *IEEE Trans. Power Syst.* 32(6), 4398–4408, 2017.
  doi:10.1109/TPWRS.2017.2669343

pandapower controller documentation used in Section 3.5:
https://pandapower.readthedocs.io/en/latest/control/controller.html
(`DiscreteTapControl`, `ContinuousTapControl`, `DiscreteShuntController`) and
https://pandapower.readthedocs.io/en/latest/control/run.html (`run_control`).
