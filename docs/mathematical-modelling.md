# Mathematical Modelling

This page describes every OPF formulation implemented in `potpourri`, from the raw pandapower network to the full set of Pyomo variables, constraints, and objectives.  All quantities are expressed in the **per-unit system** on the network base `baseMVA`.

---

## From pandapower to Pyomo

`Basemodel.__init__` extracts the following data from the pandapower network object and registers them as Pyomo Sets and Parameters:

| pandapower table | Pyomo Set | Description |
|---|---|---|
| `net.bus` | $\mathcal{B}$ | All buses |
| `net.bus` (type 3) | $\mathcal{B}_0$ | Slack / reference buses |
| `net.bus` (type 2) | $\mathcal{B}_{PV}$ | PV buses (voltage-controlled generators) |
| `net.line` (in service) | $\mathcal{L}$ | Lines |
| `net.trafo` (in service) | $\mathcal{T}$ | Transformers |
| `net.sgen` (in service) | $\mathcal{G}_s$ | Static generators (PV, wind, …) |
| `net.ext_grid` ∪ `net.gen` | $\mathcal{G}$ | External grids and synchronous generators |
| `net.ext_grid` (ref) | $\mathcal{G}_{ext}$ | Slack generators |
| `net.load` (in service) | $\mathcal{D}$ | Loads |
| `net.shunt` (in service) | $\mathcal{S}$ | Shunts |

### Admittance data

For each line $l \in \mathcal{L}$ the $\pi$-equivalent admittance parameters are extracted:

| Symbol | Meaning |
|---|---|
| $G_{ii}^{(l)}, B_{ii}^{(l)}$ | Self conductance / susceptance at the *from* end |
| $G_{ik}^{(l)}, B_{ik}^{(l)}$ | Transfer conductance / susceptance |

For each transformer $\tau \in \mathcal{T}$ the same four values are used with superscript $(\tau)$.  Additionally $\phi_\tau$ denotes the phase-shift angle in radians and $\tau_{\text{tap}}$ the complex tap ratio.

---

## 1  AC Power Flow (`AC`)

### 1.1  Variables

| Symbol | Pyomo name | Domain | Bounds | Description |
|---|---|---|---|---|
| $v_b$ | `v[b]` | $\mathbb{R}_{\ge 0}$ | $(0,\,2)$ p.u. | Voltage magnitude |
| $\delta_b$ | `delta[b]` | $\mathbb{R}$ | $(-\pi,\,\pi)$ rad | Voltage phase angle |
| $p_{g_s}$ | `psG[sg]` | $\mathbb{R}_{\ge 0}$ | — | Static generator real power |
| $q_{g_s}$ | `qsG[sg]` | $\mathbb{R}$ | — | Static generator reactive power |
| $p_g$ | `pG[g]` | $\mathbb{R}$ | — | Generator / ext-grid real power |
| $q_g$ | `qG[g]` | $\mathbb{R}$ | — | Generator / ext-grid reactive power |
| $p_d$ | `pD[d]` | $\mathbb{R}$ | — | Load real power (fixed by default) |
| $q_d$ | `qD[d]` | $\mathbb{R}$ | — | Load reactive power (fixed by default) |
| $P_{l}^{\text{from}}, P_{l}^{\text{to}}$ | `pLfrom[l]`, `pLto[l]` | $\mathbb{R}$ | — | Line real power (sending / receiving end) |
| $Q_{l}^{\text{from}}, Q_{l}^{\text{to}}$ | `qLfrom[l]`, `qLto[l]` | $\mathbb{R}$ | — | Line reactive power |
| $P_{\tau}^{\text{hv}}, P_{\tau}^{\text{lv}}$ | `pThv[τ]`, `pTlv[τ]` | $\mathbb{R}$ | — | Transformer real power (HV / LV side) |
| $Q_{\tau}^{\text{hv}}, Q_{\tau}^{\text{lv}}$ | `qThv[τ]`, `qTlv[τ]` | $\mathbb{R}$ | — | Transformer reactive power |
| $\tau_{\text{tap}}$ | `Tap[τ]` | $\mathbb{R}$ | — | Tap ratio (fixed by default) |

Reference-bus angles are fixed parameters: $\delta_b = \delta_b^{(0)}$ for all $b \in \mathcal{B}_0$.

### 1.2  Kirchhoff's Current Law (bus balance)

For every bus $b \in \mathcal{B}$:

**Real power balance (KCL\_real)**

$$
\sum_{\substack{g_s \in \mathcal{G}_s \\ (g_s,b)\in\mathcal{G}_s^\text{bus}}} p_{g_s}
+ \sum_{\substack{g \in \mathcal{G} \\ (g,b)\in\mathcal{G}^\text{bus}}} p_g
=
\sum_{\substack{d \in \mathcal{D} \\ (b,d)\in\mathcal{D}^\text{bus}}} p_d
+ \sum_{\substack{l \in \mathcal{L} \\ A_{l,1}=b}} P_l^\text{from}
+ \sum_{\substack{l \in \mathcal{L} \\ A_{l,2}=b}} P_l^\text{to}
+ \sum_{\substack{\tau \in \mathcal{T} \\ A^\tau_{\tau,1}=b}} P_\tau^\text{hv}
+ \sum_{\substack{\tau \in \mathcal{T} \\ A^\tau_{\tau,2}=b}} P_\tau^\text{lv}
+ \sum_{\substack{s \in \mathcal{S} \\ (b,s)\in\mathcal{S}^\text{bus}}} G_s^{(s)} v_b^2
$$

**Reactive power balance (KCL\_reactive)**

$$
\sum_{g_s} q_{g_s} + \sum_g q_g
=
\sum_d q_d
+ \sum_{\mathcal{L}} Q_l^\text{from/to}
+ \sum_{\mathcal{T}} Q_\tau^\text{hv/lv}
- \sum_{\substack{s \in \mathcal{S} \\ (b,s)\in\mathcal{S}^\text{bus}}} B_s^{(s)} v_b^2
$$

(The summation structure over generators, lines and transformers is identical to the real-power case.)

### 1.3  Branch flow equations (lines)

For each line $l \in \mathcal{L}$ let $i = A_{l,1}$ (from-bus) and $j = A_{l,2}$ (to-bus):

**Real power, from end (KVL\_real\_from)**

$$
P_l^\text{from} = G_{ii}^{(l)} v_i^2 + v_i v_j \!\left( B_{ik}^{(l)} \sin(\delta_i - \delta_j) + G_{ik}^{(l)} \cos(\delta_i - \delta_j) \right)
$$

**Real power, to end (KVL\_real\_to)**

$$
P_l^\text{to} = G_{ii}^{(l)} v_j^2 + v_i v_j \!\left( B_{ik}^{(l)} \sin(\delta_j - \delta_i) + G_{ik}^{(l)} \cos(\delta_j - \delta_i) \right)
$$

**Reactive power, from end (KVL\_reactive\_from)**

$$
Q_l^\text{from} = -B_{ii}^{(l)} v_i^2 + v_i v_j \!\left( G_{ik}^{(l)} \sin(\delta_i - \delta_j) - B_{ik}^{(l)} \cos(\delta_i - \delta_j) \right)
$$

**Reactive power, to end (KVL\_reactive\_to)**

$$
Q_l^\text{to} = -B_{ii}^{(l)} v_j^2 + v_i v_j \!\left( G_{ik}^{(l)} \sin(\delta_j - \delta_i) - B_{ik}^{(l)} \cos(\delta_j - \delta_i) \right)
$$

### 1.4  Transformer branch flow equations

For each transformer $\tau \in \mathcal{T}$ let $i = A^\tau_{\tau,1}$ (HV bus), $j = A^\tau_{\tau,2}$ (LV bus), $a^{hv} = $ `Tap[τ]` (ratio of the ideal transformer on the HV side), $a^{lv} = $ `Tap_lv[τ]` (ratio of the ideal transformer on the LV side), and $\phi = \phi_\tau$ (phase-shift angle, zero if no phase-shifter). The $\pi$ section sits between the two ideal transformers:

**HV real power (KVL\_real\_fromTransf)**

$$
P_\tau^\text{hv} = \frac{G_{ii}^{(\tau)}}{(a^{hv})^2} v_i^2 + \frac{v_i v_j}{a^{hv} a^{lv}} \!\left( G_{ik}^{(\tau)} \cos(\delta_i - \delta_j - \phi) + B_{ik}^{(\tau)} \sin(\delta_i - \delta_j - \phi) \right)
$$

**LV real power (KVL\_real\_toTransf)**

$$
P_\tau^\text{lv} = \frac{G_{ii}^{(\tau)}}{(a^{lv})^2} v_j^2 + \frac{v_i v_j}{a^{hv} a^{lv}} \!\left( B_{ik}^{(\tau)} \sin(\delta_j - \delta_i + \phi) + G_{ik}^{(\tau)} \cos(\delta_j - \delta_i + \phi) \right)
$$

**HV reactive power (KVL\_reactive\_fromTransf)**

$$
Q_\tau^\text{hv} = -\frac{B_{ii}^{(\tau)}}{(a^{hv})^2} v_i^2 + \frac{v_i v_j}{a^{hv} a^{lv}} \!\left( -B_{ik}^{(\tau)} \cos(\delta_i - \delta_j - \phi) + G_{ik}^{(\tau)} \sin(\delta_i - \delta_j - \phi) \right)
$$

**LV reactive power (KVL\_reactive\_toTransf)**

$$
Q_\tau^\text{lv} = -\frac{B_{ii}^{(\tau)}}{(a^{lv})^2} v_j^2 + \frac{v_i v_j}{a^{hv} a^{lv}} \!\left( -B_{ik}^{(\tau)} \cos(\delta_j - \delta_i + \phi) + G_{ik}^{(\tau)} \sin(\delta_j - \delta_i + \phi) \right)
$$

When $\phi = 0$ the phase-shift terms vanish. **By default both ratios are fixed**: $a^{hv}$ at pandapower's off-nominal ratio (ppc `TAP`, which already contains the tap position the network was built with) and $a^{lv} = 1$, so the equations are exactly the MATPOWER/PowerModels branch model with the transformer on the from side. The second ratio only comes into play when an LV-side on-load tap changer is made controllable (Section 8.1); it is how pandapower's convention of referring the short-circuit impedance to the *tapped* LV winding is expressed without changing the admittance parameters.

### 1.5  PV-bus voltage setpoint

For voltage-controlled generators ($b \in \mathcal{B}_{PV}$) the magnitude is fixed to the generator setpoint $v_b^{\text{PV}}$:

$$
v_b = v_b^{\text{PV}}, \quad b \in \mathcal{B}_{PV}
$$

---

## 2  DC Power Flow (`DC`)

The DC model linearises the AC equations by assuming $v_b \approx 1$ p.u. and $\delta_i - \delta_j \approx \sin(\delta_i - \delta_j)$ (small angle), and by neglecting resistance and reactive power.

### 2.1  Additional variables

| Symbol | Pyomo name | Description |
|---|---|---|
| $\Delta_l$ | `deltaL[l]` | Angle difference across line $l$ |
| $\Delta_\tau$ | `deltaLT[τ]` | Angle difference across transformer $\tau$ |

### 2.2  Real power balance (KCL\_const)

$$
\sum_{g_s} p_{g_s} + \sum_g p_g = \sum_d p_d + \sum_{\mathcal{L}} (P_l^\text{from/to}) + \sum_{\mathcal{T}} (P_\tau^\text{hv/lv}) + \sum_s G_s^{(s)}
$$

### 2.3  Line flow (KVL\_real\_from / to)

$$
P_l^\text{from} = -B_l \,\Delta_l, \qquad P_l^\text{to} = B_l \,\Delta_l
$$

$$
\Delta_l = \delta_i - \delta_j, \quad l \in \mathcal{L}
$$

where $B_l = 1/x_l$ is the line susceptance.

### 2.4  Transformer flow (KVL\_trans\_from / to)

$$
P_\tau^\text{hv} = -B_\tau \,\Delta_\tau, \qquad P_\tau^\text{lv} = B_\tau \,\Delta_\tau
$$

$$
\Delta_\tau = \delta_i - \delta_j - \phi_\tau, \quad \tau \in \mathcal{T}
$$

---

## 3  OPF Operational Constraints (`OPF`)

`OPF` is a mixin that adds box constraints for generators, loads, and transformer taps.  It is combined with `AC` or `DC` to form `ACOPF` or `DCOPF`.

### 3.1  Generator active power limits

$$
P_{g_s}^\text{min} \le p_{g_s} \le P_{g_s}^\text{max}, \quad g_s \in \mathcal{G}_s^c
$$

$$
P_g^\text{min} \le p_g \le P_g^\text{max}, \quad g \in \mathcal{G}
$$

where $\mathcal{G}_s^c \subseteq \mathcal{G}_s$ is the subset of *controllable* static generators.

### 3.2  Load limits

$$
P_d^\text{min} \le p_d \le P_d^\text{max}, \quad d \in \mathcal{D}^c
$$

### 3.3  Transformer tap constraints

Controllable tap changers are described in [Section 8](#8-controllable-dso-network-equipment) (`enable_oltc`). The older methods `add_tap_changer_linear()` and `add_tap_changer_discrete()` are **deprecated** but unchanged: the first bounds the HV-side ratio of *every* transformer by the ratios at `tap_min` and `tap_max`, the second adds an integer `Tap_pos` with `Tap = 1 + (Tap_pos − tap_neutral) s` for an HV-side and the reciprocal for an LV-side tap. The latter ignores the mismatch between rated and bus voltages and pandapower's referral of an LV-side tapped impedance, which is why Section 8 replaces it.

---

## 4  Full AC OPF (`ACOPF`)

`ACOPF` inherits from both `AC` and `OPF` and adds voltage limits, line thermal limits, and reactive power limits.

### 4.1  Voltage limits

$$
V_b^\text{min} \le v_b \le V_b^\text{max}, \quad b \in \mathcal{B}
$$

### 4.2  Line thermal limit

The apparent power at each end of a line must not exceed the thermal rating $S_l^\text{max}$.  The constraint is expressed as a **current-based quadratic inequality**:

$$
\left(P_l^\text{from}\right)^2 + \left(Q_l^\text{from}\right)^2 \le \left(S_l^\text{max}\right)^2 v_{A_{l,1}}^2, \quad l \in \mathcal{L}
$$

$$
\left(P_l^\text{to}\right)^2 + \left(Q_l^\text{to}\right)^2 \le \left(S_l^\text{max}\right)^2 v_{A_{l,2}}^2, \quad l \in \mathcal{L}
$$

Equivalent form: $|I_l|^2 \le (S_l^\text{max})^2$, since $|I| = \sqrt{P^2+Q^2}/V$ in p.u.

### 4.3  Transformer thermal limit

$$
\left(P_\tau^\text{hv}\right)^2 + \left(Q_\tau^\text{hv}\right)^2 \le \left(S_\tau^\text{max}\right)^2 v_{A^\tau_{\tau,1}}^2, \quad \tau \in \mathcal{T}
$$

$$
\left(P_\tau^\text{lv}\right)^2 + \left(Q_\tau^\text{lv}\right)^2 \le \left(S_\tau^\text{max}\right)^2 v_{A^\tau_{\tau,2}}^2, \quad \tau \in \mathcal{T}
$$

### 4.4  Reactive power limits

For controllable static generators:

$$
Q_{g_s}^\text{min} \le q_{g_s} \le Q_{g_s}^\text{max}, \quad g_s \in \mathcal{G}_s^c
$$

For external grids and synchronous generators:

$$
Q_g^\text{min} \le q_g \le Q_g^\text{max}, \quad g \in \mathcal{G}
$$

### 4.5  Wind generator Q-curve constraints

For wind generators $w \in \mathcal{G}_s^W$ that provide reactive support (grid-code variant), the reactive power is bounded by **linear piecewise functions of $p_w$** and **of the bus voltage $v_b$**:

**P–Q characteristic:**

$$
q_w \le b_{qp}^\text{max}(v) \cdot P_w^\text{inst} + m_{qp}^\text{max}(v) \cdot p_w
$$

$$
q_w \ge b_{qp}^\text{min}(v) \cdot P_w^\text{inst} + m_{qp}^\text{min}(v) \cdot p_w
$$

**Q–V characteristic:**

$$
q_w \ge \left( m_{qv} \cdot v_b + b_{qv}^\text{min} \right) P_w^\text{inst}
$$

$$
q_w \le \left( m_{qv} \cdot v_b + b_{qv}^\text{max} \right) P_w^\text{inst}
$$

where $P_w^\text{inst}$ is the installed capacity and the slope / intercept parameters $m, b$ encode the grid-code curve.

### 4.6  Objectives

**Voltage deviation minimisation (default)**

$$
\min \sum_{b \in \mathcal{B} \setminus \mathcal{B}_0} (v_b - 1)^2 + \sum_{b \in \mathcal{B}_0} (v_b - v_b^{(0)})^2
$$

**Reactive power minimisation**

$$
\min \sum_{g_s \in \mathcal{G}_s} q_{g_s}^2
$$

---

## 5  Hosting Capacity Analysis (`HC_ACOPF`)

`HC_ACOPF` extends `ACOPF` with **binary placement variables** to maximise the total active power that can be injected by a set of candidate wind sites $\mathcal{G}_s^{HC} \subseteq \mathcal{G}_s$.

### 5.1  Binary placement variable

$$
y_w \in \{0, 1\}, \quad w \in \mathcal{G}_s^{HC}
$$

$y_w = 1$ means generator $w$ is active (connected to the grid).

### 5.2  Apparent power constraints per generator

$$
p_w^2 + q_w^2 \le \left(S_w^\text{max}\right)^2 y_w, \quad w \in \mathcal{G}_s^{HC}
$$

$$
p_w^2 + q_w^2 \ge \left(S_w^\text{min}\right)^2 y_w, \quad w \in \mathcal{G}_s^{HC}
$$

The big-M structure forces $p_w = q_w = 0$ whenever $y_w = 0$.

### 5.3  Grid-code Q-curve (simplified)

A simplified version of the grid-code constraint is used for HC:

$$
-0.41\, p_w \le q_w \le 0.48\, p_w, \quad w \in \mathcal{G}_s^{HC}
$$

Combined with a linear Q–V characteristic identical to Section 4.5.

### 5.4  Objective

**Total hosted power minus network losses (default):**

$$
\max \sum_{w \in \mathcal{G}_s^{HC}} p_w - \sum_{l \in \mathcal{L}} \left(P_l^\text{from} + P_l^\text{to}\right) - \sum_{\tau \in \mathcal{T}} \left(P_\tau^\text{hv} + P_\tau^\text{lv}\right)
$$

**Weighted variant (with loss weighting $\varepsilon \in [0,1]$):**

$$
\max\; \varepsilon \sum_{w \in \mathcal{G}_s^{HC}} p_w + (1-\varepsilon) \left( -\sum_{l \in \mathcal{L}} \left(P_l^\text{from} + P_l^\text{to}\right) \right)
$$

---

## 6  Multi-Period Extension

All single-period models are extended by a discrete time index $t \in \mathcal{T}_\text{sim} = \{t_0, t_0{+}1, \ldots, T{-}1\}$ with step length $\Delta t$ (default 0.25 h = 15 min).  Every decision variable gains a time dimension, e.g. $v_{b,t}$, $p_{g_s,t}$, and every constraint is replicated for each $t$.

Load and generation profiles are provided by SimBench and fixed as parameters $P_d^{(t)}, Q_d^{(t)}$ unless a flexible device explicitly unfixes them.

### 6.1  ACOPF multi-period objective (voltage deviation)

$$
\min \sum_{t \in \mathcal{T}_\text{sim}} \left[ \sum_{b \in \mathcal{B} \setminus \mathcal{B}_0} (v_{b,t} - 1)^2 + \sum_{b \in \mathcal{B}_0} (v_{b,t} - v_b^{(0)})^2 \right]
$$

---

## 7  Flexible Devices

Flexible devices are Pyomo *mixin objects* — they attach their own Sets, Params, Vars, and Constraints to an existing multi-period model instance.

### 7.1  Battery (`Battery_multi_period`)

**Variables**

| Symbol | Pyomo name | Domain | Description |
|---|---|---|---|
| $p_{b,t}^\text{bat}$ | `BAT_P[bat,t]` | $\mathbb{R}$ | Battery power (+: charge, −: discharge) |
| $e_{b,t}$ | `BAT_SOC[bat,t]` | $\mathbb{R}$ | State of charge (p.u. of capacity) |

**Power bound**

$$
P_b^\text{min} \le p_{b,t}^\text{bat} \le P_b^\text{max}, \quad b \in \mathcal{B}^\text{bat},\ t \in \mathcal{T}_\text{sim}
$$

**SOC bound**

$$
E_b^\text{min} \le e_{b,t} \le E_b^\text{max}, \quad t > t_0
$$

**Initial SOC**

$$
e_{b,t_0} = 0.5\, E_b^\text{max}
$$

**SOC dynamics**

$$
e_{b,t} = e_{b,t-1} + \Delta t \cdot \frac{\eta_b\, p_{b,t}^\text{bat}}{C_b}, \quad t > t_0
$$

where $\eta_b$ is the round-trip efficiency and $C_b$ is the capacity in MWh.

### 7.2  Heat Pump (`Heatpump_multi_period`)

**Variables**

| Symbol | Pyomo name | Domain | Description |
|---|---|---|---|
| $p_{h,t}^\text{hp}$ | `hp_p[hp,t]` | $\mathbb{R}$ | Electrical power consumed |
| $\theta_{h,t}$ | `temp[hp,t]` | $\mathbb{R}$ | Indoor temperature (K) |

**Power bound**

$$
P_h^\text{min} \le p_{h,t}^\text{hp} \le P_h^\text{max}, \quad h \in \mathcal{H},\ t \in \mathcal{T}_\text{sim}
$$

**Temperature bound**

$$
\Theta_h^\text{min} \le \theta_{h,t} \le \Theta_h^\text{max}
$$

**Thermal dynamics**

$$
\theta_{h,t} = \theta_{h,t-1} + \Delta t \cdot \frac{\text{COP}_h\, p_{h,t}^\text{hp} - Q_{h,t}^\text{loss}}{C_h^\text{therm}}, \quad t > t_0
$$

where $\text{COP}_h$ is the coefficient of performance, $Q_{h,t}^\text{loss}$ the heat loss at time $t$, and $C_h^\text{therm}$ the thermal capacity of the building.

### 7.3  Photovoltaic (`PV_multi_period`)

PV units have time-varying upper bounds derived from the SimBench irradiance profile:

$$
0 \le p_{g_s,t}^\text{PV} \le P_{g_s,t}^\text{PV,max}, \quad g_s \in \mathcal{G}_s^\text{PV},\ t \in \mathcal{T}_\text{sim}
$$

Curtailment is modelled implicitly: the optimizer may dispatch below the available potential.

### 7.4  Wind Power (`Windpower_multi_period`)

Wind generators support the same grid-code Q-curve as the single-period case (Section 4.5), replicated for each time step $t$.  For hosting-capacity analyses the binary variable $y_{w,t} \in \{0,1\}$ and the apparent-power big-M constraints of Section 5 are also replicated per $t$:

$$
p_{w,t}^2 + q_{w,t}^2 \le \left(S_{w,t}^\text{max}\right)^2 y_{w,t}, \quad w \in \mathcal{G}_s^{HC},\ t \in \mathcal{T}_\text{sim}
$$

---

## 8  Controllable DSO Network Equipment

Opt-in decision variables for the equipment a distribution system operator
operates. Nothing in this section exists on a model until `enable_oltc` or
`enable_shunt_control` is called; without them the formulations above are
unchanged. Derivations, the pandapower verification and the literature are in
`docs/research/dso_controllable_equipment.md`.

### 8.1  On-Load Tap-Changing Transformers

pandapower's longitudinal (`tap_changer_type = "Ratio"`) tap changer scales the
rated voltage of the tapped winding by the **tap factor**

$$
n_\tau(k) = 1 + (k - k_\tau^\text{neutral})\, s_\tau, \qquad s_\tau = \frac{\texttt{tap\_step\_percent}}{100},
$$

and builds the off-nominal ratio $\tau_0 = r_0\, n$ (HV-side tap) or $\tau_0 = r_0 / n$ (LV-side tap), where $r_0 = (v_{n,hv}/v_{n,lv}) / (V^{bus}_{hv}/V^{bus}_{lv})$ is the nominal mismatch. For an LV-side tap it also refers the short-circuit impedance to the tapped LV voltage, so every per-unit admittance of the branch scales with $1/n^2$.

For a controlled transformer $\tau \in \mathcal{T}^\text{OLTC}$ the position becomes a variable and the ratios of Section 1.4 follow from it:

$$
k_\tau \in [k^\text{min}_\tau, k^\text{max}_\tau], \qquad
n_\tau = 1 + (k_\tau - k^\text{neutral}_\tau)\, s_\tau \quad \text{(\texttt{trafo\_tap\_factor\_def})}
$$

$$
\text{HV-side tap:}\quad a^{hv}_\tau = r_{0,\tau}\, n_\tau,\; a^{lv}_\tau = 1 \quad \text{(\texttt{trafo\_tap\_ratio\_hv\_def})}
$$

$$
\text{LV-side tap:}\quad a^{hv}_\tau = \tau_{0,\tau},\; a^{lv}_\tau = \frac{n_\tau}{n_{0,\tau}} \quad \text{(\texttt{trafo\_tap\_ratio\_lv\_def})}
$$

with $n_{0,\tau}$ the tap factor the network was built with. Substituting pandapower's $\tau(k) = r_0/n$ and $Y(k) = Y^{(1)}/n^2$ into the MATPOWER branch model gives $Y_{ff} = Y^{(1)}/r_0^2$ (tap-invariant), $Y_{ft} = -y_s^{(1)}/(r_0 n)$ and $Y_{tt} = Y^{(1)}/n^2$, which is the two-sided form above written with the stored admittances $Y_0 = Y^{(1)}/n_0^2$. At $k = k_0$ both cases reduce to the fixed model; at every other position they reproduce `pp.runpp` to solver tolerance (`tests/unit_tests/test_oltc.py`).

| Pyomo component | index | meaning |
|---|---|---|
| `TRANSF_OLTC`, `TRANSF_OLTC_HV`, `TRANSF_OLTC_LV` | — | controlled transformers, by tapped side |
| `trafo_tap_pos_min/max/neutral/init`, `trafo_tap_step`, `trafo_tap_ratio_nominal`, `trafo_tap_factor_base` | $\tau$ | $k^\text{min}, k^\text{max}, k^\text{neutral}, k_0, s, r_0, n_0$ |
| `trafo_tap_position` | $\tau$ (× $t$) | $k$ — integer or real |
| `trafo_tap_factor` | $\tau$ (× $t$) | $n$ |
| `Tap`, `Tap_lv` | $\tau$ (× $t$) | $a^{hv}$, $a^{lv}$; unfixed for the controlled side |

Eligibility (all required): in service; `tap_side` in {hv, lv}; `tap_changer_type == "Ratio"`; `tap_step_degree` 0 or NaN; no `tap_dependency_table`; no second changer; integer `tap_neutral`, `tap_min < tap_max`, `tap_step_percent != 0`, `tap_min <= tap_pos <= tap_max`; and the ratio pandapower built equals the one the formula predicts.

### 8.2  Continuous vs. Discrete Tap Control

| mode | domain of `trafo_tap_position` | model class | meaning |
|---|---|---|---|
| fixed (default) | no variable | NLP | pandapower's `tap_pos` |
| `continuous` | $\mathbb{R}$ | NLP | the LP/NLP relaxation of the discrete model: a lower bound and a warm start, **not** an implementable position |
| `discrete` | $\mathbb{Z}$ | nonconvex MINLP | the physical changer |

The tap factor is affine in the position, so one general integer per
transformer (and period) is an exact representation; no one-hot selection
and no big-M product is needed because the polar AC equations carry
$1/a^2$ natively. `solve_oltc_round_and_fix` solves the continuous model,
rounds each position to the nearest admissible integer (in time order and
within the per-step change limit), fixes it and re-solves the NLP — a
heuristic whose relaxed objective bounds the result.

### 8.3  Multi-Period Tap Scheduling

With $k_{\tau,0}$ the network's `tap_pos` (or a user reference) and $t^-$
the previous period:

$$
k_{\tau,t} - k_{\tau,t^-} = u_{\tau,t} - d_{\tau,t}, \quad u, d \ge 0 \quad \text{(\texttt{trafo\_tap\_movement\_def})}
$$

$$
u_{\tau,t} + d_{\tau,t} \le \Delta k^\text{max}_\tau \quad \text{(\texttt{trafo\_tap\_change\_limit})}, \qquad
\sum_t \left(u_{\tau,t} + d_{\tau,t}\right) \le N^\text{max}_\tau \quad \text{(\texttt{trafo\_tap\_operations\_limit})}
$$

$$
C^\text{switch} = \sum_\tau c_\tau \sum_t \left(u_{\tau,t} + d_{\tau,t}\right) \quad \text{(\texttt{trafo\_tap\_movement\_cost}, added by \texttt{penalize\_tap\_movement})}
$$

$u + d \ge |k_t - k_{t^-}|$ always holds and is tight whenever the cost is
positive or the operation limit binds; the reported operation count is
computed from the positions themselves. The single-period model carries the
same variables relative to the initial position.

### 8.4  Switchable Reactive Compensation

A pandapower shunt consumes $(p_s + jq_s)\,\text{step}\,\rho_s\, v^2$ with
$\rho_s = (V^{bus}_n/\texttt{vn\_kv}_s)^2$. For a controlled shunt
$s \in \mathcal{S}^\text{ctrl}$ the step becomes a variable:

$$
m_s \in [0, m^\text{max}_s], \qquad
p^{sh}_{s} = \frac{p_s \rho_s}{S_N}\, m_s\, v_b^2, \qquad
q^{sh}_{s} = \frac{q_s \rho_s}{S_N}\, m_s\, v_b^2,
$$

replacing the constant $G_s v_b^2$ / $-B_s v_b^2$ terms in `KCL_real` /
`KCL_reactive` for those shunts only (the balance is rebuilt by
`rebuild_kcl`). Positive $q_s$ is an inductive reactor, negative a capacitor
(pandapower's load sign). `shunt_step` is integer in `discrete` mode and
real in `continuous` mode; `shunt_step_up/down`, `shunt_step_change_limit`,
`shunt_step_operations_limit` and `shunt_switching_cost` mirror Section 8.3.

### 8.5  Solver Requirements

| formulation | class | solvers |
|---|---|---|
| fixed taps and steps | NLP | IPOPT, Gurobi 12+, NEOS |
| continuous tap / continuous step | NLP (one more bilinear term) | same |
| discrete tap and/or step | nonconvex MINLP | `gurobi_direct_minlp` (global), MindtPy (local), NEOS Bonmin/Couenne; `solve_*_round_and_fix` with IPOPT (heuristic) |
| DC OPF | LP | no tap or shunt control (no voltage magnitude) |

`solve()` raises a `ValueError` when free integer variables meet a
continuous-only solver (IPOPT, CONOPT, SNOPT, MINOS, also through NEOS),
unless `relax_integrality=True`.

### 8.6  Pandapower Mapping

| pandapower field | used as |
|---|---|
| `trafo.tap_side` | which ratio is freed (`Tap` for hv, `Tap_lv` for lv) |
| `trafo.tap_pos` | $k_0$: start value, movement reference, consistency check against ppc `TAP` |
| `trafo.tap_neutral`, `tap_step_percent` | $k^\text{neutral}$, $s$ |
| `trafo.tap_min`, `tap_max` | bounds of `trafo_tap_position` |
| `trafo.tap_changer_type` | must be `"Ratio"` |
| `trafo.tap_step_degree`, `tap_dependency_table`, `tap2_*` | must be absent/zero |
| `trafo.vn_hv_kv`, `vn_lv_kv`, `bus.vn_kv` | $r_0$ |
| `trafo.oltc` | **not used** (short-circuit flag) |
| `shunt.p_mw`, `q_mvar`, `vn_kv`, `step`, `max_step` | per-step power, $\rho_s$, start value, bound |
| `shunt.step_dependency_table` | must be False |
| results | `res_trafo.tap_pos`, `res_trafo.tap_factor`, `res_shunt.step`; `apply_tap_positions` / `apply_shunt_steps` write `trafo.tap_pos` / `shunt.step` on request |

### 8.7  Limitations

Cross regulators (`tap_step_degree != 0`), symmetrical and ideal phase
shifters, tabular characteristics, second tap changers and three-winding
transformers are detected and rejected, not approximated. Transformers
without `tap_changer_type` are rejected because pandapower ignores their
tap. No OLTC in DC; no network reconfiguration. The discrete modes
are nonconvex MINLPs: the global solver is practical for single-period
cases, the rounding heuristic is the default for horizons.

---

## 9  Symbol Reference

| Symbol | Pyomo name | Description |
|---|---|---|
| $\mathcal{B}$ | `B` | Set of all buses |
| $\mathcal{B}_0$ | `b0` | Slack / reference buses |
| $\mathcal{B}_{PV}$ | `bPV` | PV buses |
| $\mathcal{G}_s$ | `sG` | Static generators |
| $\mathcal{G}_s^c$ | `sGc` | Controllable static generators |
| $\mathcal{G}_s^{HC}$ | `WIND_HC` | Candidate hosting-capacity sites |
| $\mathcal{G}$ | `G` | External grids + synchronous generators |
| $\mathcal{G}_{ext}$ | `eG` | Slack generators (ref=True) |
| $\mathcal{D}$ | `D` | Loads |
| $\mathcal{D}^c$ | `Dc` | Controllable loads |
| $\mathcal{L}$ | `L` | Lines |
| $\mathcal{T}$ | `TRANSF` | Transformers |
| $\mathcal{S}$ | `SHUNT` | Shunts |
| $\mathcal{T}_\text{sim}$ | `T` | Time-period index |
| $A_{l,1}, A_{l,2}$ | `A[l,1]`, `A[l,2]` | From-bus / to-bus of line $l$ |
| $A^\tau_{\tau,1}, A^\tau_{\tau,2}$ | `AT[τ,1]`, `AT[τ,2]` | HV-bus / LV-bus of transformer $\tau$ |
| $\Delta t$ | `deltaT` | Time-step length in hours |
| $\text{baseMVA}$ | `baseMVA` | System base power (MVA) |
