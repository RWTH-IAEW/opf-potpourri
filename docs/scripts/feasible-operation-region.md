# Feasible Operation Region

**Script:** `scripts/feasible_operation_region.py`

Computes and plots the **Feasible Operation Region (FOR)** — the set of all
active / reactive power (P, Q) combinations that the external grid connection
points of a distribution network can realise while satisfying all network
constraints (voltage limits, line loading limits, generator capability curves).

---

## Concept

A distribution system operator (DSO) connected to a transmission grid needs
to communicate what net power exchange (P, Q) is feasible at their grid
connection point. The FOR is the boundary of that set in the P–Q plane.

```
         Q [MVar]
         ▲
         │   ╭────╮
         │  /      \
    ─────┼─/────────\──── P [MW]
         │ \        /
         │  \──────╯
         │
```

The script traces this boundary by repeatedly solving the AC-OPF in different
directions of the P–Q plane and collecting the extreme points.

---

## Algorithms

Four sampling strategies are provided, from coarse to adaptive-refined:

### `for_setpoint_based` — 8 cardinal directions

Solves 8 OPF problems, one per cardinal direction (±P, ±Q, four diagonals),
by maximising `alpha·P + beta·Q` for each `(alpha, beta)` pair.  Returns the
boundary as a flat list of P and Q values per ext_grid.

```
Directions: (1,0), (1,1), (0,1), (−1,1), (−1,0), (−1,−1), (0,−1), (1,−1)
```

### `for_setpoint_based_with_directions` — adaptive refinement

Extends the 8-direction coarse scan with automatic gap refinement:

```
coarse pass (8 directions)
  → close polygon, compute arc-length steps
  → identify gaps where step > stepsize/delta_max
  → for P-dominated gaps: sweep P linearly, optimise Q freely
  → for Q-dominated gaps: sweep Q linearly, optimise P freely
return coarse + refined boundary points
```

The `stepsize` parameter controls the target normalised arc-length step.
Smaller values produce denser sampling but require more OPF solves.

### `run_feasible_operation_region` — angle sweep

Sweeps 36 power-factor angles θ ∈ [0, 2π) and maximises P² + Q² subject to
the power-factor constraint Q = tan(θ)·P.

### `for_angle_based_sampling` — quadrant-split angle sweep

Splits the angle sweep into four quadrants and uses sign-corrected linear
objectives per quadrant to avoid degenerate solutions near the axes.

---

## Usage

### Basic FOR with 8 directions

```python
import simbench as sb
from potpourri.models.ACOPF_base import ACOPF
from scripts.feasible_operation_region import for_setpoint_based, get_pq_to_plot

net = sb.get_simbench_net("1-LV-rural1--0-sw")
factors = net.loadcases.loc["lW"]
net.load.p_mw   *= factors["pload"]
net.load.q_mvar *= factors["qload"]
net.ext_grid.vm_pu = factors["Slack_vm"]

acopf = ACOPF(net)
acopf.add_OPF()

p, q = for_setpoint_based(acopf)
pq, pq_tot = get_pq_to_plot(p, q)
print(f"P range: {pq[0][:,0].min():.4f} … {pq[0][:,0].max():.4f} pu")
print(f"Q range: {pq[0][:,1].min():.4f} … {pq[0][:,1].max():.4f} pu")
```

### Adaptive FOR with refinement

```python
from scripts.feasible_operation_region import for_setpoint_based_with_directions

p, q, v, nets = for_setpoint_based_with_directions(acopf, stepsize=0.5, solver="ipopt")
pq, pq_tot = get_pq_to_plot(p, q)
print(f"Total boundary points: {pq_tot.shape[0]}")
```

### Run the main script

```bash
conda activate potpourri_env
python scripts/feasible_operation_region.py
```

The `__main__` block operates on the HV-mixed simbench network and requires
the hosting-capacity columns (`wind_hc`, `var_q`) to be set up via
`net_augmentation/prepare_net.py` first.

---

## Example output (LV rural1, load case `lW`)

```
=== for_setpoint_based (8 cardinal directions) ===
  Boundary points: 8
  P range [pu]: [-0.11896, -0.11123]
  Q range [pu]: [-0.12576,  0.12697]

=== for_setpoint_based_with_directions (stepsize=0.5) ===
  Total boundary points after refinement: 9
  P range [pu]: [-0.11896, -0.11123]
  Q range [pu]: [-0.12576,  0.12697]
```

The P range is very narrow (≈ 0.008 pu) because the LV rural network's active
power import is nearly fixed by the loads and the small PV generation.  The Q
range is wider (≈ 0.25 pu), reflecting the reactive flexibility available from
the PV generators and the ext_grid.

---

## Plotting

The `plot_hull` function visualises the boundary points and their convex or
concave hull per ext_grid. It requires `shapely` and `geopandas`, which are
optional dependencies (uncomment the relevant imports at the top of the
script):

```python
# in feasible_operation_region.py — uncomment:
from shapely import concave_hull, MultiPoint, convex_hull
import geopandas as gpd

from scripts.feasible_operation_region import plot_hull
figs, axs = plot_hull(p, q, ratio=0.1)
figs[0].savefig("for_total.pdf", bbox_inches="tight")
```

---

## API reference

| Function | Description |
|---|---|
| `for_setpoint_based(opf)` | 8-direction boundary scan; returns `(p, q)` flat lists |
| `for_setpoint_based_with_directions(opf, stepsize, solver)` | Adaptive refinement; returns `(p, q, v, nets)` nested lists |
| `for_angle_based_sampling(opf, n)` | Angle sweep with quadrant correction; returns `(p, q, v)` |
| `run_feasible_operation_region(opf)` | 36-angle sweep with P²+Q² objective; returns `(p, q)` |
| `node_for_setpoint_based_with_directions(opf, w, stepsize)` | Per-generator node FOR; returns `(p, q, v, nets)` |
| `get_pq_to_plot(p, q)` | Flatten boundary lists to Nx2 arrays; returns `(pq, pq_tot)` |
| `plot_hull(p, q, ratio)` | Plot boundary with convex/concave hull; returns `(figs, axs)` |

---

## Notes

- All P and Q values in Pyomo are in **per-unit** (divided by `baseMVA`).
  Multiply by `acopf.net.sn_mva` to convert to MW / MVar.
- The `for_setpoint_based_with_directions` function modifies the Pyomo model
  in-place (adds `p_sp`, `q_sp` parameters and constraints). Pass a **fresh
  model** for each call or reconstruct the `ACOPF` object between runs.
- `nets` in the return value of the `*_with_directions` functions contains
  deep-copied network states at each coarse corner point, allowing
  post-hoc inspection of line loadings and voltage profiles.
