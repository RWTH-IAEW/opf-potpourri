# Utilities API

---

## Net Augmentation

```
potpourri.net_augmentation.prepare_net
```

Utility functions for pre-processing pandapower networks before building OPF models.

### `apply_loadcase_to_sb_net(net, case)`

Scales load and generation according to a named SimBench load case.

**Parameters:**

- `net` — pandapower network with `net.loadcases` DataFrame
- `case` — load case name (e.g. `'lW'` for low wind, `'hL'` for high load)

Applies scaling factors from `net.loadcases.loc[case]` to the following columns: `pload`, `qload`, `Wind_p`, `PV_p`, `RES_p`, `Slack_vm`.

**Example:**

```python
from potpourri.net_augmentation.prepare_net import apply_loadcase_to_sb_net
apply_loadcase_to_sb_net(net, case='lW')
```

### `add_regulatory_q_control_to_wind(net, variant)`

Marks wind generators as controllable and assigns a Q-V characteristic variant for reactive power control.

**Parameters:**

- `net` — pandapower network
- `variant` — Q-V curve variant (integer 0–2):
    - `0` — constant Q, no voltage dependency
    - `1` — Q-P droop (grid code variant 1)
    - `2` — Q-U droop (grid code variant 2)

Sets `sgen.controllable = True` and `sgen.var_q = variant` for all sgens of type `'Wind'`. Also copies `p_mw` to `p_inst_mw` for installed capacity tracking.

**Example:**

```python
from potpourri.net_augmentation.prepare_net import add_regulatory_q_control_to_wind
add_regulatory_q_control_to_wind(net, variant=1)
```

### `upgrade_pandapower_net(old_net)`

Migrates a pandapower network saved with an older version to the current pandapower format.

**Parameters:**

- `old_net` — network object from an older pandapower version

**Returns:** A fresh pandapower network with all DataFrames copied from `old_net`.

**Example:**

```python
from potpourri.net_augmentation.prepare_net import upgrade_pandapower_net
net = upgrade_pandapower_net(old_net)
```

---

## Plotting

```
potpourri.plotting.plot_functions
```

Visualisation utilities for OPF results and grid code characteristics. All plotting functions use Plotly or Matplotlib depending on the function.

### `set_plt_config()`

Applies standard IAEW figure configuration (text width, colour scheme, font settings) to Matplotlib.

### `plot_wind_hc_results(nets, offset=None)`

Plots wind generator locations and their optimised output across multiple HC_ACOPF solutions.

**Parameters:**

- `nets` — list of pandapower networks with HC_ACOPF results
- `offset` — optional coordinate offset for label placement

### `plot_pq_gridcodes()`

Plots Q-P grid code envelopes for all three reactive power control variants. Returns a Plotly figure.

### `plot_pq_res(nets, labels=None)`

Overlays actual Q-P operating points from `nets` onto the grid code envelope plot.

**Parameters:**

- `nets` — list of pandapower networks with `res_sgen` populated
- `labels` — optional list of legend labels

### `plot_qu_gridcodes()`

Plots Q-U grid code envelopes (reactive power vs. voltage magnitude).

### `plot_qu_res(nets, labels=None)`

Overlays actual Q-U operating points from solved networks onto the Q-U envelope.

### `plot_all_pG(hcs)`

Plots a bar chart of wind generator real power output across a list of HC_ACOPF model instances.

**Parameters:**

- `hcs` — list of solved `HC_ACOPF` model instances
