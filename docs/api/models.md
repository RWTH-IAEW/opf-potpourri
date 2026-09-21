# Single-Period Models API

Single-period model classes, in `src/potpourri/models/`.

The stack layers a network reader, a power-flow formulation and
an OPF layer: `Basemodel` -> (`AC` | `DC`) -> `OPF` ->
`ACOPF_base` / `DCOPF`. Start from `ACOPF_base` for a solvable
AC OPF; the formulation modules are the building blocks it uses.

The equations behind these classes are derived in the
[Mathematical Modelling](../mathematical-modelling.md) guide.

---

::: potpourri.models.basemodel

---

::: potpourri.models.AC

---

::: potpourri.models.DC

---

::: potpourri.models.OPF

---

::: potpourri.models.ACOPF_base

---

::: potpourri.models.DCOPF

---

::: potpourri.models.HC_ACOPF

---

::: potpourri.models.cost_objective

---

::: potpourri.models.pyo_to_net

---

::: potpourri.models.init_pyo_from_pp_res
