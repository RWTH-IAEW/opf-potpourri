# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Work out why an OPF did not solve, or why its answer looks wrong.

The entry point is a method on the model:

```python
opf = ACOPF(net)
opf.add_OPF()
opf.add_voltage_deviation_objective()
opf.solve(solver="ipopt")

report = opf.diagnose()
print(report)
```

What separates this from a generic solver diagnostic is that findings are
expressed in terms of the network the user built. Pyomo can say
`line_lim_from[23]`; this package says `net.line[7] "Cable 7"` — and keeps
the Pyomo identifier alongside, so the trail from report to constraint to
formulation stays intact.

`diagnose()` works on a model that has not been solved, one that failed,
and one that succeeded, and it never modifies the model it is given.
"""

from potpourri.diagnostics.context import DiagnosticContext
from potpourri.diagnostics.mappings import ElementRef, IndexMap
from potpourri.diagnostics.report import (
    DiagnosticCategory,
    DiagnosticIssue,
    DiagnosticReport,
    DiagnosticSeverity,
)
from potpourri.diagnostics.runner import LEVELS, diagnose

__all__ = [
    "LEVELS",
    "DiagnosticCategory",
    "DiagnosticContext",
    "DiagnosticIssue",
    "DiagnosticReport",
    "DiagnosticSeverity",
    "ElementRef",
    "IndexMap",
    "diagnose",
]
