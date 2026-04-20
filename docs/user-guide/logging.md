# Logging

POTPOURRI uses [loguru](https://loguru.readthedocs.io/) for structured logging. Following the library best-practice, **all output is suppressed by default** so that the library never pollutes your application's log stream.

## Enabling logging

```python
from loguru import logger
logger.enable("potpourri")
```

Call this once at the start of your script. All subsequent POTPOURRI calls will emit log records to loguru's default sink (stderr with colours).

## Log levels

| Level | What is logged |
|-------|---------------|
| `DEBUG` | Low-level model internals: solver options, variable fix/unfix operations, result mapping to `net.res_*` |
| `INFO` | Normal milestones: model creation timestamp, solve started (with solver name), optimal solution found, NEOS job submitted |
| `WARNING` | Non-fatal issues that may affect results: generator voltage limits overriding bus limits, missing `p_inst_mw` attribute (fallback used), non-optimal solver termination, unknown model component in `change_vals` / `fix_vars` |
| `ERROR` | Caught exceptions: solver errors (mindtpy `ValueError`), termination-condition check failures, `fix_vars` / `unfix_vars` attribute errors |

## Filtering by level

Show only `INFO` and above (suppress `DEBUG`):

```python
import sys
from loguru import logger

logger.remove()                          # remove the default sink
logger.add(sys.stderr, level="INFO")     # add it back at INFO
logger.enable("potpourri")
```

Show only `WARNING` and above (quiet normal operation):

```python
logger.remove()
logger.add(sys.stderr, level="WARNING")
logger.enable("potpourri")
```

Show everything including `DEBUG`:

```python
logger.remove()
logger.add(sys.stderr, level="DEBUG")
logger.enable("potpourri")
```

## Writing logs to a file

```python
logger.add("potpourri.log", level="DEBUG", rotation="10 MB")
logger.enable("potpourri")
```

loguru rotates the file automatically once it reaches 10 MB. See the [loguru docs](https://loguru.readthedocs.io/en/stable/api/logger.html#loguru._logger.Logger.add) for retention and compression options.

## Disabling logging again

```python
logger.disable("potpourri")
```
