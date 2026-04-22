"""potpourri — multi-period Optimal Power Flow for distribution grids.

Built on Pyomo and pandapower.  Key entry points::

    from potpourri.models.ACOPF_base import ACOPF
    from potpourri.models_multi_period.ACOPF_multi_period import (
        ACOPF_multi_period,
    )
"""

__all__ = ["__version__"]
__version__ = "0.2.0"

from loguru import logger

# Suppress loguru output from this library by default.
# Users can enable it with:
#   from loguru import logger; logger.enable("potpourri")
logger.disable("potpourri")
