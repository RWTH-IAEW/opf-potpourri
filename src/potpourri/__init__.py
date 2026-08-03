"""potpourri — multi-period Optimal Power Flow for distribution grids.

Built on Pyomo and pandapower.  Key entry points::

    from potpourri.models.ACOPF_base import ACOPF
    from potpourri.models_multi_period.ACOPF_multi_period import (
        ACOPF_multi_period,
    )
"""

from importlib.metadata import PackageNotFoundError, version

from loguru import logger

__all__ = ["__version__"]

# Read the version from the installed distribution metadata rather than
# hard-coding it, so `potpourri.__version__` cannot drift from pyproject.toml.
try:
    __version__ = version("opf-potpourri")
except PackageNotFoundError:  # running from a source tree without install
    __version__ = "0.0.0.dev0"

# Suppress loguru output from this library by default.
# Users can enable it with:
#   from loguru import logger; logger.enable("potpourri")
logger.disable("potpourri")
