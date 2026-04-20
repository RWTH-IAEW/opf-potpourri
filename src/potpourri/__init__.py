"""POTPOURRI: multi-period Optimal Power Flow tool built on Pyomo and pandapower."""

__all__ = [""]
__version__ = "0.1.0"

from loguru import logger

# Suppress loguru output from this library by default.
# Users can enable it with: from loguru import logger; logger.enable("potpourri")
logger.disable("potpourri")
