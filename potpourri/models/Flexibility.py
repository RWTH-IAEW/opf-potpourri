import pandas as pd
from pyomo.environ import *
from math import pi
import copy
import numpy as np
import pandapower as pp

from potpourri.models.pyo_to_net import pyo_sol_to_net_res
class Flexibility:
    def __init__(self, net):
        self.net = net

