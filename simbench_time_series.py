#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Small description of simbench_time_series

Copyright (c) by Institute for High Voltage Equipment and Grids,
Digitalization and Energy Economics (IAEW), RWTH Aachen University,
10.05.2024, s.kortmann. All rights reserved.
"""

import pandapower as pp
import simbench as sb

sb_code = "1-MV-rural--0-sw"
net = sb.get_simbench_net(sb_code)


