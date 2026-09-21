# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Checks that run against an installed potpourri, not the source.

Run from outside the repository so an import resolves to the
installed distribution; the source tree would otherwise shadow
it and the check would prove nothing.
"""
