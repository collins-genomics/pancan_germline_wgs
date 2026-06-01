#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2026-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Math functions
"""


from numpy import nanmin, nanmax


def is_inside(x, interval, closed_bounds=True):
	"""
	Checks whether a numeric value, `x`, is inside a numeric iterable, `interval`
	"""

	x = float(x)

	if closed_bounds:
		return x >= nanmin(interval) and x <= nanmax(interval)
	else:
		return x > nanmin(interval) and x < nanmax(interval)
