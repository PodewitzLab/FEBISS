#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright University Innsbruck, Institute for General, Inorganic, and Theoretical Chemistry, Podewitz Group
See LICENSE for details
"""

def distance_squared(a_array, b_array) -> float:
	return sum(((a-b)**2 for a, b in zip(a_array, b_array)))
