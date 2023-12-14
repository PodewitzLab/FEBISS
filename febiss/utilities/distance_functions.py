#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright University Innsbruck, Institute for General, Inorganic, and Theoretical Chemistry, Podewitz Group
See LICENSE for details
"""

import numpy as np

def distance_squared(a_array, b_array) -> float:
	return sum(((a-b)**2 for a, b in zip(a_array, b_array)))

def distance_between_quaternions(q1 : np.ndarray, q2 : np.ndarray): #maybe use quaternion python package
	NotImplementedError()