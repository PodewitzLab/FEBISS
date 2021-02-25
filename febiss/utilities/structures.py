#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright University Innsbruck, Institute for General, Inorganic, and Theoretical Chemistry, Podewitz Group
See LICENSE for details
"""

import numpy as np

from .distance_functions import distance_squared


class Solute:
    def __init__(self, polar_cutoff: float = 1.1):
        # if solute H further away from any non-C and non-H atom than this cutoff, it is apolar
        self.polar_cutoff = polar_cutoff
        self.elements = []
        self.atoms = []
        self.polars = []
        self.values = []

    def determine_polar_hydrogen_and_non_hydrogen(self):
        polars = []  # includes all solute non-H atoms and H atoms not bond to C or H
        for i, elei in enumerate(self.elements):
            polar = True
            if elei == 'H':
                polar = False
                for j, elej in enumerate(self.elements):
                    if elej not in ['C', 'H'] and \
                            distance_squared(self.atoms[i], self.atoms[j]) < self.polar_cutoff ** 2:
                        polar = True
                        break
            if polar:
                polars.append(self.atoms[i])
        self.polars = np.asarray(polars)


class Water:
    def __init__(self):
        self.elements = []
        self.atoms = []
        self.values = []
        self.all_values = []

    def sort_by_value(self):
        self.elements = [x for _, x in
                         sorted(zip(self.all_values, self.elements), key=lambda pair: pair[0], reverse=True)]
        self.atoms = [x for _, x in sorted(zip(self.all_values, self.atoms), key=lambda pair: pair[0], reverse=True)]
        self.values.sort(reverse=True)
        self.all_values.sort(reverse=True)
