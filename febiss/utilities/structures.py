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
        self.elements = [] #contains element names
        self.atoms = [] #contains xyz coordinates
        self.polars = [] #contain polar atoms of solute determined with determine_polar_... method.
        self.values = [] #contain the temperature factor. TODO:verify maybe it is the energy value

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


class Solvent:
    def __init__(self, abb : str = "WAT", solv_file : bool = False, rigid_atom_0 : int = 0, rigid_atom_1 : int = 1, rigid_atom_2 : int = 2): #before 25 September 2023: __init__(self, top, abb, size, rigid_atom_0, rigid_atom_1, rigid_atom_2):
        #self.top = top #can be None if using water
        self.abb = abb #can be None if using water
        self.solv_file = solv_file
        self.rigid_atom_0 = rigid_atom_0 #O if using water
        self.rigid_atom_1 = rigid_atom_1 #H if using water
        self.rigid_atom_2 = rigid_atom_2 #H if using water

        #placement
        self.equivalent_structures = [] #new. stores the coordinates of the symmetry equivalent structures. LM20231005
        self.reference_quats = []

        #post placement. LM20231005
        self.elements = [] #contains element names
        self.atoms = [] #contains xyz coordinates
        self.values = [] #contains the
        self.all_values = [] #contain the temperature factor for each atom 3 times. TODO:verify. maybe it is the energy value

    def find_equivalent_structures(self, file):
        pass

    def calc_quats(self):
        pass

    def compare_quats(self):
        pass

    def rotate_quats(self):
        pass

    def average_quats(self):
        pass

    def rotate_coords(self):
        pass



    def sort_by_value(self):
        self.elements = [x for _, x in
                         sorted(zip(self.all_values, self.elements), key=lambda pair: pair[0], reverse=True)]
        self.atoms = [x for _, x in sorted(zip(self.all_values, self.atoms), key=lambda pair: pair[0], reverse=True)]
        self.values.sort(reverse=True)
        self.all_values.sort(reverse=True)
