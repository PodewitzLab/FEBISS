#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright Technische Universität Wien, Institute of Materials Chemistry, Podewitz Group
See LICENSE for details
"""

from febiss.structures.quat_handling import *


class Solute:
    def __init__(self):

        # holds 11-tuples for every solute atom in solute.pdb:
        # ATOM, ordinal number, atomlabel, residuelabel, 1, x, y, z, 1.00, energy, elementlabel
        self.data = []

        self.elements = []  # contains element names
        self.coords = None  #contains xyz coordinates

    def get_coord_set(self):
        coord_list = []
        for atom in self.data:
            coord_list.append((float(atom[-6]),float(atom[-5]),float(atom[-4]))) #TODO: Prone to ValueError
        self.coords = np.asarray(coord_list)

    def get_elements(self): #TODO: merge with get_coord_set since the loop is the same
        for atom in self.data:
            self.elements.append(atom[-1])
