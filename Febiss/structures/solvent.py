#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright Technische Universität Wien, Institute of Materials Chemistry, Podewitz Group
See LICENSE for details
"""

import os
import shutil
import typing

import numpy

from Febiss.structures.mol2_to_xyz import converter
from Febiss.structures.write_xyz import write_xyz as write
from Febiss.structures.quat_handling import *
from pymatgen.core import Molecule
from pymatgen.symmetry import analyzer as ana
import quaternion as quat


class Solvent:
    """
    This class serves two purposes:
    1) A Solvent object is created which gets all the data from the datafile created by Analysis_Febiss.h/cpp.
       This object makes use of the attributes data, coords (for the COM coords) and values and the methods sort_by_energy, get_coord_set and get_energy.
    2) Solvent objects are also created for each solvent selected from the interactive barplot.
       These objects make use of the attribute coords (for the COM coord),
       quats (which are loaded from the gist-quats.dat file), elements and elem_coords (for the solvent atom coords).
    """
    def __init__(self, nvoxels = 0, ref_evv = 0):
        self.data = [] #data (i.e. voxel x y z energy) is passed as 5-tuple
        self.coords = [] #holds coords of COMs
        self.values = [] #holds energy values of solvent molecules. needed for plotting
        self.elements = [] #holds the elements of the solvent to be placed
        self.quats = self._prep_dict(nvoxels) #holds all the quats from gist-quats.dat associated with. keys are voxel numbers
        self.elem_coords = []
        self.ref_evv = ref_evv

    def sort_by_energy(self):
        self.data = sorted(self.data, key=lambda tpl: tpl[-1],reverse=True)
        with open('sorted_data.dat','w') as f:
            for i in self.data:
                f.write("{0} {1} {2} {3} {4}\n".format(i[0], i[1], i[2], i[3], i[4]))

    def get_coord_set(self):
        coord_list = []
        for com in self.data:
            coord_list.append((float(com[1]),float(com[2]),float(com[3])))
        self.coords = np.asarray(coord_list)

    def get_energy(self):
        for com in self.data:
            self.values.append(float(com[-1])+self.ref_evv)

    def _prep_dict(self, nvoxels):
        quat_dict = {}
        if nvoxels != 0:
            for voxel in range(nvoxels):
                quat_dict[voxel] = []
        return quat_dict
