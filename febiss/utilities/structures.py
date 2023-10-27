#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright University Innsbruck, Institute for General, Inorganic, and Theoretical Chemistry, Podewitz Group
See LICENSE for details
"""

import numpy as np
from utilities.mol2_to_xyz import converter
from utilities.quaternion import gist_quat,gigist_quat
from pymatgen.core import Molecule
from pymatgen.symmetry import analyzer as ana
import quaternion as quat #installed on 06 Sept 2023 in febiss_pymatgen env via "conda install -c conda-forge quaternion" according to https://quaternion.readthedocs.io/en/latest/
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
    def __init__(self, solv_file, abb : str = "WAT", rigid_atom_0 : int = 0, rigid_atom_1 : int = 1, rigid_atom_2 : int = 2): #before 25 September 2023: __init__(self, top, abb, size, rigid_atom_0, rigid_atom_1, rigid_atom_2):
        #self.top = top #can be None if using water
        self.solv_file = solv_file
        self.abb = abb  # can be None if using water
        self.mol = Molecule.from_file(converter(self.solv_file, self.abb)) #pymatgen interface. converter is from mol2_to_xyz and returns path of generated xyz-file. TODO: take care of water!
        self.rigid_atom_idx_0 = rigid_atom_0 #1 = O if using water
        self.rigid_atom_idx_1 = rigid_atom_1 #2 = H if using water
        self.rigid_atom_idx_2 = rigid_atom_2 #3 = H if using water


        #placement
        self._process_equivalent_structures(True)
        self.equivalent_structures = [] #new. stores the coordinates of the symmetry equivalent structures as dict. LM20231005
        self.symmetry_rots_as_quats = [] #new. stores the quaternions of the rotation matrices that transform the equivalent structure into the initial structure. LM20231006
        self.reference_quats = [] #new. stores the quaternions associated to the equivalent structures calculated from the same 3 rigid atoms LM20231006

        #post placement. LM20231005
        self.elements = [] #contains element names
        self.atoms = [] #contains xyz coordinates
        self.values = [] #contains the
        self.all_values = [] #contain the temperature factor for each atom 3 times. TODO:verify. maybe it is the energy value
        self.keep = False #new: cleans up xyz-files. not implemented yet as of 20231006

    def _process_equivalent_structures(self,write: bool = False):

        """
        analyzes the given solvent, gives coordinates for equivalent structures and sets:
        1) self.equivalent_structures
        2) self.symmetry_rots_as_quats
        3) self.reference_quats
        """

        pga = ana.PointGroupAnalyzer(self.mol)
        symm_ops = pga.get_symmetry_operations()
        num = 0
        coord_after_rot_dict = {}
        for symm in symm_ops:
            self.symmetry_rots_as_quats.append(quat.from_rotation_matrix(symm.rotation_matrix)) #transform the rotation matrix into quaternion using quaternion package. LM20231006
            coords_after_rot = []
            for idx in range(len(self.mol.cart_coords())):
                coords_after_rot.append(symm.apply_rotation_only(self.mol.cart_coords[idx]))
            coord_after_rot_dict[num] = coords_after_rot
            num += 1
        self.equivalent_structures = coord_after_rot_dict

        for i in range(len(coord_after_rot_dict.keys())):



    def _calc_quats(self):
        X = self.mol.cart_coords[self.rigid_atom_idx_1] - self.mol.cart_coords[self.rigid_atom_idx_0]
        V2 = self.mol.cart_coords[self.rigid_atom_idx_2] - self.mol.cart_coords[self.rigid_atom_idx_0]
        return quat.from_float_array(gigist_quat(X, V2)) #interface to quaternion package

    def _compare_quats(self):
        pass

    def _rotate_quats(self):
        pass

    def _average_quats(self):
        pass

    def _rotate_coords(self):
        pass



    def sort_by_value(self):
        self.elements = [x for _, x in
                         sorted(zip(self.all_values, self.elements), key=lambda pair: pair[0], reverse=True)]
        self.atoms = [x for _, x in sorted(zip(self.all_values, self.atoms), key=lambda pair: pair[0], reverse=True)]
        self.values.sort(reverse=True)
        self.all_values.sort(reverse=True)
