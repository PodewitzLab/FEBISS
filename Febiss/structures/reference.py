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

from Febiss.structures.mol2_to_xyz import converter
from Febiss.structures.write_xyz import write_xyz as write
from Febiss.structures.quat_handling import *
from pymatgen.core import Molecule
from pymatgen.symmetry import analyzer as ana
import quaternion as quat


class Reference:
    """
    This class contains all necessary information on the used solvent and has methods
    to deal with quaternions and equivalent structures.
    """
    def __init__(self, com, solv_abb, solv_file,
                 rigid_atom_0 : int = 0, rigid_atom_1 : int = 1, rigid_atom_2 : int = 2):

        self.verbose = False # additional output for debugging
        self.com = com
        self.solv_abb = solv_abb
        self.solv_file = solv_file

        if self.solv_file.split(".")[-1] == "xyz":
            self.xyz_path = self.solv_file
        elif self.solv_file.split(".")[-1] == "mol2":
            self.xyz_path = converter(os.path.dirname(self.solv_file), solv_abb) #returns path of generated xyz-file.
        else:
            quit("Only xyz and mol2 files can be given as solvent molecules!")

        self.refdir_path = None #Gets defined in _process_equivalent_structures(). Contains path to parent directory.
        self.mol = Molecule.from_file(self.xyz_path) #pymatgen interface
        self.cmol = self.mol.get_centered_molecule()

        self.rigid_atom_idx_0 = rigid_atom_0
        self.rigid_atom_idx_1 = rigid_atom_1
        self.rigid_atom_idx_2 = rigid_atom_2

        if self.com:
            print("Using rigidatoms {0} and {1} for quaternion determination!".format(self.rigid_atom_idx_1,
                                                                                      self.rigid_atom_idx_2))
        else:
            print("Using rigidatoms {0}, {1} and {2} for quaternion determination!".format(self.rigid_atom_idx_0,
                                                                                           self.rigid_atom_idx_1,
                                                                                           self.rigid_atom_idx_2))

        self.char_q = calc_quats(self.cmol, self.rigid_atom_idx_0, self.rigid_atom_idx_1, self.rigid_atom_idx_2)

        #placement
        # stores...
        # 1) the rotation quat,
        # 2) the equivalent structure as Molecule object and
        # 3) the characteristic quat for the orientation of the molecule as tuple per symm_op
        self.eq_dict = {}
        self._process_equivalent_structures(True)

        #post placement
        self.elements = []
        self.atoms = []
        self.values = []
        self.all_values = []

    def _makedir(self, path: str) -> str:
        """used primarily to create a folder that contains the reference structure
        of pyConSolv solvents"""
        if self.verbose: # creates a folder with the reference structures everytime a FEBISS analysis is conducted
            import datetime
            path = path + "_{0}".format(str(datetime.date.today()))
            if not os.path.exists(path):
                os.mkdir(path)
            else:
                num = 1
                while os.path.exists(path + "_{0}".format(str(num))):
                    num += 1
                path = path + "_{0}".format(str(num))
                os.mkdir(path)
            return os.path.abspath(path)

        else:
            if not os.path.exists(path):
                os.mkdir(path)

            else:
                shutil.rmtree(path)
                os.mkdir(path)

            return os.path.abspath(path)

    def _process_equivalent_structures(self, write_out: bool = True):

        """
        analyzes the given solvent, gives coordinates for equivalent structures and sets the following:
        1) self.eq_dict
        """

        # directory containing the xyz files of equivalent structures
        self.refdir_path = self._makedir("REF_{0}".format(self.solv_abb))


        #last placeholder is for enumeration of files used in write()
        path_template = self.refdir_path+"/{0}".format(self.solv_abb)+"_{0}.xyz"

        #building up self.equivalent_structures and self.symmetry_rots_as_quats
        pga = ana.PointGroupAnalyzer(self.cmol)
        symm_ops = pga.get_symmetry_operations()

        num = 0
        for symm in symm_ops:
            coord_list = []
            rot = symm.rotation_matrix
            return_path = path_template.format(num) #enumerates the xyz-files
            # adopted from analyzer.get_rotational_symmetry_number() (line 1284).
            # Filters only for those transformations that have adeterminant > 1 which means no mirror.
            if np.abs(np.linalg.det(rot) - 1) < 1e-4:
                rot_quat = quat.from_rotation_matrix(rot) #transform the rotation matrix into quaternion
            else:
                continue #rotation with mirror has been found. skipping this symm_op.

            for idx in range(len(self.mol.cart_coords)):
                coord_list.append(symm.apply_rotation_only(self.cmol.cart_coords[idx]))
            write(self.xyz_path, return_path, coord_list)
            ref = Molecule.from_file(return_path.format(num))

            char_quat = calc_quats(ref, self.rigid_atom_idx_0, self.rigid_atom_idx_1, self.rigid_atom_idx_2)


            if distance(self.char_q,char_quat) < 0.05: #0.05 rad is around 3.18°
                os.rename(return_path,self.refdir_path+"/{0}".format(self.solv_abb)+"_ori.xyz")

                # also for the identity there is an entry in eq_dict. This facilitates the quaternion cleanup
                # in _find_avg_solvent.
                self.eq_dict[num] = (rot_quat, ref, char_quat)

            else:
                self.eq_dict[num] = (rot_quat, ref, char_quat)

            num += 1

    def _find_avg_solvent(self, voxel: int, quats: list[quat.quaternion],
                          grid_point: tuple[typing.Any,typing.Any,typing.Any]):
        """
        This method finally:
        1) cleans up characteristic quaternions (quats) stored for the given voxel (voxel) by comparing them to the
        characteristic quaternions of the equivalent structures stored in self.eq_dict.
        2) calculates the average quaternion from the cleaned up list
        3) rotates the original orientation to fit the average quaternion
        4) translates the coordinates of the rotated solvent to the COM (com)

        Returns the element names of the solvent atoms, their coordinates and the associated energy value.

        :param voxel: the voxel under consideration
        :param quats: the quaternions associated with the voxel under consideration
        :param com: the coordinates where the center of mass shall be placed.
        :return: elements, coords, values
        """
        #step 1)
        key_list = list(self.eq_dict.keys())
        for i in range(len(quats)):
            distance_list = []
            for j in range(len(key_list)):
                distance_list.append(distance(quats[i], self.eq_dict[key_list[j]][2]))

            if self.verbose:
                print("\n\nquat: {0}".format(i))
                print("distance_list: {0}".format(distance_list))

            min_idx = distance_list.index(min(distance_list))
            if self.verbose:
                print("min_idx: {0}".format(min_idx))
                print("quats[{1}] before clean-up: {0}".format(quats[i],i))

            quats[i] = quats[i] * inv(self.char_q) * inv(self.eq_dict[key_list[min_idx]][0]) * self.char_q

            if self.verbose:
                print("quats[{1}] after clean-up: {0}".format(quats[i],i))
                print("distance now: {0}".format(distance(quats[i], self.eq_dict[key_list[min_idx]][2])))

        #step 2)
        q_avg = quat.from_float_array(avg(create_Q_matrix(quats)))

        #step 3) and 4)
        qt = q_avg * inv(self.char_q) #quaternion for rotation of original orientation to orientation described by q_avg
        new_coords = new_coord_gen(self.cmol, qt, np.array(grid_point))

        if not self.com: #in the nocom case, we translate the com by the coordinate of the rigid_atom_0
            if self.verbose:
                print('new coords before = {0}'.format(new_coords))
            new_coords = list(np.array(new_coords) + (np.array(grid_point)-np.array(new_coords[0])))
            if self.verbose:
                print('new coords after = {0}'.format(new_coords))
                print('grid point = {0}'.format(grid_point))


        if self.verbose:
            write(self.xyz_path, self._makedir("quats")+"/avg_at_voxel_{0}", new_coords, voxel)

        return self.cmol.labels, new_coords