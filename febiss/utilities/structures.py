#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright University Innsbruck, Institute for General, Inorganic, and Theoretical Chemistry, Podewitz Group
See LICENSE for details
"""

import os
import numpy as np
import typing
from ..utilities.mol2_to_xyz import converter
from ..utilities.structure_randomizer import write_xyz as write
from ..utilities.quat_handling import *
from pymatgen.core import Molecule
from pymatgen.symmetry import analyzer as ana
import quaternion as quat #installed on 06 Sept 2023 in febiss_pymatgen env via "conda install -c conda-forge quaternion" according to https://quaternion.readthedocs.io/en/latest/
from .distance_functions import distance_squared


class Solute:
    def __init__(self): #, polar_cutoff: float = 1.1): #not needed LM20231123
        # if solute H further away from any non-C and non-H atom than this cutoff, it is apolar
        #self.polar_cutoff = polar_cutoff #UNCLEAR: what is this polar_cutoff value based on? how to deal with inductive/mesomeric effects? LM20231027
        self.data = [] #new LM20231123: holds 11-tuples for every solute atom in solute.pdb: ATOM, ordinal number, atomlabel, residuelabel, 1, x, y, z, 1.00, energy, elementlabel.
        self.elements = []  # contains element names
        self.coords = None #renamed atoms to coords LM20231123.  #contains xyz coordinates
        #self.polars = [] #not needed since  LM20231123 #contain polar atoms of solute determined with determine_polar_... method.
        #self.values = [] #not needed LM20231123 #contains energy

    def get_coord_set(self):
        coord_list = []
        for atom in self.data:
            coord_list.append((float(atom[-6]),float(atom[-5]),float(atom[-4]))) #TODO: Prone to ValueError LM20231127
        self.coords = np.asarray(coord_list)

    def get_elements(self): #TODO: merge with get_coord_set since the loop is the same
        for atom in self.data:
            self.elements.append(atom[-1])




    # not needed right now
    # def determine_polar_hydrogen_and_non_hydrogen(self):
    #     polars = []  # includes all solute non-H atoms and H atoms not bond to C or H
    #     for i, elei in enumerate(self.elements):
    #         polar = True
    #         if elei == 'H':
    #             polar = False
    #             for j, elej in enumerate(self.elements):
    #                 if elej not in ['C', 'H'] and \
    #                         distance_squared(self.atoms[i], self.atoms[j]) < self.polar_cutoff ** 2:
    #                     polar = True
    #                     break
    #         if polar:
    #             polars.append(self.atoms[i])
    #     self.polars = np.asarray(polars)


class Reference:
    """
    This class contains all necessary information on the used solvent and has methods
    to deal with quaternions and equivalent structures.
    """
    def __init__(self, case, com, solv_file = None, abb : str = "WAT", rigid_atom_0 : int = 0, rigid_atom_1 : int = 1, rigid_atom_2 : int = 2): #3 rigid atoms again as of 28 Nov 2023 to account for case 1. only 2 rigid_atoms LM20231113. #before 25 September 2023: __init__(self, top, abb, size, rigid_atom_0, rigid_atom_1, rigid_atom_2):
        #self.top = top #can be None if using water
        self.solv_file = solv_file #contains path to mol2-file, new LM20231128: can also contain path to TP3.xyz file
        self.abb = abb  # can be None if using water

        if case == 1: #for case 1 we use TP3.xyz directly. No need for converter. LM20231128.
            self.xyz_path = solv_file
        else:
            self.xyz_path = converter(os.path.dirname(self.solv_file), self.abb) #converter is from mol2_to_xyz and returns path of generated xyz-file.

        self.refdir_path = None #new LM20231124. Gets defined in _process_equivalent_structures(). Contains path to parent directory.
        self.mol = Molecule.from_file(self.xyz_path) #pymatgen interface
        self.cmol = self.mol.get_centered_molecule()

        if case == 1: #new LM20231128
            #self.rigid_atom_idx_0 = None #1 = O if using water #used for testing if everything went well TODO: Implement com = False case. LM20231128
            self.rigid_atom_idx_1 = None #used for testing if everything went well
            self.rigid_atom_idx_2 = None #used for testing if everything went well
            if type(rigid_atom_1) == str and type(rigid_atom_2) == str: #if for some reason indices are given instead of atom names
                first = False #if H and H are chosen as rigid atoms
                for atom in range(len(self.mol.labels)):
                    if not first and rigid_atom_1 == self.mol.labels[atom]:
                        self.rigid_atom_idx_1 = atom
                        first = True
                    elif rigid_atom_2 == self.mol.labels[atom]:
                        self.rigid_atom_idx_2 = atom
                if self.rigid_atom_idx_1 == None or self.rigid_atom_idx_2 == None:
                    quit("Something is wrong with declaring the rigid atoms. Please check all-settings.yaml!")
            else:
                # self.rigid_atom_idx_0 = rigid_atom_0  # 1 = O if using water TODO: Implement com = False case. LM20231128
                self.rigid_atom_idx_1 = rigid_atom_1  # 2 = H if using water
                self.rigid_atom_idx_2 = rigid_atom_2 #3 = H if using water
        else:
            # self.rigid_atom_idx_0 = rigid_atom_0 #1 = O if using water
            self.rigid_atom_idx_1 = rigid_atom_1 #2 = H if using water
            self.rigid_atom_idx_2 = rigid_atom_2 #3 = H if using water TODO: Implement com = False case. LM20231128

        print("Using rigidatoms {0} and {1} for quaternion determination!".format(self.rigid_atom_idx_1, self.rigid_atom_idx_2))

        self.char_q = calc_quats(self.cmol, self.rigid_atom_idx_1, self.rigid_atom_idx_2)

        #placement
        self.eq_dict = {} #new. stores 1) the rotation quat, 2) the equivalent structure as Molecule object and 3) the characteristic quat for the orientation of the molecule as tuple per symm_op. LM20231113
        self._process_equivalent_structures(True)
        #self.symmetry_rots_as_quats = {} #DEPRECATED LM20231113. stores the quaternions of the rotation matrices that transform the equivalent structure into the initial structure. LM20231006
        #self.equivalent_structures = {} ##DEPRECATED LM20231113. stores symmetry equivalent structures as Molecule object
        #self.reference_quats = {} ##DEPRECATED LM20231113. stores the quaternions associated to the equivalent structures calculated from the same 3 rigid atoms LM20231006

        #post placement. LM20231005
        self.elements = [] #contains element names
        self.atoms = [] #contains xyz coordinates
        self.values = [] #contains the
        self.all_values = [] #contain the temperature factor for each atom 3 times. TODO:verify. maybe it is the energy value
        self.keep = False #cleans up xyz-files. not implemented yet as of 20231006

    def _makedir(self, path: str) -> str:  # used primarily to create a folder that contains the reference structure of pyConSolv solvents in structures.py
        import datetime
        path = path + "_{0}".format(str(datetime.date.today()))
        if not os.path.exists(path):
            os.makedirs(path)
        else:
            num = 1
            while os.path.exists(path + "_{0}".format(str(num))):
                num += 1
            path = path + "_{0}".format(str(num))
            os.makedirs(path)
        return os.path.abspath(path)

    def _process_equivalent_structures(self, write_out: bool = True):

        """
        analyzes the given solvent, gives coordinates for equivalent structures and sets the following:
        1) self.eq_dict
        """

        self.refdir_path = self._makedir("REF_{0}".format(self.abb)) #directory containing the xyz files of equivalent structures
        path_template = self.refdir_path+"/{0}".format(self.abb)+"_{0}.xyz" #last placeholder is for enumeration of files used in write()

        #building up self.equivalent_structures and self.symmetry_rots_as_quats
        pga = ana.PointGroupAnalyzer(self.cmol)
        symm_ops = pga.get_symmetry_operations()

        num = 0
        for symm in symm_ops:
            coord_list = []
            rot = symm.rotation_matrix
            return_path = path_template.format(num) #enumerates the xyz-files
            if np.abs(np.linalg.det(rot) - 1) < 1e-4: #adopted from analyzer.get_rotational_symmetry_number() (line 1284). Filters only for those transformations that have adeterminant > 1 which means no mirror. LM20231027
                rot_quat = quat.from_rotation_matrix(rot) #transform the rotation matrix into quaternion using quaternion package. LM20231006
            else:
                continue #rotation with mirror has been found. skipping this symm_op.

            for idx in range(len(self.mol.cart_coords)):
                coord_list.append(symm.apply_rotation_only(self.cmol.cart_coords[idx]))
            write(self.xyz_path, return_path, coord_list)
            ref = Molecule.from_file(return_path.format(num))

            char_quat = calc_quats(ref, self.rigid_atom_idx_1, self.rigid_atom_idx_2) #LM20231128: Changed from rigid_atom_0/1 to ..._1/2


            if distance(self.char_q,char_quat) < 0.05: #0.05 rad is around 3.18°
                os.rename(return_path,self.refdir_path+"/{0}".format(self.abb)+"_ori.xyz")
                #continue #the original orientation has been found. LM20231130: now also entry in eq_dict.
                self.eq_dict[num] = (rot_quat, ref, char_quat) #new LM20231130: now also for the identity there is an entry in eq_dict. This facilitates the quaternion cleanup in _find_avg_solvent.

            else:
                self.eq_dict[num] = (rot_quat, ref, char_quat)

            num += 1

    def _find_avg_solvent(self, voxel: int, quats: list[quat.quaternion], com: tuple[typing.Any,typing.Any,typing.Any], verbose = False): #new LM20231124
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

            if verbose:
                print("\n\nquat: {0}".format(i))
                print("distance_list: {0}".format(distance_list))

            min_idx = distance_list.index(min(distance_list)) #source: https://stackoverflow.com/questions/52294174/python-finds-the-index-of-the-smallest-element-in-the-list-a-from-index-k-onwa (accessed 24 November 2023)

            if verbose:
                print("min_idx: {0}".format(min_idx))
                print("quats[{1}] before clean-up: {0}".format(quats[i],i))

            quats[i] = quats[i] * inv(self.char_q) * inv(self.eq_dict[key_list[min_idx]][0]) * self.char_q #LM20231222: Inverse since we're not interested in the rotation from the origin to the equivalent representation but the other way round.

            if verbose:
                print("quats[{1}] after clean-up: {0}".format(quats[i],i))
                print("distance now: {0}".format(distance(quats[i], self.eq_dict[key_list[min_idx]][2])))

        #step 2)
        q_avg = quat.from_float_array(avg(create_Q_matrix(quats)))

        #step 3) and 4)
        qt = q_avg * inv(self.char_q) #quaternion for rotation of original orientation to orientation described by q_avg
        new_coords = new_coord_gen(self.cmol, qt, np.array(com))

        if verbose:
            write(self.xyz_path, self._makedir("quats")+"/avg_at_voxel_{0}", new_coords, voxel)

        return self.cmol.labels, new_coords





class Solvent: #in original Febiss: class Water
    """
    This class serves two purposes:
    1) A Solvent object is created which gets all the data from the datafile created by Analysis_Febiss.h/cpp.
       This object makes use of the attributes data, coords (for the COM coords) and values and the methods sort_by_energy, get_coord_set and get_energy.
    2) Solvent objects are also created for each solvent selected from the interactive barplot.
       These objects make use of the attribute coords (for the COM coord),
       quats (which are loaded from the gist-quats.dat file), elements and elem_coords (for the solvent atom coords).
    """
    def __init__(self, nvoxels = 0):
        self.data = [] #new LM20231123: data (i.e. voxel x y z energy) is passed as 5-tuple
        self.coords = [] #new LM20231123: holds coords of COMs
        self.values = [] #new (actually also in original Febiss) LM20231123: holds energy values of solvent molecules. needed for plotting
        self.elements = [] #new (actually also in original Febiss) LM20231123: holds the elements of the solvent to be placed
        self.quats = self._prep_dict(nvoxels) #new LM20231123: holds all the quats from gist-quats.dat associated with. keys are voxel numbers
        self.elem_coords = []

    def sort_by_energy(self): #renamed from sort_by_value. LM20231123
        self.data = sorted(self.data, key=lambda tpl: tpl[-1],reverse=True)
        with open('sorted_data.dat','w') as f:
            for i in self.data:
                f.write("{0} {1} {2} {3} {4}\n".format(i[0], i[1], i[2], i[3], i[4]))

    def get_coord_set(self): #new LM20231123:
        coord_list = []
        for com in self.data:
            coord_list.append((float(com[1]),float(com[2]),float(com[3]))) #TODO: Prone to ValueError. LM20231127
        self.coords = np.asarray(coord_list)

    def get_energy(self): #TODO: merge with get_coord_set
        for com in self.data:
            self.values.append(float(com[-1])) #TODO: Prone to ValueError LM20231127

    def _prep_dict(self, nvoxels):
        quat_dict = {}
        if nvoxels != 0:
            for voxel in range(nvoxels):
                quat_dict[voxel] = []
        return quat_dict

