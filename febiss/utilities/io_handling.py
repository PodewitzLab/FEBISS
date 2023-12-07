#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright University Innsbruck, Institute for General, Inorganic, and Theoretical Chemistry, Podewitz Group
See LICENSE for details
"""

from typing import Union
import numpy as np
import os

from ..utilities.structures import Solute, Solvent, Reference


# def read_pdb(pdb, solute: Solute, solvent: Reference):
#     with open(pdb, 'r') as f:
#         for line in f:
#             if 'HETATM' in line:
#                 row = line.split()
#                 solvent.elements.append(row[-1])
#                 solvent.atoms.append(np.array([float(r) for r in row[-6:-3]]))
#                 if row[2] == solvent.rigid_atom_idx_0: #Changed from row[-1] (which contains the ambiguous element name) to row[2] which contains the atom label. TODO: Update pdb file generation in CPPTRAJ
#                     solvent.values.append(-1 * float(row[-2]))
#                     for i in range(solvent.size): #adds the tempfactor/energy value for each solvent atom
#                         solvent.all_values.append(-1 * float(row[-2]))
#
#                 #elif row[-1] != 'H':
#                 #    raise NotImplementedError("ERROR: NON-solvent HETATM present in pdb file")
#             elif 'ATOM' in line:
#                 row = line.split()
#                 solute.elements.append(row[-1])
#                 solute.atoms.append(np.array([float(r) for r in row[-6:-3]]))
#                 solute.values.append(0.0)
#
#     solute.atoms = np.asarray(solute.atoms)
#     solute.determine_polar_hydrogen_and_non_hydrogen()
#     solvent.atoms = np.asarray(solvent.atoms)
#     solvent.sort_by_value()

# def read_febiss_file(febiss_file : str, )


def write_pdb(pdb: str, structure: Union[Solute, Solvent], abb, solute: bool = False):
    if solute:
        atomcounter = 1
        f = open(pdb, 'w')
    else:
        atomcounter = len(open(pdb, 'r').readlines()) + 1
        f = open(pdb, 'a')
        #print(structure.elements)
        #print(structure.coords)
    for count, (ele, atom) in enumerate(zip(structure.elements, structure.coords)): #changed from structure.atoms LM20231130
        j = []
        if solute:
            j.append('ATOM'.ljust(6))  # atom#6s
        else:
            j.append('HETATM'.ljust(6))  # atom#6s
        j.append(str(atomcounter + count).rjust(5))  # aomnum#5d
        j.append(ele.center(4))  # atomname$#4s
        if solute:
            j.append('SOL'.ljust(3))  # resname#1s
        else:
            j.append(abb.ljust(3))  # resname#1s abb instead of "FEB"
        j.append('A'.rjust(1))  # Astring
        if solute:
            j.append('1'.rjust(4))  # resnum
        else:
            j.append('2'.rjust(4))  # resnum
        j.append(str('%8.3f' % (float(atom[0]))).rjust(8))  # x
        j.append(str('%8.3f' % (float(atom[1]))).rjust(8))  # y
        j.append(str('%8.3f' % (float(atom[2]))).rjust(8))  # z
        j.append(str('%6.2f' % 1.0).rjust(6))  # occ
        if solute: #new LM20231130: introduced since solute does not have an attribute "values" anymore
            value = 0.0
        else:
            #print(structure.values[count])
            value = float(structure.values[count])
        if value == 0.0:
            j.append(str('%7.2f' % value).ljust(7))  # delta G
        else:
            j.append(str('%7.2f' % (-1 * value)).ljust(7))  # delta G
        j.append(ele.rjust(12))  # elname
        f.write("%s%s %s %s %s%s    %s%s%s%s%s%s\n" % (
            j[0], j[1], j[2], j[3], j[4], j[5], j[6], j[7], j[8], j[9], j[10], j[11]))
    f.close()


def write_style_file() -> str:
    filename = 'style.pml'
    with open(filename, 'w') as f:
        f.write('hide everything\n')
        f.write('show sticks\n')
        f.write('set stick_radius, .15\n')
        f.write('set sphere_scale, .18\n')
        f.write('set sphere_scale, .13, elem H\n')
        f.write('set bg_rgb=[1, 1, 1]\n')
        f.write('set stick_quality, 50\n')
        f.write('set sphere_quality, 4\n')
        f.write('color gray35, elem C\n')
        f.write('color red, elem O\n')
        f.write('color blue, elem N\n')
        f.write('color gray98, elem H\n')
        f.write('set ray_texture, 2\n')
        f.write('set antialias, 3\n')
        f.write('set ambient, 0.5\n')
        f.write('set spec_count, 5\n')
        f.write('set shininess, 50\n')
        f.write('set specular, 1\n')
        f.write('set reflect, .1\n')
        f.write('set cartoon_ring_finder, 4\n')
        f.write('set cartoon_ring_mode,1\n')
        f.write('set cartoon_ring_transparency, 0.6\n')
        f.write('set cartoon_ring_color, black\n')
        f.write('show cartoon\n')
        f.write('set h_bond_cutoff_center, 3.5\n')
        f.write('set h_bond_cutoff_edge, 3.5\n')
        f.write('set h_bond_max_angle, 135\n')
        f.write('set dash_gap, .25\n')
        f.write('set dash_length, .02\n')
        f.write('set dash_round_ends, 1\n')
        f.write('set dash_radius, .05\n')
        f.write('set opaque_background, off\n')
        f.write('set stick_h_scale, 1\n')
        f.write('set label_digits, 2\n')
        f.write('label ele o and resn FEB, b\n')
        f.write('select solute, not resn "FEB"\n')
        f.write('select solventMolecules, resn "FEB"\n')
        f.write('distance solute-solvent, solute, solventMolecules, cutoff=3.2, mode=2\n')
        f.write('set dash_color, green\n')
        f.write('spectrum b, magenta_white_yellow, ele o and resn FEB\n')
        f.write('hide labels, solute-solvent\n')
        f.write('center\n')
    return filename


class Input:
    def __init__(self, prompt, **kwargs):
        self.prompt = str(prompt)
        self.input = input(self.prompt +"\n")
        if len(kwargs) > 0:
            self.form = kwargs["type"]
            self.__assertion()

    def yn(self) -> bool:
        while self.input not in ['Y', 'N', 'y', 'n']:
            self.prompt = "Please respond with 'y' or 'n' only!\n"
            self.input = input(self.prompt + "\n")
        if self.input.lower() == "y":
            return True
        else:
            return False

    def reassure(self):
        if not Input("You provided the following input: " + str(self.input) + "\n Is that correct?").yn():
            self.input = input(self.prompt + "\n")
            self.__assertion()

    def __assertion(self):
        try:
            self.input = self.form(self.input + "\n")
            self.reassure()
        except ValueError:
            print("The input does not fit the requirements, must be " + str(self.form.__name__))
            self.input = input(self.prompt + "\n")
            self.__assertion()
        except SyntaxError:
            print("The input does not fit the requirements, must be " + str(self.form.__name__))
            self.input = input(self.prompt + "\n")
            self.__assertion()
