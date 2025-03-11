#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright Technische Universität Wien, Institute of Materials Chemistry, Podewitz Group
See LICENSE for details
"""

from typing import Union

from ..utilities.structures import Solute, Solvent

def write_pdb(pdb: str, structure: Union[Solute, Solvent], abb, solute: bool = False):
    if solute:
        atomcounter = 1
        f = open(pdb, 'w')
    else:
        atomcounter = len(open(pdb, 'r').readlines()) + 1
        f = open(pdb, 'a')

    for count, (ele, atom) in enumerate(zip(structure.elements, structure.coords)):
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
        if solute:
            value = 0.0
        else:
            value = float(structure.values[count])
        if value == 0.0:
            j.append(str('%7.2f' % value).ljust(7))  # delta G
        else:
            j.append(str('%7.2f' % (-1 * value)).ljust(7))  # delta G
        j.append(ele.rjust(12))  # elname
        f.write("%s%s %s %s %s%s    %s%s%s%s%s%s\n" % (
            j[0], j[1], j[2], j[3], j[4], j[5], j[6], j[7], j[8], j[9], j[10], j[11]))
    f.close()

def write_xyz(ori_path, return_path, new_coords, idx=0, labels = None):
    """
    :param ori_path: ori_path: path of the original xyz file that acts as a template
    :param return_path: path of the files to be created. accepts one placeholder {0} for an index that goes from 0 to num
    :param new_coords: newly created coordinates after rotation, translation, ...
    :param idx: Optional. Used for numbering the created xyz-files
    :param labels: Optional. If labels of atoms are already known they can be passed as list. Otherwise this function
           will determine them from the xyz file using the pymatgen.core package "Molecule".
    :return path of written file


    Takes an xyz-file (ori_path) as template to create an xyz-file (return_path) with the new_coords. Return_path
    accepts a placeholder {0} for numbering with idx (default: 0).
    """

    if labels is None:
        from pymatgen.core import Molecule
        mol = Molecule.from_file(ori_path)
        labels = mol.labels

    xyzlines = []
    text = open(ori_path, 'r').readlines()
    xyzlines.extend([text[0], text[1]])
    for k in range(len(new_coords)):
        line = "{0:<2} {1:>12}{2:>12}{3:>12}\n".format(labels[k], round(new_coords[k][0], 5),
                                                       round(new_coords[k][1], 5), round(new_coords[k][2], 5))
        xyzlines.append(line)
    with open(return_path.format(idx),'w') as newfile:
        newfile.writelines(xyzlines)

    return return_path.format(idx)

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
