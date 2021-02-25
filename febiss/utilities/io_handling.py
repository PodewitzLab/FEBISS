#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright University Innsbruck, Institute for General, Inorganic, and Theoretical Chemistry, Podewitz Group
See LICENSE for details
"""

from typing import Union
import numpy as np

from febiss.utilities.structures import *


def read_pdb(pdb, solute: Solute, water: Water):
    with open(pdb, 'r') as f:
        for line in f:
            if 'HETATM' in line:
                row = line.split()
                water.elements.append(row[-1])
                water.atoms.append(np.array([float(r) for r in row[-6:-3]]))
                if row[-1] == "O":
                    water.values.append(-1 * float(row[-2]))
                    water.all_values.append(-1 * float(row[-2]))
                    water.all_values.append(-1 * float(row[-2]))
                    water.all_values.append(-1 * float(row[-2]))
                elif row[-1] != 'H':
                    raise NotImplementedError("ERROR: NON-WATER HETATM present in pdb file")
            elif 'ATOM' in line:
                row = line.split()
                solute.elements.append(row[-1])
                solute.atoms.append(np.array([float(r) for r in row[-6:-3]]))
                solute.values.append(0.0)

    solute.atoms = np.asarray(solute.atoms)
    solute.determine_polar_hydrogen_and_non_hydrogen()
    water.atoms = np.asarray(water.atoms)
    water.sort_by_value()


def write_pdb(pdb: str, structure: Union[Solute, Water], solute: bool = False):
    if solute:
        atomcounter = 1
        f = open(pdb, 'w')
    else:
        atomcounter = len(open(pdb, 'r').readlines()) + 1
        f = open(pdb, 'a')
    for count, (ele, atom) in enumerate(zip(structure.elements, structure.atoms)):
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
            j.append('FEB'.ljust(3))  # resname#1s
        j.append('A'.rjust(1))  # Astring
        if solute:
            j.append('1'.rjust(4))  # resnum
        else:
            j.append('2'.rjust(4))  # resnum
        j.append(str('%8.3f' % (float(atom[0]))).rjust(8))  # x
        j.append(str('%8.3f' % (float(atom[1]))).rjust(8))  # y
        j.append(str('%8.3f' % (float(atom[2]))).rjust(8))  # z
        j.append(str('%6.2f' % 1.0).rjust(6))  # occ
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
        f.write('select waterMolecules, resn "FEB"\n')
        f.write('distance solute-water, solute, waterMolecules, cutoff=3.2, mode=2\n')
        f.write('set dash_color, green\n')
        f.write('spectrum b, magenta_white_yellow, ele o and resn FEB\n')
        f.write('hide labels, solute-water\n')
        f.write('center\n')
    return filename
