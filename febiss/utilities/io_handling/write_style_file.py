#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright Technische Universität Wien, Institute of Materials Chemistry, Podewitz Group
See LICENSE for details
"""

def write_style_file() -> str:
    """
    :return: No actual return but writes style.pml file that contains all necessary settings for displaying the
    microsolvated structure with PyMol
    """
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