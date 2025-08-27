#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright Technische Universität Wien, Institute of Materials Chemistry, Podewitz Group
See LICENSE for details
"""

def write_out_gist_grid(analyzer):
    gist_data = open(analyzer.gist_out_file, 'r').readlines()
    solute_elements = []
    solute_atoms = []

    with open('solute.pdb', 'r') as f:
        for line in f:
            row = line.split()
            if row[0] == 'ATOM':
                solute_elements.append(row[2])
                solute_atoms.append([row[5], row[6], row[7]])
    n_voxels = len(gist_data) - 2

    with open(analyzer.gist_grid_file, 'w') as f:
        f.write(str(n_voxels + len(solute_elements)) + '\n\n')
        for e, a in zip(solute_elements, solute_atoms):
            f.write(e + '\t' + a[0] + '\t' + a[1] + '\t' + a[2] + '\n')
        for line in gist_data[2:]:
            row = line.split()
            f.write('H\t' + row[1] + '\t' + row[2] + '\t' + row[3] + '\n')