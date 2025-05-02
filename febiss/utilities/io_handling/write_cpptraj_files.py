#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright Technische Universität Wien, Institute of Materials Chemistry, Podewitz Group
See LICENSE for details
"""

import os
import shutil
from ...utilities.io_handling.input import Input

# Open file and write
def _write_gist_cpptraj_file(analyzer):  # this creates the content for gist.in.
    with open(analyzer.gist_cpptraj_command_file, 'w') as f:
        f.write('parm ' + analyzer.top + '\n')
        f.write('trajin ' + analyzer.trajectory_file)
        if analyzer.frame_selection is not None and analyzer.frame_selection.lower() != 'none':
            # if not given CPPTRAJ uses all frames
            f.write(' ' + analyzer.frame_selection)
        f.write('\n')

        f.write('solvent ' + f':{analyzer.solv_abb}\n')
        f.write('center ' + analyzer.solute_residues + ' origin\n')
        f.write('image origin center familiar\n')

def _write_rdf_input_line(analyzer):
    file_list = []
    for key, val in analyzer._rdf_names.items():
        path = analyzer._rdf_cpptraj_command_file.format(key)
        with open(path, 'w') as f:
            f.write('parm ' + analyzer.top + '\n')
            f.write('trajin ' + analyzer.trajectory_file)
            if analyzer.frame_selection is not None and analyzer.frame_selection.lower() != 'none':
                # if not given CPPTRAJ uses all frames
                f.write(' ' + analyzer.frame_selection)
            f.write('\n')
            f.write('solvent ' + f':{analyzer.solv_abb}\n')
            f.write('center ' + analyzer.solute_residues + ' origin\n')
            f.write('image origin center familiar\n')
            f.write('radial ')
            f.write('spacing 0.05 10 ')
            f.write('density ' + str(analyzer.refdens) + ' ')
            if key == 'center2':
                f.write(':{0} '.format(analyzer.solv_abb) + analyzer.solute_residues + ' byres1 ' + 'byres2' + ' ')
            else:
                f.write(':{0} '.format(analyzer.solv_abb) + analyzer.solute_residues + '@/' + key + ' ')
            f.write("out 'rdf-" + val + ".dat' ")
            f.write("intrdf 'int-rdf-" + val + ".dat' ")
            f.write("rawrdf 'raw-rdf-" + val + ".dat'\n")
            f.write('run')
            f.close()
        file_list.append(path)

    return file_list


def _write_gist_input_line(analyzer):
    _write_gist_cpptraj_file(analyzer)
    with open(analyzer.gist_cpptraj_command_file, 'a') as f:
        f.write('gist ')
        if analyzer.grid_center is not None and analyzer.grid_center.lower() != 'none':
            # if not given CPPTRAJ uses origin
            f.write('gridcntr ')
            analyzer._write_variable(f, analyzer.grid_center)
        f.write('griddim ')
        analyzer._write_variable(f, analyzer.grid_lengths)
        f.write('gridspacn ' + str(analyzer.grid_spacing) + ' ')
        f.write('refdens ' + str(analyzer.refdens) + ' ')
        f.write('rigid_idx ' + str(analyzer.rigid_atom_0) + ' '
                + str(analyzer.rigid_atom_1) + ' '
                + str(analyzer.rigid_atom_2) + ' ')
        f.write('out ' + analyzer.gist_out_file + ' ')
        f.write('quat ')
        if not analyzer.com:
            f.write('nocom ')
        f.write('norm ')
        # f.write('dx\n')
        f.write('febiss\n')  # enables febiss placement in cpptraj
        f.write('run\n')
        f.write('center :1 origin\n')
        f.write('strip :{0}\n'.format(analyzer.solv_abb))
        f.write('trajout solute.pdb onlyframes 1\n')
        f.write('run')