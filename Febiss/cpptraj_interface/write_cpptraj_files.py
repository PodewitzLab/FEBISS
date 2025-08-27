#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright Technische Universität Wien, Institute of Materials Chemistry, Podewitz Group
See LICENSE for details
"""



# Open file and write
def write_gist_cpptraj_file(analyzer):  # this creates the content for cpptraj.in.
    with open(analyzer.gist_cpptraj_command_file, 'w') as f:
        f.write('parm ' + analyzer.top + '\n')
        f.write('trajin ' + analyzer.trajectory_file)
        if str(analyzer.frame_selection).lower() not in ['', 'none']:
            # if not given CPPTRAJ uses all frames
            f.write(' ' + analyzer.frame_selection)
        f.write('\n')

        f.write('solvent ' + f':{analyzer.solv_abb}\n')
        f.write('center ' + analyzer.solute_residues + ' origin\n')
        f.write('image origin center familiar\n')

        #GIST command
        f.write('gist ')
        if str(analyzer.grid_center).lower() not in ['', 'none']:
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

        #additional write_out of first frame and solute
        f.write('center :1 origin\n')
        f.write('trajout first_frame.pdb onlyframes 1\n')
        f.write('run\n')
        f.write('center :1 origin\n')
        f.write('strip :{0}\n'.format(analyzer.solv_abb))
        f.write('trajout solute.pdb onlyframes 1\n')
        f.write('run')


def write_rdf_input_files(analyzer):
    file_list = []
    for key, val in analyzer._rdf_names.items():
        path = analyzer._rdf_cpptraj_command_file.format(key)
        with open(path, 'w') as f:
            f.write('parm ' + analyzer.top + '\n')
            f.write('trajin ' + analyzer.trajectory_file)
            if str(analyzer.frame_selection).lower() not in ['', 'none']:
                # if not given CPPTRAJ uses all frames
                f.write(' ' + analyzer.frame_selection)
            f.write('\n')
            f.write('solvent ' + f':{analyzer.solv_abb}\n')
            f.write('center ' + analyzer.solute_residues + ' origin\n')
            f.write('image origin center familiar\n')
            f.write('radial ')
            f.write('spacing ' + str(analyzer.rdf_spacing) + ' ' + str(analyzer.rdf_maximum) + ' ')
            f.write('density ' + str(analyzer.refdens) + ' ')
            if key == 'center':
                f.write(':{0} '.format(analyzer.solv_abb) + analyzer.solute_residues + ' byres1 ' + 'byres2 ')
            else:
                f.write(':{0} '.format(analyzer.solv_abb) + analyzer.solute_residues + '@/' + key + ' byres1 ')
            f.write("out 'rdf-" + val + ".dat' ")
            f.write("intrdf 'int-rdf-" + val + ".dat' ")
            f.write("rawrdf 'raw-rdf-" + val + ".dat'\n")
            f.write('run')
            f.close()
        file_list.append(path)

    return file_list