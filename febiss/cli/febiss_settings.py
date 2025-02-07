#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright Technische Universität Wien, Institute of Materials Chemistry, Podewitz Group
See LICENSE for details
"""

import sys
import os
from ..utilities.io_handling import Input
from ..solvents import CASE_DICT, SOLVENT_LIST, FILE_DICT, RIGID_ATOMS_DICT, REF_DENS_DICT, REF_EWW

def help_message():
    print("\nThis program writes all possible settings of the GIST analysis "
          "and the plotting into 'all-settings.yaml'.")
    print('It requires no arguments but the topology and trajectory file can be passed (in this order).')
    sys.exit()

def pyconsolv_abb():
    """
    Function used if pyConSolvents are used (contained in this package in /febiss/solvents/). pyConSolv solvent will
    be chosen according to the user's input.

    :return: 3 letter abbreviation of the pyConSolv solvent used for this simulation to be analyzed.
    """
    solv_abb = input(
                "\n\nPlease enter the 3 character long string in the parentheses after the solvent you used. "
                + "   ".join(SOLVENT_LIST) + ":\n")
    while solv_abb not in list(FILE_DICT.keys()):
        solv_abb = input(
            "\n\nPlease enter the 3 character long string in the parentheses after the solvent you used. "
            "Write 'quit' to quit program or 'list' to see solvent list again:\n")
        if solv_abb.lower() == "quit":
            quit("\n---Quit as you wished!---")
        if solv_abb.lower() == "list":
            solv_abb = input("\n\nHere is the list of PyConSolv solvents:\n\n" + "   ".join(
                SOLVENT_LIST) + "\n\nEnter the 3 character long string of the solvent you used:\n")
    return solv_abb

def main():
    top = 'SOLVBOX_TOPOLOGY'
    traj = 'TRAJECTORY'

    if len(sys.argv) > 3 or len(sys.argv) == 2:
        help_message()
    elif len(sys.argv) == 3:
        top = sys.argv[1]
        traj = sys.argv[2]

    if Input("\nDo you use TIP3P water? [y/n]").yn():
            case = 1

    elif Input("\nDid you use one of these solvents from the PyConSolv package "
                          "(https://github.com/PodewitzLab/PyConSolv/tree/main/src/PyConSolv/solvents)? [y/n]").yn():
        case = 3
        solv_abb = pyconsolv_abb()

    else:
        case = 4

    com = False
    if case in [1, 3] and Input("\nDo you want to use the center of mass (COM)? [y/n]").yn():
        com = True

    from ..utilities.gist import GistAnalyser
        # the arguments are necessary in the init, otherwise exception
    if case == 1:
        CASE_DICT[1] = os.path.abspath(os.path.join(__file__, "../../solvents/TP3.xyz"))
        analyser = GistAnalyser(case, com, **{'top': top,
                                              'trajectory_file': traj,
                                              'solv_abb': 'WAT',
                                              'solv_file': CASE_DICT[1],
                                              'refdens':REF_DENS_DICT['TP3'],
                                              'ref_eww':REF_EWW['TP3'],
                                              'rigid_atom_0': RIGID_ATOMS_DICT['TP3'][1],
                                              'rigid_atom_1': RIGID_ATOMS_DICT['TP3'][0],
                                              'rigid_atom_2': RIGID_ATOMS_DICT['TP3'][2]
                                              })

    elif case == 3:
        CASE_DICT[3] = os.path.abspath(os.path.join(__file__, "../../solvents/{0}.mol2".format(solv_abb)))


        analyser = GistAnalyser(case, com,  **{'top': top,
                                               'trajectory_file': traj,
                                               'solv_abb': solv_abb,
                                               'solv_file': CASE_DICT[3],
                                               'refdens': REF_DENS_DICT[solv_abb],
                                               'ref_eww': REF_EWW[solv_abb],
                                               'rigid_atom_0' : RIGID_ATOMS_DICT[solv_abb][1],
                                               'rigid_atom_1': RIGID_ATOMS_DICT[solv_abb][0],
                                               'rigid_atom_2': RIGID_ATOMS_DICT[solv_abb][2]
                                               })

    else:
        analyser = None
        quit("\n\nFebiss in it's current state does not allow for user-defined solvents. We are working on it!\n"
             "---Quit program---\n")

    # write all analysis options
    with open('all-settings.yaml', 'w') as f:
        f.write("# Case 1: TIP3P water as solvent\n")
        f.write("# Case 3: pyConSolv solvent used. The reference xyz file is stored in febiss/solvents.)\n")
        f.write("# Case 4: User defined solvent used. Currently not available. Terminates program.)\n\n")
        f.write("# Two general blocks 'gist' and 'plotting' are given.\n")
        f.write('# For each block all settings and the default values are given and\n')
        f.write('# this file can be used directly after filling out the two required settings.\n')
        f.write('# Where it may be useful, a possible input is given behind the default value as a comment.\n')
        f.write("header:\n")
        f.write('  case: ' + str(case) + '\n')
        f.write('  com: ' + str(com) + '\n')
        f.write('gist:\n')

        for key in sorted(analyser.required_keys):
            if case == 3 and com and key == 'rigid_atom_0':
                f.write('  ' + str(key) + ": " + str(analyser.__dict__[key]) + '\n')
                # no check required since -1 is a necessary value when using com then

            elif case in [1, 3] and key in ['solv_abb',
                                            'solv_file',
                                            'rigid_atom_0',
                                            'rigid_atom_1',
                                            'rigid_atom_2',
                                            'refdens',
                                            'ref_eww']:
                f.write('  ' + str(key) + ": " + str(analyser.__dict__[key]) + ' # please check if this is correct\n')
                # if pyconsolv solvents are used, solv_abb is typically the 3 letter abbreviation given as input.
                # it has to be checked nonetheless just like the file path and the rigid_atom indices.

            else:
                f.write('  ' + str(key) + ": " + str(analyser.__dict__[key]) + ' # this has to be filled out\n')

        for key in sorted(analyser.allowed_keys):
            if key == 'frame_selection':
                f.write('  ' + str(key) + ": " + str(analyser.__dict__[key]) + " # example: '1 1000 5'\n")
            elif key == 'grid_center':
                f.write('  ' + str(key) + ": " + str(analyser.__dict__[key]) + " # example: (10.0, 5.0, 0.0)\n")
            elif key == 'rdf_names':
                f.write('  ' + str(key) + ": " + str(analyser.__dict__[key]).replace('center of solute',
                                                                                     'center\ of\ solute') + '\n')
            else:
                if analyser.__dict__[key] != None:
                    f.write('  ' + str(key) + ": " + str(analyser.__dict__[key]) + '\n')

    from ..plotting.display import Plot
    display = Plot()
    # write all plotting options
    with open('all-settings.yaml', 'a') as f:
        f.write('plotting:\n')
        for key in sorted(display.allowed_keys):
            f.write('  ' + str(key) + ": " + str(display.__dict__[key]) + '\n')
    print("Wrote all possible settings into 'all-settings.yaml' in the current directory.")
    sys.exit()


if __name__ == '__main__':
    main()
