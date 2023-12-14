#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright University Innsbruck, Institute for General, Inorganic, and Theoretical Chemistry, Podewitz Group
See LICENSE for details
"""

import sys
import os
from ..utilities.io_handling import Input
from ..solvents import CASE_DICT, SOLVENT_LIST, FILE_DICT, RIGID_ATOMS_DICT, REF_DENS_DICT

def help_message():
    print("\nThis program writes all possible settings of the GIST analysis "
          "and the plotting into 'all-settings.yaml'.")
    print('It requires no arguments.')
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
    if len(sys.argv) > 1:
        help_message()

    if Input("Is water your main solvent? [y/n]").yn():
        if Input("\nDo you use TIP3P water? [y/n]").yn():
            case = 1
        else:
            case = 2
    else:
        if Input("\nDid you use one of these solvents from the PyConSolv package "
                          "(https://github.com/PodewitzLab/PyConSolv/tree/main/src/PyConSolv/solvents)? [y/n]").yn():
            case = 3
            solv_abb = pyconsolv_abb()
        else:
            case = 4

    if case in [1, 2, 3] and Input("\nDo you want to use the center of mass (COM)? [y/n]").yn():
        com = True
    else:
        com = False

    from ..utilities.gist import GistAnalyser
        # the arguments are necessary in the init, otherwise exception
    if case == 1:
        CASE_DICT[1] = os.path.abspath(os.path.join(__file__, "../../solvents/TP3.xyz")) #new LM20231128 #https://stackoverflow.com/questions/27844088/python-get-directory-two-levels-up (accessed 14 Nov 2023).
        analyser = GistAnalyser(case, com, **{'top': 'SOLVBOX_TOPOLOGY',
                                              'trajectory_file': 'TRAJECTORY',
                                              'solv_abb': 'WAT',
                                              'solv_file': CASE_DICT[1] #new LM20231128
                                              }) #changed from tracectory_name to trajectory_file LM20231115
        #solv_file = False #DEPRECATED LM20231114

    elif case == 2:
        analyser = GistAnalyser(case, com, **{'top': 'SOLVBOX_TOPOLOGY',
                                              'trajectory_file': 'TRAJECTORY',  #changed from tracectory_name to trajectory_file LM20231115
                                              #'solv_top':'SOLVENT_TOPOLOGY',
                                              'solv_abb': 'ABB',
                                              'solv_file': CASE_DICT[2],
                                              'refdens':'REFDENS',
                                              'rigid_atom_0': "IDX_0",  #new: also in this case the rigid atom indices have to be defined
                                              'rigid_atom_1': "IDX_1",  #new: also in this case the rigid atom indices have to be defined
                                              'rigid_atom_2': "IDX_2",  # new: also in this case the rigid atom indices have to be defined
                                              #'char_angle':'CHAR_ANGLE' #not needed when using quaternions
                                              })
        #solv_file = False #DEPRECATED LM20231114

    elif case == 3:
        #solv_file = FILE_DICT[solv_abb] #DEPRECATED

        CASE_DICT[3] = os.path.abspath(os.path.join(__file__, "../../solvents/{0}.mol2".format(solv_abb))) #https://stackoverflow.com/questions/27844088/python-get-directory-two-levels-up (accessed 14 Nov 2023).


        analyser = GistAnalyser(case, com,  **{'top': 'SOLVBOX_TOPOLOGY',
                                               'trajectory_file': 'TRAJECTORY', #changed from tracectory_name to trajectory_file LM20231115
                                               #'solv_top':'SOLVENT_TOPOLOGY',
                                               'solv_abb': solv_abb,
                                               'solv_file': CASE_DICT[3],
                                               #'solv_size':'SOLV_SIZE',
                                               'refdens': REF_DENS_DICT[solv_abb],
                                               #'char_angle':'CHAR_ANGLE',
                                               'rigid_atom_0' : RIGID_ATOMS_DICT[solv_abb][1], #LM20231116: changed the number back again from 0 to 1 since one shall be able to choose if COM or central atom shall be used #changed number in bracket from 1 to 0 since the COM and not a central atom will be used for the characteristic quat calculation
                                               'rigid_atom_1': RIGID_ATOMS_DICT[solv_abb][0], #LM20231116: changed the number back again from 1 to 0 since one shall be able to choose if COM or central atom shall be used #changed number in bracket from 0 to 1 since the COM and not a central atom will be used for the characteristic quat calculation
                                               'rigid_atom_2': RIGID_ATOMS_DICT[solv_abb][2]
                                               })

    else:
        analyser = None
        quit("\n\nFebiss in it's current state does not allow for user-defined solvents. We are working on it!\n"
             "---Quit program---\n")

        #analyser = GistAnalyser(water, tip3p, **{'top': 'SOLVBOX_TOPOLOGY', 'trajectory_name': 'TRAJECTORY',
        #                                  'solv_top':'SOLVENT_TOPOLOGY', 'solv_abb':'ABB','solv_size':'SOLV_SIZE','refdens':'REFDENS',
        #                                  'char_angle':'CHAR_ANGLE','rigid_atom_0':'RIGID_ATOM_0',
        #                                  'rigid_atom_1':'RIGID_ATOM_1','rigid_atom_2':'RIGID_ATOM_2'})
    # write all analysis options
    with open('all-settings.yaml', 'w') as f:
        f.write("# Case 1: TIP3P water as solvent\n") # changed header LM20231114
        f.write("# Case 2: Non-TIP3P water as solvent. Path to reference xyz file of the water molecule must be passed.)\n")
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
            if case == 2 and com and key == 'rigid_atom_0':
                f.write('  ' + str(key) + ": " + str(analyser.__dict__[key]) + '\n') # no check required since -1 is a necessary value when using com then
            elif case == 3 and com and key == 'rigid_atom_0':
                f.write('  ' + str(key) + ": " + str(analyser.__dict__[key]) + '\n') # no check required since -1 is a necessary value when using com then
            elif case in [1, 3] and key in ['solv_abb','solv_file','rigid_atom_1','rigid_atom_2','refdens']: # if pyconsolv solvents are
                # used, solv_abb is typically the 3 letter abbreviation given as input. it has to be checked nonetheless
                # just like the file path and the rigid_atom indices. LM20231114
                #new: LM20231128 also case 1 now needs a reference file. Path to TP3.xyz is given.
                f.write('  ' + str(key) + ": " + str(analyser.__dict__[key]) + ' # please check if this is correct\n')
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
