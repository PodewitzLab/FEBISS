#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright University Innsbruck, Institute for General, Inorganic, and Theoretical Chemistry, Podewitz Group
See LICENSE for details
"""

import sys
from ..utilities.io_handling import Input
from ..solvents import SOLVENT_LIST, FILE_DICT, RIGID_ATOMS_DICT

def help_message():
    print("\nThis program writes all possible settings of the GIST analysis "
          "and the plotting into 'all-settings.yaml'.")
    print('It requires no arguments.')
    sys.exit()


def main():
    
    if len(sys.argv) > 1:
        help_message()
    water = Input("Is water your main solvent? [y/n]").yn()
    if water:
        tip3p = Input("\n\nDo you use TIP3P water?").yn()
        pyconsolv = False
    else:
        tip3p = False
        pyconsolv = Input("\n\nDid you use one of these solvents from the PyConSolv package "
                          "(https://github.com/PodewitzLab/PyConSolv/tree/main/src/PyConSolv/solvents)?:\n"
                          +"   ".join(SOLVENT_LIST)).yn()

    from ..utilities.gist import GistAnalyser
        # the arguments are necessary in the init, otherwise exception
    if tip3p:
        analyser = GistAnalyser(water, tip3p, pyconsolv, **{'top': 'SOLVBOX_TOPOLOGY', 'trajectory_name': 'TRAJECTORY'})
        solv_file = False
    elif water:
        analyser = GistAnalyser(water, tip3p, pyconsolv, **{'top': 'SOLVBOX_TOPOLOGY',
                                                 'trajectory_name': 'TRAJECTORY',
                                                 #'solv_top':'SOLVENT_TOPOLOGY',
                                                 'solv_abb':'ABB',
                                                 'refdens':'REFDENS',
                                                 #'char_angle':'CHAR_ANGLE' #not needed when using quaternions
                                                 })
        solv_file = False
    elif pyconsolv:
        solv_abb = input("\n\nPlease enter the 3 character long string in the parentheses after the solvent you used:\n")
        while solv_abb not in list(FILE_DICT.keys()):
            solv_abb = input("\n\nPlease enter the 3 character long string in the parentheses after the solvent you used. "
                             "Write 'quit' to quit program or 'list' to see solvent list again:\n")
            if solv_abb.lower() == "quit":
                quit("\n---Quit as you wished!---")
            if solv_abb.lower() == "list":
                solv_abb = input("\n\nHere is the list of PyConSolv solvents again:\n\n"+"   ".join(SOLVENT_LIST)+"\n\nEnter the 3 character long string of the solvent you used:\n")

        solv_file = FILE_DICT[solv_abb]

        analyser = GistAnalyser(water, tip3p, pyconsolv, **{'top': 'SOLVBOX_TOPOLOGY',
                                                 'trajectory_name': 'TRAJECTORY',
                                                 #'solv_top':'SOLVENT_TOPOLOGY',
                                                 'solv_abb':solv_abb,
                                                 #'solv_size':'SOLV_SIZE',
                                                 'refdens':'REFDENS',
                                                 #'char_angle':'CHAR_ANGLE',
                                                 'rigid_atom_0':RIGID_ATOMS_DICT[solv_abb][1], #center atom
                                                 'rigid_atom_1':RIGID_ATOMS_DICT[solv_abb][0],
                                                 'rigid_atom_2':RIGID_ATOMS_DICT[solv_abb][2]
                                                 })
    else:
        quit("\n\nFebiss in it's current state does not allow for user-defined solvents. We are working on it!\n"
                                  "---Quit program---\n")
        #analyser = GistAnalyser(water, tip3p, **{'top': 'SOLVBOX_TOPOLOGY', 'trajectory_name': 'TRAJECTORY',
        #                                  'solv_top':'SOLVENT_TOPOLOGY', 'solv_abb':'ABB','solv_size':'SOLV_SIZE','refdens':'REFDENS',
        #                                  'char_angle':'CHAR_ANGLE','rigid_atom_0':'RIGID_ATOM_0',
        #                                  'rigid_atom_1':'RIGID_ATOM_1','rigid_atom_2':'RIGID_ATOM_2'})
    # write all analysis options
    with open('all-settings.yaml', 'w') as f:
        f.write("# After a header which contains information whether or not (TIP3P) water or a PyConSolv solvent was used,\n")
        f.write("# two general blocks 'gist' and 'plotting' are given.\n")
        f.write('# For each block all settings and the default values are given and\n')
        f.write('# this file can be used directly after filling out the two required settings.\n')
        f.write('# Where it may be useful, a possible input is given behind the default value as a comment.\n')
        f.write("header:\n")
        f.write('  water: ' + str(water) + '\n')
        f.write('  tip3p: ' + str(tip3p) + '\n')
        f.write('  PyConSolv solvent: ' + str(solv_file) + '\n')
        f.write('gist:\n')
        for key in sorted(analyser.required_keys):
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
