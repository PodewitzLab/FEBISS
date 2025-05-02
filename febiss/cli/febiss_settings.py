#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright Technische Universität Wien, Institute of Materials Chemistry, Podewitz Group
See LICENSE for details
"""

import sys
import os
from ..utilities.io_handling.input import Input
from ..utilities.io_handling.write_settings_file import write_settings_file
from ..solvents import SOLVENT_DICT


def help_message():
    print("\nThis program writes all possible settings of the GIST analysis "
          "and the plotting into 'all-settings.yaml'.")
    print('It requires no arguments, however, the topology and trajectory file can be passed (in this order)'
          'for faster usability.')
    sys.exit()

def pyconsolv_abb():
    """
    Function used if PyConSolvents are used (contained in this package in /febiss/solvents/). PyConSolv solvent will
    be chosen according to the user's input.

    :return: 3 letter abbreviation of the PyConSolv solvent used for this simulation to be analyzed.
    """

    abb_list = [list(i)[0] for i in SOLVENT_DICT.values()]
    zip_list = list(zip(SOLVENT_DICT.keys(), abb_list))
    solvent_list = [str(j[0])+" ("+str(j[1])+")" for j in zip_list]

    solv_abb = input(
                "\n\nPlease enter the 3 character long string in the parentheses after the solvent you used:\n-"
                + "\n-".join(list(solvent_list)) + "\n\nYou used: ").upper()

    while solv_abb not in list(abb_list):
        solv_abb = input(
            "\n\nPlease enter the 3 character long string in the parentheses after the solvent you used. "
            "Write 'quit' to quit program or 'list' to see solvent list again:\n").upper()
        if solv_abb.lower() == "quit":
            quit("\n---Quit as you wished!---")
        if solv_abb.lower() == "list":
            solv_abb = input("\n\nHere is the list of PyConSolv solvents:\n-" + "\n-".join(
                solvent_list) + "\n\nEnter the 3 character long string of the solvent you used: ").upper()
    return solv_abb

def gui(mode: int):
    from ..utilities.io_handling.gui_behavior import Window
    from PyQt6.QtWidgets import QApplication

    # Create the application
    app = QApplication(sys.argv)

    # Create and show the application's main window. window instantiates a GistAnalyser object (analyser) and a Plot object (display)
    window = Window(mode)

    # Run the application's main loop
    app.exec()

    if window.analyser.signal: #yaml write out only if "OK" was pushed in GUI
        return window.analyser, window.display

    else:
        sys.exit()


def cli():
    top = 'SOLVBOX_TOPOLOGY'
    traj = 'TRAJECTORY'

    #lists for easier access
    abb_list = []
    solv_file_list = []
    rigid_atoms_list = []
    refdens_list = []
    ref_eww_list = []

    for val in SOLVENT_DICT.values():
        abb_list.append(list(val)[0])
        solv_file_list.append(list(val)[1])
        rigid_atoms_list.append(list(val)[2])
        refdens_list.append(list(val)[3])
        ref_eww_list.append(list(val)[4])


    if len(sys.argv) > 3 or len(sys.argv) == 2:
        help_message()
    elif len(sys.argv) == 3:
        top = sys.argv[1]
        traj = sys.argv[2]

    #defining the solvent abbreviation
    solv_abb = None #3 letter abbrevation of the solvent

    if Input("\nDo you use TIP3P water? [y/n]").yn():
        solv_abb = 'WAT'

    elif Input("\nDid you use one of these solvents from the PyConSolv package "
               "(https://github.com/PodewitzLab/PyConSolv/tree/main/src/PyConSolv/solvents)? [y/n]").yn():
        solv_abb = pyconsolv_abb()

    else:
        quit("\n\nFebiss in it's current state does not allow for user-defined solvents. We are working on it!\n"
             "---Quit program---\n")

    solv_idx = abb_list.index(solv_abb) #index of solv_abb in abb_list

    #defining if com shall be used for binning
    com = False
    if Input("\nDo you want to use the center of mass (COM) as central point of the solvent molecule? [y/n]").yn():
        com = True

    from ..utilities.gist import GistAnalyser
    from ..plotting.display import Plot
    # the arguments are necessary in the init, otherwise exception
    analyser = GistAnalyser(**{
        'top': top,
        'trajectory_file': traj,
        'solv_abb': solv_abb,
        'com' : com,
        'solv_file': os.path.abspath(os.path.join(__file__, "../../solvents/{0}".format(solv_file_list[solv_idx]))),
        'refdens': refdens_list[solv_idx],
        'ref_eww': ref_eww_list[solv_idx],
        'rigid_atom_0': rigid_atoms_list[solv_idx][1],
        'rigid_atom_1': rigid_atoms_list[solv_idx][0],
        'rigid_atom_2': rigid_atoms_list[solv_idx][2]
    })


    display = Plot()

    return analyser, display


def main():
    if Input("\nDo you want to use the GUI? [y/n]").yn():
        analyser, display = gui(0)

    else:
        analyser, display = cli()

    sys.exit(write_settings_file(analyser, display))



if __name__ == '__main__':
    main()
