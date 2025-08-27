#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright Technische Universität Wien, Institute of Materials Chemistry, Podewitz Group
See LICENSE for details
"""

import sys
import os
from ..solvents import NUMBER_DICT, SOLVENT_DICT
from ..utilities.colorgen import Color
from Febiss.utilities.input import Input
from Febiss.utilities.write_settings_file import write_settings_file


def help_message():
    print("\nThis program writes all possible settings of the GIST analysis "
          "and the plotting into 'all-settings.yaml'.")
    print('It requires no arguments, however, the topology and trajectory file can be passed (in this order)'
          'for faster usability.')
    sys.exit()

def pyconsolv_name():
    """
    Function used if PyConSolvents are used (contained in this package in /febiss/solvents/). PyConSolv solvent will
    be chosen according to the user's input.

    :return: Name of the PyConSolv solvent used for this simulation to be analyzed.
    """

    solvent_list = ["{0}: {1}".format(i, j) for i, j in NUMBER_DICT.items()]

    solv_num = input(
                "\n\nPlease enter the number of the solvent you used:\n-"
                + "\n-".join(list(solvent_list)) + "\n\nYou used: ")

    while solv_num not in list(NUMBER_DICT.keys()):
        solv_num = input(
            "\n\nPlease enter the number of the solvent the solvent you used. "
            "Write 'quit' to quit program or 'list' to see solvent list again:\n")
        if solv_num.lower() == "quit":
            quit("\n---Quit as you wished!---")
        if solv_num.lower() == "list":
            solv_num = input("\n\nHere is the list of PyConSolv solvents:\n-" + "\n-".join(
                solvent_list) + "\n\nEnter the 3 character long string of the solvent you used: ")
    return NUMBER_DICT[solv_num] #name of the solvent

def gui(mode: int):
    from Febiss.gui.gui_behavior import Window
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

    if len(sys.argv) > 3 or len(sys.argv) == 2:
        help_message()
    elif len(sys.argv) == 3:
        top = sys.argv[1]
        traj = sys.argv[2]

    #defining the solvent name
    solv_name = None

    if Input("\nDo you use TIP3P water? [y/n]").yn():
        solv_name = 'Water'

    elif Input("\nDid you use one of these solvents from the PyConSolv package "
               "(https://github.com/PodewitzLab/PyConSolv/tree/main/src/PyConSolv/solvents)? [y/n]").yn():
        solv_name = pyconsolv_name()

    else:
        print(Color.RED + "\n\tCustom solvent. Please specify the solvent you used in all-settings.yaml.\033[0m\n")
        solv_name = 'Custom'

    #defining if com shall be used for binning
    com = False
    if Input("\nDo you want to use the center of mass (COM) as central point of the solvent molecule? [y/n]").yn():
        com = True

    if solv_name == 'Custom':
        solv_file = 'Choose File ...'
    else:
        solv_file = os.path.abspath(os.path.join(__file__, "../../solvents/{0}".format(SOLVENT_DICT[solv_name][1])))

    from Febiss.cpptraj_interface.gist import GistAnalyser
    from ..plotting.display import Plot
    # the arguments are necessary in the init, otherwise exception
    analyser = GistAnalyser(**{
        'top': top,
        'trajectory_file': traj,
        'solv_abb': SOLVENT_DICT[solv_name][0],
        'com' : com,
        'solv_file': solv_file,
        'refdens': SOLVENT_DICT[solv_name][3],
        'ref_evv': SOLVENT_DICT[solv_name][4],
        'rigid_atom_0': SOLVENT_DICT[solv_name][2][1],
        'rigid_atom_1': SOLVENT_DICT[solv_name][2][0],
        'rigid_atom_2': SOLVENT_DICT[solv_name][2][2]
    })


    display = Plot()

    return analyser, display


def main():
    if Input("\nDo you want to use the GUI? [y/n]").yn():
        analyser, display = gui(0)

    else:
        analyser, display = cli()

    write_settings_file(analyser, display)
    sys.exit()



if __name__ == '__main__':
    main()
