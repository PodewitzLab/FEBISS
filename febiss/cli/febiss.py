#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright University Innsbruck, Institute for General, Inorganic, and Theoretical Chemistry, Podewitz Group
See LICENSE for details
"""

from shutil import which
import subprocess
import sys
import yaml

from ..solvents import SOLVENT_LIST, RIGID_ATOMS_DICT, RIGID_ATOMS_DICT
from ..utilities.mol2_to_xyz import converter
from ..utilities.structures import Solute, Solvent
from ..utilities.io_handling import read_pdb, write_style_file, Input
from ..plotting.display import Plot


def help_message():
    print('\nThis program takes settings from a yaml file given as command line argument.')
    print("To obtain a suitable file execute 'get_febiss_settings'.")
    print("This program can both facilitate the execution of the GIST analysis via CPPTRAJ "
          "and display and select results from the FEBISS analysis within the GIST analysis in CPPTRAJ.")
    print("To execute CPPTRAJ, it has to be first setup via 'setup_febiss'.\n")
    sys.exit()


def main():
    if len(sys.argv) != 2:
        help_message()
    arg = sys.argv[1]
    if arg.lower() in ['-h', '--help']:
        help_message()
    with open(arg, 'r') as yamlfile:
        param = yaml.safe_load(yamlfile.read().replace("\t", "  ").replace("    ", "  "))

    # instantiate solvent and solute
    #solvent = Solvent()
    solute = Solute()

    #path = ..solvents.__path__[0]  # +"/{0}".format(solv_file) #commented out because path is passed to converter
    #converter(path, solv_abb)

    # check for header
    if "header" not in param.keys():
        water = Input("\n\nCould not find a header in the yaml file. Is water your main solvent? [y/n]\n").yn()
        if water:
            tip3p = Input("\n\nDo you use TIP3P water?\n").yn()
            pyconsolv = False
        else:
            tip3p = False
            pyconsolv = Input("\n\nDid you use one of these solvents from the PyConSolv package "
                              "(https://github.com/PodewitzLab/PyConSolv/tree/main/src/PyConSolv/solvents)?:\n"
                              + "   ".join(SOLVENT_LIST)).yn()
    else:
        water = param["header"]["water"]
        tip3p = param["header"]["tip3p"]
        pyconsolv = param["header"]["pyconsolv"]

    # check for gist analysis
    if "gist" in param.keys():
        from ..utilities.gist import GistAnalyser
        analyser = GistAnalyser(water,tip3p,pyconsolv,**param["gist"])
        solvent = Solvent(analyser.solv_abb,
                          analyser.pyconsolv,
                          analyser.rigid_atom_0,
                          analyser.rigid_atom_1,
                          analyser.rigid_atom_2)
        analyser.perform_gist_analysis()
    #elif "GIST" in param.keys():
    #    from ..utilities.gist import GistAnalyser
    #    analyser = GistAnalyser(water,tip3p,**param["GIST"])
    #    solvent = Solvent(analyser.solv_top,
    #                      analyser.solv_abb,
    #                      analyser.solv_size,
    #                      analyser.rigid_atom_0,
    #                      analyser.rigid_atom_1,
    #                      analyser.rigid_atom_2)
    #    analyser.perform_gist_analysis()
    else:
        if not water:
            solvent = Solvent(Input('Could not find a "GIST" block in the yaml file. Since you are not using water as main solvent,'
                                    'please specify the following parameters: \nName of solvent topology file: ',type=str),
                              Input('3 character abbreviation of the used solvent: ', type=str),
                              Input('How many atoms does your solvent molecule contain?: ',type=int),
                              Input('Please specify "rigid_atom_0: ', type=str),
                              Input('Please specify "rigid_atom_1: ', type=str),
                              Input('Please specify "rigid_atom_2: ', type=str))

    # now check for plotting, else assume that only plotting options are given directly #

    for key in param.keys():
        if key.lower() == 'plotting':
            param = param['plotting']
            break
    #solute = Solute() #info about solute will be given in read_pdb
    #solvent = Solvent(solv_size) #info about solvent will be given in read_pdb
    febiss_file = param.get('febiss_file', 'febiss-solvents.pdb')
    read_pdb(febiss_file, solute, solvent)
    display = Plot(**param)
    filename = display.gui(solute, solvent)
    if which('pymol') is not None:
        style_file = param.get('style_file', write_style_file())
        subprocess.call(['pymol', filename, style_file])
    print('FEBISS ended successfully')
    sys.exit()

if __name__ == '__main__':
    main()
