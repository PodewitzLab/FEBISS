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

import quaternion
import yaml
import numpy as np

#import ..utilities.gist
from ..solvents import SOLVENT_LIST, RIGID_ATOMS_DICT, RIGID_ATOMS_DICT
from ..utilities.mol2_to_xyz import converter
from ..utilities.structures import Solute, Reference, Solvent
from ..utilities.io_handling import write_style_file, Input #,read_pdb
from ..utilities.gist import GistAnalyser
from ..plotting.display import Plot


def help_message():
    print('\nThis program takes settings from a yaml file given as command line argument.')
    print("To obtain a suitable file execute 'get_febiss_settings'.")
    print("This program can both facilitate the execution of the GIST analysis via CPPTRAJ "
          "and display and select results from the FEBISS analysis within the GIST analysis in CPPTRAJ.")
    print("To execute CPPTRAJ, it has to be first setup via 'setup_febiss'.\n")
    sys.exit()

def read_data(febiss_file, solute: Solute, solvent: Solvent): #TODO: pass solute.pdb as variable to avoid hardcoding problems
    """
    Reads in febiss file (default febiss.dat) and solute_file (default: solute.pdb) and assigns the values to passed Solute and
    Solvent objects

    :param febiss_file: Contains name of the febiss datafile from Analysis_FEBISS.h/cpp. Default: febiss.dat
    :param solute: Solute object
    :param solvent: Solvent object
    """
    # deal with febiss.dat
    with open(febiss_file,'r') as f: #assumed format: no header, no footer, data: voxel, x, y, z, energy
        data = f.readlines()
        for line in data:
            solvent.data.append(tuple(line.split()))
    solvent.sort_by_energy() #TODO: catch case where several solvents occupy the same voxel which should not be the case
    solvent.get_coord_set()
    solvent.get_energy()
    print("\nSolvent coords: \n")
    print(solvent.coords)

    #read solute.pdb
    with open('solute.pdb','r') as f: #assumed format of ATOM line: ATOM, ordinal number, label, residue, number, x, y, z, 1.00, energy, element
        for line in f:
            if 'ATOM' in line:
                #solute.elements.append(row[-1]) #not needed LM20231123
                #solute.atoms.append(np.array([float(r) for r in row[-6:-3]])) #not needed LM20231123
                #solute.values.append(0.0) #not needed LM20231123
                solute.data.append(line.split())
    solute.get_coord_set()
    solute.get_elements()
    print("\nSolute coords: \n")
    print(solute.coords)
    print("\nSolute elements: \n")
    print(solute.elements)

    #solute.atoms = np.asarray(solute.atoms) #not needed right now, LM20231123
    #solute.determine_polar_hydrogen_and_non_hydrogen() #not needed right now, LM20231123
    #solvent.atoms = np.asarray(solvent.atoms) #maybe needed LM20231123

    #read gist-quats.dat
    with open('gist-quats.dat','r') as f: #assumed format: header: 3 rows, data: voxel xcoord ycoord zcoord w x y z (several rows per voxel)
        data = f.readlines()[3:]
        for line in data: #solvent.quats is prepared in febiss
            array = np.array([float(line.split()[-4]),float(line.split()[-3]),float(line.split()[-2]),float(line.split()[-1])])
            solvent.quats[int(line.split()[0])].append(quaternion.from_float_array(array))








def main():
    if len(sys.argv) != 2:
        help_message()
    arg = sys.argv[1]
    if arg.lower() in ['-h', '--help']:
        help_message()
    with open(arg, 'r') as yamlfile:
        param = yaml.safe_load(yamlfile.read().replace("\t", "  ").replace("    ", "  "))

    #path = ..solvents.__path__[0]  # +"/{0}".format(solv_file) #commented out because path is passed to converter
    #converter(path, solv_abb)

    # check for header
    #if "header" not in param.keys():

        #water = Input("\n\nCould not find a header in the yaml file. Is water your main solvent? [y/n]\n").yn()
        #if water:
        #    tip3p = Input("\n\nDo you use TIP3P water?\n").yn()
        #    pyconsolv = False
        #else:
        #    tip3p = False
        #    pyconsolv = Input("\n\nDid you use one of these solvents from the PyConSolv package "
        #                      "(https://github.com/PodewitzLab/PyConSolv/tree/main/src/PyConSolv/solvents)?:\n"
        #                      + "   ".join(SOLVENT_LIST)).yn()

    if "header" in param.keys():
        case = param['header']['case']
        com = param['header']['com']
    else:
        case = 0
        com = False
        quit("Something is very wrong with the all-settings.yaml file. Before running febiss, run febiss_settings"
             "to create the all-settings.yaml file which contain crucial informations to perform a gist analysis.")

    # check for gist analysis
    if "gist" in param.keys():

        analyser = GistAnalyser(case, com, **param["gist"])

        # solvent = Solvent(analyser.solv_abb,
        #                   analyser.pyconsolv,
        #                   analyser.rigid_atom_0,
        #                   analyser.rigid_atom_1,
        #                   analyser.rigid_atom_2)
    else:
        analyser = GistAnalyser(case, com)
        quit("gist parameters are missing in the all-settings.yaml file. Before running febiss, run febiss_settings"
             "to create the all-settings.yaml file which contain crucial informations to perform a gist analysis.")

    # instantiate solvent and solute and reference

    grid_length_0 = int(analyser.grid_lengths.strip('()').split(',')[0]) #TODO: Error prone if grind_lengths is not given as '(x,y,z)'
    grid_length_1 = int(analyser.grid_lengths.strip('()').split(',')[1])
    grid_length_2 = int(analyser.grid_lengths.strip('()').split(',')[2])
    nvoxels = grid_length_0*grid_length_1*grid_length_2

    solvent = Solvent(nvoxels=nvoxels)
    solute = Solute()
    reference = Reference(case, com, analyser.solv_file, analyser.solv_abb, analyser.rigid_atom_0, analyser.rigid_atom_1, analyser.rigid_atom_2) #new: pass all 3 rigid_atoms to account for case 1.
    #analyser.perform_solute_write_out() #DEPRECATED LM20231124. #new: instead of direct calling write_solute_pdb.py
    analyser.perform_gist_analysis()
    analyser.perform_febiss_analysis()

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
    # else:
    #     if not water:
    #         solvent = Solvent(Input('Could not find a "GIST" block in the yaml file. Since you are not using water as main solvent,'
    #                                 'please specify the following parameters: \nName of solvent topology file: ',type=str),
    #                           Input('3 character abbreviation of the used solvent: ', type=str),
    #                           Input('How many atoms does your solvent molecule contain?: ',type=int),
    #                           Input('Please specify "rigid_atom_0: ', type=str),
    #                           Input('Please specify "rigid_atom_1: ', type=str),
    #                           Input('Please specify "rigid_atom_2: ', type=str))


    # now check for plotting, else assume that only plotting options are given directly #
    for key in param.keys():
        if key.lower() == 'plotting':
            param = param['plotting']
            break

    febiss_file = param.get('febiss_file', 'febiss.dat') #not needed#changed to febiss.dat
    #read_pdb(febiss_file, solute, solvent) #new reading in procedure will be used LM20231027
    read_data(febiss_file, solute, solvent)
    display = Plot(**param)
    filename = display.gui(analyser.solv_abb, solute, solvent, reference)
    if which('pymol') is not None:
        style_file = param.get('style_file', write_style_file())
        subprocess.call(['pymol', filename, style_file])
    print('FEBISS ended successfully')
    sys.exit()

if __name__ == '__main__':
    main()
