#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright Technische Universität Wien, Institute of Materials Chemistry, Podewitz Group
See LICENSE for details
"""

from shutil import which
import subprocess
import sys

import quaternion
import numpy as np

from Febiss.structures.solvent import Solvent
from Febiss.structures.reference import Reference
from Febiss.structures.solute import Solute
from Febiss.utilities.write_settings_file import write_settings_file
from Febiss.utilities.write_style_file import write_style_file
from Febiss.cpptraj_interface.read_settings import read_settings
from Febiss.utilities.input import Input


def help_message():
    print('\nThis program takes settings from a yaml file given as command line argument.')
    print("To obtain a suitable file execute 'get_febiss_settings'.")
    print("This program can both facilitate the execution of the GIST analysis via CPPTRAJ "
          "and display and select results from the FEBISS analysis within the GIST analysis in CPPTRAJ.")
    print("To execute CPPTRAJ, it has to be first setup via 'setup_febiss'.\n")
    sys.exit()

def read_data(febiss_file, solute: Solute, solvent: Solvent):
    """
    Reads in febiss file (default febiss.dat) and solute_file (default: solute.pdb) and assigns the values to passed Solute and
    Solvent objects

    :param febiss_file: Contains name of the febiss datafile from Analysis_FEBISS.h/cpp. Default: febiss.dat
    :param solute: Solute object
    :param solvent: Solvent object
    """
    # deal with febiss.dat
    with open(febiss_file, 'r') as f: #assumed format: no header, no footer, data: voxel, x, y, z, energy
        data = f.readlines()
        for line in data:
            line_entries = line.split()
            line_entries = [float(entry) for entry in line_entries]
            line_entries[-1] = -line_entries[-1]
            solvent.data.append(tuple(line_entries))
    solvent.sort_by_energy() #assuming GIST worked and one voxel is occupied not more than once
    solvent.get_coord_set()
    solvent.get_energy()

    #read solute.pdb
    with open('solute.pdb', 'r') as f: #assumed format of ATOM line: ATOM, ordinal number, label, residue, number, x, y, z, 1.00, energy, element
        for line in f:
            if 'ATOM' in line:
                solute.data.append(line.split())
    solute.get_coord_set()
    solute.get_elements()

    with open('gist-quats.dat', 'r') as f: #assumed format: header: 3 rows, data: voxel xcoord ycoord zcoord w x y z (several rows per voxel)
        data = f.readlines()[3:]
        for line in data: #solvent.quats is prepared in febiss
            array = np.array([float(line.split()[-4]),float(line.split()[-3]),float(line.split()[-2]),float(line.split()[-1])])
            solvent.quats[int(line.split()[0])].append(quaternion.from_float_array(array))

def get_settings():
    """
    Function which loads the settings for the FEBISS analysis and the subsequent plotting either from a yaml file or
    directly from the GUI.
    :return: A GistAnalyser and Plot object.
    """

    # GUI path
    if len(sys.argv) == 1:
        from Febiss.cli.febiss_settings import gui
        analyser, display = gui(1)
        febiss_file = write_settings_file(analyser, display)
        analyser, display = read_settings(febiss_file) # Reread to be consistent with the way FEBISS works via the CLI.


    # Help message in case of wrong usage
    elif len(sys.argv) > 2:
        analyser, display = None, None
        help_message()

    # CLI path
    else:
        febiss_file = sys.argv[1]
        if febiss_file.lower() in ['-h', '--help']:
            help_message()
        analyser, display = read_settings(febiss_file)

    return analyser, display


def main():
        analyser, display = get_settings()
        if (analyser or display) is None:
            sys.exit(print("Something went wrong loading the settings for the FEBISS analysis."))

        # instantiate solvent and solute and reference

        analyser._nvoxels = int(analyser.grid_lengths.strip('()').split(',')[0])*\
                            int(analyser.grid_lengths.strip('()').split(',')[1])*\
                            int(analyser.grid_lengths.strip('()').split(',')[2])

        solvent = Solvent(nvoxels=analyser._nvoxels, ref_evv=analyser.ref_evv)
        solute = Solute()
        reference = Reference(analyser.com, analyser.solv_abb, analyser.solv_file, analyser.rigid_atom_0, analyser.rigid_atom_1, analyser.rigid_atom_2)

        analyser.perform_rdf_analysis()
        analyser.perform_gist_analysis()
        analyser.perform_gist_grid_write_out()

        read_data(display.febiss_file, solute, solvent)

        filename = display.gui(analyser, solute, solvent, reference)

        if which('pymol') is not None:
            style_file = write_style_file()
            subprocess.call(['pymol', filename, style_file])

        else:
            if Input("PyMol not found. Do you still want to print a style file for PyMol? [y/n]").yn():
                write_style_file()

        print('FEBISS ended successfully')
        sys.exit()


if __name__ == '__main__':
    main()
