#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright University Innsbruck, Institute for General, Inorganic, and Theoretical Chemistry, Podewitz Group
See LICENSE for details
"""

from warnings import warn
import glob
import os
import subprocess

from febiss import SETTINGS


class MissingSettingException(Exception):
    pass


class InvalidInputException(Exception):
    pass


class MissingCpptrajException(Exception):
    pass


class UnsuccessfulAnalysisException(Exception):
    pass


class GistAnalyser:
    def __init__(self, **kwargs):
        self._set_defaults()
        self.required_keys = {'top', 'trajectory_name'}
        for key in self.required_keys:
            if key not in kwargs.keys():
                raise MissingSettingException(
                    "The setting " + str(key) + " is required for the GIST analysis, but was not given.")
        self.__dict__.update((key, kwargs[key]) for key in self.required_keys)
        self.allowed_keys = {'frame_selection', 'grid_center', 'trajectory_format', 'grid_spacing', 'grid_lengths',
                             'refdens', 'solute_residues', 'cpptraj_command_file', 'gist_out_file', 'rdf', 'rdf_names',
                             'gist_grid_file', 'water_angle'}
        self.__dict__.update((k, v) for k, v in kwargs.items() if k in self.allowed_keys)

    def perform_gist_analysis(self):
        if os.path.exists('febiss-waters.pdb'):
            warn("WARNING: 'febiss-waters.pdb' is already present in directory.")
            while True:
                inp = input("Do you want to analyze the trajectory again? [y/n] ")
                if inp == 'y' or inp == 'Y' or inp == 'yes' or inp == 'Yes' or inp == ' y':
                    break
                elif inp == 'n' or inp == 'N' or inp == 'no' or inp == 'No' or inp == ' n':
                    print("Skipping analysis and reading existing data.")
                    return
                else:
                    print("Sorry wrong input, just 'y' or 'n'.")
        if self.rdf:
            self.perform_rdf_analysis()
        self._write_cpptraj_file()
        self._write_gist_input_line()
        self._execute_cpptraj()
        if not os.path.exists('febiss-waters.pdb'):
            raise UnsuccessfulAnalysisException(
                "'febiss-waters.pdb' is not present, the CPPTRAJ analysis did not work.")
        self._write_out_gist_grid()

    def perform_rdf_analysis(self):
        for key, value in self.rdf_names.items():
            self._write_cpptraj_file()
            self._write_rdf_input_line(value, key)
            self._execute_cpptraj()

    def _set_defaults(self):
        self.top = None
        self.trajectory_name = None
        self.frame_selection = None
        self.grid_center = None
        self.trajectory_format = "cdf"
        self.grid_spacing = 0.5
        self.grid_lengths = (60, 60, 60)
        self.refdens = 0.0329  # tip3p
        self.water_angle = 104.57  # tip3p
        self.solute_residues = ':1'
        self.cpptraj_command_file = 'cpptraj.in'
        self.gist_out_file = 'gistout.dat'
        self.gist_grid_file = 'gist_grid.xyz'
        self.rdf = True
        self.rdf_names = {'center2': 'center of solute', 'C': 'carbon', 'O': 'oxygen', 'N': 'nitrogen', 'P': 'phosphor'}

    def _write_cpptraj_file(self):
        self._sanity_check()
        with open(self.cpptraj_command_file, 'w') as f:
            f.write('parm ' + self.top + '\n')
            f.write('trajin ' + self.trajectory_name + '*' + self.trajectory_format)
            if self.frame_selection is not None and self.frame_selection.lower() != 'none':
                # if not given CPPTRAJ uses all frames
                f.write(' ' + self.frame_selection)
            f.write('\n')
            f.write('center ' + self.solute_residues + ' origin\n')
            f.write('image origin center familiar\n')

    def _sanity_check(self):
        if not os.path.exists(self.top):
            raise InvalidInputException('The given top file does not exist')
        elif len(glob.glob(self.trajectory_name + '*' + self.trajectory_format)) == 0:
            raise InvalidInputException('The given trajectory name or format is invalid')

    def _write_gist_input_line(self):
        with open(self.cpptraj_command_file, 'a') as f:
            f.write('gigist ')
            if self.grid_center is not None and self.grid_center.lower() != 'none':
                # if not given CPPTRAJ uses origin
                f.write('gridcntr ')
                self._write_variable(f, self.grid_center)
            f.write('griddim ')
            self._write_variable(f, self.grid_lengths)
            f.write('gridspacn ' + str(self.grid_spacing) + ' ')
            f.write('refdens ' + str(self.refdens) + ' ')
            f.write('febiss ')# partial fix, water angle is hardcoded to TIP3P
            #f.write('febiss ' + str(self.water_angle) + ' ')  # enables febiss placement in cpptraj
            f.write('out ' + self.gist_out_file + '\n')
            f.write('run')

    def _write_variable(self, f, variable):
        # assumes same value three times if single value
        if type(variable) == int or type(variable) == float:
            f.write(str(variable) + ' ')
            f.write(str(variable) + ' ')
            f.write(str(variable) + ' ')
        elif type(variable) == str:
            # clean up string
            variable = variable.strip('[').strip(']').strip('(').strip(')').replace(',', '')
            # assumes user to give all three values as string
            if ' ' in variable:
                f.write(variable)
                f.write(' ')  # to be sure to have separation to next gist argument
            # assumes same value three times
            else:
                f.write(variable + ' ')
                f.write(variable + ' ')
                f.write(variable + ' ')
        # assumes same value three times if single value
        elif len(variable) == 1:
            f.write(str(variable[0]) + ' ')
            f.write(str(variable[0]) + ' ')
            f.write(str(variable[0]) + ' ')
        elif len(variable) == 3:
            f.write(str(variable[0]) + ' ')
            f.write(str(variable[1]) + ' ')
            f.write(str(variable[2]) + ' ')
        else:
            raise InvalidInputException('Given Setting value ' + str(variable) + ' cannot be interpreted correctly.')

    def _write_rdf_input_line(self, name, symbol):
        with open(self.cpptraj_command_file, 'a') as f:
            f.write('radial ')
            f.write('spacing 0.05 10 ')
            f.write('density ' + str(self.refdens) + ' ')
            if symbol == 'center2':
                f.write(':WAT@O ' + self.solute_residues + ' ' + symbol + ' ')
            else:
                f.write(':WAT@O ' + self.solute_residues + '@/' + symbol + ' ')
            f.write("out 'rdf-" + name + ".dat' ")
            f.write("intrdf 'int-rdf-" + name + ".dat' ")
            f.write("rawrdf 'raw-rdf-" + name + ".dat'\n")
            f.write('run')

    def _execute_cpptraj(self):
        if 'CPPTRAJ_BIN' not in SETTINGS:
            raise MissingCpptrajException(
                "The path to the CPPTRAJ binary is not given in the global settings. Execute 'febiss_setup' first")
        subprocess.call([SETTINGS['CPPTRAJ_BIN'], '-i', self.cpptraj_command_file])

    def _write_out_gist_grid(self):
        gist_data = open(self.gist_out_file, 'r').readlines()
        solute_elements = []
        solute_atoms = []
        with open('febiss-waters.pdb', 'r') as f:
            for line in f:
                row = line.split()
                if row[3] == 'SOL':
                    solute_elements.append(row[2])
                    solute_atoms.append([row[5], row[6], row[7]])
        n_voxels = len(gist_data) - 2
        with open(self.gist_grid_file, 'w') as f:
            f.write(str(n_voxels + len(solute_elements)) + '\n\n')
            for e, a in zip(solute_elements, solute_atoms):
                f.write(e + '\t' + a[0] + '\t' + a[1] + '\t' + a[2] + '\n')
            for line in gist_data[2:]:
                row = line.split()
                f.write('H\t' + row[1] + '\t' + row[2] + '\t' + row[3] + '\n')
