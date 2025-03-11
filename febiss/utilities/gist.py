#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright Technische Universität Wien, Institute of Materials Chemistry, Podewitz Group
See LICENSE for details
"""

from warnings import warn
import glob
import os
import subprocess
from .io_handling import Input
from ..solvents import CASE_DICT

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
    # case 1: "TIP3P" (in case TIP3P water is used),
    # case 3: absolute path of /febiss/solvents/[chosen solvent] (in case pyconsolv solvents are used),
    # case 4: "PATH_TO_SOLVENT_FILE" (in case none of the options above apply. this leads still to quitting the program
    # as we cannot deal with user defined solvents yet.

    def __init__(self, case = 0, com = False, **kwargs):
        self._set_defaults(case, com)

        #set required keys
        self.required_keys = {'top',
                              'trajectory_file',
                              'solv_file',
                              'refdens',
                              'ref_eww',
                              'solv_abb',
                              'rigid_atom_0',
                              'rigid_atom_1',
                              'rigid_atom_2'
                              }


        for key in self.required_keys:
            if key not in kwargs.keys():
                raise MissingSettingException(
                    "The setting " + str(key) + " is required for the GIST analysis, but was not given.")
        self.__dict__.update((key, kwargs[key]) for key in self.required_keys)

        #set allowed keys
        self.allowed_keys = {'frame_selection',
                             'grid_center',
                             'grid_spacing',
                             'grid_lengths',
                             'solute_residues',
                             'gist_cpptraj_command_file',
                             'gist_out_file',
                             'rdf',
                             'rdf_names',
                             'gist_grid_file'
                             }


        self.__dict__.update((k, v) for k, v in kwargs.items() if k in self.allowed_keys)


        if com:
            self.__dict__['rigid_atom_0'] = -1

    def perform_solute_write_out(self):
        solute_in = self._write_solute_pdb(self.top, self.trajectory_file, self.solv_abb)
        self._execute_cpptraj(solute_in)
        if not os.path.exists('solute.pdb'):
            raise UnsuccessfulAnalysisException(
                "'solute.pdb' is not present, writing out the solute from the simulation did not work.")


    def perform_gist_analysis(self):
        if os.path.exists('febiss.dat'):
            warn("WARNING: 'gist-output.dat' is already present in directory.")

            while True:
                if Input("Do you want to analyze the trajectory again? [y/n]").yn():
                    break

                else:
                    print("Skipping analysis and reading existing data.")
                    return

        if self.rdf:
            self.perform_rdf_analysis()
        self._write_gist_cpptraj_file()
        self._write_gist_input_line()
        self._execute_cpptraj(self.gist_cpptraj_command_file)

        if not os.path.exists('febiss.dat'):
            raise UnsuccessfulAnalysisException(
                "'febiss.dat' is not present, the CPPTRAJ analysis did not work.")



    def perform_rdf_analysis(self):
        for key, value in self.rdf_names.items():
            self._write_gist_cpptraj_file()
            self._write_rdf_input_line(value, key)
            self._execute_cpptraj(self.gist_cpptraj_command_file)


    def _set_defaults(self,case,com):
        self.case = case  # possible values: 1: "TIP3P", 3: abs path of pyconsolvsolvent, 4: "PATH_TO_SOLVENT_FILE".
        self.top = None
        self.trajectory_file = None
        self.com = com
        self.solv_file = CASE_DICT[case] #stores the name of the solvent file
        self.solv_abb = None #stores the abbreviation of the solvent molecule
        #self.solv_size = 3
        self.rigid_atom_0 = 0
        self.rigid_atom_1 = 1
        self.rigid_atom_2 = 2
        self.ref_eww = None
        self.frame_selection = None
        self.grid_center = None
        self.grid_spacing = 0.5
        self.grid_lengths = (60, 60, 60)
        self.refdens = None
        self.quatfile = 'gist-quats.dat'

        self.nsolvent = 0 #number of solvents used in the simulation
        self.nframes = 0 #number of simulation frames
        self.solute_residues = ':1'
        self.gist_cpptraj_command_file = 'gist.in'
        self.febiss_cpptraj_command_file = 'febiss.in'
        self.gist_out_file = 'gist-output.dat'
        self.gist_grid_file = 'gist_grid.xyz'
        self.rdf = False
        self.rdf_names = {'center2': 'center of solute', 'C': 'carbon', 'O': 'oxygen', 'N': 'nitrogen', 'P': 'phosphor'}

    def _write_gist_cpptraj_file(self): #this creates the content for gist.in.
        self._sanity_check()
        with open(self.gist_cpptraj_command_file, 'w') as f:
            f.write('parm ' + self.top + '\n')
            f.write('trajin ' + self.trajectory_file)
            if self.frame_selection is not None and self.frame_selection.lower() != 'none':

                # if not given CPPTRAJ uses all frames
                f.write(' ' + self.frame_selection)
            f.write('\n')

            if self.case == 3:
                f.write('solvent ' + f':{self.solv_abb}\n')
            f.write('center ' + self.solute_residues + ' origin\n')
            f.write('image origin center familiar\n')


    def _sanity_check(self):
        # TODO: typedict for other keys?
        if not os.path.exists(self.top):
            raise InvalidInputException('The given top file does not exist')
        elif len(glob.glob(self.trajectory_file)) == 0:
            raise InvalidInputException('The given trajectory name or format is invalid')

    def _write_gist_input_line(self):
        with open(self.gist_cpptraj_command_file, 'a') as f:
            f.write('gist ')
            if self.grid_center is not None and self.grid_center.lower() != 'none':
                # if not given CPPTRAJ uses origin
                f.write('gridcntr ')
                self._write_variable(f, self.grid_center)
            f.write('griddim ')
            self._write_variable(f, self.grid_lengths)
            f.write('gridspacn ' + str(self.grid_spacing) + ' ')
            f.write('refdens ' + str(self.refdens) + ' ')
            f.write('rigid_idx ' + str(self.rigid_atom_0) + ' '
                    + str(self.rigid_atom_1)  + ' '
                    + str(self.rigid_atom_2) + ' ')
            f.write('out ' + self.gist_out_file + ' ')
            f.write('quat ')
            if not self.com:
                f.write('nocom ')
            f.write('norm ')
            #f.write('dx\n')
            f.write('febiss\n')  # enables febiss placement in cpptraj
            f.write('run\n')
            f.write('center :1 origin\n')
            f.write('strip :{0}\n'.format(self.solv_abb))
            f.write('trajout solute.pdb onlyframes 1\n')
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
        with open(self.gist_cpptraj_command_file, 'a') as f:
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

    def _execute_cpptraj(self,file):
        if 'CPPTRAJ_BIN' not in SETTINGS:
            raise MissingCpptrajException(
                "The path to the CPPTRAJ binary is not given in the global settings. Execute 'febiss_setup' first")
        subprocess.call([SETTINGS['CPPTRAJ_BIN'], '-i', file])

    def _write_out_gist_grid(self): #TODO: Needs to be reworked.
        gist_data = open(self.gist_out_file, 'r').readlines()
        solute_elements = []
        solute_atoms = []

        with open('febiss-solvents.pdb', 'r') as f:
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
                f.write('H\t' + row[1] + '\t' + row[2] + '\t' + row[3] + '\n') #TODO: Needs to be reworked.

    def _write_solute_pdb(self, top : str, trajin : str, abb : str) -> str: #trajin is passed as trajin+format
        file_in = "solute.in" #a pdb file is not created here but only a .in file for cpptraj
        file_out = "solute.pdb"
        if os.path.isfile(file_out):
            if not Input("\n\nDo you want to overwrite solute.pdb in your directory?").yn():
                file_out = input("\n\nPlease provide a name for the solute pdb file "
                                 "('.pdb' is automatically appended to your input):\n")

        with open(file_in, 'w') as f:
            f.write('parm ' + top + '\n')
            f.write('trajin ' + trajin + '\n')
            f.write('strip :' + abb + '\n')
            f.write('trajout ' + file_out + ' onlyframes 1\n')
            f.write('run\n')
            f.write('quit\n')

        return file_in