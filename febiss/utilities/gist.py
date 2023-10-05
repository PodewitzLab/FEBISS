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
import sys
from .io_handling import Input

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
    def __init__(self, water=False, tip3p=False, pyconsolv=False, **kwargs): #TODO: change passing of 3 args to passing of 1 int arg e.g. 1 = tip3p, 2 = water, 3 = pyconsolv, 4 = random
        self._set_defaults(water,tip3p,pyconsolv)
        self.required_keys = {'top', 'trajectory_name'}
        if not self.tip3p:
            self.required_keys.update(['refdens',
                                       #'char_angle'
                                       ])
        if not self.water and not self.pyconsolv:
            self.required_keys.update(
                ['rigid_atom_0','rigid_atom_1','rigid_atom_2'])
        for key in self.required_keys:
            if key not in kwargs.keys():
                raise MissingSettingException(
                    "The setting " + str(key) + " is required for the GIST analysis, but was not given.")
        self.__dict__.update((key, kwargs[key]) for key in self.required_keys)
        self.allowed_keys = {'frame_selection', 'grid_center', 'trajectory_format', 'grid_spacing', 'grid_lengths',
                             'refdens', 'solute_residues', 'cpptraj_command_file', 'gist_out_file', 'rdf', 'rdf_names',
                             'gist_grid_file', 'char_angle', 'solv_abb', 'solv_top', 'rigid_atom_0','rigid_atom_1','rigid_atom_2'}
        self.__dict__.update((k, v) for k, v in kwargs.items() if k in self.allowed_keys)
        if not self.tip3p:
            self.allowed_keys.difference_update([
                'refdens',
                'solv_abb',
                 #'char_angle'
                 ])
        if not self.water and self.tip3p:
            self.allowed_keys.difference_update([
                #'solv_top',
                #'solv_size',
                'rigid_atom_0','rigid_atom_1','rigid_atom_2'
            ])



    def perform_gist_analysis(self):
        #changed from checking for febiss-solvents.pdb to gist-output.dat
        if os.path.exists('gist-output.dat'):
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
        # changed from febiss-solvents.pdb to gist-output.dat
        if not os.path.exists('gist-output.dat'):
            raise UnsuccessfulAnalysisException(
                "'gist-output.dat' is not present, the CPPTRAJ analysis did not work.")
        self._write_out_gist_grid()


    def perform_rdf_analysis(self):
        for key, value in self.rdf_names.items():
            self._write_gist_cpptraj_file()
            self._write_rdf_input_line(value, key)
            self._execute_cpptraj(self.gist_cpptraj_command_file)

    def perform_febiss_analysis(self):
        #following if branch adopted from perform_gist_analysis(self). LM20231005
        if os.path.exists('febiss-solvents.pdb'):
            warn("WARNING: 'febiss-solvents.pdb' is already present in directory.")
            while True:
                if Input("Do you want to analyze the trajectory again? [y/n]").yn():
                    break
                else:
                    print("Skipping analysis and reading existing data.")
                    return
        self._write_febiss_cpptraj_file()
        self._execute_cpptraj(self.febiss_cpptraj_command_file)
        #following if branch adopted from perform_gist_analysis(self). LM20231005
        if not os.path.exists('febiss-solvents.pdb'):
            raise UnsuccessfulAnalysisException(
                "'febiss-solvents.pdb' is not present, the CPPTRAJ analysis did not work.")


    def _set_defaults(self,water,tip3p,pyconsolv):
        self.top = None
        self.trajectory_name = None
        self.water = water
        self.tip3p = tip3p # set to true only if TIP3P water is used
        self.pyconsolv = pyconsolv
        self.solv_top = None #stores the name of the solvent topology
        self.solv_abb = None #stores the abbreviation of the solvent molecule
        self.solv_size = 3
        self.rigid_atom_0 = 'O'
        self.rigid_atom_1 = 'H'
        self.rigid_atom_2 = 'H'
        self.frame_selection = None
        self.grid_center = None
        self.trajectory_format = "cdf"
        self.grid_spacing = 0.5
        self.grid_lengths = (60, 60, 60)
        if self.water:
            self.refdens = 0.0329  # tip3p
            self.char_angle = 104.57  # tip3p
        else:
            self.refdens = None#Input("Provide reference density as float:\n", type=float)
            self.char_angle = None#Input("Provide characteristic angle as float:\n", type=float)
        self.quatfile = 'gist-quats.dat'
        self.occurrence = 0 #number of how often the central atom species is present in the solvent molecule
        self.nsolvent = 0 #number of solvents used in the simulation
        self.nframes = 0 #number of simulation frames
        self.solute_residues = ':1'
        self.gist_cpptraj_command_file = 'gist.in' #changes: 1) originally: cpptraj_command_file, 2) originally: cpptraj.in. LM20231005
        self.febiss_cpptraj_command_file = 'febiss.in'  # new. LM20231005
        self.gist_out_file = 'gistout.dat'
        self.gist_grid_file = 'gist_grid.xyz'
        self.rdf = True
        self.rdf_names = {'center2': 'center of solute', 'C': 'carbon', 'O': 'oxygen', 'N': 'nitrogen', 'P': 'phosphor'} #TODO: sense?

    def _write_gist_cpptraj_file(self): #this creates the content for gist.in. renamed from _write_cpptraj_file(). LM20231005
        self._sanity_check()
        with open(self.gist_cpptraj_command_file, 'w') as f:
            f.write('parm ' + self.top + '\n')
            f.write('trajin ' + self.trajectory_name + '*' + self.trajectory_format)
            if self.frame_selection is not None and self.frame_selection.lower() != 'none':
                # if not given CPPTRAJ uses all frames
                f.write(' ' + self.frame_selection)
            f.write('\n')
            #if not self.water:
            #    f.write('solvent' + self.solv_top + f':{self.solv_abb}\n')
            f.write('center ' + self.solute_residues + ' origin\n')
            f.write('image origin center familiar\n')

    def _write_febiss_cpptraj_file(self): #this creates the content for gist.in. renamed from _write_cpptraj_file(). LM20231005
        self._find_febiss_info()
        with open(self.febiss_cpptraj_command_file, 'w') as f:
            f.write('readdata population.dx')
            f.write('readdata dTSorient_dens.dx')
            f.write('readdata dTStrans_dens.dx')
            f.write('readdata Esw_dens.dx')
            f.write('readdata Eww_dens.dx')
            f.write('febiss refdens '+ str(self.refdens) + 'occurrence' + str(self.occurrence) +
                    'nsolvent' + str(self.nsolvent) + 'nframes' + str(self.nframes))

    def _sanity_check(self):
        # TODO: typedict for other keys?
        if not os.path.exists(self.top):
            raise InvalidInputException('The given top file does not exist')
        elif len(glob.glob(self.trajectory_name + '*' + self.trajectory_format)) == 0:
            raise InvalidInputException('The given trajectory name or format is invalid')

    def _find_febiss_info(self):  # assumption: occurrence, nsolvents, nframes are in self.quatfile -> 2nd row
        with open(self.quatfile, 'r') as f:
            line = f.readlines()[1].split()
            self.occurrence = line[0]
            self.nsolvents = line[1]
            self.nframes = line[2]

    def _write_gist_input_line(self):
        with open(self.cpptraj_command_file, 'a') as f:
            f.write('gist')
            if self.grid_center is not None and self.grid_center.lower() != 'none':
                # if not given CPPTRAJ uses origin
                f.write('gridcntr ')
                self._write_variable(f, self.grid_center)
            f.write('griddim ')
            self._write_variable(f, self.grid_lengths)
            f.write('gridspacn ' + str(self.grid_spacing) + ' ')
            f.write('refdens ' + str(self.refdens) + ' ')
            f.write('rigidatoms ' + str(self.rigid_atom_0) + ' ' + str(self.rigid_atom_1) + ' ' + str(self.rigid_atom_2) + ' ')
            f.write('out ' + self.gist_out_file)
            f.write('dx\n')
            #f.write('febiss ' + str(self.temp) + '\n')  # enables febiss placement in cpptraj
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

    def _execute_cpptraj(self,file): #changed from (self) to (self,file)
        if 'CPPTRAJ_BIN' not in SETTINGS:
            raise MissingCpptrajException(
                "The path to the CPPTRAJ binary is not given in the global settings. Execute 'febiss_setup' first")
        subprocess.call([SETTINGS['CPPTRAJ_BIN'], '-i', file]) #changed from self.cpptraj_command_file to file

    def _write_out_gist_grid(self):
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
                f.write('H\t' + row[1] + '\t' + row[2] + '\t' + row[3] + '\n') #grid points are represented as H nuclei TODO:change?
