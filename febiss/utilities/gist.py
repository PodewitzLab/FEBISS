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
from ..utilities.write_solute_pdb import write_solute_pdb
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
    # major change 14 Nov 2023: changed from passing 3 arguments (tip3p, water, pyconsolv) to passing of 1 arg (solvent)
    # which can have the following values which stand for a string defined in CASE_DICT:
    # 1: "TIP3P" (in case TIP3P water is used),
    # 2: "PATH_TO_WATER_FILE" (in case non-TIP3P water is used),
    # 3: absolute path of /febiss/solvents/[chosen solvent] (in case pyconsolv solvents are used),
    # 4: "PATH_TO_SOLVENT_FILE" (in case none of the options above apply. this leads still to quitting the program as
    # we cannot deal with user defined solvents yet.

    def __init__(self, case = 0, com = False, **kwargs):
        self._set_defaults(case, com)

        #set required keys
        self.required_keys = {'top', 'trajectory_file', 'solv_file'} #new LM20231128: also case 1 needs solv_file (TP3.xyz in febiss/solvents). #changed to trajectory_file since trajectory_format is not used anymore LM20231115
        if case in [2, 3]: #also in case two with a custom water file, one has to define the rigid atoms.
            self.required_keys.update(['refdens',
                                       'solv_abb',
                                       'rigid_atom_0',
                                       'rigid_atom_1',
                                       'rigid_atom_2'
                                       #'char_angle'
                                       ])

        for key in self.required_keys:
            if key not in kwargs.keys():
                raise MissingSettingException(
                    "The setting " + str(key) + " is required for the GIST analysis, but was not given.")
        self.__dict__.update((key, kwargs[key]) for key in self.required_keys)

        #set allowed keys
        self.allowed_keys = {'frame_selection', 'grid_center', 'grid_spacing', 'grid_lengths',
                             'refdens', 'solute_residues', 'gist_cpptraj_command_file', 'gist_out_file', 'rdf', 'rdf_names',
                             'gist_grid_file', 'solv_abb', 'solv_file', 'rigid_atom_0', 'rigid_atom_1', 'rigid_atom_2' #'char_angle', 'trajectory_format',
                             } # changed back to 3 possible rigid atoms since the user shall be able to choose whether to use the com or a central atom LM20231116
                               # #only 2 rigid atoms since COM will be used
        if case in [2, 3]:
            self.allowed_keys.difference_update([
                'refdens',
                'solv_abb',
                'rigid_atom_0',
                'rigid_atom_1',
                'rigid_atom_2'
                 #'char_angle'
                 ])

        self.__dict__.update((k, v) for k, v in kwargs.items() if k in self.allowed_keys)

        #if case == 1 and com:
        #    self.__dict__['rigid_atom_0'] = 'COM'

        if com: #LM20231206: No distinguish between cases.
            self.__dict__['rigid_atom_0'] = -1

    def perform_solute_write_out(self): #DEPRECATED LM20231124. Solute.pdb is now created after the gist analysis. The commands for cpptraj are now part of the gist.in file.
        #new: executes cpptraj to write out solute.pdb from solute.in which is created with utilities.write_solute_pdb.py
        solute_in = write_solute_pdb(self.top, self.trajectory_file, self.solv_abb)
        self._execute_cpptraj(solute_in)
        if not os.path.exists('solute.pdb'):
            raise UnsuccessfulAnalysisException(
                "'solute.pdb' is not present, writing out the solute from the simulation did not work.")


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
        #self._write_out_gist_grid()


    def perform_rdf_analysis(self):
        for key, value in self.rdf_names.items():
            self._write_gist_cpptraj_file()
            self._write_rdf_input_line(value, key)
            self._execute_cpptraj(self.gist_cpptraj_command_file)

    def perform_febiss_analysis(self): #TODO: take different value in all-settings.yaml for febiss-file into account
        #following if-branch adopted from perform_gist_analysis(self). LM20231005
        if os.path.exists('febiss.dat'): #changed from febiss-solvents.pdb LM20231123
            warn("WARNING: 'febiss.dat' is already present in directory.") #changed from febiss-solvents.pdb LM20231123
            while True:
                if Input("Do you want to analyze the trajectory again? [y/n]").yn():
                    break
                else:
                    print("Skipping analysis and reading existing data.")
                    return
        self._write_febiss_cpptraj_file()
        self._execute_cpptraj(self.febiss_cpptraj_command_file)
        #following if branch adopted from perform_gist_analysis(self). LM20231005
        if not os.path.exists('febiss.dat'): #new LM20231122: changed filename to febiss-solvents.pdb
            raise UnsuccessfulAnalysisException(
                "'febiss.dat' is not present, the CPPTRAJ analysis did not work.")


    def _set_defaults(self,case,com):
        self.case = case  # possible values: 1: "TIP3P", 2: "PATH_TO_WATER_FILE", 3: abs path of pyconsolvsolvent, 4: "PATH_TO_SOLVENT_FILE".
        self.top = None
        self.trajectory_file = None #changed from trajectory_name LM20231115
        self.com = com
        self.solv_file = CASE_DICT[case] #stores the name of the solvent file
        self.solv_abb = None #stores the abbreviation of the solvent molecule
        #self.solv_size = 3
        self.rigid_atom_0 = 1 #TODO: Change to 1 as default? LM20231128. Done 20231206
        self.rigid_atom_1 = 2 #TODO: Change to 2 as default? LM20231128. Done 20231206
        self.rigid_atom_2 = 3 #TODO: Change to 3 as default? LM20231128. Done 20231206
        #changed back to 3 rigid atoms since the user shall be able to choose whether COM is used or not LM20231116 #only two rigid atoms since COM will be used
        self.frame_selection = None
        self.grid_center = None
        #self.trajectory_format = "cdf" #not used anymore LM20231115
        self.grid_spacing = 0.5
        self.grid_lengths = (60, 60, 60)
        if self.case in [1, 2]:
            self.refdens = 0.0334 #LM20231207: changed from 0.0329 # tip3p
            #self.char_angle not needed anymore since tip3p water is given as xyz file. LM20231128
            #self.char_angle = 104.52  # tip3p #changed from 104.57 according to Jorgensen et al., The Journal of Chemical Physics 1983, 79 (2), 926–935. https://doi.org/10.1063/1.445869. LM20231128
        else:
            self.refdens = None#Input("Provide reference density as float:\n", type=float)
            #self.char_angle = None#Input("Provide characteristic angle as float:\n", type=float)
        self.quatfile = 'gist-quats.dat'
        #self.occurrence = 0 #number of how often the central atom species is present in the solvent molecule
        self.nsolvent = 0 #number of solvents used in the simulation
        self.nframes = 0 #number of simulation frames
        self.solute_residues = ':1'
        self.gist_cpptraj_command_file = 'gist.in' #changes: 1) originally: cpptraj_command_file, 2) originally: cpptraj.in. LM20231005
        self.febiss_cpptraj_command_file = 'febiss.in'  # new. LM20231005
        self.gist_out_file = 'gist-output.dat'
        self.gist_grid_file = 'gist_grid.xyz'
        self.rdf = False #Changed from True to False as default. LM20231128
        self.rdf_names = {'center2': 'center of solute', 'C': 'carbon', 'O': 'oxygen', 'N': 'nitrogen', 'P': 'phosphor'} #TODO: sense?

    def _write_gist_cpptraj_file(self): #this creates the content for gist.in. renamed from _write_cpptraj_file(). LM20231005
        self._sanity_check()
        with open(self.gist_cpptraj_command_file, 'w') as f:
            f.write('parm ' + self.top + '\n')
            f.write('trajin ' + self.trajectory_file) # removed self.trajectory_format LM20231115
            if self.frame_selection is not None and self.frame_selection.lower() != 'none':
                # if not given CPPTRAJ uses all frames
                f.write(' ' + self.frame_selection)
            f.write('\n')
            if self.case in [2,3]:
                f.write('solvent ' + f':{self.solv_abb}\n')
            f.write('center ' + self.solute_residues + ' origin\n')
            f.write('image origin center familiar\n')

    def _write_febiss_cpptraj_file(self): #this creates the content for febiss.in. renamed from _write_cpptraj_file(). LM20231005
        self._find_febiss_info()
        with open(self.febiss_cpptraj_command_file, 'w') as f:
            f.write('readdata gist-population.dx\n')
            f.write('readdata gist-dTSorient_norm.dx\n')
            f.write('readdata gist-dTStrans_norm.dx\n')
            f.write('readdata gist-Esw_norm.dx\n')
            f.write('readdata gist-Eww_norm.dx\n')
            f.write('febiss refdens '+ str(self.refdens) + #'occurrence' + str(self.occurrence) +
                    ' solvnum ' + str(self.nsolvent) + ' nframes ' + str(self.nframes) + '\n')
            f.write('run')

    def _sanity_check(self):
        # TODO: typedict for other keys?
        if not os.path.exists(self.top):
            raise InvalidInputException('The given top file does not exist')
        elif len(glob.glob(self.trajectory_file)) == 0:
            raise InvalidInputException('The given trajectory name or format is invalid')

    def _find_febiss_info(self):  # assumption: (occurrence, not anymore LM20231122), nsolvents, nframes are in self.quatfile -> 2nd row
        with open(self.quatfile, 'r') as f:
            line = f.readlines()[1].split()
            #self.occurrence = line[0]
            help1 = line[0].split("=")
            self.nsolvent = help1[1]
            help2 = line[1].split("=")
            self.nframes = help2[1]

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
                    + str(self.rigid_atom_2) + ' ') #LM20231116: changed back to 3 rigid atoms since the user shall be able to decide whether to use COM or a central atom #only 2 rigid atoms will be used
            f.write('out ' + self.gist_out_file + ' ')
            f.write('quat ')
            f.write('norm\n')
            #f.write('dx\n')
            #f.write('febiss ' + str(self.temp) + '\n')  # enables febiss placement in cpptraj
            f.write('run\n')
            f.write('center :1 origin\n') #new LM20231124
            f.write('strip :{0}\n'.format(self.solv_abb)) #new LM20231124
            f.write('trajout solute.pdb onlyframes 1\n') #new LM20231124
            f.write('run') #new LM20231124

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

    def _execute_cpptraj(self,file): #changed from (self) to (self,file)
        if 'CPPTRAJ_BIN' not in SETTINGS:
            raise MissingCpptrajException(
                "The path to the CPPTRAJ binary is not given in the global settings. Execute 'febiss_setup' first")
        subprocess.call([SETTINGS['CPPTRAJ_BIN'], '-i', file]) #changed from self.cpptraj_command_file to file

    def _write_out_gist_grid(self): #according to manual this file serves to take a look if the gist grid is big enough especially for big solutes. TODO: Needs to be reworked.
        gist_data = open(self.gist_out_file, 'r').readlines()
        solute_elements = [] #contains the
        solute_atoms = [] #contains the coordinates
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
