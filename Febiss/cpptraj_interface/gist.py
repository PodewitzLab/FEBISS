#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright Technische Universität Wien, Institute of Materials Chemistry, Podewitz Group
See LICENSE for details
"""

import os
import sys
import subprocess

from Febiss.utilities import ELEMENT_DICT
from collections import OrderedDict
from Febiss.utilities.check_settings import Checker
from Febiss.utilities.input import Input
from Febiss.cpptraj_interface.write_cpptraj_files import write_gist_cpptraj_file, write_rdf_input_files
from Febiss.utilities.write_grid import write_out_gist_grid

from Febiss import SETTINGS, SETTINGS_FILE


class MissingSettingException(Exception):
    pass


class InvalidInputException(Exception):
    pass


class MissingCpptrajException(Exception):
    pass


class UnsuccessfulAnalysisException(Exception):
    pass


class GistAnalyser(Checker):

    def __init__(self, **kwargs):
        Checker.__init__(self)
        self._set_defaults()

        #set required keys
        self.required_keys = OrderedDict({
            'top': ('path', 'str'),
            'trajectory_file': ('path', 'str'),
            'solv_file': ('path', 'str'),
            'com': ('type', 'bool'),
            'refdens': ('type', 'float'),
            'ref_evv': ('type', 'float'),
            'solv_abb': ('type', 'str'),
            'rigid_atom_0': ('format', '-1|[0-9]+'),
            'rigid_atom_1': ('format', '[0-9]+'),
            'rigid_atom_2': ('format', '[0-9]+'),
        })


        for key in self.required_keys.keys():
            if key not in kwargs.keys():
                raise MissingSettingException(
                    "The setting " + str(key) + " is required for the GIST analysis, but was not given.")
        self.__dict__.update((key, kwargs[key]) for key in self.required_keys.keys())

        if self.com:
            self.__dict__['rigid_atom_0'] = -1


        #set allowed keys
        self.allowed_keys = OrderedDict({
            'frame_selection': ('format', '[0-9]+ [0-9]+ [0-9]+|(?i)None|^$'),
            'grid_center': ('format', '\([0-9]+.[0-9]+, [0-9]+.[0-9]+, [0-9]+.[0-9]+\)|(?i)None|^$'),
            'grid_spacing': ('type', 'float'),
            'grid_lengths': ('format', '\([0-9]+, [0-9]+, [0-9]+\)'),
            'solute_residues': ('format', '[a-zA-Z0-9:@\-]+'),
            'gist_cpptraj_command_file': ('format', '[a-zA-Z0-9:@\._\-]+\.in'),
            'gist_out_file': ('format', '[a-zA-Z0-9:@\._\-]+\.dat'),
            'gist_grid_file': ('format', '[a-zA-Z0-9:@\._\-]+\.xyz'),
            'rdf_elements': ('format', '[a-zA-Z,]+|(?i)None|^$'),
            'rdf_maximum':('type','float'),
            'rdf_spacing':('type','float')
        })

        self.__dict__.update((k, v) for k, v in kwargs.items() if k in self.allowed_keys.keys())

        # Set up auxiliary dictionaries
        self._rdf_names = {}

        if str(self.rdf_elements).lower() not in ['', 'none']:
            not_recognized = []
            for element in self.rdf_elements.split(','):
                try:
                    if element.lower() == 'center':
                        self._rdf_names['center']='Center'
                    else:
                        self._rdf_names[element]=ELEMENT_DICT[element.lower()][1]
                except KeyError:
                    not_recognized.append(element)
            if len(not_recognized) != 0:
                print("The following elements were not recognized, thus no RDF will be performed for them:")
                for unknown in not_recognized:
                    print("\t-{0}".format(unknown))
                print("\n")

        self._all_keys.update(self.required_keys)
        self._all_keys.update(self.allowed_keys)


    def perform_gist_analysis(self):
        if self.check_path('febiss.dat'):
            print("\nWARNING: 'febiss.dat' is already present in directory. ")

            if not Input("Do you want to analyze the trajectory again? [y/n]").yn():
                print("Skipping analysis and using existing data.")
                return

        write_gist_cpptraj_file(self)
        self._execute_cpptraj(self.gist_cpptraj_command_file)

        if not os.path.exists('febiss.dat'):
            raise UnsuccessfulAnalysisException(
                "'febiss.dat' is not present, the CPPTRAJ analysis did not work.")

    def perform_rdf_analysis(self):
        if len(self._rdf_names) != 0:
            if self.check_path('rdf*dat'):
                print("\nWARNING: Found some RDF data in directory. ")

                if not Input("Do you want to calculate the RDFs again? [y/n]").yn():
                    print("Skipping analysis and using existing data.")
                    return

            file_list = write_rdf_input_files(self)
            if len(file_list) != 0:
                for file in file_list:
                    self._execute_cpptraj(file)

    def perform_gist_grid_write_out(self):
        write_out_gist_grid(self)


    def _set_defaults(self):
        # command signal for the gui()
        self.signal = False

        # required keys
        self.top = None
        self.trajectory_file = None
        self.solv_file = "PATH_TO_SOLVENT_FILE"
        self.com = False
        self.refdens = None
        self.ref_evv = None
        self.solv_abb = None
        self.rigid_atom_0 = 0
        self.rigid_atom_1 = 1
        self.rigid_atom_2 = 2

        # allowed keys
        self.frame_selection = None
        self.grid_center = None
        self.grid_spacing = 0.5
        self.grid_lengths = (60, 60, 60)

        self.solute_residues = ':1'
        self.gist_cpptraj_command_file = 'cpptraj.in'
        self.gist_out_file = 'gist-output.dat'
        self.gist_grid_file = 'gist_grid.xyz'
        self.rdf_elements = None
        self.rdf_maximum = 10.0
        self.rdf_spacing = 0.05

        #hidden attributes
        self._quatfile = 'gist-quats.dat'
        self._rdf_names = {} # Will be updated after initializing required_keys and allowed_keys.
        self._nvoxels = 0 # Will be calculated from grid lengths after a successful sanity check.
        self._all_keys = OrderedDict({}) # Will be updated after initializing required_keys and allowed_keys.
        self._rdf_cpptraj_command_file = 'rdf_{0}.in'


    def _sanity_check(self):
        for key, val in self._all_keys.items():
            if val[0] == 'path':
                if self.check_path(str(self.__dict__[key])):
                    continue

                else:
                    self.err_string += '\n\t-{0}. Path not found.'.format(key)
                    continue

            elif val[0] == 'type':
                if self.check_type(str(self.__dict__[key]), eval(val[1])):
                    continue

                else:
                    self.err_string += '\n\t-{0}. Required type: {1}'.format(key, val[1])
                    continue

            elif val[0] == 'format':
                if self.check_format(str(self.__dict__[key]), val[1]): # Conversion of __dict__[val] to string as required by check_format
                    continue

                else:
                    self.err_string += '\n\t-{0}. Required format: {1}'.format(key, val[1])
                    continue

            else:
                print("You've just discovered a bug (required formats). Please reach out to us!")
                sys.exit()

        if len(set([self.rigid_atom_0, self.rigid_atom_1, self.rigid_atom_2])) < 3:
            self.err_string += '\n\t-The rigid atoms must be unique!'

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

    def _execute_cpptraj(self, file):
        if 'CPPTRAJ_BIN' not in SETTINGS:
            raise MissingCpptrajException(
                "The path to the CPPTRAJ binary is not given in the global settings. Execute 'febiss_setup' first")

        elif not os.path.exists(SETTINGS['CPPTRAJ_BIN']):
            sys.exit(print('The CPPTRAJ_BIN ' + SETTINGS['CPPTRAJ_BIN'] + ' as declared in ' + SETTINGS_FILE + ' could '
                                                                                                       'not be found!'))
        subprocess.call([SETTINGS['CPPTRAJ_BIN'], '-i', file])

