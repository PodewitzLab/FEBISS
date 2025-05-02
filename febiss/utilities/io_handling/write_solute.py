#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright Technische Universität Wien, Institute of Materials Chemistry, Podewitz Group
See LICENSE for details
"""

import os
from .input import Input

def _write_solute_pdb(analyzer) -> str:  # trajin is passed as trajin+format
    file_in = "solute.in"  # a pdb file is not created here but only a .in file for cpptraj
    file_out = "solute.pdb"
    if os.path.isfile(file_out):
        if not Input("\n\nDo you want to overwrite solute.pdb in your directory?").yn():
            file_out = input("\n\nPlease provide a name for the solute pdb file "
                             "('.pdb' is automatically appended to your input):\n")

    with open(file_in, 'w') as f:
        f.write('parm ' + analyzer.top + '\n')
        f.write('trajin ' + analyzer.trajectory_file + '\n')
        f.write('strip :' + analyzer.solv_abb + '\n')
        f.write('trajout ' + file_out + ' onlyframes 1\n')
        f.write('run\n')
        f.write('quit\n')

    return file_in