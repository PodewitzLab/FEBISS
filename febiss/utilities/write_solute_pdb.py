#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright University Innsbruck, Institute for General, Inorganic, and Theoretical Chemistry, Podewitz Group
See LICENSE for details
"""

import os
from io_handling import Input

"""adopted from ..utilities/gist/_write_cpptraj_file"""
def write_solute_pdb(top : str, trajin : str, abb : str) -> str: #trajin is passed as trajin+format
    file = "solute.pdb"
    if os.path.isfile(file):
        if not Input("\n\nDo you want to overwrite solute.pdb in your directory?\n").yn():
            file = input("\n\nPlease provide a name for the solute pdb file ('.pdb' is automatically appended to your input):\n")
    with open(file, 'w') as f:
        f.write('parm ' + top + '\n')
        f.write('trajin ' + trajin + '\n')
        f.write('strip :' + abb + '\n')
        f.write('trajout ' + file + 'onlyframes 1\n')
        f.write('run\n')
        f.write('quit\n')
    return file