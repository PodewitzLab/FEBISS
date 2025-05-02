#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright Technische Universität Wien, Institute of Materials Chemistry, Podewitz Group
See LICENSE for details
"""

import os
from ..gist import GistAnalyser
from ...plotting.display import Plot
def write_settings_file(analyser:GistAnalyser, display:Plot):
    # write all analysis options
    settings_file = 'all-settings.yaml'
    with open(settings_file, 'w') as f:
        f.write("# Two general blocks 'gist' and 'plotting' are given.\n")
        f.write('# For each block all settings and the default values are given and\n')
        f.write('# this file can be used directly after filling out the two required settings.\n')
        f.write('# Where it may be useful, a possible input is given behind the default value as a comment.\n')
        f.write('gist:\n')

        for key in analyser._all_keys.keys():
                f.write('  ' + str(key) + ": " + str(analyser.__dict__[key])+'\n')

    # write all plotting options
        f.write('plotting:\n')
        for key in display.allowed_keys.keys():
            f.write('  ' + str(key) + ": " + str(display.__dict__[key]) + '\n')

    print("Wrote all possible settings into 'all-settings.yaml' in the current directory.")
    return os.path.abspath(settings_file)