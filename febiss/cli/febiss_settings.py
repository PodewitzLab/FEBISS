#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright University Innsbruck, Institute for General, Inorganic, and Theoretical Chemistry, Podewitz Group
See LICENSE for details
"""

import sys


def help_message():
    print("\nThis program writes all possible settings of the GIST analysis "
          "and the plotting into 'all-settings.yaml'.")
    print('It requires no arguments.')
    sys.exit()


def main():
    if len(sys.argv) > 1:
        help_message()
    from ..utilities.gist import GistAnalyser
    # the arguments are necessary in the init, otherwise exception
    analyser = GistAnalyser(**{'top': 'molecule.parm7', 'trajectory_name': 'molecule'})
    # write all analysis options
    with open('all-settings.yaml', 'w') as f:
        f.write("# Two general blocks 'gist' and 'plotting' are given.\n")
        f.write('# For each block all settings and the default values are given and\n')
        f.write('# this file can be used directly after filling out the two required settings.\n')
        f.write('# Where it may be useful, a possible input is given behind the default value as a comment.\n')
        f.write('gist:\n')
        for key in analyser.required_keys:
            f.write('  ' + str(key) + ": " + str(analyser.__dict__[key]) + ' # this has to be filled out\n')
        for key in sorted(analyser.allowed_keys):
            if key == 'frame_selection':
                f.write('  ' + str(key) + ": " + str(analyser.__dict__[key]) + " # example: '1 1000 5'\n")
            elif key == 'grid_center':
                f.write('  ' + str(key) + ": " + str(analyser.__dict__[key]) + " # example: (10.0, 5.0, 0.0)\n")
            elif key == 'rdf_names':
                f.write('  ' + str(key) + ": " + str(analyser.__dict__[key]).replace('center of solute',
                                                                                     'center\ of\ solute') + '\n')
            else:
                f.write('  ' + str(key) + ": " + str(analyser.__dict__[key]) + '\n')
    from ..plotting.display import Plot
    display = Plot()
    # write all plotting options
    with open('all-settings.yaml', 'a') as f:
        f.write('plotting:\n')
        for key in sorted(display.allowed_keys):
            f.write('  ' + str(key) + ": " + str(display.__dict__[key]) + '\n')
    print("Wrote all possible settings into 'all-settings.yaml' in the current directory.")
    sys.exit()


if __name__ == '__main__':
    main()
