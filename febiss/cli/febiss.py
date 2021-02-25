#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright University Innsbruck, Institute for General, Inorganic, and Theoretical Chemistry, Podewitz Group
See LICENSE for details
"""

from shutil import which
import subprocess
import sys
import yaml

from ..utilities.structures import Solute, Water
from ..utilities.io_handling import read_pdb, write_style_file
from ..plotting.display import Plot


def help_message():
    print('\nThis program takes settings from a yaml file given as command line argument.')
    print("To obtain a suitable file execute 'get_febiss_settings'.")
    print("This program can both facilitate the execution of the GIST analysis via CPPTRAJ "
          "and display and select results from the FEBISS analysis within the GIST analysis in CPPTRAJ.")
    print("To execute CPPTRAJ, it has to be first setup via 'setup_febiss'.\n")
    sys.exit()


def main():
    if len(sys.argv) != 2:
        help_message()
    arg = sys.argv[1]
    if arg.lower() in ['-h', '--help']:
        help_message()
    with open(arg, 'r') as yamlfile:
        param = yaml.safe_load(yamlfile.read().replace("\t", "  ").replace("    ", "  "))
    # first check for gist analysis
    if "gist" in param.keys():
        from ..utilities.gist import GistAnalyser
        analyser = GistAnalyser(**param["gist"])
        analyser.perform_gist_analysis()
    elif "GIST" in param.keys():
        from ..utilities.gist import GistAnalyser
        analyser = GistAnalyser(**param["GIST"])
        analyser.perform_gist_analysis()
    # now check for plotting, else assume that only plotting options are given directly
    for key in param.keys():
        if key.lower() == 'plotting':
            param = param['plotting']
            break
    solute = Solute()
    water = Water()
    febiss_file = param.get('febiss_file', 'febiss-waters.pdb')
    read_pdb(febiss_file, solute, water)
    display = Plot(**param)
    filename = display.gui(solute, water)
    if which('pymol') is not None:
        style_file = param.get('style_file', write_style_file())
        subprocess.call(['pymol', filename, style_file])
    print('FEBISS ended successfully')
    sys.exit()


if __name__ == '__main__':
    main()
