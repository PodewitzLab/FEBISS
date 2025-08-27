#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright Technische Universität Wien, Institute of Materials Chemistry, Podewitz Group
See LICENSE for details
"""

import yaml
import sys
from Febiss.cpptraj_interface.gist import GistAnalyser
from Febiss.plotting.display import Plot

def read_settings(febiss_file):
    with open(febiss_file, 'r') as yamlfile:
        param = yaml.safe_load(yamlfile)

        #debug
        #print(param)
        #for key in param.keys():
        #    for key2 in param[key]:
        #        print("The key is {0}, the value is {1} and the type is {2}".format(key2,param[key][key2],type(param[key][key2])))

    # check for GIST analysis
    if "gist" and "plotting" in param.keys():
        analyser = GistAnalyser(**param["gist"])
        display = Plot(**param["plotting"])

        # Sanity check
        analyser._sanity_check()
        display.sanity_checks()

        # debug
        # for key in param["gist"].keys():
        #        print("The key is {0}, the value is {1} and the type is {2}".format(key,analyser.__dict__[key],type(analyser.__dict__[key])))
        # for key in param["plotting"].keys():
        #        print("The key is {0}, the value is {1} and the type is {2}".format(key,display.__dict__[key],type(display.__dict__[key])))

        if (analyser.err_string or display.err_string) != '':
            print('The following settings did not pass the checks:' + analyser.err_string + display.err_string)
            sys.exit()

    else:
        analyser, display = None, None
        quit("GIST parameters are missing in the all-settings.yaml file. Before running febiss, run febiss_settings"
             "to create the all-settings.yaml file which contain crucial informations to perform a GIST analysis "
             "or run febiss without arguments which opens a GUI window.")

    #debug
    #for key in param["gist"].keys():
    #        print("The key is {0}, the value is {1} and the type is {2}".format(key,analyser.__dict__[key],type(analyser.__dict__[key])))
    #for key in param["plotting"].keys():
    #        print("The key is {0}, the value is {1} and the type is {2}".format(key,display.__dict__[key],type(display.__dict__[key])))

    return analyser, display