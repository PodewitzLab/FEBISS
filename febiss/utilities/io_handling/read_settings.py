#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright Technische Universität Wien, Institute of Materials Chemistry, Podewitz Group
See LICENSE for details
"""

import os
import yaml
from ..gist import GistAnalyser
from ...plotting.display import Plot

def read_settings(febiss_file):
    with open(febiss_file, 'r') as yamlfile:
        param = yaml.safe_load(yamlfile.read().replace("\t", "  ").replace("    ", "  "))

    # check for gist analysis
    if "gist" and "plotting" in param.keys():

        analyser = GistAnalyser(**param["gist"])
        display = Plot(**param["plotting"])

    else:
        analyser, display = None, None
        quit("gist parameters are missing in the all-settings.yaml file. Before running febiss, run febiss_settings"
             "to create the all-settings.yaml file which contain crucial informations to perform a gist analysis "
             "or run febiss without arguments which opens a GUI window.")

    return analyser, display