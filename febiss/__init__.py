#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright University Innsbruck, Institute for General, Inorganic, and Theoretical Chemistry, Podewitz Group
See LICENSE for details
"""

from typing import Dict
import os
import yaml

__author__ = "Miguel Steiner"
__email__ = "steiner.mig@gmail.com"
__maintainer__ = "Maren Podewitz"
__maintainer_email__ = "maren.podewitz@uibk.ac.at"

SETTINGS_FILE = os.path.join(os.path.expanduser("~"), ".febissrc.yaml")


def _load_febiss_settings() -> Dict[str, str]:
    try:
        with open(SETTINGS_FILE, "rt") as f:
            d = yaml.safe_load(f)
    except IOError:
        # If there are any errors, default to using environment variables
        # if present.
        d = {}
        for k, v in os.environ.items():
            if k in ['CPPTRAJ_BIN', 'CPPTRAJ_HOME, ''GIGIST_HOME', 'FEBISS_HOME']:
                d[k] = v
    return dict(d)


SETTINGS = _load_febiss_settings()
