#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright Technische Universität Wien, Institute of Materials Chemistry, Podewitz Group
See LICENSE for details
"""

from typing import Dict
import os
import yaml

__author__ = "Miguel Steiner, Lukas Magenheim"
__email__ = "steiner.mig@gmail.com, l.magenheim@protonmail.com"
__maintainer__ = "Maren Podewitz"
__maintainer_email__ = "maren.podewitz@tuwien.ac.at"

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
            if k in ['CPPTRAJ_BIN', 'CPPTRAJ_HOME', 'FEBISS_HOME']: #deleted GIGIST_HOME LM20231214
                d[k] = v
    return dict(d)


SETTINGS = _load_febiss_settings()
