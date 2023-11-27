#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright University Innsbruck, Institute for General, Inorganic, and Theoretical Chemistry, Podewitz Group
See LICENSE for details
"""


CASE_DICT = {1: "TIP3P",
             2: "PATH_TO_WATER_FILE",
             3: "PATH_TO_PYCONSOLV_FILE",
             4: "PATH_TO_SOLVENT_FILE"}
"""
This case dict is used in febiss_settings to create the proper yaml-file.
"""

SOLVENT_LIST = ["Acetonitrile (ACN)",
                "Acetone (ACT)",
                "Benzene (BNZ)",
                # "Cyclohexane (CHX)",
                "Chloroform (CL3)",
                "Carbon tetrachloride (CL4)",
                "Dichloromethane (DCM)",
                # "Dimethylformamide (DMF)",
                "Dimethyl sulfoxide (DMS)",
                "Ethanol (ETL)",
                # "n-Hexane (HEX)",
                "Methanol (MTL)",
                "Ammonia (NH3)",
                # "n-Octanol (OCT)",
                "Pyridine (PYR)",
                "Tetrahydrofurane (THF)",
                "Toluene (TOL)"]
"""
These solvents are taken from PyConSolv (https://github.com/PodewitzLab/PyConSolv/tree/e893eb0da780cdd7205f5f94433a599119f1b5b5/src/PyConSolv/solvents) 
"""

FILE_DICT = {}
for solvent in SOLVENT_LIST:
    FILE_DICT[solvent[-4:-1]] = "{0}.mol2".format(solvent[-4:-1])

NAME_DICT = {
    "ACN": "Acetonitrile",
    "ACT": "Acetone",
    "BNZ": "Benzene",
    # "CHX": "Cyclohexane",
    "CL3": "Chloroform",
    "CL4": "Carbon tetrachloride",
    "DCM": "Dichloromethane",
    # "DMF": "Dimethyl formamide",
    "DMS": "Dimethyl sulfoxide",
    "ETL": "Ethanol",
    # "HEX":(),
    "MTL": "Methanol",
    "NH3": "Ammonia",
    # "OCT":"n-Octane",
    "PYR": "Pyridine",
    "THF": "Tetrahydrofurane",
    "TOL": "Toluene"
}

RIGID_ATOMS_DICT = {
    "ACN": (4, 2, 1),  # H1-C-C
    "ACT": (10, 2, 1),  # O-C-C
    "BNZ": (1, 2, 3),  # C-C-C
    # "CHX":(?),
    "CL3": (2, 1, 3),  # Cl-C-H
    "CL4": (2, 1, 3),  # Cl-C-Cl
    "DCM": (2, 1, 3),  # Cl-C-H
    # "DMF":(?),
    "DMS": (10, 1, 2),  # O-S-C
    "ETL": (1, 2, 8),  # C-C-O
    # "HEX":(),
    "MTL": (2, 5, 6),  # C-O-H
    "NH3": (2, 1, 3),  # H-N-H
    # "OCT":(?),
    "PYR": (2, 1, 9),  # C-N-C
    "THF": (4, 9, 6),  # C-O-C
    "TOL": (6, 7, 9)  # C-C-C (all ring)
    }
