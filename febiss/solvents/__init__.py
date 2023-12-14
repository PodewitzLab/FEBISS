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
FILE_DICT["TP3"] = "TP3.xyz"

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
    "TOL": "Toluene",
    "TP3" : "TIP3P Water"
}

RIGID_ATOMS_DICT = { #LM20231207: added -1 on every index since the first atom has an index of 0
    "ACN": (3, 1, 0),  # H1-C-C
    "ACT": (9, 1, 0),  # O-C-C
    "BNZ": (0, 1, 2),  # C-C-C
    # "CHX":(?),
    "CL3": (1, 0, 2),  # Cl-C-H
    "CL4": (1, 0, 2),  # Cl-C-Cl
    "DCM": (1, 0, 2),  # Cl-C-H
    # "DMF":(?),
    "DMS": (9, 0, 1),  # O-S-C
    "ETL": (0, 1, 7),  # C-C-O
    # "HEX":(),
    "MTL": (1, 4, 5),  # C-O-H
    "NH3": (1, 0, 2),  # H-N-H
    # "OCT":(?),
    "PYR": (1, 0, 8),  # C-N-C
    "THF": (3, 8, 5),  # C-O-C
    "TOL": (5, 6, 8),  # C-C-C (all ring)
    "TP3": (1, 0, 2)   # H-O-H
    }

REF_DENS_DICT = {
    #densities (at temperatures 20/25 °C) obtained from CRC Handbook of Chemistry and Physics, 97th ed.;
    # Haynes, W. M., Lide, D. R., Bruno, T. J., Eds.; CRC Press, Taylor & Francis Group: Boca Raton, FL, 2017.
    # https://doi.org/10.1201/9781315380476. p. 15-13 ff.,
    # NH3: NIST database at 1.013 bar and 25 C
    # TP3: CPPTRAJ manual (December 16, 2022), p. 123
    "ACN" : 0.0115,
    "ACT" : 0.0082,
    "BNZ" : 0.0068,
    "CL3" : 0.0075,
    "CL4" : 0.0062,
    "DCM" : 0.0094,
    "DMS" : 0.0085,
    "ETL" : 0.0103,
    "MTL" : 0.0149,
    "NH3" : 0.000025,
    "PYR" : 0.0075,
    "THF" : 0.0074,
    "TOL" : 0.0057,
    "TP3" : 0.0334,
}

