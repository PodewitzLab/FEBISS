#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright Technische Universität Wien, Institute of Materials Chemistry, Podewitz Group
See LICENSE for details
"""

from collections import OrderedDict

NUMBER_DICT = OrderedDict({
    '1':'Acetonitrile',
    '2':'Acetone',
    '3':'Benzene',
    '4':'Chloroform',
    '5':'Tetrachloro methane',
    '6':'Dichloro methane',
    '7':'Dimethyl sulfoxide',
    '8':'Ethanol',
    '9':'Methanol',
    '10':'Ammonia',
    '11':'Pyridine',
    '12':'Tetrahydrofuran',
    '13':'Toluene',
    #'14':'Water',
    #'15':'Custom'
})

SOLVENT_DICT = {
    #Format: Solvent name: (Abbreviation, Solvent file, Rigid Atoms, Reference Density, Reference Evv)
    #These solvents are taken from PyConSolv (https://github.com/PodewitzLab/PyConSolv/tree/e893eb0da780cdd7205f5f94433a599119f1b5b5/src/PyConSolv/solvents)

    'Acetonitrile': ('ACN', 'ACN.mol2', (3, 1, 0), 0.0121, -11.51),
    'Acetone': ('ACT', 'ACT.mol2', (9, 1, 0), 0.0085, -9.33),
    'Benzene': ('BNZ', 'BNZ.mol2', (0, 1, 2), 0.0066, -6.95),
    # "Cyclohexane (CHX)"
    'Chloroform': ('CL3', 'CL3.mol2', (1, 0, 2), 0.0073, -6.49),
    'Tetrachloro methane': ('CL4', 'CL4.mol2', (1, 0, 2), 0.0061, -7.25),
    'Dichloro methane': ('DCM', 'DCM.mol2', (1, 0, 2), 0.0094, -6.02),
    # "Dimethylformamide (DMF)",
    'Dimethyl sulfoxide': ('DMS', 'DMS.mol2', (9, 0, 1), 0.0089, -14.77),
    'Ethanol': ('ETL', 'ETL.mol2', (0, 1, 7), 0.0106, -12.8),
    # "n-Hexane (HEX)",
    'Methanol': ('MTL', 'MTL.mol2', (1, 4, 5), 0.0154, -10.17),
    'Ammonia': ('NH3', 'NH3.mol2', (1, 0, 2), 2.5e-05, -8.073),
    # "n-Octanol (OCT)",
    'Pyridine': ('PYR', 'PYR.mol2', (1, 0, 8), 0.0075, 0.0),
    'Tetrahydrofuran': ('THF', 'THF.mol2', (3, 8, 5), 0.0074, 0.0),
    'Toluene': ('TOL', 'TOL.mol2', (5, 6, 8), 0.0055, 0.0),
    'Water': ('WAT', 'TP3.xyz', (1, 0, 2), 0.0329, -9.544),
    'Custom':('UNL', '', (0, 1, 2), 0.0000, 0.0)
}