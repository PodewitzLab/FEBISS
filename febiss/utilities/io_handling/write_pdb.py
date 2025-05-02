#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright Technische Universität Wien, Institute of Materials Chemistry, Podewitz Group
See LICENSE for details
"""

from typing import Union
from ...utilities.structures import Solute, Solvent



def write_pdb(pdb: str, structure: Union[Solute, Solvent], abb, solute: bool = False):
    """
    Function that writes the microsolvated structure to a PDB file.
    :param pdb: Name of the pdb-file to be written into.
    :param structure: Either a Solute or a Solvent object.
    :param abb: Holds the three letter abbreviation of the solvent (e.g. CL3 for chloroform)
    :param solute: Bool variable to distinguish whether the solute or solvent molecules get written.
    :return: No actual return but writes PDB file that holds the microsolvated structure.
    """
    if solute:
        atomcounter = 1
        f = open(pdb, 'w')
    else:
        atomcounter = len(open(pdb, 'r').readlines()) + 1
        f = open(pdb, 'a')

    for count, (ele, atom) in enumerate(zip(structure.elements, structure.coords)):
        j = []
        if solute:
            j.append('ATOM'.ljust(6))  # atom#6s
        else:
            j.append('HETATM'.ljust(6))  # atom#6s
        j.append(str(atomcounter + count).rjust(5))  # aomnum#5d
        j.append(ele.center(4))  # atomname$#4s
        if solute:
            j.append('SOL'.ljust(3))  # resname#1s
        else:
            j.append(abb.ljust(3))  # resname#1s abb instead of "FEB"
        j.append('A'.rjust(1))  # Astring
        if solute:
            j.append('1'.rjust(4))  # resnum
        else:
            j.append('2'.rjust(4))  # resnum
        j.append(str('%8.3f' % (float(atom[0]))).rjust(8))  # x
        j.append(str('%8.3f' % (float(atom[1]))).rjust(8))  # y
        j.append(str('%8.3f' % (float(atom[2]))).rjust(8))  # z
        j.append(str('%6.2f' % 1.0).rjust(6))  # occ
        if solute:
            value = 0.0
        else:
            value = float(structure.values[count])
        if value == 0.0:
            j.append(str('%7.2f' % value).ljust(7))  # delta G
        else:
            j.append(str('%7.2f' % (-1 * value)).ljust(7))  # delta G
        j.append(ele.rjust(12))  # elname
        f.write("%s%s %s %s %s%s    %s%s%s%s%s%s\n" % (
            j[0], j[1], j[2], j[3], j[4], j[5], j[6], j[7], j[8], j[9], j[10], j[11]))
    f.close()