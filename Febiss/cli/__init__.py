#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright Technische Universität Wien, Institute of Materials Chemistry, Podewitz Group
See LICENSE for details
"""

from ..utilities.colorgen import Color

version = "1.9.9.4"
print(Color.GREEN + r"""
        ______ ______ ____   ____ _____ _____
       / ____// ____// __ ) /  _// ___// ___/
      / /_   / __/  / __  | / /  \__ \ \__ \ 
     / __/  / /___ / /_/ /_/ /  ___/ /___/ / 
    /_/    /_____//_____//___/ /____//____/
                v. {}
    """.format(version) + Color.END)