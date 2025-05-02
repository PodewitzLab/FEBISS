#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright Technische Universität Wien, Institute of Materials Chemistry, Podewitz Group
See LICENSE for details
"""

import os
import glob
import re
import sys


class Checker:
    def __init__(self):
        self.err_string = ""

    @staticmethod
    def check_type(val, type0: type, *args):
        """
        Checks if the type of the given val can be converted into type0
        and returns True if this can be done and False otherwise.

        :param val: Value to be checked.
        :param type0: Desired type of val.
        :param args: Optional argument for type conversion. For example to specify the base in int(...)
        :return: True or False
        """
        try:
            type0(val, *args) # check if val can be converted into type0
            return True
        except ValueError:
            return False

    @staticmethod
    def check_format(val: str, pattern: str):
        if type(re.fullmatch(pattern, val)) == re.Match:
            return True
        else:
            return False

    @staticmethod
    def check_path(path: str):
        if "*" in path: # path contains "*" as wildcard. checks if any matching file exists
            for p in glob.glob(path):
                if os.path.exists(p):
                    return True
                else:
                    return False

        else: # checks if given path exists
            if os.path.exists(path):
                return True
            else:
                return False

    @staticmethod
    def check_uniques(vals: list):
        if len(vals) == len(set(vals)):
            return True
        else:
            return False
