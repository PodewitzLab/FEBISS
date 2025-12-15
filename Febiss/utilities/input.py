#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright Technische Universität Wien, Institute of Materials Chemistry, Podewitz Group
See LICENSE for details
"""
class Input:
    """
    Input class that is used for user input requests in the CLI. Also holds additional input checks.
    """
    def __init__(self, prompt, **kwargs):
        self.prompt = str(prompt)
        self.input = input(self.prompt +": ")
        if len(kwargs) > 0:
            self.form = kwargs["type"]
            self.__assertion()

    def yn(self) -> bool:
        while self.input not in ['Y', 'N', 'y', 'n']:
            self.prompt = "Please respond with 'y' or 'n' only: "
            self.input = input(self.prompt)
        if self.input.lower() == "y":
            return True
        else:
            return False

    def reassure(self):
        if not Input("You provided the following input: " + str(self.input) + "\n Is that correct? ").yn():
            self.input = input(self.prompt)
            self.__assertion()

    def __assertion(self):
        try:
            self.input = self.form(self.input)
            self.reassure()
        except ValueError:
            print("The input does not fit the requirements, must be " + str(self.form.__name__))
            self.input = input(self.prompt + ": ")
            self.__assertion()
        except SyntaxError:
            print("The input does not fit the requirements, must be " + str(self.form.__name__))
            self.input = input(self.prompt + ": ")
            self.__assertion()
