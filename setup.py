#!/usr/bin/env python3
# -*- coding: utf-8 -*-
__copyright__ = """
This code is licensed under the MIT license.
Copyright Technische Universität Wien, Institute of Materials Chemistry, Podewitz Group
See LICENSE for details
"""

from os import path
from setuptools import setup, find_packages
import sys

min_version = (3, 6)
if sys.version_info < min_version:
    error = """
Puffin does not support Python {0}.{1}.
Python {2}.{3} and above is required. Check your Python version like so:

python3 --version

This may be due to an out-of-date pip. Make sure you have pip >= 9.0.1.
Upgrade pip like so:

pip install --upgrade pip
""".format(*(sys.version_info[:2] + min_version))
    sys.exit(error)

here = path.dirname(path.realpath(__file__))

with open(path.join(here, 'README.rst'), encoding='utf-8') as readme_file:
    readme = readme_file.read()

with open(path.join(here, 'requirements.txt')) as requirements_file:
    requirements = [line for line in requirements_file.read().splitlines()
                    if not line.startswith('#')]
setup(
    name="Febiss",
    version="1.9.0",
    description="Tool to ease GIST analysis and display and select FEBISS waters",
    long_description=readme,
    author="Miguel Steiner, Lukas Magenheim",
    author_email="steiner.mig@gmail.com, l.magenheim@protonmail.com",
    python_requires='>={}'.format('.'.join(str(n) for n in min_version)),
    install_requires=requirements,
    zip_safe=False,
    license="MIT",
    classifiers=[
        "Programming Language :: Python",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry"
    ],
    packages=['febiss',
              'febiss.cli',
              'febiss.plotting',
              'febiss.solvents',
              'febiss.utilities'],
    package_data={
        'febiss': [
            'manual/*',
            'solvents/*.mol2',
            'solvents/TP3.xyz'

        ]
    },
    py_modules=[
        "febiss.plotting",
        "febiss.utilities",
        "febiss.cli",
        "febiss.gist",
    ],
    entry_points={
        'console_scripts': [
            'febiss = febiss.cli.febiss:main',
            'febiss_setup = febiss.cli.febiss_setup:main',
            'febiss_settings = febiss.cli.febiss_settings:main',
        ],
    },
)