#! /usr/bin/env python

import sys
import re
import os

extra = {}

try:
    from setuptools import setup
    if sys.version_info < (2, 7):
        extra['install_requires'] = ['argparse']
    if sys.version_info >= (3,):
        extra['use_2to3'] = True
except ImportError:
    from distutils.core import setup
    if sys.version_info < (2, 7):
        extra['dependencies'] = ['argparse']

version = re.search(
    '^__version__\s*=\s*"(.*)"',
    open('pipelines/pipelines.py').read(),
    re.M
    ).group(1)

setup(
    name='pypiper',
    packages=['pypiper']
    version=version,
    description='A lightweight python toolkit for gluing together restartable, robust command line pipelines',
    author='Nathan Sheffield, Johanna Klughammer, ',
    author_email='nsheffield@cemm.oeaw.ac.at, jklughammer@cemm.oeaw.ac.at',
    url='https://github.com/ComputationalEpigenetics/pypiper/',
)
