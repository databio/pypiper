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

setup(
    name='pypiper',
    packages=['pypiper'],
    version='0.1',
    description='A lightweight python toolkit for gluing together restartable, robust command line pipelines',
    author='Nathan Sheffield, Johanna Klughammer, ',
    author_email='nsheffield@cemm.oeaw.ac.at, jklughammer@cemm.oeaw.ac.at',
    url='https://github.com/ComputationalEpigenetics/pypiper/'
)
