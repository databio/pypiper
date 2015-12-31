#! /usr/bin/env python

import sys
import pypiper

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
    version=pypiper.__version__,
    description='A lightweight python toolkit for gluing together restartable, robust command line pipelines',
    author='Nathan Sheffield, Johanna Klughammer, Andre Rendeiro',
    author_email='nathan@code.databio.org, jklughammer@cemm.oeaw.ac.at',
    url='https://github.com/epigen/pypiper/'
)
