#! /usr/bin/env python

import sys

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


with open("pypiper/_version.py", 'r') as versionfile:
    version = versionfile.readline().split()[-1].strip("\"'\n")

setup(
    name='pypiper',
    packages=['pypiper'],
    version=version,
    description='A lightweight python toolkit for gluing together restartable, robust command line pipelines',
    author='Nathan Sheffield, Johanna Klughammer, Andre Rendeiro',
    author_email='nathan@code.databio.org, jklughammer@cemm.oeaw.ac.at, arendeiro@cemm.oeaw.ac.at',
    url='https://github.com/epigen/pypiper/',
    **extra
)
