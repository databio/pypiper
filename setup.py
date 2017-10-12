#! /usr/bin/env python

import os
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


with open(os.path.join("pypiper", "_version.py"), 'r') as versionfile:
    version = versionfile.readline().split()[-1].strip("\"'\n")


with open("requirements-test.txt", 'r') as test_reqs_file:
    test_reqs = [l for l in test_reqs_file.readlines()
                 if l and not l.startswith("#")]


setup(
    name='pypiper',
    packages=['pypiper'],
    version=version,
    description='A lightweight python toolkit for gluing together restartable, robust command line pipelines',
    author='Nathan Sheffield, Johanna Klughammer, Andre Rendeiro',
    author_email='nathan@code.databio.org, jklughammer@cemm.oeaw.ac.at, arendeiro@cemm.oeaw.ac.at',
    url='https://github.com/epigen/pypiper/',
    test_suite="tests",         # python setup.py test
    tests_require=test_reqs,    # Test-specific package dependencies
    # Extra package if doing `python setup.py test`
    setup_requires=(["pytest-runner"] if {"test", "pytest", "ptr"} & set(sys.argv) else []),
    # Version-specific items
    **extra
)
