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


def read_reqs_file(reqs_name):
    """ Read requirements file for given requirements group. """
    path_reqs_file = os.path.join(
            "requirements", "reqs-{}.txt".format(reqs_name))
    with open(path_reqs_file, 'r') as reqs_file:
        return [pkg.rstrip() for pkg in reqs_file.readlines()
                if not pkg.startswith("#")]


with open(os.path.join("pypiper", "_version.py"), 'r') as versionfile:
    version = versionfile.readline().split()[-1].strip("\"'\n")


basic_reqs = read_reqs_file("pypiper")

# Requirements for tests
test_reqs = read_reqs_file("test")

# Allow specification of desired features, which implies dependencies.
addl_reqs = {bundle_name: read_reqs_file(bundle_name)
             for bundle_name in ["ngstk", "plot"]}

# Complete collection of user requirements.
addl_reqs["all"] = list({pkg for bundle in addl_reqs.values() for pkg in bundle})

# Dev installation is full user + test.
addl_reqs["dev"] = list(set(test_reqs + addl_reqs["all"]))

# Handle the pypi README formatting.
try:
    import pypandoc
    long_description = pypandoc.convert_file('README.md', 'rst')
except(IOError, ImportError, OSError):
    long_description = open('README.md').read()

setup(
    name='piper',
    packages=['pypiper'],
    install_requires=basic_reqs,
    version=version,
    description='A lightweight python toolkit for gluing together restartable, robust command line pipelines',
    long_description=long_description,
    long_description_content_type='text/markdown',
    classifiers=[
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 2.7",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],    
    author='Nathan Sheffield, Johanna Klughammer, Andre Rendeiro',
    author_email='nathan@code.databio.org, jklughammer@cemm.oeaw.ac.at, arendeiro@cemm.oeaw.ac.at',
    url='https://github.com/databio/pypiper/',
    license="BSD2",
    test_suite="tests",         # python setup.py test
    tests_require=test_reqs,    # Test-specific package dependencies
    # Extra package if doing `python setup.py test`
    setup_requires=(["pytest-runner"] if {"test", "pytest", "ptr"} & set(sys.argv) else []),
    extras_require=addl_reqs,
    # Version-specific items
    **extra
)
