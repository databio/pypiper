""" Helpers for tests """

from functools import partial
import pytest

__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"



def named_param(argnames, argvalues):
    """
    Improve pytest's native labeling of test case parameterization.

    This function thinly wraps the 'parametrize' mark from pytest, adding
    clearer labeling of each individual parameterized test case, overriding
    the index-based labeling that pytest uses by default.

    :param str argnames: Single parameter name, named in the plural only for
        concordance with the native pytest name.
    :param Iterable argvalues: Arguments for the parameter, what define the
        distinct test cases.
    :return functools.partial: Parameterize version of parametrize, with
        values and ids fixed.
    """
    return partial(pytest.mark.parametrize(
                   argnames=argnames, argvalues=argvalues,
                   ids=lambda val: "{}={}".format(argnames, val)))
