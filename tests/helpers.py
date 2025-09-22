"""Helpers for tests"""

import glob
import os
from functools import partial

import pytest

from pypiper import Pipeline
from pypiper.utils import checkpoint_filepath

__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


def assert_equal_dirpath(p1, p2):
    """
    Assert that a pair of folder paths has two equal members.

    :param str p1: One path to compare.
    :param str p2: Other path to compare.
    """
    assert p1.rstrip(os.sep) == p2.rstrip(os.sep)


def fetch_checkpoint_files(pm):
    """
    Fetch all of a manager's checkpoint file paths.

    :param pyiper.PipelineManager pm: manager for which checkpoint files'
        paths are of interest.
    :return Iterable[str]: collection of all of given manager's checkpoint
        files' paths.
    """
    pattern = checkpoint_filepath("*", pm)
    return glob.glob(pattern)


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
    return partial(
        pytest.mark.parametrize(
            argnames=argnames,
            argvalues=argvalues,
            ids=lambda val: "{}={}".format(argnames, val),
        )
    )


class SafeTestPipeline(Pipeline):
    """Pipeline for tests that protects against bad file descriptor."""

    def __init__(self, *args, **kwargs):
        kwd_args = {"multi": True}  # Like interactive mode.
        kwd_args.update(kwargs)
        super(SafeTestPipeline, self).__init__(*args, **kwd_args)


# originally imported from veracitools v0.1.4 (now deprecated)
class ExpectContext(object):
    """Pytest validation context, a framework for varied kinds of expectations."""

    def __init__(self, expected, test_func):
        """
        Create the test context by specifying expectation and function.

        :param object | type expected: expected result or exception
        :param callable test_func: the callable object to test
        """
        self._f = test_func
        self._exp = expected

    def __enter__(self):
        """Return the instance for use as a callable."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass

    def __call__(self, *args, **kwargs):
        """Execute the instance's function, passing given args/kwargs."""
        if isinstance(self._exp, type) and issubclass(self._exp, Exception):
            with pytest.raises(self._exp):
                self._f(*args, **kwargs)
        else:
            assert self._exp == self._f(*args, **kwargs)
