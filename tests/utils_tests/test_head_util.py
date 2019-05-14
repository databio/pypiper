""" Tests for the head() utility function """

import random
import string
import pytest
from hypothesis import given, strategies as st
from pypiper.utils import head


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


NUMBERS_AND_LETTERS = list(string.ascii_letters) + list(range(-9, 10))

# Strategy for generating a pretty arbitrary atomic
ATOMICS = st.deferred(lambda: st.booleans() | st.characters() | st.integers() |
                              st.floats(allow_nan=False) | st.text())


def pytest_generate_tests(metafunc):
    """ Test case generation/parameterization for this module. """
    if "seqtype" in metafunc.fixturenames:
        metafunc.parametrize("seqtype", [tuple, list])
    if "iter_cast" in metafunc.fixturenames:
        metafunc.parametrize("iter_cast", [lambda c: c, lambda c: iter(c)])
    if "h" in metafunc.fixturenames and "xs" in metafunc.fixturenames:
        metafunc.parametrize(
            ["h", "xs"],
            [(random.choice(NUMBERS_AND_LETTERS),
              [random.choice(NUMBERS_AND_LETTERS)
               for _ in range(random.randint(5, 10))]) for _ in range(10)])


@given(obj=ATOMICS)
def test_head_atomic(obj):
    """ head() of an atomic object is the object itself. """
    assert obj == head(obj)


def test_head_empty_string():
    """ Empty string is exception to exceptional-ness of empty collection. """
    assert "" == head("")


@pytest.mark.parametrize("coll", [dict(), set(), tuple(), list()])
def test_head_empty_collection(coll):
    """ Request for first element from an empty Iterable is exceptional. """
    with pytest.raises(ValueError):
        head(coll)


def test_head_nonempty_sequential_collection(h, xs, seqtype, iter_cast):
    """ Verify accuracy of request for first element from nonempty Iterable. """
    c = seqtype([h]) + seqtype(xs)
    assert h == head(iter_cast(c))


def test_head_nonempty_set():
    """ Verify that head of nonempty set is non-exceptional. """
    head({-1, 0, 1})


def test_head_nonempty_dict():
    """ Verify that head of nonempty dictionary is non-exceptional. """
    head({"a": 1, "b": 2})
