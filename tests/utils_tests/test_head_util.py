""" Tests for the head() utility function """

import pytest
from hypothesis import given, strategies as st
from pypiper.utils import head


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


# Strategy for generating a pretty arbitrary atomic
#ATOMICS = st.deferred(lambda: )


#@st.composite
#@pytest.mark.parametrize("coll", [st.lists, st.tuples])
#def head_coll_pair(draw, coll, elements=ATOMICS):
#    """ Strategy to generate a pair of head element and rest of Iterable. """
#    return draw(coll(elements)), draw(elements)


#@given("obj", st.booleans() | st.characters() | st.integers() |
#                              st.floats() | st.text())
#def test_head_atomic(obj):
#    """ head() of an atomic object is the object itself. """
#    assert obj == head(obj)


def test_head_empty_string():
    """ Empty string is exception to exceptional-ness of empty collection. """
    assert "" == head("")


@pytest.mark.skip("not implemented")
@pytest.mark.parametrize("coll", [dict(), set(), tuple(), list()])
def test_head_empty_collection(coll):
    """ Request for first element from an empty Iterable is exceptional. """
    with pytest.raises(ValueError):
        head(coll)


"""
@pytest.mark.skip("not implemented")
@pytest.mark.parametrize("coll_type", [list, tuple])
@pytest.mark.parametrize("transform", [lambda c: c, lambda c: iter(c)])
@given(st.data())
def test_head_nonempty_sequential_collection(data, coll_type, transform):
    #" Verify accuracy of request for first element from nonempty Iterable. "
    h = data.draw(ATOMICS)
    c = coll_type([h]) + coll_type(data.draw(st.lists(ATOMICS)))
    assert h == head(transform(c))
"""

def test_head_nonempty_set():
    """ Verify that head of nonempty set is non-exceptional. """
    head({-1, 0, 1})


def test_head_nonempty_dict():
    """ Verify that head of nonempty dictionary is non-exceptional. """
    head({"a": 1, "b": 2})


#@pytest.mark.skip("not implemented")
#@given(st.data())
#@pytest.mark.parametrize("coll_type", [list, tuple])
#def test_head_nonempty_sequential_collection(data, coll_type):
#    """ Verify accuracy of request for first element from nonempty Iterable. """
#    h = data.draw(ATOMICS)
#    c = coll_type([h]) + coll_type(data.draw(st.lists(ATOMICS)))
#    assert h == head(iter(c))
