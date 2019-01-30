""" Tests for pypiper implementation of AttributeDict """

import pytest
from pypiper import AttributeDict as AD


@pytest.mark.parametrize("base", ["random", "irrelevant", "arbitrary"])
@pytest.mark.parametrize("protect", [False, True])
def test_echo_respects_protected(base, protect):
    """ Regardless of AttributeDict behavior, protected member isn't echoed. """
    ad = AD({}, default=True)
    if protect:
        with pytest.raises(AttributeError):
            ad.__getattr__("__{}__".format(base))
    else:
        assert base == ad.__getattr__(base)
