""" Validate what's available directly on the top-level import. """

from inspect import isfunction

import pytest

__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


@pytest.mark.parametrize(
    ["obj_name", "typecheck"],
    [
        ("add_logging_options", isfunction),
        ("check_all_commands", isfunction),
        ("determine_uncallable", isfunction),
        ("logger_via_cli", isfunction),
    ],
)
def test_top_level_exports(obj_name, typecheck):
    """At package level, validate object availability and type."""
    import pypiper

    try:
        obj = getattr(pypiper, obj_name)
    except AttributeError:
        pytest.fail("Unavailable on {}: {}".format(pypiper.__name__, obj_name))
    else:
        assert typecheck(obj)
