"""Tests for checking a collection of commands for callability"""

import mock
import pytest
from ubiquerg import powerset

from pypiper import utils as piper_utils
from tests.helpers import ExpectContext

__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


def pytest_generate_tests(metafunc):
    """Dynamic test case generation and parameterization for this module"""
    if "str_list_monad" in metafunc.fixturenames:
        metafunc.parametrize("str_list_monad", [lambda s: s, lambda s: [s]])


@pytest.mark.parametrize(
    ["uncall_result", "expectation"], [([], True), ([("noncmd", "noncmd")], TypeError)]
)
@pytest.mark.parametrize("handler", [lambda: True, "not-a-function"])
def test_check_all_bad_handler_is_type_error_iff_uncallability_exists(
    uncall_result, str_list_monad, handler, expectation
):
    """Invalid handler evaluation is conditional having >= 1 uncallable command."""
    cmd = "noncmd"
    with (
        mock.patch.object(piper_utils, "determine_uncallable", return_value=uncall_result),
        ExpectContext(expectation, piper_utils.check_all_commands) as check,
    ):
        check(cmds=str_list_monad(cmd), handle=handler)


@pytest.mark.parametrize(
    ["create_result", "expected"],
    [
        (
            lambda bads: Exception("{} bad commands: {}".format(len(bads), bads)),
            Exception,
        ),
        (lambda bads: "{} bad commands: {}".format(len(bads), bads), False),
    ],
)
def test_check_all_result_is_conjunctive(create_result, expected, str_list_monad):
    """Even one uncallable means result is False or an Exception occurs."""
    cmd = "noncmd"
    with (
        mock.patch.object(piper_utils, "determine_uncallable", return_value=[(cmd, cmd)]),
        ExpectContext(expected, piper_utils.check_all_commands) as check,
    ):
        check(cmds=str_list_monad(cmd), get_bad_result=create_result)


@pytest.mark.parametrize("commands", ["man", "ls", ["man", "ls"]])
@pytest.mark.parametrize(
    ["transforms", "expectation"],
    [(arg, lambda res: isinstance(res, list)) for arg in [None, []]]
    + [(arg, TypeError) for arg in [1, "a"]],
)
def test_check_all_requires_iterable_transformations_argument(commands, transforms, expectation):
    """If transformations arg is non-null, it must be iterable."""

    def call():
        return piper_utils.determine_uncallable(commands, transformations=transforms)

    if isinstance(expectation, type) and issubclass(expectation, Exception):
        with pytest.raises(expectation):
            call()
    else:
        assert expectation(call())


@pytest.mark.parametrize("commands", powerset(["ls", "picard.jar", "$ENVVAR"], nonempty=True))
def test_transformation_accumulation(commands):
    """Accumulation of transformations works as expected"""
    mapjar = lambda c: "java -jar {}".format(c)
    envjar = "env.jar"
    transforms = [
        (lambda c: c == "$ENVVAR", lambda _: envjar),
        (lambda c: c.endswith(".jar"), mapjar),
    ]
    exps = {"ls": "ls", "picard.jar": mapjar("picard.jar"), "$ENVVAR": mapjar(envjar)}
    with mock.patch.object(piper_utils, "is_command_callable", return_value=False):
        res = piper_utils.determine_uncallable(
            commands, transformations=transforms, accumulate=True
        )
    expectation = [(c, exps[c]) for c in commands]
    print("EXPECTED: {}".format(expectation))
    print("OBSERVED: {}".format(res))
    assert expectation == res


@pytest.mark.parametrize(
    "transforms",
    [
        {(lambda _: True, lambda c: c), (lambda _: False, lambda c: c)},
        {
            "id": (lambda _: True, lambda c: c),
            "java": (lambda c: c.endswith(".jar"), lambda c: "java -jar {}".format(c)),
        },
    ],
)
def test_non_accumulative_but_unordered_transformation_is_exceptional(transforms):
    with pytest.raises(Exception) as err_ctx:
        piper_utils.determine_uncallable("ls", transformations=transforms)
    exp_msg = (
        "If transformations are unordered, non-accumulation of "
        "effects may lead to nondeterministic behavior."
    )
    assert str(err_ctx.value) == exp_msg
