""" Tests for checking a collection of commands for callability """

import mock
import os
import pytest
from pypiper import utils as piper_utils
from veracitools import ExpectContext

__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"


EXTENSIONS = [".py", ".rb", ".sh", ".java", ".jar", ".pl", ".o", ".R", ".r",
              ".cpp", ".c", ".hs", ".scala", ".class"]


def _touch(f):
    """ 'touch' the given file.

    :param str f: filepath to create
    """
    with open(f, 'w'):
        print("touch: {}".format(f))


def _make_exec(f):
    """
    'touch' a file and set exec bit.

    :param str f: path to create
    """
    import subprocess
    _touch(f)
    subprocess.check_call(["chmod", "+x", f])


def pytest_generate_tests(metafunc):
    """ Dynamic test case generation and parameterization for this module """
    if "str_list_monad" in metafunc.fixturenames:
        metafunc.parametrize("str_list_monad", [lambda s: s, lambda s: [s]])


@pytest.mark.parametrize("filename", ["testfile" + x for x in EXTENSIONS])
@pytest.mark.parametrize(["setup", "pretest", "exp_miss"], [
    (lambda _: None,
     lambda f: not os.path.exists(f),
     lambda _: True),
    (_touch,
     lambda f: os.path.isfile(f) and not os.access(f, os.X_OK),
     lambda f: not f.endswith(".jar")),
    (_make_exec,
     lambda f: os.path.isfile(f) and os.access(f, os.X_OK),
     lambda _: False)
])
def test_callability_checker_defaults(tmpdir, filename, setup, pretest, exp_miss):
    cmd = os.path.join(tmpdir.strpath, filename)
    setup(cmd)
    assert pretest(cmd)
    extra_commands = ["this-is-not-a-program", "man", "ls"]
    expected = ["this-is-not-a-program"]
    if exp_miss(cmd):
        expected.append(cmd)
    observed = [c for c, _ in piper_utils.determine_uncallable([cmd] + extra_commands)]
    print("expected: {}".format(expected))
    print("observed: {}".format(observed))
    assert len(expected) == len(observed)
    assert set(expected) == set(observed)


@pytest.mark.parametrize(
    ["uncall_result", "expectation"],
    [([], True), ([("noncmd", "noncmd")], TypeError)])
@pytest.mark.parametrize("handler", [lambda: True, "not-a-function"])
def test_check_all_bad_handler_is_type_error_iff_uncallability_exists(
        uncall_result, str_list_monad, handler, expectation):
    """ Invalid handler evaluation is conditional having >= 1 uncallable command. """
    cmd = "noncmd"
    with mock.patch.object(piper_utils, "determine_uncallable",
                           return_value=uncall_result), \
         ExpectContext(expectation, piper_utils.check_all_commands) as check:
        check(cmds=str_list_monad(cmd), handle=handler)


@pytest.mark.parametrize(["create_result", "expected"], [
    (lambda bads: Exception("{} bad commands: {}".format(len(bads), bads)), Exception),
    (lambda bads: "{} bad commands: {}".format(len(bads), bads), False)
])
def test_check_all_result_is_conjunction(create_result, expected, str_list_monad):
    """ Even one uncallable means result is False or an Exception occurs. """
    cmd = "noncmd"
    with mock.patch.object(piper_utils, "determine_uncallable",
                           return_value=[(cmd, cmd)]), \
        ExpectContext(expected, piper_utils.check_all_commands) as check:
        check(cmds=str_list_monad(cmd), get_bad_result=create_result)
