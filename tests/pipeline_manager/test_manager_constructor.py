""" Test effects of construction of a pipeline manager. """

import pytest


__author__ = "Vince Reuter"
__email__ = "vreuter@virginia.edu"



@pytest.mark.skip("Not implemented")
class ManagerConstructorCheckpointSpecificationTests:
    """ Tests for manager's constructor's ability to parse and set
    checkpoint specifications, which can determine aspects of control flow. """


    def test_no_checkpoint_specifications(self):
        """ A manager may be constructed without any checkpoint provision. """
        pass


    def test_just_start(self):
        """ Starting point may be set from command-line or ctor keyword. """
        pass

    
    def test_just_stop(self):
        """ Particular stopping type is set correctly. """
        pass


    def test_start_and_stop(self):
        """ Specifying both start and stop works just fine. """
        pass


    def test_both_stop_modes_is_prohibited(self):
        """ Provision of both prospective and retrospective stop is bad. """
        pass


    def test_complementary_specification_modes(self):
        """ Command-line and keyword specifciations can harmonize. """
        pass


    def test_command_line_beats_constructor_keyword(self):
        """ Command-line specification is favored over constructor keyword. """
        pass
