"""Tests for EchoDict class"""

import pytest

from pypiper.echo_dict import EchoDict


class TestEchoDict:
    """Test the EchoDict class functionality"""

    def test_attribute_access_for_existing_keys(self):
        """Test that we can access existing keys as attributes"""
        d = EchoDict({"foo": "bar", "num": 42})
        assert d.foo == "bar"
        assert d.num == 42

    def test_echo_behavior_for_missing_keys(self):
        """Test that missing keys return the key name (echo behavior)"""
        d = EchoDict({"foo": "bar"})
        assert d.missing == "missing"
        assert d.another_missing == "another_missing"

    def test_nested_dict_wrapping(self):
        """Test that nested dicts are automatically wrapped as EchoDict"""
        d = EchoDict({"outer": {"inner": "value"}})
        assert isinstance(d.outer, EchoDict)
        assert d.outer.inner == "value"
        # Echo behavior should work on nested dicts too
        assert d.outer.missing == "missing"

    def test_setting_attributes(self):
        """Test that we can set attributes"""
        d = EchoDict()
        d.foo = "bar"
        assert d.foo == "bar"
        assert d["foo"] == "bar"

    def test_dict_methods_work(self):
        """Test that standard dict methods still work"""
        d = EchoDict({"a": 1, "b": 2, "c": 3})
        assert set(d.keys()) == {"a", "b", "c"}
        assert set(d.values()) == {1, 2, 3}
        assert ("a", 1) in d.items()

    def test_construction_from_yaml_dict(self):
        """Test that EchoDict works with dict from YAML loading"""
        # Simulate what load_yaml returns
        yaml_data = {
            "tools": {"samtools": "/usr/bin/samtools", "bwa": "/usr/bin/bwa"},
            "parameters": {"threads": 4},
        }
        d = EchoDict(yaml_data)
        assert d.tools.samtools == "/usr/bin/samtools"
        assert d.tools.bwa == "/usr/bin/bwa"
        assert d.parameters.threads == 4

    def test_echo_behavior_in_nested_dict(self):
        """Test that echo behavior works in nested structures"""
        d = EchoDict({"tools": {}})
        # Even though "samtools" is not defined, it should echo
        assert d.tools.samtools == "samtools"

    def test_private_attributes_raise_error(self):
        """Test that private attributes (starting with _) raise AttributeError"""
        d = EchoDict({"foo": "bar"})
        with pytest.raises(AttributeError):
            _ = d._private

    def test_delattr(self):
        """Test that deleting attributes works"""
        d = EchoDict({"foo": "bar"})
        assert d.foo == "bar"
        del d.foo
        # After deletion, should echo
        assert d.foo == "foo"

    def test_delattr_missing_key(self):
        """Test that deleting non-existent attribute raises AttributeError"""
        d = EchoDict()
        with pytest.raises(AttributeError):
            del d.missing

    def test_mixed_access_patterns(self):
        """Test that dict-style and attribute-style access can be mixed"""
        d = EchoDict()
        d["key1"] = "value1"
        d.key2 = "value2"
        assert d.key1 == "value1"
        assert d["key2"] == "value2"

    def test_nested_echo_with_multiple_levels(self):
        """Test echo behavior with multiple nesting levels"""
        d = EchoDict({"level1": {"level2": {}}})
        assert d.level1.level2.level3 == "level3"
