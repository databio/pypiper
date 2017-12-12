
import os

class AttributeDict(object):
    """
    A class to convert a nested Dictionary into an object with key-values
    accessibly using attribute notation (AttributeDict.attribute) instead of
    key notation (Dict["key"]). This class recursively sets Dicts to objects,
    allowing you to recurse down nested dicts (like: AttributeDict.attr.attr)
    """
    def __init__(self, entries, default=False):
        """
        :param entries: A dictionary (key-value pairs) to add as attributes.
        :param default: Should this AttributeDict object return Default
        values for attributes that it does not have? If True, then
        AttributeDict.java would return "java" instead of raising an error,
        if no .java attribute were found.
        """
        self.add_entries(entries, default)
        self.return_defaults = default

    def add_entries(self, entries, default=False):
        for key, value in entries.items():
            if type(value) is dict:
                self.__dict__[key] = AttributeDict(value, default)
            else:
                try:
                    # try expandvars() to allow the yaml to use
                    # shell variables.
                    self.__dict__[key] = os.path.expandvars(value)  # value
                except TypeError:
                    # if value is an int, expandvars() will fail; if 
                    # expandvars() fails, just use the raw value
                    self.__dict__[key] = value

    def __getitem__(self, key):
        """
        Provides dict-style access to attributes
        """
        return getattr(self, key)

    def __repr__(self):
        return str(self.__dict__)

    def __getattr__(self, name):
        if name in self.__dict__.keys():
            return self.name
        else:
            if self.return_defaults:
                # If this object has default mode on, then we should
                # simply return the name of the requested attribute as
                # a default, if no attribute with that name exists.
                return name
            else:
                raise AttributeError("No attribute " + name)
