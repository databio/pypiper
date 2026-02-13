"""Minimal AttMapEcho replacement for backwards compatibility."""

__all__ = ["EchoDict"]


class EchoDict(dict):
    """Dict with attribute access and echo-on-missing-key behavior.

    This is a minimal shim to maintain compatibility with pipelines
    that use pm.config.tools.samtools style access.

    For new code, prefer explicit config classes.
    """

    def __getattr__(self, key):
        if key.startswith("_"):
            raise AttributeError(key)
        try:
            val = self[key]
        except KeyError:
            return key  # echo: return the key name
        if isinstance(val, dict) and not isinstance(val, EchoDict):
            val = EchoDict(val)
            self[key] = val
        return val

    def __setattr__(self, key, value):
        if key.startswith("_"):
            super().__setattr__(key, value)
        else:
            self[key] = value

    def __delattr__(self, key):
        try:
            del self[key]
        except KeyError:
            raise AttributeError(key)
