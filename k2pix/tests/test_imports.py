"""Simple import tests to guard against syntax errors."""
import pytest

def test_imports():
    """A simple test to make sure k2pix's key functions can be imported."""
    from ..main import k2pix
    from ..tpf import TargetPixelFile
    from ..figure import K2Fig


def test_main_function():
    """The `k2pix` command-line tool calls the function `k2pix.main.k2pix()`,
    let's make sure it runs!"""
    from ..main import k2pix
    with pytest.raises(SystemExit):
        k2pix()
