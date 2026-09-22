# tests/test_package.py
"""Smoke test: ensures the package installs and imports under the project environment."""

import raymondlab


def test_package_imports_and_has_a_version():
    assert isinstance(raymondlab.__version__, str)
    assert raymondlab.__version__


def test_subpackages_import():
    import raymondlab.apps
    import raymondlab.core
    import raymondlab.gui
    import raymondlab.rig
    import raymondlab.session
