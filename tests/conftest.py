# tests/conftest.py
"""Shared pytest fixtures. Pytest loads this file by itself; nothing imports it.

Every test that touches the rig folder asks for `rig_root`. It points RAYMONDLAB_RIG_ROOT at a
temp folder, so no test ever reads or writes the real rig folder on C:\\.
"""

from pathlib import Path
import pytest


@pytest.fixture
def rig_root(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """A fresh, empty rig root for one test.

    The folder is not created here; the code under test must create it, and one test checks
    that it does. `tmp_path` is deleted and the env change is undone after the test.
    """
    root = tmp_path / "RaymondLab"
    monkeypatch.setenv("RAYMONDLAB_RIG_ROOT", str(root))
    return root
