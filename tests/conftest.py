"""Shared pytest fixtures.

pytest loads this file by itself, before the tests, so every test can use
what it holds.
"""

from pathlib import Path

import pytest


@pytest.fixture
def rig_root(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Give the test an empty temp folder as the rig root.

    Points the RAYMONDLAB_RIG_ROOT environment variable at that folder, so
    no test touches the real C:\\RaymondLab\\ folder.
    """
    root = tmp_path / "RaymondLab"
    root.mkdir()
    monkeypatch.setenv("RAYMONDLAB_RIG_ROOT", str(root))
    return root
