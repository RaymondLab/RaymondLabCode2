"""Shared pytest fixtures.

pytest loads this file by itself, before the tests, so every test can use
what it holds.
"""

import os
import shutil
import tempfile
from pathlib import Path

import pytest


@pytest.fixture
def rig_root(monkeypatch: pytest.MonkeyPatch) -> Path:
    """Give the test an empty temp folder as the rig root.

    Points the RAYMONDLAB_RIG_ROOT environment variable at that folder, so
    no test touches the real C:\\RaymondLab\\ folder.

    The folder is made under /tmp when that exists, because a staging_root
    path may be at most 80 characters and pytest's own tmp_path on macOS is
    longer than that.
    """
    short_tmp = "/tmp" if os.path.isdir("/tmp") else None
    base = Path(tempfile.mkdtemp(prefix="rl-", dir=short_tmp))
    root = base / "RaymondLab"
    root.mkdir()
    monkeypatch.setenv("RAYMONDLAB_RIG_ROOT", str(root))
    yield root
    shutil.rmtree(base, ignore_errors=True)
