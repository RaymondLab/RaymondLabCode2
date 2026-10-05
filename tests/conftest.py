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


@pytest.fixture
def seeded_rig(rig_root: Path) -> Path:
    """The rig_root folder, set up as rig D241A with one project (okr2026). Returns the root."""
    from raymondlab.rig.setup import apply_setup

    apply_setup(rig_root, {
        "rig_id": "D241A",
        "mint_slot": "7k",
        "staging_root": str(rig_root / "staging"),
        "default_experimenter": "bangeles",
        "projects": [{"project_id": "okr2026", "dataset_dir": "RaymondLab_OKR2026",
                      "protocol_id": "APLAC-12345", "lab": "Raymond Lab",
                      "institution": "Stanford University"}],
    })
    return rig_root
