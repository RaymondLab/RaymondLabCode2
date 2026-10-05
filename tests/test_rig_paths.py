from pathlib import Path

from raymondlab.rig import paths


def test_rig_root_uses_the_environment_override(rig_root: Path):
    assert paths.rig_root() == rig_root


def test_rig_root_default_is_the_windows_lab_folder(monkeypatch):
    monkeypatch.delenv("RAYMONDLAB_RIG_ROOT", raising=False)
    assert paths.rig_root() == Path(r"C:\RaymondLab")


def test_file_names():
    assert paths.RIG_JSON == "rig.json"
    assert paths.PROJECTS_JSON == "projects.json"
    assert paths.CURRENT_SESSION_JSON == "current-session.json"
