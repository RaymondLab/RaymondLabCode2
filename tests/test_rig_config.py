from dataclasses import replace
from pathlib import Path

import pytest

from raymondlab.core import jsonio
from raymondlab.rig import config
from raymondlab.rig.config import Project, RigConfig


def good_rig(staging_root: Path) -> RigConfig:
    return RigConfig(schema_version="1.0.0", rig_id="D241A", mint_slot="7k",
                     staging_root=str(staging_root), default_experimenter="bangeles",
                     created_at="2026-09-21T10:00:00-07:00")


def good_project(**over) -> Project:
    base = dict(project_id="eyecal-2026", dataset_dir="eyecal_2026", protocol_id="APLAC-12345",
                lab="Raymond", institution="Stanford")
    base.update(over)
    return Project(**base)


def test_load_rig_returns_none_when_missing(rig_root: Path):
    assert config.load_rig(rig_root) is None


def test_save_load_rig_round_trip(rig_root: Path):
    cfg = good_rig(rig_root / "staging")
    config.save_rig(rig_root, cfg)
    assert config.load_rig(rig_root) == cfg
    on_disk = jsonio.read_json(rig_root / "rig.json")
    assert on_disk["mint_slot"] == "7k"
    assert on_disk["schema_version"] == "1.0.0"


def test_save_load_projects_round_trip(rig_root: Path):
    rig_root.mkdir(exist_ok=True)
    ps = [good_project(), good_project(project_id="okr-2026", dataset_dir="okr_2026")]
    config.save_projects(rig_root, ps)
    assert config.load_projects(rig_root) == ps
    on_disk = jsonio.read_json(rig_root / "projects.json")
    assert on_disk["schema_version"] == "1.0.0"
    assert len(on_disk["projects"]) == 2


def test_validate_rig_good(rig_root: Path):
    assert config.validate_rig(good_rig(rig_root / "staging")) == []


@pytest.mark.parametrize("field,value,word", [
    ("mint_slot", "7kk", "slot"),
    ("mint_slot", "o1", "slot"),
    ("staging_root", "C:\\" + "x" * 80, "staging_root"),
    ("rig_id", "", "rig_id"),
    ("staging_root", "staging", "staging_root"),
    ("staging_root", "RaymondLab\\staging", "staging_root"),
])
def test_validate_rig_catches_bad_values(rig_root: Path, field, value, word):
    cfg = replace(good_rig(rig_root / "staging"), **{field: value})
    errors = config.validate_rig(cfg)
    assert errors, f"{field}={value!r} should be rejected"
    assert any(word in e for e in errors)


def test_validate_rig_accepts_a_windows_absolute_staging_root(rig_root: Path):
    assert config.validate_rig(replace(good_rig(rig_root), staging_root=r"C:\RaymondLab\staging")) == []


def test_validate_project_good():
    assert config.validate_project(good_project(), []) == []


@pytest.mark.parametrize("over,word", [
    ({"dataset_dir": "d" * 41}, "dataset_dir"),
    ({"dataset_dir": "has space"}, "dataset_dir"),
    ({"dataset_dir": "eyecal/2026"}, "dataset_dir"),
    ({"dataset_dir": "a\\b"}, "dataset_dir"),
    ({"dataset_dir": ".."}, "dataset_dir"),
    ({"project_id": "Bad_Id"}, "project_id"),
    ({"protocol_id": ""}, "protocol_id"),
])
def test_validate_project_catches_bad_values(over, word):
    errors = config.validate_project(good_project(**over), [])
    assert errors
    assert any(word in e for e in errors)


def test_validate_project_rejects_duplicate_id():
    errors = config.validate_project(good_project(), [good_project()])
    assert any("project_id" in e for e in errors)


def test_slot_is_locked_false_on_empty_staging(rig_root: Path):
    staging = rig_root / "staging"
    staging.mkdir(parents=True)
    assert config.slot_is_locked(good_rig(staging)) is False


def test_slot_is_locked_true_after_a_minted_id(rig_root: Path):
    staging = rig_root / "staging"
    (staging / "X" / "raw" / "sub-m7k222").mkdir(parents=True)
    assert config.slot_is_locked(good_rig(staging)) is True


def test_slot_is_locked_ignores_other_slots(rig_root: Path):
    staging = rig_root / "staging"
    (staging / "X" / "raw" / "sub-m9q222").mkdir(parents=True)
    assert config.slot_is_locked(good_rig(staging)) is False
