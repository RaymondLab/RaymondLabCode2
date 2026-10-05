from datetime import datetime
from pathlib import Path

import pytest

from raymondlab.rig import config, setup


def answers(staging_root: Path, **over) -> dict:
    a = dict(rig_id="D241A", mint_slot="7k", staging_root=str(staging_root),
             default_experimenter="bangeles",
             projects=[dict(project_id="eyecal-2026", dataset_dir="eyecal_2026",
                            protocol_id="APLAC-12345", lab="Raymond", institution="Stanford")])
    a.update(over)
    return a


def test_first_run_creates_root_staging_and_both_files(rig_root: Path):
    root = rig_root / "fresh"          # does not exist yet
    staging = root / "staging"
    cfg = setup.apply_setup(root, answers(staging))
    assert root.is_dir()
    assert staging.is_dir()
    assert (root / "rig.json").is_file()
    assert (root / "projects.json").is_file()
    assert cfg == config.load_rig(root)
    assert cfg.mint_slot == "7k"
    assert cfg.schema_version == "1.0.0"
    assert datetime.fromisoformat(cfg.created_at).utcoffset() is not None
    assert [p.project_id for p in config.load_projects(root)] == ["eyecal-2026"]


def test_rerun_with_new_slot_on_unlocked_rig_succeeds(rig_root: Path):
    staging = rig_root / "staging"
    first = setup.apply_setup(rig_root, answers(staging))
    second = setup.apply_setup(rig_root, answers(staging, mint_slot="9q"))
    assert second.mint_slot == "9q"
    assert second.created_at == first.created_at
    assert config.load_rig(rig_root).mint_slot == "9q"


def test_rerun_with_new_slot_on_locked_rig_raises(rig_root: Path):
    staging = rig_root / "staging"
    setup.apply_setup(rig_root, answers(staging))
    (staging / "eyecal_2026" / "raw" / "sub-m7k222").mkdir(parents=True)
    with pytest.raises(setup.SlotLocked):
        setup.apply_setup(rig_root, answers(staging, mint_slot="9q"))
    assert config.load_rig(rig_root).mint_slot == "7k"


def test_rerun_same_slot_on_locked_rig_is_fine(rig_root: Path):
    staging = rig_root / "staging"
    setup.apply_setup(rig_root, answers(staging))
    (staging / "eyecal_2026" / "raw" / "sub-m7k222").mkdir(parents=True)
    cfg = setup.apply_setup(rig_root, answers(staging, default_experimenter="someone"))
    assert cfg.default_experimenter == "someone"


def test_rerun_adding_a_project_leaves_rig_json_bytes_unchanged(rig_root: Path):
    staging = rig_root / "staging"
    a = answers(staging)
    setup.apply_setup(rig_root, a)
    before = (rig_root / "rig.json").read_bytes()
    a["projects"].append(dict(project_id="okr-2026", dataset_dir="okr_2026",
                              protocol_id="APLAC-12345", lab="Raymond", institution="Stanford"))
    setup.apply_setup(rig_root, a)
    assert (rig_root / "rig.json").read_bytes() == before
    assert [p.project_id for p in config.load_projects(rig_root)] == ["eyecal-2026", "okr-2026"]


def test_bad_answers_raise_with_every_error_and_write_nothing(rig_root: Path):
    root = rig_root / "fresh"
    a = answers(root / "staging", rig_id="", mint_slot="o1")
    a["projects"].append(dict(a["projects"][0]))   # duplicate project_id
    with pytest.raises(ValueError) as exc:
        setup.apply_setup(root, a)
    msg = str(exc.value)
    assert "rig_id" in msg
    assert "mint_slot" in msg
    assert "project_id" in msg
    assert msg.count("\n") >= 2
    assert not root.exists()
