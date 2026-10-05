import copy
import json
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pytest

from raymondlab.core import facet_registry
from raymondlab.core.jsonio import read_json, write_json
from raymondlab.session import store

FIXTURE = Path(__file__).parent / "fixtures" / "calibration_answers.json"
NOW = datetime(2026, 4, 10, 15, 5, 30, tzinfo=timezone(timedelta(hours=-7)))
CURRENT_KEYS = {"session_dir", "entity_stem", "task", "next_run", "duration_s", "written_at"}


def leave_open_session(root: Path) -> None:
    """Drop current-session.json, so that the next create is not refused for an open session."""
    (root / "current-session.json").unlink()


@pytest.fixture
def answers() -> dict:
    return json.loads(FIXTURE.read_text(encoding="utf-8"))


@pytest.fixture
def session(seeded_rig, answers) -> Path:
    return store.create_session(seeded_rig, answers, now=NOW)


def test_create_makes_folder_and_valid_sidecar(session, answers):
    assert session.is_dir()
    assert session.parts[-3:] == ("raw", "sub-m7k222", "ses-20260410T1505")
    sidecar = read_json(session / "_metadata.json")
    assert facet_registry.check_required(
        sidecar, facet_registry.required_keys("common")
        + facet_registry.required_keys("eye-calibration")) == []
    assert sidecar["subject_id"] == "m7k222"
    assert sidecar["session_start_time"] == "2026-04-10T15:05:30-07:00"
    assert sidecar["runs"] == []
    assert sidecar["rig_id"] == "D241A"


def test_create_writes_flat_current_session(seeded_rig, session, answers):
    current = read_json(seeded_rig / "current-session.json")
    assert set(current) == CURRENT_KEYS
    assert all(isinstance(v, (str, int)) and not isinstance(v, bool) for v in current.values())
    assert "\\" not in current["session_dir"]
    assert current["session_dir"] == session.resolve().as_posix()
    assert current["entity_stem"] == "sub-m7k222_ses-20260410T1505_task-eyecal"
    assert current["task"] == "eyecal"
    assert current["next_run"] == 1
    assert current["duration_s"] == answers["session"]["duration_s"] == 120
    assert datetime.fromisoformat(current["written_at"]).utcoffset() is not None


def test_create_uses_given_subject_id_without_minting(seeded_rig, answers):
    answers["subject_id"] = "m7kabc"
    path = store.create_session(seeded_rig, answers, now=NOW)
    assert path.parent.name == "sub-m7kabc"
    assert read_json(path / "_metadata.json")["subject_id"] == "m7kabc"


def test_create_mints_next_id_for_a_new_mouse(seeded_rig, answers):
    first = store.create_session(seeded_rig, answers, now=NOW)
    leave_open_session(seeded_rig)
    answers["mouse"]["system_id"] = "777"   # another mouse: a known system_id is not minted again
    second = store.create_session(seeded_rig, answers, now=NOW)
    assert first.parent.name == "sub-m7k222"
    assert second.parent.name == "sub-m7k223"


def test_naive_now_gets_an_offset(seeded_rig, answers):
    path = store.create_session(seeded_rig, answers, now=datetime(2026, 4, 10, 15, 5, 30))
    start = read_json(path / "_metadata.json")["session_start_time"]
    assert datetime.fromisoformat(start).utcoffset() is not None


def test_default_now_is_current_time(seeded_rig, answers):
    path = store.create_session(seeded_rig, answers)
    start = datetime.fromisoformat(read_json(path / "_metadata.json")["session_start_time"])
    assert start.utcoffset() is not None
    assert abs(datetime.now().astimezone() - start) < timedelta(minutes=1)


def test_same_minute_twice_refuses(seeded_rig, answers):
    answers["subject_id"] = "m7k222"
    path = store.create_session(seeded_rig, answers, now=NOW)
    leave_open_session(seeded_rig)
    with pytest.raises(FileExistsError, match="ses-20260410T1505"):
        store.create_session(seeded_rig, answers, now=NOW)
    assert (path / "_metadata.json").exists()


def test_unknown_project_names_it(seeded_rig, answers):
    answers["project_id"] = "nope2026"
    with pytest.raises(ValueError, match="nope2026"):
        store.create_session(seeded_rig, answers, now=NOW)


def test_bad_answers_write_nothing(seeded_rig, answers):
    bad = copy.deepcopy(answers)
    bad["mouse"]["sex"] = "X"
    with pytest.raises(ValueError):
        store.create_session(seeded_rig, bad, now=NOW)
    assert not (seeded_rig / "current-session.json").exists()
    assert list((seeded_rig / "staging").glob("*/raw/*")) == []


def test_read_current_none_when_no_session(seeded_rig):
    assert store.read_current(seeded_rig) is None


def test_read_current_returns_the_file(seeded_rig, session):
    assert store.read_current(seeded_rig) == read_json(seeded_rig / "current-session.json")


def test_add_run_twice(seeded_rig, session):
    assert store.add_run(seeded_rig) == 1
    assert store.add_run(seeded_rig) == 2
    sidecar = read_json(session / "_metadata.json")
    assert sidecar["runs"] == [{"run_index": 1}, {"run_index": 2}]
    current = read_json(seeded_rig / "current-session.json")
    assert current["next_run"] == 3
    assert current["duration_s"] == 120
    assert current["task"] == "eyecal"
    assert set(current) == CURRENT_KEYS


def test_add_run_without_session_raises(seeded_rig):
    with pytest.raises(RuntimeError, match="no session is open"):
        store.add_run(seeded_rig)


def test_add_run_on_closed_session_raises(seeded_rig, session):
    sidecar = read_json(session / "_metadata.json")
    sidecar["closed_at"] = "2026-04-10T16:00:00-07:00"
    write_json(session / "_metadata.json", sidecar)
    with pytest.raises(RuntimeError, match="sub-m7k222_ses-20260410T1505_task-eyecal"):
        store.add_run(seeded_rig)
    assert read_json(seeded_rig / "current-session.json")["next_run"] == 1


def test_write_current_casts_duration_to_int(seeded_rig, tmp_path):
    store.write_current(seeded_rig, tmp_path / "ses", "stem", "eyecal", 4, duration_s=90.0)
    current = read_json(seeded_rig / "current-session.json")
    assert current["duration_s"] == 90
    assert isinstance(current["duration_s"], int)
    assert current["next_run"] == 4


@pytest.mark.parametrize("duration", ["abc", None, 90.5, True])
def test_bad_duration_writes_nothing(seeded_rig, answers, duration):
    answers["session"]["duration_s"] = duration
    with pytest.raises(ValueError, match="duration_s"):
        store.create_session(seeded_rig, answers, now=NOW)
    assert not (seeded_rig / "current-session.json").exists()
    assert list((seeded_rig / "staging").glob("*/raw/*")) == []


def test_missing_duration_writes_nothing(seeded_rig, answers):
    del answers["session"]["duration_s"]
    with pytest.raises(ValueError, match="duration_s"):
        store.create_session(seeded_rig, answers, now=NOW)
    assert not (seeded_rig / "current-session.json").exists()
    assert list((seeded_rig / "staging").glob("*/raw/*")) == []


def test_find_mouse_by_padded_system_id(seeded_rig, answers):
    answers["mouse"]["system_id"] = "1218"
    store.create_session(seeded_rig, answers, now=NOW)
    found = store.find_mouse(seeded_rig / "staging", "00000000001218")
    assert found["subject_id"] == "m7k222"
    assert found["colony_mouse_id"] == "RL-0421"
    assert store.find_mouse(seeded_rig / "staging", "00000000009999") is None


def test_last_sidecar_is_the_most_recent_session(seeded_rig, answers):
    assert store.last_sidecar(seeded_rig / "staging") is None
    later = copy.deepcopy(answers)
    later["mouse"]["system_id"] = "777"
    later["session"]["preamp_gain"] = 80
    # created first but started later: the order comes from session_start_time, not the folder
    store.create_session(seeded_rig, later, now=NOW + timedelta(days=1))
    leave_open_session(seeded_rig)
    store.create_session(seeded_rig, answers, now=NOW)
    last = store.last_sidecar(seeded_rig / "staging")
    assert last["system_id"] == "00000000000777"
    assert last["preamp_gain"] == 80


def test_create_refuses_while_a_session_is_open(seeded_rig, session, answers):
    with pytest.raises(RuntimeError, match="a session is still open: .*ses-20260410T1505"):
        store.create_session(seeded_rig, answers, now=NOW + timedelta(minutes=1))
    assert store.read_current(seeded_rig)["session_dir"] == session.resolve().as_posix()
    assert [p.name for p in (seeded_rig / "staging").glob("*/raw/*")] == ["sub-m7k222"]


def test_create_ignores_a_pointer_to_a_closed_or_missing_session(seeded_rig, session, answers):
    sidecar = read_json(session / "_metadata.json")
    sidecar["closed_at"] = "2026-04-10T16:00:00-07:00"
    write_json(session / "_metadata.json", sidecar)
    answers["subject_id"] = "m7k222"   # the same mouse again: its system_id is known
    second = store.create_session(seeded_rig, answers, now=NOW + timedelta(minutes=1))
    store.write_current(seeded_rig, seeded_rig / "staging" / "gone", "stem", "eyecal", 1)
    third = store.create_session(seeded_rig, answers, now=NOW + timedelta(minutes=2))
    assert second.is_dir() and third.is_dir()


def test_known_system_id_without_subject_id_is_refused(seeded_rig, answers):
    answers["mouse"]["system_id"] = "1218"
    store.create_session(seeded_rig, answers, now=NOW)
    leave_open_session(seeded_rig)
    answers["mouse"]["system_id"] = "0001218"   # the same mouse once padded
    with pytest.raises(ValueError,
                       match="system_id 00000000001218 already belongs to subject m7k222; pick that mouse"):
        store.create_session(seeded_rig, answers, now=NOW + timedelta(minutes=1))
    assert not (seeded_rig / "current-session.json").exists()
    assert [p.name for p in (seeded_rig / "staging").glob("*/raw/sub-*/ses-*")] == ["ses-20260410T1505"]


def test_known_system_id_with_another_subject_id_is_refused(seeded_rig, answers):
    store.create_session(seeded_rig, answers, now=NOW)
    leave_open_session(seeded_rig)
    answers["subject_id"] = "m7kabc"
    with pytest.raises(ValueError, match="system_id 20262280001218 belongs to subject m7k222, "
                                         "not to subject m7kabc"):
        store.create_session(seeded_rig, answers, now=NOW + timedelta(minutes=1))
    assert not (seeded_rig / "current-session.json").exists()
    assert [p.name for p in (seeded_rig / "staging").glob("*/raw/*")] == ["sub-m7k222"]
    assert len(list((seeded_rig / "staging").glob("*/raw/sub-*/ses-*"))) == 1


def test_known_system_id_with_its_own_subject_id_creates(seeded_rig, answers):
    first = store.create_session(seeded_rig, answers, now=NOW)
    leave_open_session(seeded_rig)
    answers["subject_id"] = "m7k222"
    second = store.create_session(seeded_rig, answers, now=NOW + timedelta(minutes=1))
    assert second.parent == first.parent
    assert read_json(second / "_metadata.json")["subject_id"] == "m7k222"


def test_add_run_takes_the_index_from_the_sidecar_when_the_pointer_is_stale(seeded_rig, session):
    store.add_run(seeded_rig)
    store.add_run(seeded_rig)
    # as after a crash between the sidecar write and the pointer write of the second add_run
    current = store.read_current(seeded_rig)
    store.write_current(seeded_rig, session, current["entity_stem"], current["task"], 2,
                        current["duration_s"])
    assert store.add_run(seeded_rig) == 3
    assert read_json(session / "_metadata.json")["runs"] == [
        {"run_index": 1}, {"run_index": 2}, {"run_index": 3}]
    assert store.read_current(seeded_rig)["next_run"] == 4


@pytest.mark.parametrize("duration", [0, -1, -120])
def test_zero_or_negative_duration_writes_nothing(seeded_rig, answers, duration):
    answers["session"]["duration_s"] = duration
    with pytest.raises(ValueError, match="duration_s must be more than 0 seconds"):
        store.create_session(seeded_rig, answers, now=NOW)
    assert not (seeded_rig / "current-session.json").exists()
    assert list((seeded_rig / "staging").glob("*/raw/*")) == []
