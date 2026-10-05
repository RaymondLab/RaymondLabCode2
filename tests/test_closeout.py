import hashlib
import json
from datetime import datetime, timedelta, timezone
from pathlib import Path

import pytest

from raymondlab.core.jsonio import read_json
from raymondlab.core.naming import run_stem
from raymondlab.rig.paths import CURRENT_SESSION_JSON
from raymondlab.session import closeout, model, store

FIXTURE = Path(__file__).parent / "fixtures" / "calibration_answers.json"
NOW = datetime(2026, 4, 10, 15, 5, 30, tzinfo=timezone(timedelta(hours=-7)))
CLOSE = datetime(2026, 4, 10, 16, 20, 0, tzinfo=timezone(timedelta(hours=-7)))
QC = {1: ("usable", None), 2: ("excluded", "hardware-failure")}
SMRX_ONE = b"fake spike2 data one"
SMRX_TWO = b"fake spike2 data two, a little longer"
BIN = bytes(range(256)) * 5


@pytest.fixture
def answers() -> dict:
    return json.loads(FIXTURE.read_text(encoding="utf-8"))


@pytest.fixture
def session(seeded_rig, answers) -> Path:
    """An open session with two runs and three fake data files."""
    directory = store.create_session(seeded_rig, answers, now=NOW)
    store.add_run(seeded_rig)
    store.add_run(seeded_rig)
    stem = store.read_current(seeded_rig)["entity_stem"]
    (directory / (run_stem(stem, 1) + ".smrx")).write_bytes(SMRX_ONE)
    (directory / (run_stem(stem, 2) + ".smrx")).write_bytes(SMRX_TWO)
    frames = directory / "frames" / "run-01"
    frames.mkdir(parents=True)
    (frames / (run_stem(stem, 1) + "_cam-1.bin")).write_bytes(BIN)
    return directory


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def test_close_lists_every_file_with_hash(seeded_rig, session):
    stem = store.read_current(seeded_rig)["entity_stem"]
    final = closeout.close_session(session, "bangeles", QC, "production", now=CLOSE)
    expected = {
        run_stem(stem, 1) + ".smrx": SMRX_ONE,
        run_stem(stem, 2) + ".smrx": SMRX_TWO,
        "frames/run-01/" + run_stem(stem, 1) + "_cam-1.bin": BIN,
    }
    files = final["closeout_files"]
    assert [f["path"] for f in files] == sorted(expected)
    for entry in files:
        data = expected[entry["path"]]
        assert "\\" not in entry["path"]
        assert entry["bytes"] == len(data)
        assert entry["sha256"] == _sha(data)
    assert set(files[0]) == {"path", "bytes", "sha256"}


def test_close_writes_closed_fields_and_history(seeded_rig, session):
    final = closeout.close_session(session, "bangeles", QC, "production", now=CLOSE)
    on_disk = read_json(session / "_metadata.json")
    assert on_disk == final
    assert final["closed_at"] == "2026-04-10T16:20:00-07:00"
    assert final["closed_by"] == "bangeles"
    assert final["runs"] == [{"run_index": 1}, {"run_index": 2}]
    assert final["state"]["qc"] == [
        {"run_index": 1, "verdict": "usable", "reason": None},
        {"run_index": 2, "verdict": "excluded", "reason": "hardware-failure"},
    ]
    assert len(final["state"]["lifecycle"]) == 2
    assert final["state"]["lifecycle"][1] == {
        "value": "production", "at": "2026-04-10T16:20:00-07:00", "by": "bangeles"}


def test_close_removes_current_session(seeded_rig, session):
    assert (seeded_rig / CURRENT_SESSION_JSON).exists()
    closeout.close_session(session, "bangeles", QC, "production", now=CLOSE)
    assert not (seeded_rig / CURRENT_SESSION_JSON).exists()


def test_close_keeps_pointer_to_another_session(seeded_rig, session, tmp_path):
    other = tmp_path / "other-session"
    store.write_current(seeded_rig, other, "sub-x_ses-y_task-eyecal", task="eyecal", next_run=1)
    closeout.close_session(session, "bangeles", QC, "production", now=CLOSE)
    assert (seeded_rig / CURRENT_SESSION_JSON).exists()


def test_close_without_pointer_file_is_fine(seeded_rig, session):
    (seeded_rig / CURRENT_SESSION_JSON).unlink()
    final = closeout.close_session(session, "bangeles", QC, "production", now=CLOSE)
    assert "closed_at" in final


def test_close_twice_raises(session):
    closeout.close_session(session, "bangeles", QC, "production", now=CLOSE)
    with pytest.raises(RuntimeError, match=str(session)):
        closeout.close_session(session, "bangeles", QC, "production", now=CLOSE)


def test_closed_sidecar_validates(session):
    final = closeout.close_session(session, "bangeles", QC, "production", now=CLOSE)
    sidecar = read_json(session / "_metadata.json")
    assert sidecar == final
    assert model.validate_sidecar(sidecar) == []


def test_metadata_and_tmp_files_are_not_hashed(session):
    (session / "stray.tmp").write_bytes(b"half written")
    final = closeout.close_session(session, "bangeles", QC, "production", now=CLOSE)
    paths = [f["path"] for f in final["closeout_files"]]
    assert len(paths) == 3
    assert "_metadata.json" not in paths
    assert "stray.tmp" not in paths


def test_large_file_is_hashed_in_chunks(session):
    data = b"abcdefghij" * (300 * 1024)   # a little under 3 MiB: more than one chunk
    (session / "big.bin").write_bytes(data)
    final = closeout.close_session(session, "bangeles", QC, "production", now=CLOSE)
    big = [f for f in final["closeout_files"] if f["path"] == "big.bin"][0]
    assert big["bytes"] == len(data)
    assert big["sha256"] == _sha(data)


@pytest.mark.parametrize("qc, word", [
    ({1: ("usable", None)}, "2"),                                              # run 2 missing
    ({1: ("usable", None), 2: ("usable", None), 3: ("usable", None)}, "3"),    # run 3 is extra
    ({1: ("fine", None), 2: ("usable", None)}, "fine"),                        # bad verdict
    ({1: ("usable", "hardware-failure"), 2: ("usable", None)}, "reason"),      # usable with a reason
    ({1: ("excluded", "not-a-reason"), 2: ("usable", None)}, "not-a-reason"),  # bad reason
])
def test_bad_qc_raises_and_writes_nothing(seeded_rig, session, qc, word):
    before = (session / "_metadata.json").read_bytes()
    with pytest.raises(ValueError, match=word):
        closeout.close_session(session, "bangeles", qc, "production", now=CLOSE)
    assert (session / "_metadata.json").read_bytes() == before
    assert (seeded_rig / CURRENT_SESSION_JSON).exists()


def test_bad_lifecycle_raises_and_writes_nothing(seeded_rig, session):
    before = (session / "_metadata.json").read_bytes()
    with pytest.raises(ValueError, match="lifecycle"):
        closeout.close_session(session, "bangeles", QC, "nonsense", now=CLOSE)
    assert (session / "_metadata.json").read_bytes() == before
    assert (seeded_rig / CURRENT_SESSION_JSON).exists()


def test_naive_now_is_given_an_offset(session):
    final = closeout.close_session(session, "bangeles", QC, "production",
                                   now=datetime(2026, 4, 10, 16, 20, 0))
    assert final["closed_at"][19:] != ""   # has a UTC offset after the seconds


def test_excluded_without_reason_names_the_run(seeded_rig, session):
    before = (session / "_metadata.json").read_bytes()
    with pytest.raises(ValueError, match="run 2 is excluded but has no reason"):
        closeout.close_session(session, "bangeles", {1: ("usable", None), 2: ("excluded", None)},
                               "production", now=CLOSE)
    assert (session / "_metadata.json").read_bytes() == before
    assert (seeded_rig / CURRENT_SESSION_JSON).exists()


def test_file_of_an_unregistered_run_is_refused(seeded_rig, answers):
    directory = store.create_session(seeded_rig, answers, now=NOW)
    store.add_run(seeded_rig)
    stem = store.read_current(seeded_rig)["entity_stem"]
    (directory / (run_stem(stem, 1) + ".smrx")).write_bytes(SMRX_ONE)
    stray = run_stem(stem, 2) + ".smrx"   # recorded without Add run
    (directory / stray).write_bytes(SMRX_TWO)
    before = (directory / "_metadata.json").read_bytes()
    with pytest.raises(ValueError, match=f"{stray} is for run 2, which the session does not have"):
        closeout.close_session(directory, "bangeles", {1: ("usable", None)}, "production", now=CLOSE)
    assert (directory / "_metadata.json").read_bytes() == before
    assert (seeded_rig / CURRENT_SESSION_JSON).exists()
