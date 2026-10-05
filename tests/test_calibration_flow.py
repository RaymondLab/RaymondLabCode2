import json
import os
from datetime import datetime, timedelta, timezone
from pathlib import Path

from raymondlab.core import facet_registry
from raymondlab.core.jsonio import read_json
from raymondlab.core.naming import MAX_DIRNAME, MAX_FILENAME, run_stem
from raymondlab.session import closeout, model, store

FIXTURE = Path(__file__).parent / "fixtures" / "calibration_answers.json"
NOW = datetime(2026, 4, 10, 15, 5, 30, tzinfo=timezone(timedelta(hours=-7)))
CLOSE = datetime(2026, 4, 10, 16, 20, 0, tzinfo=timezone(timedelta(hours=-7)))


def test_calibration_session_from_create_to_close(seeded_rig):
    answers = json.loads(FIXTURE.read_text(encoding="utf-8"))
    directory = store.create_session(seeded_rig, answers, now=NOW)
    assert store.add_run(seeded_rig) == 1
    assert store.add_run(seeded_rig) == 2

    stem = store.read_current(seeded_rig)["entity_stem"]
    (directory / (run_stem(stem, 1) + ".smrx")).write_bytes(b"one")
    (directory / (run_stem(stem, 2) + ".smrx")).write_bytes(b"two")
    frames = directory / "frames" / "run-01"
    frames.mkdir(parents=True)
    (frames / (run_stem(stem, 1) + "_cam-1.bin")).write_bytes(b"frames")

    closeout.close_session(directory, "bangeles",
                           {1: ("usable", None), 2: ("excluded", "hardware-failure")}, "production",
                           now=CLOSE)
    assert store.read_current(seeded_rig) is None

    # every file name is at most 64 characters, every directory name at most 40
    for _, dirnames, filenames in os.walk(directory):
        for name in dirnames:
            assert len(name) <= MAX_DIRNAME, name
        for name in filenames:
            assert len(name) <= MAX_FILENAME, name

    sidecar = read_json(directory / "_metadata.json")
    required = facet_registry.required_keys("common") + facet_registry.required_keys("eye-calibration")
    assert facet_registry.check_required(sidecar, required) == []
    assert model.validate_sidecar(sidecar) == []
