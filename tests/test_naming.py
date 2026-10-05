from datetime import datetime
from pathlib import Path

import pytest

from raymondlab.core import naming

STAMP = "20260410T1505"


def test_limits():
    assert naming.MAX_FILENAME == 64
    assert naming.MAX_DIRNAME == 40


def test_ses_stamp_drops_seconds():
    assert naming.ses_stamp(datetime(2026, 4, 10, 15, 5)) == STAMP
    assert naming.ses_stamp(datetime(2026, 4, 10, 15, 5, 42)) == STAMP
    # seconds are dropped, never rounded up
    assert naming.ses_stamp(datetime(2026, 4, 10, 15, 5, 59, 999999)) == STAMP


def test_entity_stem_exact_and_length():
    stem = naming.entity_stem("m7k3q9", STAMP, "eyecal")
    assert stem == "sub-m7k3q9_ses-20260410T1505_task-eyecal"
    assert len(stem) == 40


def test_run_stem_exact_and_length():
    stem = naming.entity_stem("m7k3q9", STAMP, "eyecal")
    rs = naming.run_stem(stem, 1)
    assert rs == "sub-m7k3q9_ses-20260410T1505_task-eyecal_run-01"
    assert len(rs) == 47
    assert naming.run_stem(stem, 12).endswith("_run-12")


def test_check_filename_accepts_63_and_64():
    rs = naming.run_stem(naming.entity_stem("m7k3q9", STAMP, "eyecal"), 1)
    name = rs + "_cam-1_first.tif"
    assert len(name) == 63
    naming.check_filename(name)
    naming.check_filename("a" * 64)


def test_check_filename_rejects_65_and_names_the_file():
    name = "a" * 65
    with pytest.raises(ValueError) as exc:
        naming.check_filename(name)
    assert name in str(exc.value)


def test_check_dirname_limits():
    naming.check_dirname("d" * 40)
    name = "d" * 41
    with pytest.raises(ValueError) as exc:
        naming.check_dirname(name)
    assert name in str(exc.value)


def test_session_dir_layout():
    root = Path("/stage")
    got = naming.session_dir(root, "okr2026", "m7k3q9", STAMP)
    assert got == root / "okr2026" / "raw" / "sub-m7k3q9" / "ses-20260410T1505"


def test_session_dir_checks_each_new_segment():
    root = Path("/stage")
    with pytest.raises(ValueError):
        naming.session_dir(root, "d" * 41, "m7k3q9", STAMP)
    with pytest.raises(ValueError):
        naming.session_dir(root, "okr2026", "x" * 40, STAMP)
    with pytest.raises(ValueError):
        naming.session_dir(root, "okr2026", "m7k3q9", "t" * 40)
