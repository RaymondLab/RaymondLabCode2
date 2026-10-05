import pytest

from raymondlab.registry import mint

SLOT = "7k"


def make_subject(staging_root, name, dataset="okr2026"):
    d = staging_root / dataset / "raw" / name
    d.mkdir(parents=True)
    return d


def test_empty_staging_gives_first_id(tmp_path):
    assert mint.next_subject_id(SLOT, tmp_path) == "m7k222"


def test_missing_staging_root_gives_first_id(tmp_path):
    assert mint.next_subject_id(SLOT, tmp_path / "nope") == "m7k222"


def test_next_is_highest_plus_one(tmp_path):
    for tail in ("222", "223", "22a"):
        make_subject(tmp_path, f"sub-m{SLOT}{tail}")
    assert mint.next_subject_id(SLOT, tmp_path) == "m7k22b"


def test_scans_every_dataset(tmp_path):
    make_subject(tmp_path, "sub-m7k223", dataset="okr2026")
    make_subject(tmp_path, "sub-m7k225", dataset="vor2026")
    assert mint.next_subject_id(SLOT, tmp_path) == "m7k226"


def test_foreign_slot_is_ignored(tmp_path):
    make_subject(tmp_path, "sub-m9q222")
    make_subject(tmp_path, "sub-m9qzzz")
    assert mint.next_subject_id(SLOT, tmp_path) == "m7k222"


def test_files_and_odd_names_are_ignored(tmp_path):
    raw = tmp_path / "okr2026" / "raw"
    raw.mkdir(parents=True)
    (raw / "sub-m7kzzz").write_text("a file, not a directory")
    (raw / "sub-m7k22").mkdir()       # tail too short
    (raw / "sub-m7k2222").mkdir()     # tail too long
    (raw / "sub-m7k2o2").mkdir()      # "o" is not in the alphabet
    (raw / "notes").mkdir()
    make_subject(tmp_path, "sub-m7k224")
    assert mint.next_subject_id(SLOT, tmp_path) == "m7k225"


def test_full_slot_raises(tmp_path):
    make_subject(tmp_path, "sub-m7kzzz")
    with pytest.raises(mint.SlotFull):
        mint.next_subject_id(SLOT, tmp_path)


def test_last_free_value_is_still_minted(tmp_path):
    make_subject(tmp_path, "sub-m7kzzy")
    assert mint.next_subject_id(SLOT, tmp_path) == "m7kzzz"
