"""Every name the pipeline writes, in one place.

Session folders and files follow one pattern:
sub-<subject_id>_ses-<YYYYMMDDThhmm>_task-<task>. Run files add _run-NN.
File names are capped at 64 characters and directory names at 40, so the full
path stays inside the limits of the archive tools. Nothing else in the code
base builds a name by hand.
"""

from datetime import datetime
from pathlib import Path

MAX_FILENAME = 64  # limit set by the lab data standard
MAX_DIRNAME = 40


def ses_stamp(t: datetime) -> str:
    """Session time stamp to the minute, no seconds: 2026-04-10 15:05 -> "20260410T1505".

    Minute resolution is enough. A second attempt in the same minute is a new
    run, not a new session.
    """
    return t.strftime("%Y%m%dT%H%M")


def entity_stem(subject_id: str, stamp: str, task: str) -> str:
    """The shared prefix of every file in a session.

    ("m7k3q9", "20260410T1505", "eyecal") -> "sub-m7k3q9_ses-20260410T1505_task-eyecal"
    """
    return f"sub-{subject_id}_ses-{stamp}_task-{task}"


def run_stem(stem: str, run: int) -> str:
    """The stem plus "_run-01" (two digits)."""
    return f"{stem}_run-{run:02d}"


def check_filename(name: str) -> None:
    """Raise ValueError, naming the file, if the name is over 64 characters."""
    if len(name) > MAX_FILENAME:
        raise ValueError(
            f"file name {name!r} is {len(name)} characters; the limit is {MAX_FILENAME}"
        )


def check_dirname(name: str) -> None:
    """Raise ValueError, naming the directory, if the name is over 40 characters."""
    if len(name) > MAX_DIRNAME:
        raise ValueError(
            f"directory name {name!r} is {len(name)} characters; the limit is {MAX_DIRNAME}"
        )


def session_dir(staging_root: Path, dataset_dir: str, subject_id: str, stamp: str) -> Path:
    """staging_root / dataset_dir / "raw" / "sub-<id>" / "ses-<stamp>"."""
    sub = f"sub-{subject_id}"
    ses = f"ses-{stamp}"
    for name in (dataset_dir, sub, ses):
        check_dirname(name)
    return staging_root / dataset_dir / "raw" / sub / ses
