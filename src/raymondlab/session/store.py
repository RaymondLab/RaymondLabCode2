"""The session on disk, and the file Spike2 reads.

This module writes the session folder and keeps current-session.json up to
date. That file tells the Spike2 recording script which session is open, what
the stem is, and which run comes next. Spike2's script language has no JSON
parser; it finds values by string search. So the file stays flat: one level,
strings and integers only, forward slashes in paths, no free text.
"""

from datetime import datetime
from pathlib import Path

from raymondlab.core.jsonio import read_json, write_json
from raymondlab.core.naming import entity_stem, ses_stamp, session_dir as build_session_dir
from raymondlab.registry.mint import next_subject_id
from raymondlab.rig.config import load_projects, load_rig
from raymondlab.rig.paths import CURRENT_SESSION_JSON
from raymondlab.session.model import build_calibration_sidecar, pad_system_id

SIDECAR_NAME = "_metadata.json"
TASK = "eyecal"


def _duration_s(value: object) -> int:
    """The typed duration as a whole number of seconds, more than 0. Raise ValueError naming
    duration_s otherwise."""
    if isinstance(value, bool) or (isinstance(value, float) and not value.is_integer()):
        raise ValueError(f"duration_s must be a whole number of seconds: {value!r}")
    try:
        seconds = int(value)
    except (TypeError, ValueError):
        raise ValueError(f"duration_s must be a whole number of seconds: {value!r}") from None
    if seconds <= 0:   # 0 is kept for experiment sessions, which have no fixed run length
        raise ValueError(f"duration_s must be more than 0 seconds: {value!r}")
    return seconds


def _check_known_mouse(staging_root: Path, answers: dict) -> None:
    """Raise ValueError when the typed system_id already belongs to a mouse on this rig and the
    answers do not name that mouse's subject_id. One mouse keeps one subject_id.

    A system_id that cannot be padded is not checked here; the model reports it with the other errors.
    """
    try:
        padded = pad_system_id(answers["mouse"].get("system_id"))
    except ValueError:
        return
    found = find_mouse(staging_root, padded)
    if found is None:
        return
    given = answers.get("subject_id")
    if not given:
        raise ValueError(f"system_id {padded} already belongs to subject {found['subject_id']}; "
                         "pick that mouse")
    if given != found["subject_id"]:
        raise ValueError(f"system_id {padded} belongs to subject {found['subject_id']}, "
                         f"not to subject {given}; pick that mouse")


def create_session(root: Path, answers: dict, now: datetime | None = None) -> Path:
    """Create a calibration session folder and return its path.

    root is the rig root: it holds rig.json, projects.json, current-session.json and staging/.
    answers holds the form's values: project_id, mouse, session, and optionally subject_id
    (an existing mouse, which keeps its own id; a new mouse gets the next free id).

    Raise ValueError for an unknown project, for any bad answer, and for a system_id that already
    belongs to a mouse on this rig when subject_id is missing or is not that mouse's subject_id.
    Raise FileExistsError when the session folder is already there. Raise RuntimeError when the
    rig is not set up or another session is still open. Nothing is written in any of these cases.
    """
    root = Path(root)
    cfg = load_rig(root)
    if cfg is None:
        raise RuntimeError(f"the rig is not set up: no rig.json in {root}")
    still_open = open_session(root)
    if still_open is not None:   # a new session would overwrite current-session.json
        raise RuntimeError(f"a session is still open: {still_open}; close it first")
    project_id = answers.get("project_id")
    matches = [p for p in load_projects(root) if p.project_id == project_id]
    if not matches:
        raise ValueError(f"unknown project_id {project_id!r}")
    project = matches[0]

    staging_root = Path(cfg.staging_root)
    _check_known_mouse(staging_root, answers)   # before minting: a known mouse gets no new id
    subject_id = answers.get("subject_id") or next_subject_id(cfg.mint_slot, staging_root)

    duration_s = _duration_s(answers["session"].get("duration_s"))   # check before any write

    now = (now or datetime.now()).astimezone()   # a naive time is taken as local time
    stamp = ses_stamp(now)
    stem = entity_stem(subject_id, stamp, TASK)
    directory = build_session_dir(staging_root, project.dataset_dir, subject_id, stamp)
    if directory.exists():
        raise FileExistsError(f"session folder already exists: {directory}")

    session = dict(answers["session"], duration_s=duration_s,
                   session_start_time=now.isoformat(timespec="seconds"))
    sidecar = build_calibration_sidecar(cfg, project, answers["mouse"], session, subject_id)

    directory.mkdir(parents=True)
    write_json(directory / SIDECAR_NAME, sidecar)
    write_current(root, directory, stem, task=TASK, next_run=1, duration_s=duration_s)
    return directory


def add_run(root: Path) -> int:
    """Register the next run of the open session and return its index.

    The sidecar is the authority: the new index is one more than the highest run_index in its runs
    list (1 when it has none). Add {"run_index": index} to that list, write next_run = index + 1 to
    current-session.json, and rewrite both files. So a crash between the two writes never records
    the same run twice. Raise RuntimeError when no session is open or the open session is closed.
    """
    root = Path(root)
    current = read_current(root)
    if current is None:
        raise RuntimeError("no session is open")
    directory = Path(current["session_dir"])
    sidecar = read_json(directory / SIDECAR_NAME)
    if "closed_at" in sidecar:
        raise RuntimeError(f"session {current['entity_stem']} is closed")

    run_index = max((run["run_index"] for run in sidecar["runs"]), default=0) + 1
    sidecar["runs"].append({"run_index": run_index})
    write_json(directory / SIDECAR_NAME, sidecar)
    write_current(root, directory, current["entity_stem"], current["task"], run_index + 1,
                  current["duration_s"])
    return run_index


def write_current(root: Path, session_dir: Path, stem: str, task: str, next_run: int,
                  duration_s: int = 0) -> None:
    """Write root / "current-session.json".

    Keys: session_dir (absolute, forward slashes), entity_stem, task, next_run, duration_s,
    written_at (ISO-8601 with the UTC offset). duration_s is the typed recording time of one
    calibration run; Spike2 sets its block length from it. An experiment session writes 0.

    Keep it flat: one level, strings and integers only. Spike2 reads this file by string search,
    so nesting or a backslash breaks the reader.
    """
    write_json(Path(root) / CURRENT_SESSION_JSON, {
        "session_dir": Path(session_dir).resolve().as_posix(),
        "entity_stem": stem,
        "task": task,
        "next_run": int(next_run),
        "duration_s": int(duration_s),
        "written_at": datetime.now().astimezone().isoformat(timespec="seconds"),
    })


def read_current(root: Path) -> dict | None:
    """The parsed current-session.json, or None when no session is open."""
    path = Path(root) / CURRENT_SESSION_JSON
    if not path.exists():
        return None
    return read_json(path)


def open_session(root: Path) -> Path | None:
    """The folder of the session that current-session.json names, while that session is open.

    None when there is no pointer file, when its folder has no sidecar any more, or when the
    sidecar has closed_at.
    """
    current = read_current(root)
    if current is None:
        return None
    directory = Path(current["session_dir"])
    sidecar = directory / SIDECAR_NAME
    if not sidecar.exists() or "closed_at" in read_json(sidecar):
        return None
    return directory


def _sidecars_by_start(staging_root: Path) -> list[dict]:
    """Every session sidecar under staging_root, oldest first by session_start_time."""
    found = [read_json(path) for path in
             Path(staging_root).glob(f"*/raw/sub-*/ses-*/{SIDECAR_NAME}")]
    # compare as times, not as text: two sidecars may carry different UTC offsets
    return sorted(found, key=lambda sidecar: datetime.fromisoformat(sidecar["session_start_time"]))


def find_mouse(staging_root: Path, system_id: str) -> dict | None:
    """The most recent sidecar under staging_root whose system_id is this one, or None.

    system_id must already be padded to 14 digits; the sidecars store it padded.
    """
    matches = [s for s in _sidecars_by_start(staging_root) if s.get("system_id") == system_id]
    return matches[-1] if matches else None


def last_sidecar(staging_root: Path) -> dict | None:
    """The sidecar of the most recent session under staging_root (any mouse), or None.

    Each rig has its own staging_root, so this is the rig's last session.
    """
    sidecars = _sidecars_by_start(staging_root)
    return sidecars[-1] if sidecars else None
