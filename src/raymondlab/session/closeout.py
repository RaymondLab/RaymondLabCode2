"""Close a session once.

Closing records the QC verdict of each run, hashes every data file so that
later analysis can prove the raw data is unchanged, stamps closed_at, and
removes current-session.json so that Spike2 sees no open session. A closed
session is never reopened.
"""

import hashlib
import os
import re
from datetime import datetime
from pathlib import Path

from raymondlab.core.jsonio import read_json, write_json
from raymondlab.rig.paths import CURRENT_SESSION_JSON, rig_root
from raymondlab.session.model import validate_sidecar
from raymondlab.session.store import SIDECAR_NAME, read_current

CHUNK_BYTES = 1024 * 1024   # the camera .bin files are gigabytes: hash a megabyte at a time
VERDICTS = ("usable", "excluded")
_RUN_IN_NAME = re.compile(r"_run-([0-9]{2})")   # a data file name carries its run as "_run-01"


def _check_qc(run_indexes: list[int], qc: dict[int, tuple[str, str | None]]) -> None:
    """Raise one ValueError naming every problem: a missing or extra run, a bad verdict, a bad reason.

    An excluded run must have a reason; a usable run must have none.
    """
    errors = []
    errors += [f"qc is missing run {i}" for i in run_indexes if i not in qc]
    errors += [f"qc has run {i}, which the session does not have" for i in sorted(qc)
               if i not in run_indexes]
    for i in sorted(qc):
        verdict, reason = qc[i]
        if verdict not in VERDICTS:
            errors.append(f"run {i}: verdict must be one of {', '.join(VERDICTS)}: {verdict!r}")
        elif verdict == "usable" and reason is not None:
            errors.append(f"run {i}: reason must be None when the verdict is usable: {reason!r}")
        elif verdict == "excluded" and reason is None:
            errors.append(f"run {i} is excluded but has no reason")
    if errors:
        raise ValueError("\n".join(errors))


def _check_run_files(session_dir: Path, run_indexes: list[int]) -> None:
    """Raise one ValueError naming every file under session_dir whose name holds "_run-NN" for a run
    the sidecar does not list. Such a file was recorded without Add run, so no run owns it.
    """
    errors = []
    for folder, _, names in os.walk(session_dir):
        for name in sorted(names):
            match = _RUN_IN_NAME.search(name)
            if match and int(match.group(1)) not in run_indexes:
                relative = (Path(folder) / name).relative_to(session_dir).as_posix()
                errors.append(f"{relative} is for run {int(match.group(1))}, which the session "
                              "does not have; press Add run for it before closing")
    if errors:
        raise ValueError("\n".join(errors))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(CHUNK_BYTES):
            digest.update(chunk)
    return digest.hexdigest()


def _closeout_files(session_dir: Path) -> list[dict]:
    """Path, size and sha256 of every file under session_dir, sorted by path.

    Paths are relative to session_dir with forward slashes. The sidecar itself and
    any *.tmp file are left out.
    """
    entries = []
    for folder, _, names in os.walk(session_dir):
        for name in names:
            path = Path(folder) / name
            relative = path.relative_to(session_dir).as_posix()
            if relative == SIDECAR_NAME or name.endswith(".tmp"):
                continue
            entries.append({"path": relative, "bytes": path.stat().st_size, "sha256": _sha256(path)})
    return sorted(entries, key=lambda entry: entry["path"])


def close_session(session_dir: Path, closed_by: str, qc: dict[int, tuple[str, str | None]],
                  lifecycle: str, now: datetime | None = None) -> dict:
    """Close a session and return the final sidecar dict.

    qc maps run_index -> ("usable" | "excluded", reason or None). Every run in the sidecar must be
    present, and no other. lifecycle is the lifecycle value recorded at close.

    Appends one qc entry per run and one lifecycle entry, lists every file with its size and
    sha256, and writes closeout_files, closed_at and closed_by in one write so that a crash never
    leaves a session marked closed without its hashes. Then deletes current-session.json in the rig
    root, but only when it points at this session.

    Raise RuntimeError when the session is already closed. Raise ValueError for a bad qc, for a
    data file of a run the sidecar does not list, or for a sidecar that fails the registry checks;
    nothing is written then.
    """
    session_dir = Path(session_dir)
    sidecar = read_json(session_dir / SIDECAR_NAME)
    if "closed_at" in sidecar:
        raise RuntimeError(f"session is already closed: {session_dir}")
    run_indexes = [run["run_index"] for run in sidecar["runs"]]
    _check_qc(run_indexes, qc)
    _check_run_files(session_dir, run_indexes)   # before hashing: hashing a large file takes time

    stamp = (now or datetime.now()).astimezone().isoformat(timespec="seconds")   # a naive time is local
    state = sidecar["state"]
    state["qc"] += [{"run_index": i, "verdict": qc[i][0], "reason": qc[i][1]} for i in run_indexes]
    state["lifecycle"].append({"value": lifecycle, "at": stamp, "by": closed_by})
    sidecar["closeout_files"] = _closeout_files(session_dir)
    sidecar["closed_at"] = stamp
    sidecar["closed_by"] = closed_by

    errors = validate_sidecar(sidecar)
    if errors:
        raise ValueError("\n".join(errors))
    write_json(session_dir / SIDECAR_NAME, sidecar)

    root = rig_root()
    current = read_current(root)
    if current is not None and Path(current["session_dir"]).resolve() == session_dir.resolve():
        (root / CURRENT_SESSION_JSON).unlink()
    return sidecar
