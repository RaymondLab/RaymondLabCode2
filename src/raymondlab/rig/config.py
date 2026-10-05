"""What a rig and a project are: the shape of rig.json and projects.json.

Reads and writes both files, checks typed answers before they reach disk, and
answers one question the setup dialog needs: may this rig's mint slot still
change?
"""

import re
from dataclasses import asdict, dataclass, fields
from pathlib import Path, PurePosixPath, PureWindowsPath

from raymondlab.core import alphabet
from raymondlab.core.jsonio import read_json, write_json
from raymondlab.rig.paths import PROJECTS_JSON, RIG_JSON

SCHEMA_VERSION = "1.0.0"

MAX_DIRNAME = 40        # a directory name segment is at most 40 characters
MAX_STAGING_ROOT = 80   # so that the full path of a session file stays under the archive tools' limit

_PROJECT_ID = re.compile(r"^[a-z0-9-]+$")


@dataclass(frozen=True)
class Project:
    """One project this rig records for: one entry of projects.json."""
    project_id: str        # slug: lowercase letters, digits and dashes; unique in projects.json
    dataset_dir: str       # directory name, at most 40 characters, no spaces
    protocol_id: str       # IACUC id, free text, non-empty
    lab: str
    institution: str


@dataclass(frozen=True)
class RigConfig:
    """The settings of one rig: the contents of rig.json."""
    schema_version: str    # "1.0.0"
    rig_id: str            # e.g. "D241A"; provisional, the final rig naming scheme is not decided
    mint_slot: str         # the 2 characters this rig owns in every id it mints
    staging_root: str      # absolute path, at most 80 characters
    default_experimenter: str
    created_at: str        # ISO-8601 with UTC offset


def _from_dict(cls, data: dict):
    """Build a dataclass from a dict, taking only the fields the class declares."""
    names = {f.name for f in fields(cls)}
    return cls(**{k: v for k, v in data.items() if k in names})


def load_rig(root: Path) -> RigConfig | None:
    """Read root/rig.json as a RigConfig, or None when the file is not there."""
    path = Path(root) / RIG_JSON
    if not path.exists():
        return None
    return _from_dict(RigConfig, read_json(path))


def save_rig(root: Path, cfg: RigConfig) -> None:
    """Write cfg to root/rig.json with write_json, so the bytes are stable and the write is atomic."""
    write_json(Path(root) / RIG_JSON, asdict(cfg))


def load_projects(root: Path) -> list[Project]:
    """Read root/projects.json and return its projects."""
    data = read_json(Path(root) / PROJECTS_JSON)
    return [_from_dict(Project, p) for p in data["projects"]]


def save_projects(root: Path, projects: list[Project]) -> None:
    """Write the projects to root/projects.json with write_json."""
    write_json(Path(root) / PROJECTS_JSON,
               {"schema_version": SCHEMA_VERSION, "projects": [asdict(p) for p in projects]})


def slot_is_locked(cfg: RigConfig) -> bool:
    """True when this rig has already minted an id, so its slot may no longer change.

    True if any directory matches <staging_root>/*/raw/sub-m<slot>*.
    """
    staging = Path(cfg.staging_root)
    return any(p.is_dir() for p in staging.glob(f"*/raw/sub-m{cfg.mint_slot}*"))


def validate_rig(cfg: RigConfig) -> list[str]:
    """Check a rig configuration. Return one message per problem; [] means good."""
    errors = []
    if not cfg.rig_id.strip():
        errors.append("rig_id must not be empty")
    if not alphabet.is_slot(cfg.mint_slot):
        errors.append(f"mint_slot {cfg.mint_slot!r} must be exactly 2 characters from the id alphabet "
                      f"({alphabet.ALPHABET})")
    if not cfg.staging_root.strip():
        errors.append("staging_root must not be empty")
    else:
        if len(cfg.staging_root) > MAX_STAGING_ROOT:
            errors.append(f"staging_root is {len(cfg.staging_root)} characters; "
                          f"at most {MAX_STAGING_ROOT} are allowed")
        if not _is_absolute(cfg.staging_root):
            errors.append(f"staging_root {cfg.staging_root!r} must be an absolute path, "
                          f"for example C:\\RaymondLab\\staging")
    return errors


def _is_absolute(path: str) -> bool:
    """True for an absolute path in either Windows or POSIX form, on any platform."""
    return PureWindowsPath(path).is_absolute() or PurePosixPath(path).is_absolute()


def validate_project(p: Project, existing: list[Project]) -> list[str]:
    """Check one project against the projects already saved. Return one message per problem."""
    errors = []
    if not _PROJECT_ID.match(p.project_id):
        errors.append(f"project_id {p.project_id!r} must use only lowercase letters, digits and dashes")
    if any(e.project_id == p.project_id for e in existing):
        errors.append(f"project_id {p.project_id!r} already exists")
    if not p.dataset_dir:
        errors.append("dataset_dir must not be empty")
    else:
        if len(p.dataset_dir) > MAX_DIRNAME:
            errors.append(f"dataset_dir {p.dataset_dir!r} is {len(p.dataset_dir)} characters; "
                          f"at most {MAX_DIRNAME} are allowed")
        if " " in p.dataset_dir:
            errors.append(f"dataset_dir {p.dataset_dir!r} must not contain spaces")
        if "/" in p.dataset_dir or "\\" in p.dataset_dir:
            errors.append(f"dataset_dir {p.dataset_dir!r} must be one folder name, not a path")
        if p.dataset_dir in (".", ".."):
            errors.append(f"dataset_dir {p.dataset_dir!r} is not a folder name")
    if not p.protocol_id.strip():
        errors.append("protocol_id must not be empty")
    if not p.lab.strip():
        errors.append("lab must not be empty")
    if not p.institution.strip():
        errors.append("institution must not be empty")
    return errors
