"""Build a rig's configuration on a bare disk, or edit it later.

Two layers. apply_setup does all the work and opens no window. run_setup is
the PySide6 dialog that only collects answers and hands them to apply_setup.
"""

from datetime import datetime
from pathlib import Path

from raymondlab.rig.config import (SCHEMA_VERSION, Project, RigConfig, load_projects, load_rig,
                                   save_projects, save_rig, slot_is_locked, validate_project,
                                   validate_rig)


class SlotLocked(Exception):
    """The rig's mint slot cannot change: at least one subject_id has already been minted in it."""


def _now_iso() -> str:
    """The current local time as ISO-8601 with the UTC offset, to the second."""
    return datetime.now().astimezone().isoformat(timespec="seconds")


def apply_setup(root: Path, answers: dict) -> RigConfig:
    """Create or update a rig's configuration from typed answers, and return the saved RigConfig.

    answers = {rig_id, mint_slot, staging_root, default_experimenter, projects: [dict, ...]}.
    projects is the complete list the rig records for; it replaces what projects.json held.

    First run: validate, make root and staging_root, save rig.json and projects.json.
    Re-run: keep created_at; if mint_slot differs and the slot is locked, raise SlotLocked,
    because a slot may not change once an id has been minted in it.

    Raises ValueError with every error joined by newlines. Nothing is written when any
    answer is bad.
    """
    root = Path(root)
    existing = load_rig(root)

    cfg = RigConfig(
        schema_version=SCHEMA_VERSION,
        rig_id=str(answers.get("rig_id", "")),
        mint_slot=str(answers.get("mint_slot", "")),
        staging_root=str(answers.get("staging_root", "")),
        default_experimenter=str(answers.get("default_experimenter", "")),
        created_at=existing.created_at if existing else _now_iso(),
    )

    errors = validate_rig(cfg)
    projects: list[Project] = []
    for raw in answers.get("projects", []):
        p = Project(project_id=str(raw.get("project_id", "")),
                    dataset_dir=str(raw.get("dataset_dir", "")),
                    protocol_id=str(raw.get("protocol_id", "")),
                    lab=str(raw.get("lab", "")),
                    institution=str(raw.get("institution", "")))
        errors += validate_project(p, projects)
        projects.append(p)
    if errors:
        raise ValueError("\n".join(errors))

    if existing and cfg.mint_slot != existing.mint_slot and slot_is_locked(existing):
        raise SlotLocked(f"mint_slot {existing.mint_slot!r} is locked: an id has already been "
                         f"minted in it, so it cannot change to {cfg.mint_slot!r}")

    root.mkdir(parents=True, exist_ok=True)
    Path(cfg.staging_root).mkdir(parents=True, exist_ok=True)
    if cfg != existing:
        save_rig(root, cfg)
    save_projects(root, projects)
    return cfg


def run_setup(parent=None) -> RigConfig | None:
    """Ask a person for the rig's settings, apply them, and return the saved configuration.

    A PySide6 QDialog. Pre-fills from load_rig() when the rig already has a configuration.
    The slot field is disabled when the slot is locked. Returns None when the person cancels.
    Qt is imported here, not at module level, so apply_setup and the tests need no Qt.
    """
    from raymondlab.gui.rig_setup_dialog import RigSetupDialog
    from raymondlab.rig.paths import rig_root

    root = rig_root()
    existing = load_rig(root)
    projects = load_projects(root) if existing and (root / "projects.json").exists() else []
    dialog = RigSetupDialog(root, existing, projects, parent=parent)
    if dialog.exec():
        return dialog.result_config
    return None
