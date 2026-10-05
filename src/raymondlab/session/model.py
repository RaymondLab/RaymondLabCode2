"""The shape of the calibration sidecar (_metadata.json).

Builds the sidecar from the rig settings, the project, the facts about the
mouse and the answers of the session form. The list of keys the sidecar must
carry comes from facets.json through facet_registry.required_keys, so this
module holds no list of its own. A typo in a key name or a missing key fails
here, before anything reaches disk.
"""

import re

from raymondlab.core import facet_registry
from raymondlab.rig.config import Project, RigConfig

SCHEMA_TYPE = "eye-calibration"
SCHEMA_VERSION = "1.0.0"   # bump when the sidecar shape changes

# Which keys of the sidecar describe the mouse, and are therefore registered at "subject" scope.
# The session form uses this list to choose which facts of an existing mouse to copy into its
# fields. The model does not use it: it does not drop or check keys by this list.
SUBJECT_KEYS = ["subject_id", "system_id", "colony_mouse_id", "species", "strain", "sex",
                "date_of_birth", "genotype", "alleles", "litter_id", "colony_origin", "surgery",
                "interval_log", "schema_version"]
SIDECAR_SCOPES = {"session", "block", "signal", "device"}   # where the rest of the keys are registered

_SYSTEM_ID = re.compile(r"[0-9]{1,14}")   # ASCII digits only: str.isdigit() also accepts other scripts
_SYSTEM_ID_LENGTH = 14


def pad_system_id(s: str) -> str:
    """Normalise a Transnetyx system id: digits only, at most 14 of them, left-padded with zeros to 14.

    Raise ValueError for anything else. Keep the result a string; a number would drop the leading zeros.
    """
    if not isinstance(s, str) or not _SYSTEM_ID.fullmatch(s):
        raise ValueError(f"system_id must be 1 to {_SYSTEM_ID_LENGTH} digits: {s!r}")
    return s.zfill(_SYSTEM_ID_LENGTH)


def _experimenters(value: object) -> list[str]:
    """The experimenter answer as a list of names. Raise ValueError when there is no name."""
    names = [value] if isinstance(value, str) else value
    if not isinstance(names, list) or not names or not all(isinstance(n, str) and n for n in names):
        raise ValueError(f"experimenter must be a name or a list of names, none empty: {value!r}")
    return list(names)


def _eye_measured(surgery: object) -> str:
    """The eye of the first surgery entry that has one. Raise ValueError when no entry has an eye."""
    if isinstance(surgery, list):
        for entry in surgery:
            if isinstance(entry, dict) and entry.get("eye"):
                return entry["eye"]
    raise ValueError("eye_measured: no surgery entry has an eye")


def _split_errors(call, *args) -> tuple[object, list[str]]:
    """Run call(*args). Return (result, []) or (None, [message]) when it raises ValueError."""
    try:
        return call(*args), []
    except ValueError as err:
        return None, [str(err)]


def _subject_key_names(registry: dict) -> set[str]:
    """The top-level key names that the registry lists at "subject" scope."""
    names = set()
    for key, entries in registry.items():
        if any("subject" in f.scope for f in entries):
            names.add(re.split(r"[.\[]", key, maxsplit=1)[0])
    return names


def validate_sidecar(sidecar: dict) -> list[str]:
    """Check a sidecar that was read back from disk. Return the errors; [] means valid.

    The top-level keys that the registry lists at "subject" scope are checked at "subject" scope,
    the other keys at the sidecar scopes. The keys that every calibration sidecar must carry are
    checked too.
    """
    registry = facet_registry.load_registry()
    subject_names = _subject_key_names(registry)
    subject_part = {k: v for k, v in sidecar.items() if k in subject_names}
    rest = {k: v for k, v in sidecar.items() if k not in subject_names}
    required = facet_registry.required_keys("common") + facet_registry.required_keys(SCHEMA_TYPE)
    return (facet_registry.check_required(sidecar, required)
            + facet_registry.validate(subject_part, {"subject"}, registry)
            + facet_registry.validate(rest, SIDECAR_SCOPES, registry))


def build_calibration_sidecar(rig: RigConfig, project: Project, mouse: dict, session: dict,
                              subject_id: str) -> dict:
    """Assemble the _metadata.json dict for a new calibration session, validate it, and return it.

    mouse holds the typed facts about the animal; session holds the calibration constants,
    session_description, experimenter and session_start_time. subject_id is the id the store minted.
    Every key of mouse and of session is copied into the sidecar, so the registry decides what is
    allowed. A key given in both mouse and session is an error.

    Fixed values: schema_type, schema_version, registered_via = "rig", task = "eyecal",
    light_state = "lit", runs = [] and
    state = {"lifecycle": [{"value": "production", "at": session_start_time, "by": <first experimenter>}],
             "qc": []}.

    system_id is padded to 14 digits. experimenter is stored as a list. eye_measured is copied from
    the mouse's surgery record, never typed, so the form cannot disagree with it.

    Raise one ValueError that holds every problem, one per line: a bad system_id, experimenter or
    eye_measured, a key given twice, then the registry checks in this order: required keys, the mouse
    keys (plus subject_id and schema_version) at "subject" scope, the rest at the sidecar scopes.
    """
    errors = []
    system_id, id_errors = _split_errors(pad_system_id, mouse.get("system_id"))
    experimenter, who_errors = _split_errors(_experimenters, session.get("experimenter"))
    eye, eye_errors = _split_errors(_eye_measured, mouse.get("surgery"))
    errors += id_errors + who_errors + eye_errors
    errors += [f"key {key} is given in both mouse and session" for key in mouse if key in session]

    sidecar = {**session, **mouse}
    sidecar.update({
        "subject_id": subject_id,
        "schema_type": SCHEMA_TYPE,
        "schema_version": SCHEMA_VERSION,
        "registered_via": "rig",
        "task": "eyecal",
        "light_state": "lit",
        "rig_id": rig.rig_id,
        "project_id": project.project_id,
        "protocol_id": project.protocol_id,
        "lab": project.lab,
        "institution": project.institution,
        "runs": [],
        "state": {
            "lifecycle": [{"value": "production", "at": session.get("session_start_time"),
                           "by": experimenter[0] if experimenter else None}],
            "qc": [],
        },
    })
    # a value that could not be made is left out (or kept as typed, for system_id) so that the
    # registry checks report it too
    if system_id is not None:
        sidecar["system_id"] = system_id
    if experimenter is not None:
        sidecar["experimenter"] = experimenter
    else:
        sidecar.pop("experimenter", None)
    if eye is not None:
        sidecar["eye_measured"] = eye
    else:
        sidecar.pop("eye_measured", None)

    required = facet_registry.required_keys("common") + facet_registry.required_keys(SCHEMA_TYPE)
    subject_names = set(mouse) | {"subject_id", "schema_version"}
    subject_part = {k: v for k, v in sidecar.items() if k in subject_names}
    rest = {k: v for k, v in sidecar.items() if k not in subject_names}
    errors += (facet_registry.check_required(sidecar, required)
               + facet_registry.validate(subject_part, {"subject"})
               + facet_registry.validate(rest, SIDECAR_SCOPES))
    if errors:
        raise ValueError("\n".join(errors))
    return sidecar
