import copy
import json
from pathlib import Path

import pytest

from raymondlab.core import facet_registry
from raymondlab.rig.config import Project, RigConfig
from raymondlab.session import model

FIXTURE = Path(__file__).parent / "fixtures" / "calibration_answers.json"
START = "2026-04-10T15:05:00-07:00"

RIG = RigConfig(schema_version="1.0.0", rig_id="D241A", mint_slot="m7", staging_root="/tmp/stage",
                default_experimenter="bangeles", created_at="2026-04-01T09:00:00-07:00")
PROJECT = Project(project_id="okr2026", dataset_dir="okr2026", protocol_id="IACUC-123",
                  lab="Raymond Lab", institution="Stanford University")


@pytest.fixture
def answers() -> dict:
    return json.loads(FIXTURE.read_text(encoding="utf-8"))


@pytest.fixture
def mouse(answers) -> dict:
    return answers["mouse"]


@pytest.fixture
def session(answers) -> dict:
    return {**answers["session"], "session_start_time": START}


def build(mouse, session, subject_id="m7k3q9"):
    return model.build_calibration_sidecar(RIG, PROJECT, mouse, session, subject_id)


# pad_system_id

def test_pad_short_id_to_14_characters():
    padded = model.pad_system_id("123")
    assert padded == "00000000000123"
    assert len(padded) == 14
    assert isinstance(padded, str)


def test_pad_keeps_a_full_id():
    assert model.pad_system_id("20262280001218") == "20262280001218"


@pytest.mark.parametrize("bad", ["", "12a4", " 123", "12.5", "-1", "123456789012345", "١٢٣", 123, None])
def test_pad_rejects_bad_ids(bad):
    with pytest.raises(ValueError):
        model.pad_system_id(bad)


# validate_sidecar

def test_validate_sidecar_accepts_a_built_sidecar(mouse, session):
    assert model.validate_sidecar(build(mouse, session)) == []


def test_validate_sidecar_names_a_misspelled_key(mouse, session):
    sidecar = build(mouse, session)
    sidecar["wieght"] = sidecar.pop("weight")
    errors = model.validate_sidecar(sidecar)
    assert any("wieght" in e for e in errors)


# build_calibration_sidecar: the good case

def test_built_sidecar_has_every_required_key(mouse, session):
    sidecar = build(mouse, session)
    keys = facet_registry.required_keys("common") + facet_registry.required_keys(model.SCHEMA_TYPE)
    assert facet_registry.check_required(sidecar, keys) == []


def test_built_sidecar_validates_in_both_parts(mouse, session):
    sidecar = build(mouse, session)
    subject_part = {k: v for k, v in sidecar.items() if k in model.SUBJECT_KEYS}
    rest = {k: v for k, v in sidecar.items() if k not in model.SUBJECT_KEYS}
    assert facet_registry.validate(subject_part, {"subject"}) == []
    assert facet_registry.validate(rest, model.SIDECAR_SCOPES) == []


def test_fixed_values(mouse, session):
    sidecar = build(mouse, session)
    assert sidecar["schema_type"] == "eye-calibration" == model.SCHEMA_TYPE
    assert sidecar["schema_version"] == model.SCHEMA_VERSION
    assert sidecar["registered_via"] == "rig"
    assert sidecar["task"] == "eyecal"
    assert sidecar["light_state"] == "lit"
    assert sidecar["runs"] == []
    assert sidecar["state"] == {
        "lifecycle": [{"value": "production", "at": START, "by": "bangeles"}],
        "qc": [],
    }


def test_values_copied_from_rig_and_project_and_arguments(mouse, session):
    sidecar = build(mouse, session, subject_id="abc123")
    assert sidecar["subject_id"] == "abc123"
    assert sidecar["rig_id"] == "D241A"
    assert sidecar["project_id"] == "okr2026"
    assert sidecar["protocol_id"] == "IACUC-123"
    assert sidecar["lab"] == "Raymond Lab"
    assert sidecar["institution"] == "Stanford University"
    assert sidecar["session_start_time"] == START


def test_mouse_and_session_values_pass_through(mouse, session):
    sidecar = build(mouse, session)
    assert sidecar["sex"] == "M"
    assert sidecar["strain"] == "C57BL/6J"
    assert sidecar["surgery"] == mouse["surgery"]
    assert sidecar["weight"] == 24.5
    assert sidecar["stimulus_frequency_hz"] == [0.5, 1.0]
    assert sidecar["ir_led_wavelength_nm"] == 875


def test_system_id_is_padded_string(mouse, session):
    mouse["system_id"] = "1218"
    assert build(mouse, session)["system_id"] == "00000000001218"


def test_misspelled_mouse_key_raises(mouse, session):
    mouse["litter_idd"] = "L-77"
    with pytest.raises(ValueError, match="litter_idd"):
        build(mouse, session)


def test_key_given_in_both_mouse_and_session_raises(mouse, session):
    session["sex"] = "F"
    with pytest.raises(ValueError, match="key sex is given in both mouse and session"):
        build(mouse, session)


def test_missing_system_id_and_weight_are_reported_together(mouse, session):
    del mouse["system_id"]
    del session["weight"]
    with pytest.raises(ValueError) as err:
        build(mouse, session)
    assert "system_id" in str(err.value)
    assert "weight" in str(err.value)


def test_inputs_are_not_changed(mouse, session):
    mouse_before, session_before = copy.deepcopy(mouse), copy.deepcopy(session)
    build(mouse, session)
    assert mouse == mouse_before
    assert session == session_before


# eye_measured

def test_eye_measured_equals_implant_eye(mouse, session):
    mouse["surgery"][0]["eye"] = "right"
    assert build(mouse, session)["eye_measured"] == "right"


def test_eye_measured_comes_from_the_first_surgery_entry(mouse, session):
    mouse["surgery"] = [
        {"date": "2026-03-20", "performed_by": "bangeles", "procedure": "sensor-implant",
         "description": "x", "eye": "right"},
        {"date": "2026-04-01", "performed_by": "bangeles", "procedure": "sensor-revision",
         "description": "y", "eye": "left"},
    ]
    assert build(mouse, session)["eye_measured"] == "right"


def test_surgery_entry_without_an_eye_fails_the_required_key_check(mouse, session):
    # the registry requires an eye on every surgery entry, so a headpost-only entry is refused
    mouse["surgery"].append({"date": "2026-03-01", "performed_by": "bangeles", "procedure": "headpost",
                             "description": "x"})
    with pytest.raises(ValueError, match=r"surgery\[\]\.eye"):
        build(mouse, session)


@pytest.mark.parametrize("surgery", [[], [{"date": "2026-03-01", "procedure": "headpost"}]])
def test_no_surgery_eye_raises(mouse, session, surgery):
    mouse["surgery"] = surgery
    with pytest.raises(ValueError) as err:
        build(mouse, session)
    assert "eye_measured" in str(err.value)
    assert "surgery" in str(err.value)


def test_missing_surgery_key_raises(mouse, session):
    del mouse["surgery"]
    with pytest.raises(ValueError, match="eye_measured"):
        build(mouse, session)


# system_id

def test_missing_system_id_raises(mouse, session):
    del mouse["system_id"]
    with pytest.raises(ValueError, match="system_id"):
        build(mouse, session)


@pytest.mark.parametrize("bad", ["ab12", "", "123456789012345"])
def test_bad_system_id_raises(mouse, session, bad):
    mouse["system_id"] = bad
    with pytest.raises(ValueError, match="system_id"):
        build(mouse, session)


# experimenter

def test_experimenter_string_becomes_list(mouse, session):
    assert build(mouse, session)["experimenter"] == ["bangeles"]


def test_experimenter_list_is_kept_and_first_one_is_lifecycle_author(mouse, session):
    session["experimenter"] = ["jsmith", "bangeles"]
    sidecar = build(mouse, session)
    assert sidecar["experimenter"] == ["jsmith", "bangeles"]
    assert sidecar["state"]["lifecycle"][0]["by"] == "jsmith"


@pytest.mark.parametrize("bad", ["", [], [""], None])
def test_empty_experimenter_raises(mouse, session, bad):
    session["experimenter"] = bad
    with pytest.raises(ValueError, match="experimenter"):
        build(mouse, session)


def test_missing_experimenter_raises(mouse, session):
    del session["experimenter"]
    with pytest.raises(ValueError, match="experimenter"):
        build(mouse, session)


# validation errors

def test_missing_weight_raises_with_weight_in_message(mouse, session):
    del session["weight"]
    with pytest.raises(ValueError, match="weight"):
        build(mouse, session)


def test_missing_session_start_time_raises(mouse, session):
    del session["session_start_time"]
    with pytest.raises(ValueError, match="session_start_time"):
        build(mouse, session)


def test_unregistered_session_key_raises(mouse, session):
    session["wieght"] = 24.5
    with pytest.raises(ValueError, match="wieght"):
        build(mouse, session)


def test_bad_enum_value_raises(mouse, session):
    mouse["sex"] = "male"
    with pytest.raises(ValueError, match="sex"):
        build(mouse, session)


def test_errors_from_all_checks_are_joined_with_newlines(mouse, session):
    del session["weight"]          # caught by the required-key check
    mouse["sex"] = "male"          # caught by the subject check
    session["wieght"] = 24.5       # caught by the session check
    with pytest.raises(ValueError) as err:
        build(mouse, session)
    lines = str(err.value).split("\n")
    assert len(lines) == 3
    assert any("weight" in line and "missing" in line for line in lines)
    assert any("sex" in line for line in lines)
    assert any("wieght" in line for line in lines)


# rig_reconfigured

def test_rig_reconfigured_is_copied_as_given(mouse, session):
    session["rig_reconfigured"] = {"value": True, "note": "moved camera"}
    assert build(mouse, session)["rig_reconfigured"] == {"value": True, "note": "moved camera"}


def test_rig_reconfigured_is_absent_when_not_given(mouse, session):
    assert "rig_reconfigured" not in build(mouse, session)
