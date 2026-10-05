import json
from importlib.resources import files
from pathlib import Path

import pytest

from raymondlab.core import facet_registry as fr

FIXTURE = Path(__file__).parent / "fixtures" / "calibration_answers.json"
SESSION_SCOPES = {"session", "block", "signal", "device"}


def _facets_file() -> dict:
    return json.loads((files("raymondlab.standard") / "facets.json").read_text(encoding="utf-8"))


def _answers() -> dict:
    with open(FIXTURE, encoding="utf-8") as f:
        return json.load(f)


# --- load_registry ---

def test_load_registry_holds_every_facet_of_the_file():
    registry = fr.load_registry()
    assert sum(len(entries) for entries in registry.values()) == len(_facets_file()["facets"])


def test_every_facet_has_key_scope_and_type():
    for key, entries in fr.load_registry().items():
        for facet in entries:
            assert facet.key == key and key
            assert facet.scope
            assert facet.type


def test_file_names_its_source_version():
    assert _facets_file()["source_version"]


def test_nested_keys_are_registered():
    registry = fr.load_registry()
    assert "surgery[].eye" in registry
    assert "closeout_files[].sha256" in registry


def test_enum_facet_carries_allowed_values_as_tuple():
    (sex,) = fr.load_registry()["sex"]
    assert sex.allowed == ("M", "F", "U", "O")
    assert isinstance(sex.scope, tuple)


def test_load_registry_reads_an_explicit_path(tmp_path):
    path = tmp_path / "facets.json"
    path.write_text(json.dumps({"facets": [
        {"key": "a", "scope": ["x"], "origin": "captured", "type": "string", "status": "confirmed", "note": "n"},
        {"key": "a", "scope": ["y"], "origin": "derived", "type": "enum", "allowed": ["p"], "status": "open"},
    ], "required": {}}), encoding="utf-8")
    entries = fr.load_registry(path)["a"]
    assert [f.scope for f in entries] == [("x",), ("y",)]
    assert entries[0].allowed == ()
    assert entries[1].allowed == ("p",)


# --- required_keys ---

def test_required_keys_for_eye_calibration():
    keys = fr.required_keys("eye-calibration")
    assert "eye_measured" in keys
    assert "runs[].run_index" in keys


def test_required_keys_common_matches_the_file():
    assert fr.required_keys("common") == _facets_file()["required"]["common"]["keys"]


def test_required_keys_unknown_kind_raises_naming_it():
    with pytest.raises(KeyError, match="nope"):
        fr.required_keys("nope")


# --- flatten ---

def test_flatten_nested_dict():
    assert fr.flatten({"a": {"b": 1}}) == [("a.b", 1)]


def test_flatten_list_of_dicts():
    assert fr.flatten({"s": [{"d": 1}]}) == [("s[].d", 1)]


def test_flatten_list_of_dicts_gives_one_pair_per_entry_key():
    assert fr.flatten({"s": [{"d": 1}, {"d": 2}]}) == [("s[].d", 1), ("s[].d", 2)]


def test_flatten_list_of_plain_values_is_one_leaf():
    assert fr.flatten({"x": [1, 2]}) == [("x", [1, 2])]


def test_flatten_empty_list_is_one_leaf_with_brackets():
    assert fr.flatten({"s": []}) == [("s[]", [])]


def test_flatten_nested_empty_list():
    assert fr.flatten({"state": {"qc": []}}) == [("state.qc[]", [])]


# --- validate ---

def test_validate_rejects_unregistered_key_naming_it():
    errors = fr.validate({"weigth": 20}, {"session"})
    assert len(errors) == 1
    assert "unregistered key: weigth" in errors[0]


def test_validate_accepts_registered_key():
    assert fr.validate({"weight": 20}, {"session"}) == []


def test_validate_rejects_enum_value_not_allowed():
    errors = fr.validate({"sex": "male"}, {"subject"})
    assert len(errors) == 1
    assert "sex" in errors[0]


def test_validate_accepts_none_for_enum_but_still_checks_a_given_value():
    usable = {"state": {"qc": [{"run_index": 1, "verdict": "usable", "reason": None}]}}
    assert fr.validate(usable, {"session"}) == []
    bogus = {"state": {"qc": [{"run_index": 1, "verdict": "excluded", "reason": "bogus"}]}}
    errors = fr.validate(bogus, {"session"})
    assert len(errors) == 1
    assert "state.qc[].reason" in errors[0]


def test_validate_accepts_enum_value_allowed():
    assert fr.validate({"sex": "M"}, {"subject"}) == []


def test_validate_enum_compares_str_of_value():
    # a channel number is stored as the integer 1 but registered as the string "1"
    registry = {"ch": [fr.Facet("ch", ("signal",), "captured", "enum", ("1", "2"), "confirmed")]}
    assert fr.validate({"ch": 1}, {"signal"}, registry) == []
    assert len(fr.validate({"ch": 3}, {"signal"}, registry)) == 1


def test_validate_enum_with_no_listed_values_accepts_anything():
    registry = {"r": [fr.Facet("r", ("session",), "captured", "enum", (), "open")]}
    assert fr.validate({"r": "whatever"}, {"session"}, registry) == []


def test_validate_scope_must_match():
    errors = fr.validate({"preamp_gain": 60}, {"session"})
    assert len(errors) == 1
    assert "preamp_gain" in errors[0]
    assert fr.validate({"preamp_gain": 60}, {"signal"}) == []


def test_validate_checks_nested_list_keys():
    good = {"surgery": [{"eye": "left"}]}
    bad = {"surgery": [{"eye": "middle"}]}
    assert fr.validate(good, {"subject"}) == []
    errors = fr.validate(bad, {"subject"})
    assert len(errors) == 1
    assert "surgery[].eye" in errors[0]


def test_validate_accepts_inner_keys_of_a_block_ancestor():
    registry = {"blk": [fr.Facet("blk", ("session",), "captured", "typed block", (), "confirmed")]}
    assert fr.validate({"blk": {"anything": {"deep": 1}}}, {"session"}, registry) == []


def test_validate_accepts_inner_keys_of_a_bool_plus_note_ancestor():
    obj = {"rig_reconfigured": {"value": True, "note": "moved the camera"}}
    assert fr.validate(obj, {"session"}) == []


def test_validate_empty_list_accepted_when_inner_keys_are_in_scope():
    assert fr.validate({"runs": []}, {"session"}) == []
    assert fr.validate({"state": {"qc": []}}, {"session"}) == []
    assert fr.validate({"alleles": [], "interval_log": []}, {"subject"}) == []


def test_validate_empty_list_rejected_when_inner_keys_are_out_of_scope():
    errors = fr.validate({"runs": []}, {"subject"})
    assert len(errors) == 1
    assert "runs[]" in errors[0]


def test_validate_empty_list_of_unknown_name_rejected():
    assert len(fr.validate({"nothing": []}, {"session", "subject"})) == 1


def test_validate_reports_every_problem():
    errors = fr.validate({"weigth": 20, "sex": "male"}, {"subject", "session"})
    assert len(errors) == 2


def test_validate_accepts_sidecar_built_from_fixture():
    answers = _answers()
    mouse = answers["mouse"]
    constants = dict(answers["session"])
    constants.pop("experimenter")
    constants.pop("weight")
    constants.pop("session_description")
    session_part = {
        "schema_type": "eye-calibration",
        "task": "eyecal",
        "rig_id": "D241A",
        "experimenter": [answers["session"]["experimenter"]],
        "weight": answers["session"]["weight"],
        "session_description": answers["session"]["session_description"],
        "state": {"lifecycle": [{"value": "pilot", "at": "2026-10-05T10:00:00-07:00", "by": "bangeles"}],
                  "qc": []},
        "runs": [],
        "light_state": "lit",
        **constants,
    }
    assert fr.validate(mouse, {"subject"}) == []
    assert fr.validate(session_part, SESSION_SCOPES) == []


# --- check_required ---

def test_check_required_names_missing_key():
    errors = fr.check_required({"a": 1}, ["a", "weight"])
    assert errors == ["missing required key: weight"]


def test_check_required_none_counts_as_missing():
    assert fr.check_required({"weight": None}, ["weight"]) == ["missing required key: weight"]


def test_check_required_all_present():
    assert fr.check_required({"a": 1, "b": 0}, ["a", "b"]) == []


def test_check_required_accepts_empty_runs_for_run_index():
    assert fr.check_required({"runs": []}, ["runs[].run_index"]) == []


def test_check_required_runs_must_be_a_list():
    assert fr.check_required({"runs": 3}, ["runs[]"]) == ["missing required key: runs[]"]
    assert fr.check_required({"runs": []}, ["runs[]"]) == []


def test_check_required_every_run_needs_run_index():
    ok = {"runs": [{"run_index": 1}, {"run_index": 2}]}
    bad = {"runs": [{"run_index": 1}, {"other": 2}]}
    assert fr.check_required(ok, ["runs[].run_index"]) == []
    assert fr.check_required(bad, ["runs[].run_index"]) == ["missing required key: runs[].run_index"]


def test_check_required_walks_dotted_paths():
    obj = {"state": {"lifecycle": []}}
    assert fr.check_required(obj, ["state.lifecycle[]"]) == []
    assert fr.check_required({"state": {}}, ["state.lifecycle[]"]) == ["missing required key: state.lifecycle[]"]


def test_check_required_missing_parent_counts_as_missing():
    assert fr.check_required({}, ["state.lifecycle[]"]) == ["missing required key: state.lifecycle[]"]
    assert fr.check_required({}, ["runs[].run_index"]) == ["missing required key: runs[].run_index"]
    assert fr.check_required({"state": 5}, ["state.qc[]"]) == ["missing required key: state.qc[]"]
