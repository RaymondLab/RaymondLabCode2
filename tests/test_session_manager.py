"""Headless tests of the session form and the session-manager entry point.

Qt runs on the "offscreen" platform, so no window opens. No test calls exec().
"""

import json
from datetime import datetime
from pathlib import Path

import pytest

from raymondlab.core.jsonio import read_json
from raymondlab.session import store

FIXTURE = Path(__file__).parent / "fixtures" / "calibration_answers.json"


@pytest.fixture(scope="module")
def qapp():
    """One QApplication for this module, on the offscreen platform."""
    with pytest.MonkeyPatch.context() as mp:
        mp.setenv("QT_QPA_PLATFORM", "offscreen")   # must be set before Qt starts
        from PySide6.QtWidgets import QApplication
        yield QApplication.instance() or QApplication([])


@pytest.fixture
def fixture_answers() -> dict:
    return json.loads(FIXTURE.read_text(encoding="utf-8"))


@pytest.fixture
def form(qapp, seeded_rig):
    from raymondlab.gui.session_form import SessionForm
    return SessionForm(seeded_rig)


@pytest.fixture
def messages(monkeypatch) -> list[tuple[str, str]]:
    """Record every QMessageBox.critical and .information as (kind, text) instead of showing it.

    Every QMessageBox.question is answered Yes without being shown (a test may patch it again).
    """
    from PySide6.QtWidgets import QMessageBox
    shown = []
    monkeypatch.setattr(QMessageBox, "question",
                        lambda parent, title, text, *a: QMessageBox.StandardButton.Yes)
    monkeypatch.setattr(QMessageBox, "critical",
                        lambda parent, title, text, *a: shown.append(("critical", text)))
    monkeypatch.setattr(QMessageBox, "information",
                        lambda parent, title, text, *a: shown.append(("information", text)))
    return shown


def fill(form, answers: dict) -> None:
    """Type the fixture's values into the form's widgets, as a person would."""
    mouse, session = answers["mouse"], answers["session"]
    form.project.setCurrentText(answers["project_id"])
    form.system_id.setText(mouse["system_id"])
    for key, field in form.mouse_fields.items():
        field.setText(mouse[key])
    form.sex.setCurrentText({"M": "Male", "F": "Female", "U": "Unknown", "O": "Other"}[mouse["sex"]])
    implant = mouse["surgery"][0]
    for key, field in form.surgery_fields.items():
        field.setText(implant[key])
    form.eye.setCurrentText(implant["eye"])
    form.eye_tracking_method.setCurrentText(session["eye_tracking_method"])
    for key, field in form.calibration_fields.items():
        value = session[key]
        field.setText(", ".join(str(v) for v in value) if isinstance(value, list) else str(value))
    for key, field in form.session_fields.items():
        field.setText(str(session[key]))


def test_answers_have_the_fixture_shape(form, fixture_answers):
    fill(form, fixture_answers)
    answers = form.answers()
    assert set(answers) == set(fixture_answers)
    assert set(answers["mouse"]) == set(fixture_answers["mouse"])
    assert set(answers["session"]) == set(fixture_answers["session"])
    assert answers["mouse"]["sex"] == "M"
    assert answers["session"]["stimulus_frequency_hz"] == [0.5, 1.0]
    assert all(isinstance(f, float) for f in answers["session"]["stimulus_frequency_hz"])
    assert "rig_reconfigured" not in answers["session"]
    assert answers == fixture_answers


def test_rig_reconfigured_is_a_value_and_a_note_when_checked(form, fixture_answers):
    fill(form, fixture_answers)
    form.rig_reconfigured.setChecked(True)
    form.reconfigured_note.setText("moved camera")
    assert form.answers()["session"]["rig_reconfigured"] == {"value": True, "note": "moved camera"}


def test_create_with_empty_system_id_shows_the_error_and_creates_nothing(form, seeded_rig, messages):
    form.create_button.click()
    assert len(messages) == 1
    kind, text = messages[0]
    assert kind == "critical"
    assert "system_id" in text
    assert list((seeded_rig / "staging").glob("*/raw/sub-*/ses-*")) == []
    assert store.read_current(seeded_rig) is None
    assert form.stack.currentWidget() is not form.status_page


def test_calibration_defaults_without_an_earlier_session(form):
    assert form.calibration_fields["stimulus_frequency_hz"].text() == "0.5, 1.0"
    assert form.calibration_fields["duration_s"].text() == "120"
    assert form.calibration_fields["camera_model"].text() == ""
    assert form.session_fields["experimenter"].text() == "bangeles"


def test_calibration_is_prefilled_from_the_last_session(qapp, seeded_rig, fixture_answers):
    from raymondlab.gui.session_form import SessionForm
    fixture_answers["session"]["preamp_gain"] = 80
    store.create_session(seeded_rig, fixture_answers)
    form = SessionForm(seeded_rig)
    assert form.calibration_fields["preamp_gain"].text() == "80"
    assert form.calibration_fields["camera_model"].text() == "OV2311"
    assert form.eye_tracking_method.currentText() == "magnetic-sensor"


def test_lookup_of_a_known_mouse_copies_its_facts(qapp, seeded_rig, fixture_answers, messages):
    from raymondlab.gui.session_form import SessionForm
    fixture_answers["mouse"]["system_id"] = "1218"
    store.create_session(seeded_rig, fixture_answers)
    form = SessionForm(seeded_rig)
    form.system_id.setText("1218")
    form.look_up_button.click()
    assert "00000000001218" in form.padded_label.text()
    assert "RL-0421" in form.found_label.text()
    assert "sex: Male" in form.found_label.text()
    assert not form.yes_button.isHidden()
    form.yes_button.click()
    answers = form.answers()
    assert answers["subject_id"] == "m7k222"
    assert answers["mouse"]["strain"] == "C57BL/6J"
    assert answers["mouse"]["sex"] == "M"
    assert answers["mouse"]["surgery"] == fixture_answers["mouse"]["surgery"]
    assert messages == []


def test_lookup_of_an_unknown_mouse_keeps_the_typed_facts(form, messages):
    form.system_id.setText("42")
    form.look_up_button.click()
    assert "00000000000042" in form.padded_label.text()
    assert form.yes_button.isHidden()
    assert "subject_id" not in form.answers()


def test_create_add_run_and_close(form, seeded_rig, fixture_answers, messages, monkeypatch):
    from raymondlab.gui.session_form import QcDialog
    fill(form, fixture_answers)
    form.create_button.click()
    assert messages == []
    assert form.stack.currentWidget() is form.status_page
    current = store.read_current(seeded_rig)
    assert current["entity_stem"] in form.status_label.text()

    form.add_run_button.click()
    assert store.read_current(seeded_rig)["next_run"] == 2

    form.session_fields["experimenter"].setText("someone-else")   # closed_by is fixed at Create
    monkeypatch.setattr(QcDialog, "exec", lambda self: True)   # accept with the default verdicts
    form.close_button.click()
    sidecar = read_json(Path(current["session_dir"]) / "_metadata.json")
    assert sidecar["closed_by"] == "bangeles"
    assert sidecar["state"]["qc"] == [{"run_index": 1, "verdict": "usable", "reason": None}]
    assert store.read_current(seeded_rig) is None
    assert [kind for kind, _ in messages] == ["information"]
    assert form.result() == form.DialogCode.Accepted


def test_a_number_field_with_text_names_the_key_and_creates_nothing(form, seeded_rig,
                                                                    fixture_answers, messages):
    fill(form, fixture_answers)
    form.session_fields["weight"].setText("abc")
    form.create_button.click()
    assert messages == [("critical", "weight is not a number: 'abc'")]
    assert list((seeded_rig / "staging").glob("*/raw/sub-*/ses-*")) == []
    assert form.stack.currentWidget() is not form.status_page


def test_padded_label_follows_every_edit_and_matches_what_is_saved(form, fixture_answers, messages):
    from PySide6.QtTest import QTest
    fill(form, fixture_answers)
    form.system_id.setText("1218")
    form.look_up_button.click()
    assert form.padded_label.text() == "system_id: 00000000001218"
    QTest.keyClicks(form.system_id, "9")   # typed by a person, after the look-up
    assert form.padded_label.text() == "system_id: 00000000012189"
    form.create_button.click()
    assert messages == []
    sidecar = read_json(form.session_dir / "_metadata.json")
    assert form.padded_label.text() == f"system_id: {sidecar['system_id']}"


def test_form_opens_on_the_status_page_when_a_session_is_open(qapp, seeded_rig, fixture_answers,
                                                             messages, monkeypatch):
    from raymondlab.gui.session_form import QcDialog, SessionForm
    fixture_answers["session"]["experimenter"] = "other-person"   # not the rig's default
    directory = store.create_session(seeded_rig, fixture_answers)
    form = SessionForm(seeded_rig)
    assert form.stack.currentWidget() is form.status_page
    assert form.session_dir.resolve() == directory.resolve()
    assert store.read_current(seeded_rig)["entity_stem"] in form.status_label.text()
    assert form.session_fields["experimenter"].text() == "other-person"

    monkeypatch.setattr(QcDialog, "exec", lambda self: True)
    form.close_button.click()
    assert read_json(directory / "_metadata.json")["closed_by"] == "other-person"


def test_leaving_with_an_open_session_asks_first(form, fixture_answers, messages, monkeypatch):
    from PySide6.QtWidgets import QMessageBox
    fill(form, fixture_answers)
    form.create_button.click()
    asked, rejected = [], []
    form.rejected.connect(lambda: rejected.append(True))
    answer = QMessageBox.StandardButton.No
    monkeypatch.setattr(QMessageBox, "question", lambda parent, title, text, *a: asked.append(text) or answer)
    form.reject()
    assert asked == ["A session is still open. Leave without closing it?"]
    assert rejected == []
    answer = QMessageBox.StandardButton.Yes
    form.reject()
    assert rejected == [True]


def test_qc_dialog_reasons_come_from_the_registry(qapp):
    from raymondlab.core.facet_registry import load_registry
    from raymondlab.gui.session_form import QcDialog
    dialog = QcDialog([1, 2])
    reasons = [dialog.reasons[2].itemText(i) for i in range(dialog.reasons[2].count())]
    assert reasons[1:] == list(load_registry()["state.qc[].reason"][0].allowed)
    dialog.verdicts[2].setCurrentText("excluded")
    dialog.reasons[2].setCurrentText("hardware-failure")
    assert dialog.qc() == {1: ("usable", None), 2: ("excluded", "hardware-failure")}


def test_main_returns_1_when_rig_setup_is_cancelled(qapp, rig_root, monkeypatch):
    from raymondlab.apps import session_manager
    from raymondlab.rig import setup
    monkeypatch.setattr(setup, "run_setup", lambda parent=None: None)
    assert session_manager.main([]) == 1


@pytest.mark.parametrize("accepted, code", [(True, 0), (False, 1)])
def test_main_returns_the_form_result(qapp, seeded_rig, monkeypatch, accepted, code):
    from raymondlab.apps import session_manager
    from raymondlab.gui.session_form import SessionForm
    monkeypatch.setattr(SessionForm, "exec", lambda self: accepted)
    assert session_manager.main([]) == code


def test_a_new_mouse_system_id_is_confirmed_before_create(form, seeded_rig, fixture_answers,
                                                          messages, monkeypatch):
    from PySide6.QtWidgets import QMessageBox
    fill(form, fixture_answers)
    form.system_id.setText("1218")
    asked = []
    answer = QMessageBox.StandardButton.No
    monkeypatch.setattr(QMessageBox, "question", lambda parent, title, text, *a: asked.append(text) or answer)
    form.create_button.click()
    assert asked == ["Save system_id 00000000001218 for a new mouse?"]
    assert list((seeded_rig / "staging").glob("*/raw/sub-*/ses-*")) == []
    assert store.read_current(seeded_rig) is None
    assert form.stack.currentWidget() is not form.status_page

    answer = QMessageBox.StandardButton.Yes
    form.create_button.click()
    assert len(asked) == 2
    assert messages == []
    assert form.stack.currentWidget() is form.status_page
    assert read_json(form.session_dir / "_metadata.json")["system_id"] == "00000000001218"


def test_a_confirmed_existing_mouse_is_not_asked_again(qapp, seeded_rig, fixture_answers, messages,
                                                       monkeypatch):
    from PySide6.QtWidgets import QMessageBox
    from raymondlab.gui.session_form import SessionForm
    from raymondlab.session.closeout import close_session
    first = store.create_session(seeded_rig, fixture_answers, now=datetime(2026, 4, 10, 15, 5))
    close_session(first, "bangeles", {}, "production")
    form = SessionForm(seeded_rig)
    fill(form, fixture_answers)
    form.look_up_button.click()
    form.yes_button.click()
    asked = []
    monkeypatch.setattr(QMessageBox, "question", lambda parent, title, text, *a: asked.append(text))
    form.create_button.click()
    assert asked == []
    assert messages == []
    assert form.session_dir.parent == first.parent


def test_status_page_names_the_next_recording(form, fixture_answers, messages):
    fill(form, fixture_answers)
    form.create_button.click()
    assert "The next recording is run 01. Press Add run after it finishes." in form.status_label.text()
    form.add_run_button.click()
    assert "The next recording is run 02. Press Add run after it finishes." in form.status_label.text()
