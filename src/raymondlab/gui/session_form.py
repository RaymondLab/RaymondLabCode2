"""The session form: collects a calibration session's answers and hands them to the store.

No checks live here. The store, the model and the closeout do every check, and
every error they raise is shown to the person with its full text, so they can
fix the answer and try again. The form only turns typed text into numbers and
lists (and says which key holds text that is not a number), and words into the
letters the standard stores.
"""

from datetime import datetime
from pathlib import Path

from PySide6.QtWidgets import (QCheckBox, QComboBox, QDialog, QDialogButtonBox, QFormLayout,
                               QGroupBox, QHBoxLayout, QLabel, QLineEdit, QMessageBox, QPushButton,
                               QStackedWidget, QVBoxLayout, QWidget)

from raymondlab.core.facet_registry import load_registry
from raymondlab.core.jsonio import read_json
from raymondlab.registry.mint import SlotFull
from raymondlab.rig.config import load_projects, load_rig
from raymondlab.session.closeout import close_session
from raymondlab.session.model import SUBJECT_KEYS, pad_system_id
from raymondlab.session.store import (SIDECAR_NAME, add_run, create_session, find_mouse,
                                      last_sidecar, open_session, read_current)

TITLE = "Session manager"
# every error the store and closeout raise; OSError covers a file that cannot be read while hashing
ERRORS = (ValueError, RuntimeError, FileExistsError, SlotFull, OSError)
READ_ERRORS = ERRORS + (KeyError,)   # reading sidecars: a hand-edited file may lack a key
LEAVE_QUESTION = "A session is still open. Leave without closing it?"
LIFECYCLE_AT_CLOSE = "production"

SEX_WORDS = [("Male", "M"), ("Female", "F"), ("Unknown", "U"), ("Other", "O")]   # shown word, stored letter
SEX_WORD_OF = {letter: word for word, letter in SEX_WORDS}
MOUSE_TEXT_KEYS = ["colony_mouse_id", "species", "strain", "date_of_birth", "genotype", "litter_id",
                   "colony_origin"]
SURGERY_TEXT_KEYS = ["date", "performed_by", "procedure", "description"]
DEFAULT_PROCEDURE = "sensor-implant"
SESSION_KEYS = ["session_description", "experimenter", "weight"]

# The rig's calibration constants, and their values on a rig that has no earlier session.
# eye_tracking_method is a separate list; it has no default.
CALIBRATION_DEFAULTS = {
    "stimulus_frequency_hz": [0.5, 1.0],
    "stimulus_amplitude_deg_per_s": 10,
    "duration_s": 120,
    "align_led_pulse_width_ms": 30,
    "align_led_pulse_rate_hz": 1,
    "preamp_gain": 60,
    "lowpass_filter_hz": 100,
    "camera_model": "",
    "camera_separation_deg": 40,
    "camera_equidistance_cm": 5,
    "video_frame_rate_hz": 30,
    "ir_led_wavelength_nm": 875,
}
TEXT_KEYS = {"camera_model"}               # kept as typed text, not made a number
LIST_KEYS = {"stimulus_frequency_hz"}      # typed as numbers split by commas

PROJECT, MOUSE, CALIBRATION, SESSION, STATUS = range(5)
PAGE_TITLES = ["Project", "Mouse", "Calibration", "Session", "Session open"]


def _allowed(key: str) -> list[str]:
    """The allowed values of an enum key, read from the registry."""
    return next((list(f.allowed) for f in load_registry().get(key, []) if f.allowed), [])


def _text(field: QLineEdit) -> str | None:
    """The typed text, or None when the field is empty, so a missing answer is reported as missing."""
    return field.text().strip() or None


def _number(key: str, text: str) -> int | float | None:
    """Typed text as an int or a float. Empty text is None. Raise ValueError naming key otherwise."""
    text = text.strip()
    if not text:
        return None
    try:
        return int(text)
    except ValueError:
        pass
    try:
        return float(text)
    except ValueError:
        raise ValueError(f"{key} is not a number: {text!r}") from None


def _numbers(key: str, text: str) -> list[float] | None:
    """Comma-separated numbers as a list of floats. Empty text is None. Raise ValueError naming key
    when a part is not a number."""
    parts = [part.strip() for part in text.split(",") if part.strip()]
    if not parts:
        return None
    values = []
    for part in parts:
        try:
            values.append(float(part))
        except ValueError:
            raise ValueError(f"{key} is not a number: {part!r}") from None
    return values


def _shown(value: object) -> str:
    """A stored value as text for a field: None is empty, a list is joined with commas."""
    if value is None:
        return ""
    if isinstance(value, list):
        return ", ".join(str(v) for v in value)
    return str(value)


def _combo(items: list[tuple[str, str]]) -> QComboBox:
    """A list of (shown text, stored value) with an empty first entry that stores None."""
    combo = QComboBox()
    combo.addItem("", None)
    for shown, stored in items:
        combo.addItem(shown, stored)
    return combo


def _select(combo: QComboBox, stored: object) -> None:
    """Select the entry that stores this value; the empty entry when no entry does."""
    combo.setCurrentIndex(max(combo.findData(stored), 0))


def _button(text: str, slot) -> QPushButton:
    """A push button that the Enter key does not press: Enter in system_id must only look up."""
    button = QPushButton(text)
    button.setAutoDefault(False)
    button.clicked.connect(slot)
    return button


class QcDialog(QDialog):
    """Ask the QC verdict of every run, and the reason for an excluded run."""

    def __init__(self, run_indexes: list[int], parent=None):
        super().__init__(parent)
        self.setWindowTitle("Close session")
        self.verdicts: dict[int, QComboBox] = {}
        self.reasons: dict[int, QComboBox] = {}
        verdicts = _allowed("state.qc[].verdict")
        reasons = [(r, r) for r in _allowed("state.qc[].reason")]

        form = QFormLayout()
        for i in run_indexes:
            self.verdicts[i] = QComboBox()
            self.verdicts[i].addItems(verdicts)
            self.reasons[i] = _combo(reasons)
            row = QHBoxLayout()
            row.addWidget(self.verdicts[i])
            row.addWidget(self.reasons[i])
            form.addRow(f"run {i}", row)
        if not run_indexes:
            form.addRow(QLabel("This session has no runs."))

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout = QVBoxLayout()
        layout.addLayout(form)
        layout.addWidget(buttons)
        self.setLayout(layout)

    def qc(self) -> dict[int, tuple[str, str | None]]:
        """{run_index: (verdict, reason or None)}, the shape close_session takes."""
        return {i: (self.verdicts[i].currentText(), self.reasons[i].currentData())
                for i in self.verdicts}


class SessionForm(QDialog):
    """Pages Project, Mouse, Calibration, Session; Create makes the session; then Add run and Close."""

    def __init__(self, root: Path, parent=None):
        super().__init__(parent)
        self.root = Path(root)
        self.setWindowTitle(TITLE)
        cfg = load_rig(self.root)
        self.staging_root = Path(cfg.staging_root) if cfg else None
        self.session_dir: Path | None = None
        self._session_open = False             # created or found open, and not closed yet
        self._closed_by: str | None = None     # the experimenter when the session was created
        self._subject_id: str | None = None    # set when the person confirms an existing mouse
        self._found: dict | None = None        # the sidecar the last look-up found
        self._forget_lists()

        self.stack = QStackedWidget()
        self.stack.addWidget(self._project_page(cfg))
        self.stack.addWidget(self._mouse_page())
        self.stack.addWidget(self._calibration_page())
        self.stack.addWidget(self._session_page(cfg))
        self.status_page = self._status_page()
        self.stack.addWidget(self.status_page)

        self.title = QLabel()
        self.back_button = _button("Back", lambda: self._show_page(self.stack.currentIndex() - 1))
        self.next_button = _button("Next", lambda: self._show_page(self.stack.currentIndex() + 1))
        self.create_button = _button("Create", self._create)
        nav = QHBoxLayout()
        nav.addWidget(self.back_button)
        nav.addStretch()
        nav.addWidget(self.next_button)
        nav.addWidget(self.create_button)

        layout = QVBoxLayout()
        layout.addWidget(self.title)
        layout.addWidget(self.stack)
        layout.addLayout(nav)
        self.setLayout(layout)

        still_open = self._read(open_session, self.root)
        if still_open is None:
            self._show_page(PROJECT)
        else:   # a session left open earlier: go straight to it, so it can be closed
            self.session_dir = still_open
            self._session_open = True
            # the experimenter who created it closes it: the first name in its sidecar
            sidecar = self._read(read_json, still_open / SIDECAR_NAME) or {}
            self._closed_by = (sidecar.get("experimenter") or [None])[0]
            self.session_fields["experimenter"].setText(_shown(self._closed_by))
            self._refresh_status()
            self._show_page(STATUS)

    def _read(self, call, *args):
        """Run a look-up that reads sidecars. Show any error it raises and return None then."""
        try:
            return call(*args)
        except READ_ERRORS as e:
            QMessageBox.critical(self, TITLE, str(e))
            return None

    # pages

    def _project_page(self, cfg) -> QWidget:
        self.project = QComboBox()
        if cfg is not None:
            self.project.addItems([p.project_id for p in load_projects(self.root)])
        form = QFormLayout()
        form.addRow("project_id", self.project)
        return self._page(form)

    def _mouse_page(self) -> QWidget:
        self.system_id = QLineEdit()
        self.system_id.returnPressed.connect(self._look_up)
        self.system_id.textEdited.connect(self._on_system_id_edited)
        self.look_up_button = _button("Look up", self._look_up)
        self.padded_label = QLabel()
        self.found_label = QLabel()
        self.found_label.setWordWrap(True)
        self.yes_button = _button("Yes, this is the mouse", self._use_found_mouse)
        self.no_button = _button("No", self._forget_mouse)
        self.yes_button.hide()
        self.no_button.hide()

        id_row = QHBoxLayout()
        id_row.addWidget(self.system_id)
        id_row.addWidget(self.look_up_button)
        answer_row = QHBoxLayout()
        answer_row.addWidget(self.yes_button)
        answer_row.addWidget(self.no_button)
        answer_row.addStretch()

        self.mouse_fields = {key: QLineEdit() for key in MOUSE_TEXT_KEYS}
        self.mouse_fields["date_of_birth"].setPlaceholderText("YYYY-MM-DD")
        self.sex = _combo(SEX_WORDS)
        form = QFormLayout()
        form.addRow("system_id", id_row)
        form.addRow("", self.padded_label)
        form.addRow("", self.found_label)
        form.addRow("", answer_row)
        for key, field in self.mouse_fields.items():
            form.addRow(key, field)
            if key == "strain":
                form.addRow("sex", self.sex)

        self.surgery_fields = {key: QLineEdit() for key in SURGERY_TEXT_KEYS}
        self.surgery_fields["date"].setPlaceholderText("YYYY-MM-DD")
        self.surgery_fields["procedure"].setText(DEFAULT_PROCEDURE)
        self.eye = _combo([(e, e) for e in _allowed("surgery[].eye")])
        surgery_form = QFormLayout()
        for key, field in self.surgery_fields.items():
            surgery_form.addRow(key, field)
        surgery_form.addRow("eye", self.eye)
        surgery_box = QGroupBox("surgery (the implant)")
        surgery_box.setLayout(surgery_form)

        layout = QVBoxLayout()
        layout.addLayout(form)
        layout.addWidget(surgery_box)
        return self._page(layout)

    def _calibration_page(self) -> QWidget:
        last = self._read(last_sidecar, self.staging_root) if self.staging_root else None
        source = last or CALIBRATION_DEFAULTS
        if last:
            note = ("These are the rig's constants from its last session "
                    f"({last.get('session_start_time')}). Change them only if the rig changed.")
        else:
            note = "This rig has no earlier session. These are the default constants."
        source_label = QLabel(note)
        source_label.setWordWrap(True)

        self.eye_tracking_method = _combo([(m, m) for m in _allowed("eye_tracking_method")])
        _select(self.eye_tracking_method, source.get("eye_tracking_method"))
        self.calibration_fields = {}
        for key, default in CALIBRATION_DEFAULTS.items():
            self.calibration_fields[key] = QLineEdit(_shown(source.get(key, default)))
        self.calibration_fields["stimulus_frequency_hz"].setPlaceholderText("numbers split by commas")

        self.rig_reconfigured = QCheckBox("rig_reconfigured")
        self.reconfigured_note = QLineEdit()
        self.reconfigured_note.setPlaceholderText("what changed")

        form = QFormLayout()
        form.addRow(source_label)
        form.addRow("eye_tracking_method", self.eye_tracking_method)
        for key, field in self.calibration_fields.items():
            form.addRow(key, field)
        form.addRow(self.rig_reconfigured, self.reconfigured_note)
        return self._page(form)

    def _session_page(self, cfg) -> QWidget:
        self.session_fields = {key: QLineEdit() for key in SESSION_KEYS}
        if cfg is not None:
            self.session_fields["experimenter"].setText(cfg.default_experimenter)
        form = QFormLayout()
        for key, field in self.session_fields.items():
            form.addRow(key, field)
        return self._page(form)

    def _status_page(self) -> QWidget:
        self.status_label = QLabel()
        self.status_label.setWordWrap(True)
        self.add_run_button = _button("Add run", self._add_run)
        self.close_button = _button("Close", self._close)
        buttons = QHBoxLayout()
        buttons.addWidget(self.add_run_button)
        buttons.addWidget(self.close_button)
        buttons.addStretch()
        layout = QVBoxLayout()
        layout.addWidget(self.status_label)
        layout.addLayout(buttons)
        layout.addStretch()
        return self._page(layout)

    @staticmethod
    def _page(layout) -> QWidget:
        page = QWidget()
        page.setLayout(layout)
        return page

    def _show_page(self, index: int) -> None:
        self.stack.setCurrentIndex(index)
        self.title.setText(PAGE_TITLES[index])
        self.back_button.setVisible(PROJECT < index < STATUS)
        self.next_button.setVisible(index < SESSION)
        self.create_button.setVisible(index == SESSION)

    # the mouse

    def _forget_lists(self) -> None:
        """The parts of a mouse record that have no field: they come only from a confirmed mouse."""
        self._alleles: list = []
        self._interval_log: list = []
        self._more_surgery: list = []    # surgery entries after the first (the implant row)

    def _forget_mouse(self) -> None:
        """Drop the confirmed mouse: the typed facts make a new mouse."""
        self._subject_id = None
        self._forget_lists()
        self.yes_button.hide()
        self.no_button.hide()
        if self._found is not None:
            self.found_label.setText("Type the mouse's facts below.")
        self._found = None

    def _on_system_id_edited(self, _text: str) -> None:
        self._forget_mouse()
        self._show_padded()

    def _show_padded(self) -> None:
        """Echo the system_id as it will be saved (padded to 14 digits), or why it cannot be."""
        text = self.system_id.text().strip()
        if not text:
            self.padded_label.setText("")
            return
        try:
            self.padded_label.setText(f"system_id: {pad_system_id(text)}")
        except ValueError as e:
            self.padded_label.setText(str(e))

    def _look_up(self) -> None:
        """Pad the typed system_id, show it, and look for a session of that mouse on this rig."""
        try:
            padded = pad_system_id(self.system_id.text().strip())
        except ValueError as e:
            QMessageBox.critical(self, TITLE, str(e))
            return
        self._forget_mouse()
        self.padded_label.setText(f"system_id: {padded}")
        found = self._read(find_mouse, self.staging_root, padded) if self.staging_root else None
        if found is None:
            self.found_label.setText("No session on this rig has this system_id. "
                                     "Type the mouse's facts below.")
            return
        self._found = found
        self.found_label.setText(
            f"colony_mouse_id: {found.get('colony_mouse_id')}\n"
            f"sex: {SEX_WORD_OF.get(found.get('sex'), found.get('sex'))}\n"
            f"date_of_birth: {found.get('date_of_birth')}\n"
            f"strain: {found.get('strain')}\n"
            "Is this the mouse?")
        self.yes_button.show()
        self.no_button.show()

    def _use_found_mouse(self) -> None:
        """Copy the found mouse's facts into the fields; the session keeps its subject_id."""
        facts = {k: self._found[k] for k in SUBJECT_KEYS if k in self._found and k != "schema_version"}
        subject_id = facts.pop("subject_id", None)
        self._forget_mouse()
        self._subject_id = subject_id
        self.system_id.setText(_shown(facts.get("system_id")))
        self._show_padded()
        for key, field in self.mouse_fields.items():
            field.setText(_shown(facts.get(key)))
        _select(self.sex, facts.get("sex"))
        surgery = facts.get("surgery") or [{}]
        for key, field in self.surgery_fields.items():
            field.setText(_shown(surgery[0].get(key)))
        _select(self.eye, surgery[0].get("eye"))
        self._more_surgery = list(surgery[1:])
        self._alleles = list(facts.get("alleles") or [])
        self._interval_log = list(facts.get("interval_log") or [])
        self.found_label.setText(f"Using mouse {subject_id}.")

    # answers and actions

    def answers(self) -> dict:
        """The typed answers, in the shape create_session takes.

        Raise ValueError naming the key when a number field holds text that is not a number.
        """
        implant = {key: _text(field) for key, field in self.surgery_fields.items()}
        implant["eye"] = self.eye.currentData()
        mouse = {"system_id": _text(self.system_id)}
        mouse.update({key: _text(field) for key, field in self.mouse_fields.items()})
        mouse.update(sex=self.sex.currentData(), alleles=list(self._alleles),
                     surgery=[implant] + self._more_surgery, interval_log=list(self._interval_log))

        session = {"eye_tracking_method": self.eye_tracking_method.currentData()}
        for key, field in self.calibration_fields.items():
            if key in TEXT_KEYS:
                session[key] = _text(field)
            elif key in LIST_KEYS:
                session[key] = _numbers(key, field.text())
            else:
                session[key] = _number(key, field.text())
        session["weight"] = _number("weight", self.session_fields["weight"].text())
        session["session_description"] = _text(self.session_fields["session_description"])
        session["experimenter"] = _text(self.session_fields["experimenter"])
        if self.rig_reconfigured.isChecked():
            session["rig_reconfigured"] = {"value": True, "note": self.reconfigured_note.text().strip()}

        answers = {"project_id": self.project.currentText() or None, "mouse": mouse, "session": session}
        if self._subject_id:
            answers["subject_id"] = self._subject_id
        return answers

    def _create(self) -> None:
        now = datetime.now().astimezone()   # the session starts the moment Create is pressed
        self._show_padded()                  # echo the system_id that is about to be saved
        try:
            answers = self.answers()
            if "subject_id" not in answers and not self._confirm_new_mouse(answers["mouse"]["system_id"]):
                return   # stay on the form: nothing is created
            self.session_dir = create_session(self.root, answers, now=now)
        except ERRORS as e:
            QMessageBox.critical(self, TITLE, str(e))
            return
        self._session_open = True
        self._closed_by = answers["session"]["experimenter"]
        self._refresh_status()
        self._show_page(STATUS)

    def _confirm_new_mouse(self, system_id: str | None) -> bool:
        """Ask whether to save this system_id for a new mouse. True when the person says Yes.

        A system_id that cannot be padded is not asked about: the store names the problem.
        """
        try:
            padded = pad_system_id(system_id)
        except ValueError:
            return True
        answer = QMessageBox.question(self, TITLE, f"Save system_id {padded} for a new mouse?")
        return answer == QMessageBox.StandardButton.Yes

    def _refresh_status(self) -> None:
        current = read_current(self.root)
        if current is None:
            self.status_label.setText(f"session folder: {self.session_dir}\nNo session is open.")
            return
        self.status_label.setText(f"session folder: {current['session_dir']}\n"
                                  f"entity_stem: {current['entity_stem']}\n"
                                  f"The next recording is run {current['next_run']:02d}. "
                                  "Press Add run after it finishes.")

    def _add_run(self) -> None:
        try:
            add_run(self.root)
        except ERRORS as e:
            QMessageBox.critical(self, TITLE, str(e))
            return
        self._refresh_status()

    def _close(self) -> None:
        """Ask a verdict for every run, close the session, show the result, and close the form."""
        current = read_current(self.root)
        runs = list(range(1, current["next_run"])) if current else []
        dialog = QcDialog(runs, parent=self)
        if not dialog.exec():
            return
        try:
            sidecar = close_session(self.session_dir, closed_by=self._closed_by, qc=dialog.qc(),
                                    lifecycle=LIFECYCLE_AT_CLOSE)
        except ERRORS as e:
            QMessageBox.critical(self, TITLE, str(e))
            return
        QMessageBox.information(self, TITLE,
                                f"Session closed at {sidecar['closed_at']}.\n"
                                f"{len(sidecar['closeout_files'])} files hashed.")
        self._session_open = False
        self.accept()

    def reject(self) -> None:
        """Esc and the window's close button come here. While a session is open, ask first."""
        if self._session_open:
            answer = QMessageBox.question(self, TITLE, LEAVE_QUESTION)
            if answer != QMessageBox.StandardButton.Yes:
                return
        super().reject()
