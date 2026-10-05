"""The rig-setup dialog: collects answers and hands them to apply_setup.

No logic lives here. Every check is done by apply_setup, and every error it
raises is shown to the person, who can then fix the answer and try again.
"""

from pathlib import Path

from PySide6.QtWidgets import (QDialog, QDialogButtonBox, QFormLayout, QGroupBox, QHBoxLayout,
                               QLabel, QLineEdit, QMessageBox, QPushButton, QTableWidget,
                               QTableWidgetItem, QVBoxLayout)

from raymondlab.rig.config import Project, RigConfig, slot_is_locked
from raymondlab.rig.setup import SlotLocked, apply_setup

PROJECT_COLUMNS = ["project_id", "dataset_dir", "protocol_id", "lab", "institution"]

SLOT_HELP = ("Two characters from the admin's slot table. Two rigs must never share a slot; "
             "nothing can check that offline, so take the value from the table only.")


class RigSetupDialog(QDialog):
    """Ask for the rig's settings and its projects. OK applies them; Cancel changes nothing."""

    def __init__(self, root: Path, existing: RigConfig | None, projects: list[Project],
                 parent=None):
        super().__init__(parent)
        self.root = Path(root)
        self.existing = existing
        self.result_config: RigConfig | None = None
        self.setWindowTitle("Rig setup" if existing is None else "Rig setup (edit)")

        self.rig_id = QLineEdit(existing.rig_id if existing else "")
        self.mint_slot = QLineEdit(existing.mint_slot if existing else "")
        self.mint_slot.setMaxLength(2)
        self.staging_root = QLineEdit(existing.staging_root if existing
                                      else str(self.root / "staging"))
        self.experimenter = QLineEdit(existing.default_experimenter if existing else "")

        form = QFormLayout()
        form.addRow("Rig root", QLabel(str(self.root)))
        form.addRow("rig_id", self.rig_id)
        form.addRow("mint_slot", self.mint_slot)
        slot_note = QLabel(SLOT_HELP)
        slot_note.setWordWrap(True)
        form.addRow("", slot_note)
        form.addRow("staging_root", self.staging_root)
        form.addRow("default_experimenter", self.experimenter)

        if existing is not None and slot_is_locked(existing):
            self.mint_slot.setEnabled(False)
            self.mint_slot.setToolTip("Locked: an id has already been minted in this slot.")

        self.table = QTableWidget(0, len(PROJECT_COLUMNS))
        self.table.setHorizontalHeaderLabels(PROJECT_COLUMNS)
        for p in projects:
            self._add_project_row([getattr(p, c) for c in PROJECT_COLUMNS])

        add = QPushButton("Add project")
        add.clicked.connect(lambda: self._add_project_row([""] * len(PROJECT_COLUMNS)))
        remove = QPushButton("Remove selected")
        remove.clicked.connect(self._remove_selected_row)
        row_buttons = QHBoxLayout()
        row_buttons.addWidget(add)
        row_buttons.addWidget(remove)
        row_buttons.addStretch()

        projects_box = QGroupBox("Projects this rig records for")
        box_layout = QVBoxLayout()
        box_layout.addWidget(self.table)
        box_layout.addLayout(row_buttons)
        projects_box.setLayout(box_layout)

        buttons = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        buttons.accepted.connect(self._on_ok)
        buttons.rejected.connect(self.reject)

        layout = QVBoxLayout()
        layout.addLayout(form)
        layout.addWidget(projects_box)
        layout.addWidget(buttons)
        self.setLayout(layout)

    def _add_project_row(self, values: list[str]) -> None:
        r = self.table.rowCount()
        self.table.insertRow(r)
        for c, v in enumerate(values):
            self.table.setItem(r, c, QTableWidgetItem(v))

    def _remove_selected_row(self) -> None:
        r = self.table.currentRow()
        if r >= 0:
            self.table.removeRow(r)

    def answers(self) -> dict:
        """The typed answers, in the shape apply_setup takes."""
        projects = []
        for r in range(self.table.rowCount()):
            row = {}
            for c, name in enumerate(PROJECT_COLUMNS):
                item = self.table.item(r, c)
                row[name] = item.text().strip() if item else ""
            projects.append(row)
        return dict(rig_id=self.rig_id.text().strip(),
                    mint_slot=self.mint_slot.text().strip(),
                    staging_root=self.staging_root.text().strip(),
                    default_experimenter=self.experimenter.text().strip(),
                    projects=projects)

    def _on_ok(self) -> None:
        try:
            self.result_config = apply_setup(self.root, self.answers())
        except (ValueError, SlotLocked) as e:
            QMessageBox.critical(self, "Rig setup", str(e))
            return
        self.accept()
