"""Entry point behind the session-manager command: set the rig up if needed, then open the form."""

import sys


def main(argv=None) -> int:
    """Console-script entry point for `session-manager`.

    If rig.json is missing, run the rig-setup dialog in this same process (not a subprocess),
    then open the session form. Return 0 when the form closes a session, 1 when the person
    cancels the rig setup or the form.
    """
    from PySide6.QtWidgets import QApplication

    from raymondlab.gui.session_form import SessionForm
    from raymondlab.rig.config import load_rig
    from raymondlab.rig.paths import rig_root
    from raymondlab.rig.setup import run_setup

    app = QApplication.instance() or QApplication(argv if argv is not None else sys.argv)
    root = rig_root()
    if load_rig(root) is None:
        run_setup()
        if load_rig(root) is None:
            print("session-manager: rig setup cancelled", file=sys.stderr)
            return 1
    return 0 if SessionForm(root).exec() else 1


if __name__ == "__main__":
    sys.exit(main())
