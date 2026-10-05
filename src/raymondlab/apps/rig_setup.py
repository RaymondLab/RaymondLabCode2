"""Entry point behind the rig-setup command: open the setup dialog, print what it made."""

import sys


def main(argv=None) -> int:
    """Console-script entry point for `rig-setup`.

    Start a QApplication, run the setup dialog, print the resulting paths and return 0.
    Return 1 when the person cancels the dialog.
    """
    from PySide6.QtWidgets import QApplication

    from raymondlab.rig.paths import PROJECTS_JSON, RIG_JSON, rig_root
    from raymondlab.rig.setup import run_setup

    app = QApplication.instance() or QApplication(argv if argv is not None else sys.argv)
    root = rig_root()
    cfg = run_setup()
    if cfg is None:
        print("rig-setup: cancelled", file=sys.stderr)
        return 1
    print(f"rig root:      {root}")
    print(f"rig config:    {root / RIG_JSON}")
    print(f"projects:      {root / PROJECTS_JSON}")
    print(f"staging root:  {cfg.staging_root}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
