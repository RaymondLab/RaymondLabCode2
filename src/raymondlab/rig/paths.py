"""The one fixed path on the rig, and the names of the files that live in it.

This is the only place in the code base that knows where the rig's root
folder is. An environment variable overrides it, so tests and development on
macOS use a temp folder instead.
"""

import os
from pathlib import Path

RIG_JSON = "rig.json"
PROJECTS_JSON = "projects.json"
CURRENT_SESSION_JSON = "current-session.json"


def rig_root() -> Path:
    """The rig's root folder. RAYMONDLAB_RIG_ROOT overrides the fixed Windows path.

    Call this inside functions, never at import time, or the override in tests
    has no effect.
    """
    return Path(os.environ.get("RAYMONDLAB_RIG_ROOT", r"C:\RaymondLab"))
