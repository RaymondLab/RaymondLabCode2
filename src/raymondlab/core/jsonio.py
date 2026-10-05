"""One way in and out of every JSON file.

Every JSON file the pipeline writes goes through write_json. It sorts the keys
and fixes the layout, so the same content always gives the same bytes and two
files can be compared byte for byte. It writes to a temporary file and then
replaces the target in one step, so a crash never leaves a half-written file.
"""

import json
import os
from pathlib import Path


def write_json(path: Path, obj) -> None:
    """Write obj to path as JSON, in the one layout the whole code base uses.

    Sorted keys, 2-space indent, non-ASCII kept as typed, '\\n' line ends,
    trailing newline, UTF-8. Sorted keys are what makes identical content give
    identical bytes.

    The text goes to <path>.tmp first, then os.replace() moves it onto path in
    one step, so a crash never leaves a half-written file behind.
    """
    path = Path(path)
    text = json.dumps(obj, sort_keys=True, indent=2, ensure_ascii=False) + "\n"
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w", encoding="utf-8", newline="\n") as f:
        f.write(text)
    os.replace(tmp, path)


def read_json(path: Path):
    """Return the parsed contents of a JSON file. Errors propagate; none is swallowed."""
    return json.loads(Path(path).read_text(encoding="utf-8"))
