"""Hand out the next subject_id.

A subject_id is "m" + a 2-character rig slot + a 3-character sequence, for
example m7k222. All characters come from the 31-letter id alphabet. Each rig
owns one slot, so two rigs can mint offline and never collide. The next id is
the highest one already on disk, plus one.
"""

from pathlib import Path

from raymondlab.core import alphabet


class SlotFull(Exception):
    """Every one of the 31**3 sequence values in this rig's slot is already used."""


def next_subject_id(slot: str, staging_root: Path) -> str:
    """Return the next unused subject_id for this rig, e.g. "m7k222".

    Scan staging_root/*/raw/ for directories named sub-m<slot>???. Read each
    3-character tail as a number. The next sequence is the highest + 1, or 0
    when no directory exists yet. Raise SlotFull when no sequence is left.
    """
    prefix = f"sub-m{slot}"
    highest = -1
    for path in Path(staging_root).glob("*/raw/*"):
        name = path.name
        if not path.is_dir() or not name.startswith(prefix) or len(name) != len(prefix) + 3:
            continue
        try:
            highest = max(highest, alphabet.decode3(name[len(prefix):]))
        except ValueError:
            continue  # tail has a character outside the id alphabet
    seq = highest + 1
    if seq >= alphabet.BASE**3:
        raise SlotFull(f"every sequence in slot {slot!r} is already used")
    return "m" + slot + alphabet.encode3(seq)
