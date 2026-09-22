"""The ready flag: the app tells Spike2 that alignment is done and recording is starting.

Spike2's ProgStatus() reports only 1 = running and 0 = terminated -- it cannot carry the child's
exit code -- so nothing about the outcome can be returned that way. Three facts carry it all:

    this app CREATES the flag     alignment accepted; recording begins now. LINE 1 of it is the
                                  session directory, for Spike2 to save its own file alongside
    ProgStatus() drops to 0
      while the flag is absent    the operator cancelled, or start-up failed

    session.json exists           the recording completed (written last, see session.py)

THE APP NEVER DELETES THE FLAG. Spike2 owns its lifetime and must delete it before ProgRun, not
after: the app takes a couple of seconds just to import cv2, so anything it cleared at start-up
would be cleared long after a Spike2 poll loop had already seen the stale file and charged ahead.

Cancel needs no signal of its own -- the flag is simply never created.

Optional gating: if Spike2 DELETES the flag once it has started sampling, the app can be told to
wait for that (handshake_wait_s) before storing its first frame. That keeps camera start-up and
alignment strobes out of the Spike2 file entirely. With handshake_wait_s at 0 the app just
announces itself and records. Both are in README.md.
"""

import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path


def write_flag(path, session_dir, note=""):
    """Create the ready flag atomically, with the session directory on LINE 1.

    LINE 1 IS A CONTRACT. It is the full path of the directory the frames will be written to,
    alone on the line and nothing else, so Spike2 can lift it with one Read() straight after
    FileOpen and save its own .smrx alongside them:

        var fh% := FileOpen(readyPath$, 8, 0);
        var sessionDir$;
        Read(sessionDir$);
        FileClose();
        ... FileSaveAs(sessionDir$ + "/eye_calibration.smrx");

    Everything after line 1 is prose for whoever opens the file by hand. Anything machine-read
    that gets added later belongs on its own numbered line BELOW the path, never above it.

    The directory does not exist yet at this point -- it is created when storage begins, a
    moment later -- so Spike2 should save at the END of its script, by which time the recording
    has run. session.json in the same directory is the signal that it completed.

    Written to a temporary file in the same directory and then renamed, because os.replace is
    atomic on Windows: Spike2 polls this path several times a second and must never open a
    half-written file, nor see one exist before its contents are complete.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    body = (f"{session_dir}\n"
            f"eye-calibration-recording, {stamp}: alignment accepted, cameras streaming, "
            f"waiting for Spike2.\n"
            f"LINE 1 above is the session directory -- Spike2 reads it to save alongside the "
            f"frames.\nDelete this file to start the recording.\n{note}\n")
    fd, tmp = tempfile.mkstemp(dir=str(path.parent), prefix=".eyecal", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            fh.write(body)
        os.replace(tmp, path)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


def flag_exists(path):
    return Path(path).exists()
