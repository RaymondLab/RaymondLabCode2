"""End-to-end tests of the whole app with fake cameras and a fake Spike2.

    python tests/test_endtoend.py

Runs main() for real -- config, camera open, alignment, ready flag, recording, session.json,
exit code -- with cv2.VideoCapture replaced by a fake that honours whatever format it is asked
for, and the preview window replaced by one that presses keys on cue.

The flag protocol is asserted from Spike2's side: the app must CREATE the flag and never delete
it, must not hang when nobody answers, and must finish in about the requested duration.
"""

import importlib.util
import json
import sys
import tempfile
import threading
import time
from pathlib import Path

import cv2
import numpy as np

APP = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(APP))

from eyecal import display

FAILS = []


def check(name, cond, extra=""):
    print(("  ok   " if cond else "  FAIL ") + name + (f"  {extra}" if extra else ""))
    if not cond:
        FAILS.append(name)


def load_entry_point():
    """The entry point has hyphens in its name, so it cannot be imported normally."""
    spec = importlib.util.spec_from_file_location("entry", APP / "eye-calibration-recording.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class FakeCap:
    """Honours whatever geometry and codec it is asked for; brightness follows the exposure."""

    DARK, BRIGHT = 30, 200

    def __init__(self, index, backend=None):
        self.index = index
        self.props = {cv2.CAP_PROP_FRAME_WIDTH: 640, cv2.CAP_PROP_FRAME_HEIGHT: 480,
                      cv2.CAP_PROP_FOURCC: float(cv2.VideoWriter_fourcc(*"YUY2")),
                      cv2.CAP_PROP_FPS: 30.0, cv2.CAP_PROP_EXPOSURE: -10.0}
        self.released = False
        self.n = 0
        self._lock = threading.Lock()

    def isOpened(self):
        return not self.released

    def set(self, prop, value):
        with self._lock:
            self.props[prop] = float(value)
        return True

    def get(self, prop):
        with self._lock:
            return self.props.get(prop, 0.0)

    def read(self):
        time.sleep(0.004)
        with self._lock:
            self.n += 1
            n = self.n
            # Brightness doubles per stop and clamps at -6, as the real OV2311 does. It has to
            # be a continuous function now, not a threshold: open_camera proves manual control
            # by driving the exposure four stops DARKER and checking the image followed.
            stops = max(min(self.props[cv2.CAP_PROP_EXPOSURE], -6.0), -13.0) + 13.0
            level = min(4.0 * (2.0 ** stops), 255.0)
            h = int(self.props[cv2.CAP_PROP_FRAME_HEIGHT])
            w = int(self.props[cv2.CAP_PROP_FRAME_WIDTH])
        f = np.full((h, w, 3), int(level), np.uint8)   # BGR, as a real MJPG decode returns
        f[0, 0, 0] = n % 251
        return True, f

    def release(self):
        self.released = True


class FakeWindow:
    """Presses a scripted key after a few pumps; otherwise reports 'no key'."""

    key_after = {}          # window-name prefix -> (n_pumps, key)
    created = []            # (name, screen_fraction) in construction order

    def __init__(self, name, fullscreen=False, screen_fraction=0.92):
        self.name, self.pumps, self.shown = name, 0, 0
        self.screen_fraction = screen_fraction
        FakeWindow.created.append((name, screen_fraction))
        self._trigger = (10 ** 9, 255)
        for prefix, spec in FakeWindow.key_after.items():
            if name.startswith(prefix):
                self._trigger = spec

    def show(self, canvas):
        self.shown += 1

    def pump(self):
        self.pumps += 1
        return self._trigger[1] if self.pumps >= self._trigger[0] else 255

    def closed(self):
        return False

    def close(self):
        pass


def fake_spike2(flag, delay=0.3, delete=True, timeout=20.0):
    """Poll for the ready flag exactly as the .s2s loop does, then optionally delete it."""
    seen = {"flag": False, "at": None}

    def run():
        deadline = time.perf_counter() + timeout
        while time.perf_counter() < deadline:
            if Path(flag).exists():
                seen["flag"], seen["at"] = True, time.perf_counter()
                if delete:
                    time.sleep(delay)
                    try:
                        Path(flag).unlink()
                    except FileNotFoundError:
                        pass
                return
            time.sleep(0.01)

    t = threading.Thread(target=run, daemon=True)
    t.start()
    return seen, t


def run_app(argv, accept_key=ord("a")):
    FakeWindow.created = []
    FakeWindow.key_after = {"camera alignment": (5, accept_key)}
    real_capture, real_window = cv2.VideoCapture, display.Window
    cv2.VideoCapture = FakeCap
    display.Window = FakeWindow
    try:
        return load_entry_point().main(argv)
    finally:
        cv2.VideoCapture, display.Window = real_capture, real_window


BASE = ["--camera", "ov2311", "--devices", "1", "2", "--seconds", "1.5",
        "--anchor-hold-s", "0.15", "--downsample", "2"]

# --- the live protocol: Spike2 samples first, app signals, app records ------------------------
print("accept -> flag -> record (no gating, as the .s2s script runs it)")
out = tempfile.mkdtemp(prefix="eyecal_e2e_")
flag = str(Path(out) / "ready.flag")
spy, thread = fake_spike2(flag, delete=False)

t0 = time.perf_counter()
code = run_app(BASE + ["--out", out, "--session-id", "run01", "--ready-flag", flag])
elapsed = time.perf_counter() - t0
thread.join(timeout=5)

# 0 or 3: the fake camera has no hardware pacing, so synthesising 1.9 MB frames on a busy
# machine trips the inter-frame-gap guard and returns 3. That guard firing is correct behaviour;
# what this asserts is that the run completed and nothing worse happened.
check("run completed (0, or 3 with timing warnings)", code in (0, 3), str(code))
check("Spike2 saw the ready flag", spy["flag"])
check("the app did NOT delete the flag -- Spike2 owns it", Path(flag).exists())
check("recording stopped at the requested duration, did not hang", elapsed < 6.0,
      f"{elapsed:.1f} s for a 1.5 s recording")

d = Path(out) / "run01"
check("session.json exists (the Spike2 success signal)", (d / "session.json").exists())
# Line 1 of the flag is what Spike2 lifts with one Read() to save its .smrx beside the frames.
flag_line1 = Path(flag).read_text(encoding="utf-8").splitlines()[0]
check("flag line 1 is the session directory", Path(flag_line1) == d.resolve(), flag_line1)
check("and that is really where session.json landed",
      (Path(flag_line1) / "session.json").exists())
info = json.loads((d / "session.json").read_text())
check("stopped on duration", info["stopReason"] == "duration reached", info["stopReason"])
check("both cameras have frames", all(c["nFrames"] > 30 for c in info["perCamera"]),
      str([c["nFrames"] for c in info["perCamera"]]))
check("format was negotiated to the preset", info["requested"]["format"] == "MJPG_1600x1200")
check("no queue-cap warnings (storage path kept up)",
      not [w for w in info["warnings"] if "queue" in w], str(info["warnings"]))
check("frames stored as mono", all(c["height"] == 1200 and c["width"] == 1600
                                   for c in info["perCamera"]))
check("anchors found at both ends",
      all(info["strobeAnchors"]["brightFrames"][c]["start"] and
          info["strobeAnchors"]["brightFrames"][c]["end"] for c in ("cam1", "cam2")),
      str(info["strobeAnchors"]["brightFrames"]["cam1"]))
check("both windows were created (alignment, then recording)",
      len(FakeWindow.created) == 2 and FakeWindow.created[0][0].startswith("camera alignment")
      and FakeWindow.created[1][0].startswith("RECORDING"), str(FakeWindow.created))
check("the recording window is smaller than the alignment one",
      FakeWindow.created[1][1] < FakeWindow.created[0][1],
      f"alignment {FakeWindow.created[0][1]}, recording {FakeWindow.created[1][1]}")

# --- a stale flag is reported, never removed ---------------------------------------------------
print("stale flag from a previous run")
out2 = tempfile.mkdtemp(prefix="eyecal_e2e_")
flag2 = str(Path(out2) / "ready.flag")
Path(flag2).write_text("stale flag from a crashed run", encoding="utf-8")
code = run_app(BASE + ["--out", out2, "--session-id", "run02", "--ready-flag", flag2,
                       "--no-preview", "--seconds", "1.0"])
check("runs anyway with a stale flag present", code in (0, 3), str(code))
check("flag still there (only Spike2 may delete it)", Path(flag2).exists())

# --- gating: Spike2 deletes the flag to say go --------------------------------------------------
print("gated start (--handshake-wait-s, Spike2 deletes the flag)")
out3 = tempfile.mkdtemp(prefix="eyecal_e2e_")
flag3 = str(Path(out3) / "ready.flag")
spy3, thread3 = fake_spike2(flag3, delay=0.5, delete=True)
t0 = time.perf_counter()
code = run_app(BASE + ["--out", out3, "--session-id", "run03", "--ready-flag", flag3,
                       "--handshake-wait-s", "5.0", "--no-preview", "--seconds", "1.0"])
elapsed = time.perf_counter() - t0
thread3.join(timeout=5)
check("gated run completed", code in (0, 3), str(code))
check("Spike2 deleted the flag", not Path(flag3).exists())
check("proceeded on the delete, not on the timeout", elapsed < 4.0, f"{elapsed:.1f} s")

# --- gating with no answer must not hang ---------------------------------------------------------
print("gated start with no answer")
out4 = tempfile.mkdtemp(prefix="eyecal_e2e_")
flag4 = str(Path(out4) / "ready.flag")
t0 = time.perf_counter()
code = run_app(BASE + ["--out", out4, "--session-id", "run04", "--ready-flag", flag4,
                       "--handshake-wait-s", "1.0", "--no-preview", "--seconds", "1.0"])
elapsed = time.perf_counter() - t0
check("records anyway when Spike2 never answers", code in (0, 3), str(code))
check("session.json still written", (Path(out4) / "run04" / "session.json").exists())
check("waited about the stated time, then carried on", 1.5 < elapsed < 5.0, f"{elapsed:.1f} s")

# --- cancel ---------------------------------------------------------------------------------
print("cancel at the alignment stage")
out5 = tempfile.mkdtemp(prefix="eyecal_e2e_")
flag5 = str(Path(out5) / "ready.flag")
code = run_app(BASE + ["--out", out5, "--session-id", "run05", "--ready-flag", flag5],
               accept_key=27)          # Esc

check("exit code 4 (cancelled)", code == 4, str(code))
check("no ready flag was ever created", not Path(flag5).exists())
check("no session directory", not (Path(out5) / "run05").exists())
check("only the alignment window was created", len(FakeWindow.created) == 1,
      str(FakeWindow.created))

# --- manual exposure has to be PROVEN, not assumed ---------------------------------------------
class AutoCap(FakeCap):
    """Starts in auto-exposure and only 0.75 or 1.0 hands over control -- the real OV2311.

    Auto pins the image near its own target no matter what exposure is written, which is exactly
    what made a whole session look plausible and be untrustworthy.
    """

    AUTO_TARGET = 120
    escapes = (0.75, 1.0)

    def __init__(self, index, backend=None):
        super().__init__(index, backend)
        self.auto = True

    def set(self, prop, value):
        if prop == cv2.CAP_PROP_AUTO_EXPOSURE:
            self.auto = float(value) not in self.escapes
        return super().set(prop, value)     # always True, as the real driver reports

    def read(self):
        ok, frame = super().read()
        if ok and self.auto:
            frame[:] = self.AUTO_TARGET
            frame[0, 0, 0] = self.n % 251
        return ok, frame


class StuckAutoCap(AutoCap):
    escapes = ()                            # nothing hands over control


def run_with(cap_class, argv):
    real_capture, real_window = cv2.VideoCapture, display.Window
    FakeWindow.created = []
    FakeWindow.key_after = {"camera alignment": (5, ord("a"))}
    cv2.VideoCapture = cap_class
    display.Window = FakeWindow
    try:
        return load_entry_point().main(argv)
    finally:
        cv2.VideoCapture, display.Window = real_capture, real_window


print("a camera that starts in auto-exposure")
out6 = tempfile.mkdtemp(prefix="eyecal_e2e_")
code = run_with(AutoCap, BASE + ["--out", out6, "--session-id", "run06", "--no-handshake",
                                 "--no-preview", "--seconds", "1.0"])
check("recovers and records", code in (0, 3), str(code))
info6 = json.loads((Path(out6) / "run06" / "session.json").read_text())
ec = info6["perCamera"][0]["exposureControl"]
check("session.json records how manual was secured", ec is not None and "manualVia" in ec)
check("it found the value that actually works", ec["manualVia"] == 0.75, str(ec["manualVia"]))
check("and recorded the measured response", ec["responseRatio"] >= ec["minRatioRequired"],
      f"{ec['responseRatio']:.2f}x")
check("the failed candidate is in the record too",
      any(a["autoExposure"] == "as applied" and a["ratio"] < 1.5 for a in ec["attempts"]),
      str([(a["autoExposure"], round(a["ratio"], 2)) for a in ec["attempts"]]))
check("both cameras carry the evidence",
      all(c["exposureControl"] for c in info6["perCamera"]))

print("\na camera that will not leave auto-exposure")
out7 = tempfile.mkdtemp(prefix="eyecal_e2e_")
try:
    run_with(StuckAutoCap, BASE + ["--out", out7, "--session-id", "run07", "--no-handshake",
                                   "--no-preview", "--seconds", "1.0"])
    check("refuses to record", False, "no exception raised")
except RuntimeError as exc:
    check("refuses to record", "not under manual control" in str(exc))
    check("names every candidate it tried", "0.75" in str(exc) and "0.25" in str(exc))
check("no session directory was left behind", not (Path(out7) / "run07").exists())

print("\n" + ("ALL TESTS PASSED" if not FAILS else f"{len(FAILS)} FAILED: {FAILS}"))
sys.exit(1 if FAILS else 0)
