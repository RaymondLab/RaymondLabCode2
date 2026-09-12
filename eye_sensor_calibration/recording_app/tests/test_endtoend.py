"""End-to-end tests of the whole app with fake cameras and a fake Spike2.

    python tests/test_endtoend.py

Runs main() for real -- config, camera identification, camera open, alignment, ready flag,
recording, session.json, exit code -- with cv2.VideoCapture replaced by a fake that offers one
camera model's mode list, and the preview window replaced by one that presses keys on cue.

The flag protocol is asserted from Spike2's side: the app must CREATE the flag and never delete
it, must not hang when nobody answers, and must finish in about the requested duration.

THE FSIN TRIGGER IS FAKED AT ITS OWN SEAM, dshow.Controls, not switched off. FakeControls is the
control write, FakeControls.pulses is the pulse train, and FakeCap reads one frame per pulse
while it runs and the driver's all-zero timeout frame while it does not -- so the trigger path,
both fallbacks and the restore on every exit are exercised end to end. Nothing here may reach a
real camera: the seam is patched in run_app, which every scenario goes through.
"""

import contextlib
import importlib.util
import io
import json
import sys
import tempfile
import threading
import time
from pathlib import Path

import comtypes
import cv2
import numpy as np

APP = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(APP))

from eyecal import display, dshow, session

# The fake pulse rate, in the fake cameras and on the command line of the scenarios that give
# --pulse-hz. 100 Hz is a 10 ms period: fast enough that a second of it is a real recording,
# slow enough that a fake camera synthesising 1.9 MB frames can keep up.
PULSE_HZ = 100.0

# What a camera in trigger mode with no pulses does: answers ok=True with an ALL-ZERO frame after
# about a second, which is the driver's timeout and not an exposure. Slept in slices so the fake
# follows a train that starts, or a restore that puts it back to free-run, without sitting out
# the whole second first. See eyecal/trigger.py for the measurement.
NO_PULSE_S = 1.0
NO_PULSE_SLICE_S = 0.05

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


class FakePulseTrain(threading.Event):
    """The fake FSIN train: an Event to start and stop, plus the ONE CLOCK both cameras read.

    The clock is the point. Timing each fake camera off its own sleep lets the two drift a frame
    apart, which is exactly what cannot happen on the rig -- one train drives both FSIN pins, so
    both sensors expose on the SAME pulses, and equal frame counts is the invariant
    session._trigger_warnings checks for. Here a pulse therefore has a NUMBER and a DUE TIME, and
    a camera that has fallen behind (synthesising a 1.9 MB frame is not free) catches up rather
    than missing pulses the other one got.
    """

    t0 = 0.0                    # perf_counter when the train started
    t_stop = 0.0                # and when it stopped: a pulse due after this never fired

    def set(self):
        self.t0, self.t_stop = time.perf_counter(), 0.0
        super().set()

    def clear(self):
        self.t_stop = time.perf_counter()
        super().clear()


class FakeControls:
    """Stands in for dshow.Controls: the ONE seam through which the app writes the FSIN trigger.

    Patched in over eyecal.dshow.Controls for every run, so nothing in this file can reach a real
    camera's controls -- and so the writes themselves can be asserted, which is the only evidence
    that the app put both cameras into trigger mode and took them out again on the way out.

    The state is class-level because the app binds a fresh Controls per write (trigger.py does,
    deliberately) and the camera is what remembers the setting, exactly as the real one does
    across a release and a reopen.
    """

    ae = {}                             # OpenCV index -> AE priority: 1 trigger, 0 free-running
    writes = []                         # (index, prop_id, value) in order, the whole write log
    fail_indices = set()                # indices whose set() raises, as a camera without prop 19
    pulses = FakePulseTrain()           # the fake FSIN train; set = running

    def __init__(self, index):
        self.index = index
        self.friendly_name = f"Fake Camera {index}"

    def set(self, interface, prop_id, value, flags=2):
        if self.index in FakeControls.fail_indices:
            # E_FAIL, not the teardown HRESULT: trigger.write_ae_priority must fail FAST on this
            # one rather than sleeping between three attempts. See trigger._is_teardown_error.
            raise comtypes.COMError(-2147467259, "fake: no such camera control", None)
        FakeControls.writes.append((self.index, prop_id, int(value)))
        FakeControls.ae[self.index] = int(value)

    def get(self, interface, prop_id):
        return FakeControls.ae.get(self.index, 0), 2

    def read(self):
        return []

    def close(self):
        pass


def writes_for(index):
    """Just the values written to one device, in order. The pre-open/accept/exit sequence."""
    return [value for i, _prop, value in FakeControls.writes if i == index]


class FakeCap:
    """Offers one camera model's mode list; brightness follows the exposure.

    SUPPORTED is the OV2311's real list, so the shipped ov2311 preset works and everything else
    does not. A geometry the camera does not have is never synthesised: DirectShow snaps the
    request to the nearest mode it does have and reports THAT back, which is the whole basis of
    cameras.identify_camera. SUPPORTED_BY_INDEX overrides it per device, for a mixed pair.
    """

    DARK, BRIGHT = 30, 200
    SUPPORTED = {(160, 120), (320, 240), (640, 480), (800, 600), (1280, 720), (1280, 960),
                 (1600, 1200)}
    SUPPORTED_BY_INDEX = {}             # OpenCV index -> geometries, when the pair is not alike

    def __init__(self, index, backend=None):
        self.index = index
        self.supported = self.SUPPORTED_BY_INDEX.get(index, self.SUPPORTED)
        self.props = {cv2.CAP_PROP_FRAME_WIDTH: 640, cv2.CAP_PROP_FRAME_HEIGHT: 480,
                      cv2.CAP_PROP_FOURCC: float(cv2.VideoWriter_fourcc(*"YUY2")),
                      cv2.CAP_PROP_FPS: 30.0, cv2.CAP_PROP_EXPOSURE: -10.0}
        self.asked = [640, 480]         # width and height arrive as two separate writes
        self.released = False
        self.n = 0
        self.pulse_t0 = None            # the train this camera is counting pulses off
        self.pulse_k = 0                # and how many of them it has exposed
        self._lock = threading.Lock()

    def isOpened(self):
        return not self.released

    def set(self, prop, value):
        with self._lock:
            if prop in (cv2.CAP_PROP_FRAME_WIDTH, cv2.CAP_PROP_FRAME_HEIGHT):
                # The driver negotiates on the PAIR, so each write is judged against the other
                # value as it stands -- which is why asking for 1600 wide lands somewhere else
                # until the matching height arrives.
                self.asked[0 if prop == cv2.CAP_PROP_FRAME_WIDTH else 1] = int(value)
                w, h = self.asked
                got = ((w, h) if (w, h) in self.supported
                       else min(self.supported, key=lambda m: (m[0] - w) ** 2 + (m[1] - h) ** 2))
                self.props[cv2.CAP_PROP_FRAME_WIDTH] = float(got[0])
                self.props[cv2.CAP_PROP_FRAME_HEIGHT] = float(got[1])
            else:
                self.props[prop] = float(value)
        return True                     # True whether or not it was honoured, as the real one is

    def get(self, prop):
        with self._lock:
            return self.props.get(prop, 0.0)

    def read(self):
        """One frame: free-running, or one per FSIN pulse while this camera is under the trigger.

        MEASURED ON THE RIG and reproduced here because the app depends on all three (see
        eyecal/trigger.py): at AE priority 1 the camera exposes only on a pulse; with no pulses it
        does NOT fail, it answers ok=True with an ALL-ZERO frame after about a second; and written
        back to 0 it free-runs again immediately, which is the fallback path's whole premise.
        """
        while FakeControls.ae.get(self.index, 0) == 1:
            train = FakeControls.pulses
            if train.is_set() or self.pulse_t0 == train.t0:
                frame = self._next_pulse(train)
                if frame is not None:
                    return True, frame
            for _ in range(int(NO_PULSE_S / NO_PULSE_SLICE_S)):
                if FakeControls.pulses.is_set() or FakeControls.ae.get(self.index, 0) != 1:
                    break                       # the train started, or free-run was restored
                time.sleep(NO_PULSE_SLICE_S)
            else:
                return True, np.zeros(self._shape() + (3,), np.uint8)
        time.sleep(0.004)                       # free-running, at the camera's own rate
        return True, self._frame()

    def _next_pulse(self, train):
        """Wait for this camera's next pulse and expose it. None once there is no next pulse.

        The pulse is identified by NUMBER off the train's own clock, not by sleeping a period, so
        both cameras deliver the same ones however unevenly the fakes are scheduled. A camera
        behind the train does not sleep at all -- it works through its backlog and stops at the
        last pulse that really fired, which is where the driver's timeout frame then comes from.
        """
        if self.pulse_t0 != train.t0:
            self.pulse_t0, self.pulse_k = train.t0, 0       # a new train; count from its first
        due = train.t0 + self.pulse_k / PULSE_HZ
        while time.perf_counter() < due:
            if not train.is_set():
                return None                                 # stopped before this one was due
            time.sleep(0.001)
        if not train.is_set() and due > train.t_stop:
            return None                                     # caught up with a stopped train
        self.pulse_k += 1
        return self._frame()

    def _shape(self):
        with self._lock:
            return (int(self.props[cv2.CAP_PROP_FRAME_HEIGHT]),
                    int(self.props[cv2.CAP_PROP_FRAME_WIDTH]))

    def _frame(self):
        """A real exposure. Never all zero, so the app's timeout test can tell the two apart."""
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
        return f

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


def fake_spike2(flag, delay=0.3, delete=True, timeout=20.0, pulses=None):
    """Poll for the ready flag exactly as the .s2s loop does, then optionally delete it.

    `pulses` is the fake FSIN train: given an Event, it is SET `delay` after the flag appears,
    which is the operator's key on the real rig -- the flag is what tells Spike2 to start pulsing.
    """
    seen = {"flag": False, "at": None}

    def run():
        deadline = time.perf_counter() + timeout
        while time.perf_counter() < deadline:
            if Path(flag).exists():
                seen["flag"], seen["at"] = True, time.perf_counter()
                if delete or pulses is not None:
                    time.sleep(delay)
                if pulses is not None:
                    pulses.set()
                if delete:
                    try:
                        Path(flag).unlink()
                    except FileNotFoundError:
                        pass
                return
            time.sleep(0.01)

    t = threading.Thread(target=run, daemon=True)
    t.start()
    return seen, t


def stop_pulses_after(seconds, event=None, timeout=20.0):
    """Stop the fake pulse train `seconds` after it STARTS -- Spike2 reaching the end of its run.

    Timed from the set rather than from now, because the train does not start until the app has
    written the flag, and how long the app takes to get there is exactly what is not fixed.
    """
    event = FakeControls.pulses if event is None else event

    def run():
        if event.wait(timeout=timeout):
            time.sleep(seconds)
            event.clear()

    t = threading.Thread(target=run, daemon=True)
    t.start()
    return t


def run_app(argv, accept_key=ord("a"), cap_class=None, fail_indices=()):
    """main() against fake cameras, a fake preview window and a fake trigger control.

    dshow.Controls is patched alongside the other two and put back in the finally. That seam is
    not optional: unpatched, every run would write AE priority on whatever real cameras are
    plugged into this machine.
    """
    FakeWindow.created = []
    FakeWindow.key_after = {"camera alignment": (5, accept_key)}
    FakeControls.ae, FakeControls.writes = {}, []
    FakeControls.fail_indices = set(fail_indices)
    FakeControls.pulses.clear()
    real_capture, real_window = cv2.VideoCapture, display.Window
    real_controls = dshow.Controls
    cv2.VideoCapture = cap_class or FakeCap
    display.Window = FakeWindow
    dshow.Controls = FakeControls
    try:
        return load_entry_point().main(argv)
    finally:
        cv2.VideoCapture, display.Window = real_capture, real_window
        dshow.Controls = real_controls
        FakeControls.pulses.clear()         # never leave a train running into the next scenario


BASE = ["--camera", "ov2311", "--devices", "1", "2", "--seconds", "1.5",
        "--anchor-hold-s", "0.15", "--downsample", "2"]

# --- the live protocol: Spike2 samples first, app signals, app records ------------------------
# --no-trigger throughout this first group: these scenarios are about the flag protocol and the
# strobe anchor, which is the FREE-RUNNING path. The trigger has its own scenarios below.
print("accept -> flag -> record (no gating, as the .s2s script runs it)")
out = tempfile.mkdtemp(prefix="eyecal_e2e_")
flag = str(Path(out) / "ready.flag")
spy, thread = fake_spike2(flag, delete=False)

t0 = time.perf_counter()
code = run_app(BASE + ["--no-trigger", "--out", out, "--session-id", "run01",
                       "--ready-flag", flag])
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
# brightFrames only exists when the anchor ran, and under the trigger it deliberately does not
# (strobeAnchors is {"enabled": False, "reason": ...} there), so the enabled flag is read first.
anchors = info["strobeAnchors"]
check("anchors found at both ends",
      anchors["enabled"] and all(anchors["brightFrames"][c]["start"] and
                                 anchors["brightFrames"][c]["end"] for c in ("cam1", "cam2")),
      str(anchors.get("brightFrames", {}).get("cam1")))
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
code = run_app(BASE + ["--no-trigger", "--out", out2, "--session-id", "run02",
                       "--ready-flag", flag2, "--no-preview", "--seconds", "1.0"])
check("runs anyway with a stale flag present", code in (0, 3), str(code))
check("flag still there (only Spike2 may delete it)", Path(flag2).exists())

# --- gating: Spike2 deletes the flag to say go --------------------------------------------------
print("gated start (--handshake-wait-s, Spike2 deletes the flag)")
out3 = tempfile.mkdtemp(prefix="eyecal_e2e_")
flag3 = str(Path(out3) / "ready.flag")
spy3, thread3 = fake_spike2(flag3, delay=0.5, delete=True)
t0 = time.perf_counter()
code = run_app(BASE + ["--no-trigger", "--out", out3, "--session-id", "run03",
                       "--ready-flag", flag3, "--handshake-wait-s", "5.0", "--no-preview",
                       "--seconds", "1.0"])
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
code = run_app(BASE + ["--no-trigger", "--out", out4, "--session-id", "run04",
                       "--ready-flag", flag4, "--handshake-wait-s", "1.0", "--no-preview",
                       "--seconds", "1.0"])
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
# The trigger is left ENABLED here, unlike the scenarios above: a cancel must still write free-run
# on the way out, because the control persists in the camera and the next run would open onto it.
check("cancel never asked for trigger mode, and still restored free-run on the way out",
      writes_for(0) == [0, 0] and writes_for(1) == [0, 0], str(FakeControls.writes))

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
    """run_app with a different camera class. One seam, patched in one place."""
    return run_app(argv, cap_class=cap_class)


print("a camera that starts in auto-exposure")
out6 = tempfile.mkdtemp(prefix="eyecal_e2e_")
code = run_with(AutoCap, BASE + ["--no-trigger", "--out", out6, "--session-id", "run06",
                                 "--no-handshake", "--no-preview", "--seconds", "1.0"])
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

# --- identifying the camera instead of being told which one it is -----------------------------
ELP_MODES = {(320, 240), (640, 480), (800, 600), (1024, 768), (1280, 720), (1920, 1080)}


def one_line(text, limit=110):
    """Captured stderr, squashed to something a check line can carry."""
    return " ".join(text.split())[:limit]


def run_app_exit(argv, cap_class=FakeCap):
    """Run main() the way __main__ does: any exception is exit 1, with the reason on stderr."""
    err = io.StringIO()
    try:
        with contextlib.redirect_stderr(err):
            return run_with(cap_class, argv), err.getvalue()
    except Exception as exc:
        return 1, err.getvalue() + f"\nFAILED: {exc}"


print("\nno --camera given at all; config.json says auto")
out8 = tempfile.mkdtemp(prefix="eyecal_e2e_")
code = run_with(FakeCap, ["--devices", "1", "2", "--seconds", "1.0", "--anchor-hold-s", "0.15",
                          "--out", out8, "--session-id", "run08", "--no-handshake",
                          "--no-preview", "--no-trigger"])
check("run completed (0, or 3 with timing warnings)", code in (0, 3), str(code))
req = json.loads((Path(out8) / "run08" / "session.json").read_text())["requested"]
check("the preset used is the one the cameras answered to", req["camera"] == "ov2311",
      str(req["camera"]))
check("what was actually asked for is recorded beside it", req["cameraRequested"] == "auto",
      str(req["cameraRequested"]))
check("every device is in the detection evidence",
      [e["deviceId"] for e in req["cameraDetection"]] == [1, 2],
      str([(e["deviceId"], e["detected"]) for e in req["cameraDetection"]]))
check("both devices identified as the same camera",
      all(e["detected"] == "ov2311" for e in req["cameraDetection"]))
# The evidence is a readback, not a claim: the two presets this camera does not have were asked
# for and answered with something else, which is what makes the third one an identification.
probes = {pr["preset"]: pr for pr in req["cameraDetection"][0]["probes"]}
check("the matching preset got exactly what it asked for",
      probes["ov2311"]["asked"] == probes["ov2311"]["got"] == [1600, 1200], str(probes["ov2311"]))
check("the other presets were snapped to a mode this camera does have",
      probes["elp"]["got"] != probes["elp"]["asked"]
      and probes["ov9281"]["got"] != probes["ov9281"]["asked"],
      f"elp {probes['elp']['got']}, ov9281 {probes['ov9281']['got']}")

print("\ntwo cameras of different models")
out9 = tempfile.mkdtemp(prefix="eyecal_e2e_")
flag9 = str(Path(out9) / "ready.flag")
FakeCap.SUPPORTED_BY_INDEX = {1: ELP_MODES}         # device 2 is an ELP, device 1 is not
try:
    code, err = run_app_exit(["--devices", "1", "2", "--seconds", "1.0", "--out", out9,
                              "--session-id", "run09", "--ready-flag", flag9, "--no-preview"])
finally:
    FakeCap.SUPPORTED_BY_INDEX = {}
check("exit code 1 (failed)", code == 1, str(code))
check("it names what each device is", "device 1: ov2311" in err and "device 2: elp" in err,
      one_line(err))
check("and says how to override it", "--camera" in err)
check("no ready flag was ever created", not Path(flag9).exists())
check("no session directory", not (Path(out9) / "run09").exists())

print("\nan explicit --camera still overrides, and still fails loudly when it is wrong")
out10 = tempfile.mkdtemp(prefix="eyecal_e2e_")
flag10 = str(Path(out10) / "ready.flag")
code, err = run_app_exit(["--camera", "elp", "--devices", "1", "2", "--seconds", "1.0",
                          "--out", out10, "--session-id", "run10", "--ready-flag", flag10,
                          "--no-preview"])
check("exit code 1 (failed)", code == 1, str(code))
check("the preset was taken as given, and the format checked against it",
      "Driver would not deliver 1920x1080" in err, one_line(err))
check("the wrong-preset hint still names the right one",
      "re-run with --camera ov2311" in err, one_line(err[-140:]))
check("no ready flag was ever created", not Path(flag10).exists())

# --- the FSIN trigger, end to end -------------------------------------------------------------
# TRIG is BASE without --seconds and without the anchor hold: every scenario below sets its own,
# and the anchor only runs on the ones that fall back.
TRIG = ["--camera", "ov2311", "--devices", "1", "2", "--downsample", "2", "--no-preview"]

print("\ntriggered recording, stopped by the pulse train ending")
out11 = tempfile.mkdtemp(prefix="eyecal_e2e_")
flag11 = str(Path(out11) / "ready.flag")
# The flag is what starts the train, exactly as on the rig: Spike2 sees it and begins pulsing.
spy11, thread11 = fake_spike2(flag11, delay=0.3, delete=False, pulses=FakeControls.pulses)
stop_pulses_after(1.2)
t0 = time.perf_counter()
code = run_app(TRIG + ["--seconds", "3", "--pulse-hz", str(PULSE_HZ), "--out", out11,
                       "--session-id", "run11", "--ready-flag", flag11])
elapsed = time.perf_counter() - t0
thread11.join(timeout=5)

check("triggered run completed (0, or 3 with timing warnings)", code in (0, 3), str(code))
check("Spike2 saw the ready flag, which is what starts the train", spy11["flag"])
# Measured 6.2 s: about 2.4 s of opening and proving manual exposure, the 0.25 s mode settle, a
# 0.3 s wait for the train, its 1.2 s, the fakes' catch-up and the 1 s of silence that ends it.
# The 8 s bound is still decisive, because the ceiling is 8 s from mark_start and that would land
# near 11 s -- endReason is what says WHICH stop this was; this says it did not hang.
check("it stopped when the train did, not on the clock", elapsed < 8.0,
      f"{elapsed:.1f} s for a 1.2 s train")
info11 = json.loads((Path(out11) / "run11" / "session.json").read_text())
trig11 = info11["trigger"]
check("recorded under the trigger", trig11["mode"] == "trigger",
      f"{trig11['mode']} -- {trig11['reason']}")
check("the pulse train ending is what ended it",
      str(trig11["endReason"]).startswith("pulse train ended"), str(trig11["endReason"]))
check("no strobe anchor under the trigger: frame k is pulse k",
      info11["strobeAnchors"]["enabled"] is False, str(info11["strobeAnchors"]))
counts11 = [c["nFrames"] for c in info11["perCamera"]]
check("about one frame per pulse for 1.2 s of train", all(60 <= n <= 160 for n in counts11),
      str(counts11))
# Every camera sees the SAME pulses, which is the invariant the whole trigger path rests on.
check("both cameras stored the same number of frames", counts11[0] == counts11[1], str(counts11))
check("one timeout frame each -- the second of silence after the last pulse",
      all(n == 1 for n in trig11["timeouts"].values()), str(trig11["timeouts"]))
frames11 = session.read_frames(Path(out11) / "run11", 1)
check("no all-zero timeout frame was ever stored",
      all(f[::64, ::64].any() for f in frames11), f"{frames11.shape[0]} frames checked")
check("the control was written 0 before the open, 1 at accept, 0 on the way out",
      writes_for(0) == [0, 1, 0] and writes_for(1) == [0, 1, 0], str(FakeControls.writes))

print("\ntrigger requested, no pulses ever arrive")
out12 = tempfile.mkdtemp(prefix="eyecal_e2e_")
code = run_app(TRIG + ["--anchor-hold-s", "0.15", "--seconds", "1.5", "--trigger-wait-s", "1.0",
                       "--no-handshake", "--out", out12, "--session-id", "run12"])
info12 = json.loads((Path(out12) / "run12" / "session.json").read_text())
trig12 = info12["trigger"]
check("fell back to free-run", trig12["mode"] == "free-run", trig12["mode"])
check("and the reason says how long it waited", "no pulses within 1 s" in trig12["reason"],
      trig12["reason"])
check("a fallback is a warning, not a failure", code == 3, str(code))
check("the fallback put every camera back to free-run, and recorded the readback",
      all(e.get("value") == 0 for e in (trig12.get("aePriorityRestored") or {}).values()),
      str(trig12.get("aePriorityRestored")))
anchors12 = info12["strobeAnchors"]
check("the strobe anchor is back on, because that is now the mapping",
      anchors12["enabled"] and all(anchors12["brightFrames"][c]["start"] and
                                   anchors12["brightFrames"][c]["end"] for c in ("cam1", "cam2")),
      str(anchors12.get("brightFrames", {}).get("cam1")))
check("written 0 pre-open, 1 at accept, 0 on the fallback and 0 again on the way out",
      writes_for(0) == [0, 1, 0, 0] and writes_for(1) == [0, 1, 0, 0], str(FakeControls.writes))

print("\na camera that does not have the control at all (the ELP case)")
out13 = tempfile.mkdtemp(prefix="eyecal_e2e_")
t0 = time.perf_counter()
code = run_app(TRIG + ["--anchor-hold-s", "0.15", "--seconds", "1.0", "--no-handshake",
                       "--out", out13, "--session-id", "run13"],
               fail_indices={1})            # device 2 raises on every write
elapsed = time.perf_counter() - t0
info13 = json.loads((Path(out13) / "run13" / "session.json").read_text())
trig13 = info13["trigger"]
check("a half-written trigger falls back rather than recording two different things",
      trig13["mode"] == "free-run", trig13["mode"])
check("naming the device that refused",
      trig13["reason"].startswith("trigger write failed on device 2"), trig13["reason"])
check("a failed write is a warning, not a failure", code == 3, str(code))
check("and it never waited for a pulse it had no reason to expect", elapsed < 6.0,
      f"{elapsed:.1f} s")
check("the readback shows which device took it and which did not",
      trig13["aePriority"]["1"]["value"] == 1 and "error" in trig13["aePriority"]["2"],
      str(trig13["aePriority"]))

print("\ntrigger switched off entirely")
out14 = tempfile.mkdtemp(prefix="eyecal_e2e_")
code = run_app(TRIG + ["--anchor-hold-s", "0.15", "--seconds", "1.0", "--no-handshake",
                       "--no-trigger", "--out", out14, "--session-id", "run14"])
check("run completed (0, or 3 with timing warnings)", code in (0, 3), str(code))
trig14 = json.loads((Path(out14) / "run14" / "session.json").read_text())["trigger"]
check("session.json says it was never asked for", trig14["requested"] is False, str(trig14))
check("and it recorded free-running", trig14["mode"] == "free-run", trig14["mode"])
# --no-trigger must reach the camera controls too: not one COM write, on any device, anywhere.
check("no camera control was written at all", not FakeControls.writes, str(FakeControls.writes))

print("\na pulse train that never stops hits the ceiling")
out15 = tempfile.mkdtemp(prefix="eyecal_e2e_")
flag15 = str(Path(out15) / "ready.flag")
spy15, thread15 = fake_spike2(flag15, delay=0.3, delete=False, pulses=FakeControls.pulses)
code = run_app(TRIG + ["--seconds", "0.6", "--trigger-end-margin-s", "0.6", "--out", out15,
                       "--session-id", "run15", "--ready-flag", flag15])
thread15.join(timeout=5)
check("ceiling run completed (0, or 3 with timing warnings)", code in (0, 3), str(code))
trig15 = json.loads((Path(out15) / "run15" / "session.json").read_text())["trigger"]
check("a train that never ends is stopped by seconds + margin",
      str(trig15["endReason"]).startswith("ceiling reached"), str(trig15["endReason"]))

print("\n" + ("ALL TESTS PASSED" if not FAILS else f"{len(FAILS)} FAILED: {FAILS}"))
sys.exit(1 if FAILS else 0)
