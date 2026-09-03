"""Offline tests: everything that does not need a real camera.

    python tests/test_offline.py

Run this after any change. It covers config merging, the crosshair geometry pixel by pixel, the
Spike2 handshake files, the anchor schedule, and a full end-to-end recording against a fake
capture whose brightness follows the exposure -- so the strobe anchor is exercised for real,
right through to the frame indices written into session.json.

What it CANNOT cover is anything about the real cameras. See README.md, "Bench checks".
"""
import json, os, sys, tempfile, threading, time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import cv2

APP = str(Path(__file__).resolve().parent.parent)
sys.path.insert(0, APP)

from eyecal import cameras, capture, config, display, record, session, spike2

FAILS = []


def check(name, cond, extra=""):
    print(("  ok   " if cond else "  FAIL ") + name + (f"  {extra}" if extra else ""))
    if not cond:
        FAILS.append(name)


# --- config ------------------------------------------------------------------------------
print("config")
cfg = config.load(None)
check("defaults load", cfg["camera"] == "auto" and cfg["seconds"] == 30.0)
real = config.load(os.path.join(APP, "config.json"), required=True)
check("config.json parses", real["camera"] == "auto" and real["out"] == "C:/Temp/test")
# The Spike2 script passes no --camera, so the shipped default is what actually runs.
check("the shipped default identifies the camera rather than assuming one",
      config.DEFAULTS["camera"] == cameras.AUTO == "auto" and real["camera"] == cameras.AUTO)
check("and auto is not itself a preset", "auto" not in cameras.CAMERA_PRESETS)
check("config.json has no unknown keys", set(real) <= set(config.DEFAULTS))

with tempfile.TemporaryDirectory() as d:
    bad = Path(d) / "bad.json"
    bad.write_text('{"camrea": "elp"}', encoding="utf-8")
    try:
        config.load(str(bad))
        check("typo in config raises", False)
    except KeyError as exc:
        check("typo in config raises", "camrea" in str(exc))


class Args:
    pass


a = Args()
a.camera, a.seconds, a.devices = "ov9281", None, [2, 1]
merged = config.apply_cli(config.load(None), a)
check("CLI overrides config", merged["camera"] == "ov9281" and merged["devices"] == [2, 1])
check("CLI None does not override", merged["seconds"] == 30.0)

# --- the recording preview's two settings are coupled ---------------------------------------
print("recording preview sizing")


def canvas_fits_window(cfg, camera, n_cams=2, screen=(1920, 1080)):
    """Does the recording canvas land INSIDE its window, or overflow it?

    Overflowing means every refresh downscales, measured 7x dearer per output pixel than the
    upscale path -- so a smaller window with an unchanged downsample is SLOWER than the full-size
    window it replaced. This is the check that keeps the two settings honest.
    """
    _, w, h = cameras.parse_format(cameras.CAMERA_PRESETS[camera]["format"])
    cw, ch = (w // cfg["rec_downsample"]) * n_cams, h // cfg["rec_downsample"]
    frac = cfg["rec_window_scale"]
    scale = min(screen[0] * frac / cw, screen[1] * frac / ch)
    return scale >= 1.0, (cw, ch), scale


shipped = config.load(os.path.join(APP, "config.json"), required=True)
for cam in sorted(cameras.CAMERA_PRESETS):
    ok, canvas, scale = canvas_fits_window(shipped, cam)
    check(f"shipped settings keep the {cam} canvas inside its window", ok,
          f"canvas {canvas[0]}x{canvas[1]}, letterbox scale {scale:.2f} "
          f"({'upscale' if scale >= 1 else 'DOWNSCALE'})")

# the mistake the coupling exists to prevent: shrink the window, forget the downsample
bad = dict(shipped, rec_downsample=shipped["downsample"])
ok, canvas, scale = canvas_fits_window(bad, "ov2311")
check("a shrunken window with the alignment downsample is caught", not ok,
      f"canvas {canvas[0]}x{canvas[1]}, scale {scale:.2f}")
check("Window accepts a screen fraction",
      "screen_fraction" in display.Window.__init__.__code__.co_varnames)

# --- format parsing / presets -------------------------------------------------------------
print("cameras")
check("parse_format", cameras.parse_format("MJPG_1920x1080") == ("MJPG", 1920, 1080))
try:
    cameras.parse_format("nope")
    check("bad format raises", False)
except ValueError:
    check("bad format raises", True)
check("all presets parse",
      all(cameras.parse_format(p["format"]) for p in cameras.CAMERA_PRESETS.values()))
# The identification rests entirely on this: if two presets shared a native full frame, asking
# for it would not say which camera answered.
check("no two presets share a native geometry",
      len({cameras.parse_format(p["format"])[1:] for p in cameras.CAMERA_PRESETS.values()})
      == len(cameras.CAMERA_PRESETS))

# --- identifying a camera from the geometries it offers ---------------------------------------
print("camera identification")


class GeometryCap:
    """Just enough camera to be identified: a mode list, and DirectShow's snapping.

    A geometry the camera does not have is never synthesised -- the driver snaps to the nearest
    mode it does have and reports that back. These mode lists are the real ones; the ELP numbers
    reproduce what was measured on the rig (ask 1280x800, get 1280x720; ask 1600x1200, get
    1920x1080).
    """

    ELP = {(320, 240), (640, 480), (800, 600), (1024, 768), (1280, 720), (1920, 1080)}
    OV2311 = {(160, 120), (320, 240), (640, 480), (800, 600), (1280, 720), (1280, 960),
              (1600, 1200)}
    OV9281 = {(320, 240), (640, 400), (640, 480), (1280, 720), (1280, 800)}
    VGA_ONLY = {(640, 480)}

    def __init__(self, supported, opened=True, width=640, height=480):
        self.supported, self.opened = supported, opened
        self.geometry = [width, height]
        self.asked = [width, height]
        self.released = False

    def isOpened(self):
        return self.opened

    def set(self, prop, value):
        if prop in (cv2.CAP_PROP_FRAME_WIDTH, cv2.CAP_PROP_FRAME_HEIGHT):
            self.asked[0 if prop == cv2.CAP_PROP_FRAME_WIDTH else 1] = int(value)
            w, h = self.asked
            self.geometry = list((w, h) if (w, h) in self.supported else
                                 min(self.supported,
                                     key=lambda m: (m[0] - w) ** 2 + (m[1] - h) ** 2))
        return True

    def get(self, prop):
        if prop == cv2.CAP_PROP_FRAME_WIDTH:
            return float(self.geometry[0])
        if prop == cv2.CAP_PROP_FRAME_HEIGHT:
            return float(self.geometry[1])
        return 0.0

    def release(self):
        self.released = True


def with_caps(caps, fn):
    """Run fn with cv2.VideoCapture handing out these caps, device 1 first."""
    real = cv2.VideoCapture
    cv2.VideoCapture = lambda index, backend=None: caps[index]
    try:
        return fn()
    finally:
        cv2.VideoCapture = real


probed = {}
for name, modes in (("elp", GeometryCap.ELP), ("ov2311", GeometryCap.OV2311),
                    ("ov9281", GeometryCap.OV9281)):
    cap = GeometryCap(modes)
    got, probes = with_caps([cap], lambda: cameras.identify_camera(1))
    probed[name] = probes
    check(f"identifies the {name}", got == name, got)
    check(f"{name}: exactly one preset got what it asked for",
          [pr["preset"] for pr in probes if pr["asked"] == pr["got"]] == [name],
          str([(pr["preset"], pr["got"]) for pr in probes]))
    check(f"{name}: the capture was released", cap.released)

snapped = {pr["preset"]: tuple(pr["got"]) for pr in probed["elp"]}
check("the ELP snaps 1280x800 to 1280x720, as measured on the rig",
      snapped["ov9281"] == (1280, 720), str(snapped["ov9281"]))
check("and 1600x1200 to 1920x1080", snapped["ov2311"] == (1920, 1080), str(snapped["ov2311"]))

unknown = GeometryCap(GeometryCap.VGA_ONLY)
try:
    with_caps([unknown], lambda: cameras.identify_camera(1))
    check("a camera matching no preset raises", False, "no exception")
except RuntimeError as exc:
    check("a camera matching no preset raises", "Could not identify" in str(exc))
    check("and shows every ask against what came back",
          all(f"{w}x{h}" in str(exc) for w, h in ((1920, 1080), (1600, 1200), (1280, 800))))
check("the capture is released on that path too", unknown.released)

shut = GeometryCap(GeometryCap.OV2311, opened=False)
try:
    with_caps([shut], lambda: cameras.identify_camera(1))
    check("a device that will not open is busy-class", False, "no exception")
except cameras.CameraBusyError:
    check("a device that will not open is busy-class", True)

dead = GeometryCap(GeometryCap.OV2311, width=-1, height=-1)
try:
    with_caps([dead], lambda: cameras.identify_camera(1))
    check("a device answering -1 is busy-class", False, "no exception")
except cameras.CameraBusyError:
    check("a device answering -1 is busy-class", True)
check("and that capture is released", dead.released)

pair = [GeometryCap(GeometryCap.OV2311), GeometryCap(GeometryCap.OV2311)]
name, evidence = with_caps(pair, lambda: cameras.detect_preset([1, 2]))
check("a matched pair identifies as one preset", name == "ov2311", name)
check("the evidence names every device and what it was",
      [(e["deviceId"], e["detected"]) for e in evidence] == [(1, "ov2311"), (2, "ov2311")],
      str([(e["deviceId"], e["detected"]) for e in evidence]))
check("the evidence carries the readbacks it was decided on",
      all(len(e["probes"]) == len(cameras.CAMERA_PRESETS) for e in evidence))

mixed = [GeometryCap(GeometryCap.OV2311), GeometryCap(GeometryCap.ELP)]
try:
    with_caps(mixed, lambda: cameras.detect_preset([1, 2]))
    check("a mixed pair raises", False, "no exception")
except RuntimeError as exc:
    check("a mixed pair raises", "not the same model" in str(exc))
    check("naming which device is which",
          "device 1: ov2311" in str(exc) and "device 2: elp" in str(exc), str(exc))
    check("and offering the way out", "--camera" in str(exc))
check("both captures released even so", all(c.released for c in mixed))

# --- display -------------------------------------------------------------------------------
print("display")
GREEN = display.COLOURS[0]
check("green is first", GREEN[0] == "green" and GREEN[1] == (0, 255, 0))
img = np.zeros((1080, 1920), np.uint8)
disp = display.display_copy(img, downsample=2)
check("display_copy downsamples and stays contiguous",
      disp.shape == (540, 960) and disp.flags["C_CONTIGUOUS"])
check("display_copy rotates", display.display_copy(img, 2, True).shape == (540, 960))

pane = display.alignment_pane(disp, 1, 7, GREEN[1], True, 49.4, 25.0, img.shape)
h, w = pane.shape[:2]
cx, cy = w // 2, h // 2
check("crosshair is one green column at the exact centre",
      tuple(pane[300, cx]) == (0, 255, 0) and tuple(pane[300, cx - 1]) == (0, 0, 0)
      and tuple(pane[300, cx + 1]) == (0, 0, 0), f"cx={cx}")
check("crosshair is one green row at the exact centre",
      tuple(pane[cy, 400]) == (0, 255, 0) and tuple(pane[cy - 1, 400]) == (0, 0, 0)
      and tuple(pane[cy + 1, 400]) == (0, 0, 0), f"cy={cy}")
bare = display.alignment_pane(disp, 1, 7, GREEN[1], False, 49.4, 25.0, img.shape)
check("hide overlay removes the reticle but keeps the border",
      tuple(bare[300, cx]) == (0, 0, 0) and tuple(bare[0, 400]) == (0, 255, 0))

two = display.tile([pane, display.alignment_pane(disp, 2, 8, display.COLOURS[1][1], True,
                                                 49.4, 25.0, img.shape)])
check("tile keeps each pane's own centre marked",
      two.shape == (540, 1920, 3) and tuple(two[300, 480]) == (0, 255, 0)
      and tuple(two[300, 960 + 480]) == (255, 0, 255) and tuple(two[300, 960]) != (0, 255, 0))
check("recording pane draws", display.recording_pane(disp, 1, 7, 123, 4).shape == (540, 960, 3))

# --- spike2 handshake ------------------------------------------------------------------------
print("spike2 handshake")
with tempfile.TemporaryDirectory() as d:
    flag = os.path.join(d, "sub", "ready.flag")
    check("no flag initially", not spike2.flag_exists(flag))
    session_dir = str(Path("C:/Temp/test/test01"))
    spike2.write_flag(flag, session_dir, note="unit test")
    check("flag created (parent dir made)", spike2.flag_exists(flag))
    lines = Path(flag).read_text(encoding="utf-8").splitlines()
    # Line 1 is a contract with the Spike2 script: one Read() after FileOpen must yield the
    # session directory and nothing else, so Spike2 can save its .smrx alongside the frames.
    check("line 1 is the session directory, alone", lines[0] == session_dir, repr(lines[0]))
    check("the path survives verbatim, backslashes and all",
          "\\" in lines[0] and lines[0].endswith("test01"), repr(lines[0]))
    check("flag is still human readable",
          any("alignment accepted" in ln for ln in lines[1:]))
    check("no temp files left behind",
          [p.name for p in Path(d, "sub").iterdir()] == ["ready.flag"])
    # The app must have no way to delete the flag: Spike2 owns its lifetime, and a stale flag
    # cleared here would be cleared seconds after Spike2 had already acted on it.
    check("the app cannot delete a flag", not hasattr(spike2, "clear_flag"))
    check("rewriting over an existing flag is fine",
          spike2.write_flag(flag, session_dir, note="again") is None
          and spike2.flag_exists(flag))

# --- anchor scheduling logic -------------------------------------------------------------
print("anchor logic")
anc = record.make_anchor(True, -10)
check("bright defaults to +4 stops", anc.bright == -6 and anc.enabled)
check("anchor off when disabled", not record.make_anchor(False, -10).enabled)
check("anchor off with no exposure", not record.make_anchor(True, None).enabled)
sched = record._anchor_schedule(anc, 30.0)
check("four scheduled events at the right times",
      [round(t, 3) for t, _, _ in sched] == [0.0, 0.25, 29.5, 29.75], str(sched))
check("start window", record._sample_window(anc, 30.0, 0.1) == "start")
check("middle is not sampled", record._sample_window(anc, 30.0, 15.0) is None)
check("end window", record._sample_window(anc, 30.0, 29.9) == "end")


# --- end to end against a fake camera --------------------------------------------------------
print("end-to-end recording (fake cameras)")

DARK, BRIGHT = 30, 200


class FakeCap:
    """Brightness follows the exposure, exactly as the anchor assumes the real one does."""

    def __init__(self, h=48, w=64, period=0.005):
        self.h, self.w, self.period = h, w, period
        self.exposure = -10.0
        self.n = 0
        self.released = False
        self._lock = threading.Lock()

    def read(self):
        time.sleep(self.period)
        with self._lock:
            self.n += 1
            level = BRIGHT if self.exposure > -8 else DARK
        f = np.full((self.h, self.w), level, np.uint8)
        f[0, 0] = self.n % 251          # a per-frame fingerprint for ordering checks
        return True, f

    def set(self, prop, value):
        if prop == cv2.CAP_PROP_EXPOSURE:
            with self._lock:
                self.exposure = float(value)
        return True

    def get(self, prop):
        return 0.0

    def release(self):
        self.released = True


def make_readers(n=2, devices=(1, 2)):
    rs = []
    for i in range(n):
        r = capture.CameraReader(FakeCap(), devices[i], want_gray=True, settle_s=0.05)
        rs.append(r)
        r.start()
    for r in rs:
        assert r.warmed.wait(5.0), "reader never warmed"
    return rs


out = tempfile.mkdtemp(prefix="eyecal_test_")
devices = [1, 2]
order = [1, 0]                       # deliberately swapped: cam1 must be device 2
readers = make_readers(2, devices)
sess = session.Session(out, "unit_test_session", 2)
anchor = record.make_anchor(True, -10, hold_s=0.1)
started = datetime.now(timezone.utc)

t_begin = time.perf_counter()
stop_reason, anchors = record.run_recording(readers, devices, order, sess, seconds=2.0,
                                            rotate180=True, anchor=anchor, preview_hz=0)
elapsed = time.perf_counter() - t_begin
check("ran for about the requested duration", 1.9 < elapsed < 3.0, f"{elapsed:.2f} s")
check("stop reason is the duration", stop_reason == "duration reached", stop_reason)
check("both cameras recorded frames", min(sess.recorded) > 50, str(sess.recorded))

request = {"camera": "test", "format": "MJPG_64x48", "deviceIds": devices, "seconds": 2.0,
           "queueCap": capture.QUEUE_MAX, "source": {"Exposure": -10}}
sess.write_sidecars()
tiffs = sess.write_first_frame_tiffs()
info = session.summarise(sess, started, request, readers, devices, order, anchors,
                         stop_reason, tiffs, True)
capture.close_all(readers)
check("close_all released every capture", all(r.cap.released for r in readers))
check("close_all joined every thread", not any(r.is_alive() for r in readers))

d = Path(info["sessionDir"])
check("session.json written last and present", (d / "session.json").exists())
check("ts.csv present", (d / "ts.csv").exists())
check("first-frame tiffs present", len(tiffs) == 2 and (d / tiffs[0]).exists())
check("device order recorded as accepted", info["acceptedDeviceOrder"] == [2, 1],
      str(info["acceptedDeviceOrder"]))
check("rotation recorded", info["rotate180Stored"] is True)

for k in (1, 2):
    side = json.loads((d / "frames" / f"c{k}.json").read_text())
    size = (d / "frames" / f"c{k}.bin").stat().st_size
    check(f"c{k}.bin size matches the sidecar",
          size == side["count"] * side["height"] * side["width"], f"{size} bytes")
    frames = session.read_frames(d, k)
    check(f"c{k} memmaps to the declared shape", frames.shape[0] == side["count"])

rows = (d / "ts.csv").read_text().strip().splitlines()
check("ts.csv has a row per frame", len(rows) - 1 == sum(sess.recorded), str(len(rows) - 1))
check("ts.csv header unchanged", rows[0] == ",".join(session.TS_COLUMNS))

anc_rep = info["strobeAnchors"]
check("anchor enabled in the report", anc_rep["enabled"] and anc_rep["brightExposure"] == -6)
for cam in ("cam1", "cam2"):
    marks = anc_rep["brightFrames"][cam]
    check(f"{cam} has bright frames at the start", len(marks["start"]) > 0, str(marks["start"]))
    check(f"{cam} has bright frames at the end", len(marks["end"]) > 0, str(marks["end"]))
    # Not frame 0: the reader applies a requested exposure after its NEXT read, so the first
    # stored frame is still dark. What matters is that the bright run is early and unbroken.
    check(f"{cam} start anchor lands within a few frames", min(marks["start"]) <= 3,
          str(marks["start"][:4]))
    for wnd in ("start", "end"):
        run = marks[wnd]
        check(f"{cam} {wnd} anchor is one contiguous run",
              run == list(range(min(run), max(run) + 1)), f"{len(run)} frames")
    check(f"{cam} start and end anchors do not overlap",
          max(marks["start"]) < min(marks["end"]))
    n_frames = info["perCamera"][int(cam[-1]) - 1]["nFrames"]
    check(f"{cam} end anchor is near the last frame",
          n_frames - max(marks["end"]) < 25, f"{max(marks['end'])} of {n_frames}")
check("no warnings on a clean run", not info["warnings"], str(info["warnings"]))

# the bright frames really are bright in the stored data
f1 = session.read_frames(d, 1)
b = anc_rep["brightFrames"]["cam1"]["start"]
mid = f1.shape[0] // 2
check("stored bright frames are saturated, middle frames are not",
      float(f1[b[0]].mean()) > 150 and float(f1[mid].mean()) < 60,
      f"anchor {float(f1[b[0]].mean()):.0f} vs middle {float(f1[mid].mean()):.0f}")

session.report(info)

# --- short recording disables the anchor rather than mangling it ----------------------------
print("short recording")
readers = make_readers(2, devices)
sess2 = session.Session(out, "short_session", 2)
stop2, anc2 = record.run_recording(readers, devices, [0, 1], sess2, seconds=0.4,
                                   anchor=record.make_anchor(True, -10, hold_s=0.25),
                                   preview_hz=0)
capture.close_all(readers)
check("anchor disabled when there is no room for it", anc2["enabled"] is False)
check("short run still recorded", min(sess2.recorded) > 5, str(sess2.recorded))

print("\n" + ("ALL TESTS PASSED" if not FAILS else f"{len(FAILS)} FAILED: {FAILS}"))
sys.exit(1 if FAILS else 0)
