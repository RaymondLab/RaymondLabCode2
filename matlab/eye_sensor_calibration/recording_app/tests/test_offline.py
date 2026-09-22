"""Offline tests: everything that does not need a real camera.

    python tests/test_offline.py

Run this after any change. It covers config merging, the crosshair geometry pixel by pixel, the
Spike2 handshake files, the anchor schedule, how manual exposure is proven (including the camera
that will not leave its own auto-exposure, and the one reopen it gets), the whole FSIN trigger
mechanism against a fake dshow.Controls and fake reader queues, and a full end-to-end recording
against a fake capture whose brightness follows the exposure -- so the strobe anchor is exercised
for real, right through to the frame indices written into session.json.

What it CANNOT cover is anything about the real cameras. See README.md, "Bench checks".
"""
import json, os, queue, sys, tempfile, threading, time
from datetime import datetime, timezone
from pathlib import Path

import comtypes
import numpy as np
import cv2

APP = str(Path(__file__).resolve().parent.parent)
sys.path.insert(0, APP)

from eyecal import (cameras, capture, config, display, dshow, record, session, spike2,
                    trigger)

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

# --- proving manual exposure, and the camera that will not leave auto ------------------------
print("manual exposure")

OV9281_FMT = cameras.CAMERA_PRESETS["ov9281"]["format"]
OV9281_SRC = dict(cameras.CAMERA_PRESETS["ov9281"]["source"])


class ExposureCap:
    """A camera whose image follows the exposure -- unless its own auto-exposure has hold of it.

    The mean is the measured shape of these sensors rather than an arbitrary curve: a black-level
    pedestal of about 28 counts that does not scale with integration time, a doubling per stop
    above it, and saturation at 255. Over the probe range that is 32 -> 255, about 8x, which is
    what manual control looks like. In auto the mean sits at 120 whatever is written, which is
    what a camera holding its own target looks like -- 1.00x, the same as a lens cap gives.

    `flips_manual` says which AUTO_EXPOSURE write hands control over: one value, a (previous
    write, this write) pair for a camera that only answers the transition, or None for a camera
    that never lets go. Once it has let go it stays manual, as a real one does.
    """

    AUTO_MEAN = 120

    def __init__(self, manual=False, flips_manual=None, h=40, w=64):
        self.manual, self.flips_manual = manual, flips_manual
        self.h, self.w = h, w
        self.exposure = -11.0
        self.geometry = [640, 480]      # positive before anything is set: _open_live checks
        self.fourcc, self.fps = 0, 0.0
        self.auto_writes = []
        self.released = False

    def isOpened(self):
        return True

    def read(self):
        level = min(255.0, 2 ** (self.exposure + 13) * 4 + 28) if self.manual else self.AUTO_MEAN
        return True, np.full((self.h, self.w), int(round(level)), np.uint8)

    def set(self, prop, value):
        if prop == cv2.CAP_PROP_EXPOSURE:
            self.exposure = float(value)
        elif prop == cv2.CAP_PROP_AUTO_EXPOSURE:
            previous = self.auto_writes[-1] if self.auto_writes else None
            self.auto_writes.append(float(value))
            want = self.flips_manual
            if isinstance(want, tuple):
                self.manual = self.manual or (previous, float(value)) == want
            elif want is not None:
                self.manual = self.manual or float(value) == float(want)
        elif prop == cv2.CAP_PROP_FRAME_WIDTH:
            self.geometry[0] = int(value)
        elif prop == cv2.CAP_PROP_FRAME_HEIGHT:
            self.geometry[1] = int(value)
        elif prop == cv2.CAP_PROP_FOURCC:
            self.fourcc = int(value)
        elif prop == cv2.CAP_PROP_FPS:
            self.fps = float(value)
        return True

    def get(self, prop):
        if prop == cv2.CAP_PROP_FRAME_WIDTH:
            return float(self.geometry[0])
        if prop == cv2.CAP_PROP_FRAME_HEIGHT:
            return float(self.geometry[1])
        if prop == cv2.CAP_PROP_FOURCC:
            return float(self.fourcc)       # what was set, so verify_format sees MJPG come back
        if prop == cv2.CAP_PROP_FPS:
            return self.fps
        return 0.0

    def release(self):
        self.released = True


def with_popped_caps(caps, fn):
    """Run fn with cv2.VideoCapture handing out a NEW cap on every call, popped off `caps`.

    The list is the caller's own, so what is left in it afterwards counts the opens -- which is
    the only way to count them when the call under test raises instead of returning.
    """
    real = cv2.VideoCapture
    cv2.VideoCapture = lambda index, backend=None: caps.pop(0)
    try:
        return fn()
    finally:
        cv2.VideoCapture = real


def with_fast_probe(fn):
    """Run fn with the exposure probe's settle time and the reopen wait at zero.

    _mean_after binds _PROBE_SETTLE_S as a DEFAULT ARGUMENT, so setting the module attribute is
    not enough on its own -- the default was captured when the function was defined. Both are
    set, and everything is put back afterwards, so nothing later runs against a patched module.
    Without this the cases below wait 0.3 s per measurement for a camera that is a numpy array.
    """
    settle, defaults = cameras._PROBE_SETTLE_S, cameras._mean_after.__defaults__
    wait = cameras.REOPEN_WAIT_S
    cameras._PROBE_SETTLE_S = 0.0
    cameras._mean_after.__defaults__ = (0.0,) + defaults[1:]
    cameras.REOPEN_WAIT_S = 0.0
    try:
        return fn()
    finally:
        cameras._PROBE_SETTLE_S, cameras._mean_after.__defaults__ = settle, defaults
        cameras.REOPEN_WAIT_S = wait


def open_one(caps):
    """cameras.open_camera as device 1, against these fake captures, without the settle waits."""
    return with_fast_probe(lambda: with_popped_caps(
        caps, lambda: cameras.open_camera(1, OV9281_FMT, OV9281_SRC)))


already = ExposureCap(manual=True)
_, description, ev = open_one([already])
check("a camera already under manual control needs no AUTO_EXPOSURE write at all",
      ev["manualVia"] == "as applied" and already.auto_writes == [], str(already.auto_writes))
check("it is proven by measurement, not by a readback",
      ev["responseRatio"] >= cameras.MIN_RESPONSE_RATIO, f"{ev['responseRatio']:.2f}x")
check("nothing was reopened", ev["reopened"] is False)
check("the exposure goes back to the preset value", already.exposure == -11.0,
      str(already.exposure))
check("the description says how manual control was secured", "as applied" in description,
      description)

escapes = ExposureCap(flips_manual=0.75)
_, _, ev = open_one([escapes])
check("a camera in auto is taken out of it by 0.75, the value measured on this rig",
      ev["manualVia"] == 0.75, str(ev["manualVia"]))
check("and the attempt that measured nothing stays in the evidence",
      [a["autoExposure"] for a in ev["attempts"]] == ["as applied", 0.75],
      str([a["autoExposure"] for a in ev["attempts"]]))
check("no reopen was needed", ev["reopened"] is False)

paired = ExposureCap(flips_manual=(0.25, 0.75))
_, _, ev = open_one([paired])
check("a camera that only answers the manual-then-auto transition is caught by the pairs",
      ev["manualVia"] == "0.25 then 0.75", str(ev["manualVia"]))
check("after every single value has been tried and measured first",
      [a["autoExposure"] for a in ev["attempts"]]
      == ["as applied", 0.75, 1.0, 0.25, 0.0, "0.25 then 0.75"],
      str([a["autoExposure"] for a in ev["attempts"]]))
check("the pair really is written in order, both values",
      paired.auto_writes == [0.75, 1.0, 0.25, 0.0, 0.25, 0.75], str(paired.auto_writes))

stuck = [ExposureCap(), ExposureCap()]
pool = list(stuck)
try:
    open_one(pool)
    check("a camera that never leaves auto refuses the run", False, "no exception")
except cameras.AutoExposureStuckError as exc:
    msg = str(exc)
    check("a camera that never leaves auto refuses the run", True)
    per_round = 1 + len(cameras.MANUAL_AE_CANDIDATES) + len(cameras.MANUAL_AE_PAIRS)
    check("every attempt's measurements are in the message, not just the last ratio",
          msg.count("at exposure ") == 2 * per_round,
          f"{msg.count('at exposure ')} of {2 * per_round} lines")
    check("and in the exception, for whoever catches it", len(exc.attempts) == 2 * per_round,
          str(len(exc.attempts)))
    check("the means are printed, because a pinned mean and a lens cap give the same ratio",
          "mean 120.0" in msg and "1.00x" in msg)
    check("both rounds are labelled", "first open" in msg and "after close and reopen" in msg)
    check("the Windows Camera app is named, and the way out of it",
          "Camera app" in msg and "Device Manager" in msg)
check("the camera was closed and reopened exactly once", pool == [], f"{len(pool)} caps unused")
check("both captures were released", all(c.released for c in stuck))

recovered = [ExposureCap(), ExposureCap(manual=True)]
pool = list(recovered)
_, description, ev = open_one(pool)
check("a reopen recovers a camera that comes up clean the second time",
      ev["manualVia"] == "as applied" and ev["reopened"] is True,
      f"{ev['manualVia']}, reopened {ev['reopened']}")
check("and the reopen is in what the operator is shown", "close and reopen" in description,
      description)
check("the stuck capture was released and the working one kept",
      recovered[0].released and not recovered[1].released)
check("the patched probe settings were put back",
      cameras._PROBE_SETTLE_S == 0.3 and cameras._mean_after.__defaults__[0] == 0.3,
      str(cameras._PROBE_SETTLE_S))

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

# --- camera_check: the parts that do not need a camera ---------------------------------------
print("camera_check")

# Imported with cv2.VideoCapture booby-trapped. A tool that claimed a device at import time
# would take the rig's cameras away from a running session every time anything so much as looked
# at it -- including this test run. The stub records the attempt rather than raising, so a
# regression is reported by check() instead of by a traceback that stops the whole file.
_opens_at_import = []
_real_videocapture = cv2.VideoCapture
cv2.VideoCapture = lambda *a, **k: _opens_at_import.append(a)
sys.path.insert(0, os.path.join(APP, "tools"))
try:
    import camera_check
finally:
    cv2.VideoCapture = _real_videocapture
check("importing camera_check opens no camera", not _opens_at_import, str(_opens_at_import))

# What dshow.Controls.read() gives back, with the properties an OV9281 does NOT expose already
# dropped -- that dropping is the module's job, so the tool must simply not invent lines for
# them. WhiteBalance really does read AUTO on both rig cameras out of the box.
ROWS = [
    {"interface": "CameraControl", "name": "Exposure", "value": -11, "min": -13, "max": -1,
     "step": 1, "default": -6, "auto": False, "flags": 2},
    {"interface": "VideoProcAmp", "name": "WhiteBalance", "value": 4600, "min": 2800,
     "max": 6500, "step": 1, "default": 4600, "auto": True, "flags": 1},
    {"interface": "VideoProcAmp", "name": "Gain", "value": 0, "min": 0, "max": 100,
     "step": 1, "default": 0, "auto": None, "flags": 0},
]

lines = camera_check.overlay_lines(2, 7, "Arducam OV9281 USB Camera", "ov9281",
                                   "Arducam OV9281", "1280x800 MJPG", 240.0, 49.2, ROWS)
texts = [t for t, _ in lines]
autos = [t for t, a in lines if a]
check("the first line is the pane, the device and the Windows name",
      texts[0] == "cam 2  device 7   Arducam OV9281 USB Camera", texts[0])
check("the second carries the preset, the media type and BOTH rates",
      texts[1] == "preset ov9281 (Arducam OV9281)   1280x800 MJPG   driver 240 fps   "
                  "measured 49.2 fps", texts[1])
check("one line per property the driver exposes", len(lines) == 2 + len(ROWS), str(len(lines)))
check("a property line carries value, range, step, default and mode",
      texts[2] == "Exposure                 -11   [-13..-1 step 1, default -6]   MANUAL",
      repr(texts[2]))
check("an AUTO property is flagged, so it draws in red", autos == [texts[3]], str(autos))
check("and says AUTO in the text too, which is what survives a screenshot",
      texts[3].endswith("AUTO") and "WhiteBalance" in texts[3], texts[3])
check("a property with neither flag is shown as unknown, not guessed as manual",
      texts[4].endswith("mode unknown"), texts[4])
# The one thing an overlay must never do is invent a control. dshow drops what GetRange refuses.
check("a property the driver does not expose gets no line at all",
      not any(name in "\n".join(texts) for name in ("Pan", "Tilt", "Zoom", "Focus", "Iris")))
check("the note from the last m is appended when there is one",
      camera_check.overlay_lines(1, 1, "n", "ov9281", "t", "1280x800 MJPG", 240.0, 49.2, ROWS,
                                 note="manual proven via AUTO_EXPOSURE=0.75 (3.02x)")[-1][0]
      == "manual proven via AUTO_EXPOSURE=0.75 (3.02x)")

no_com = camera_check.overlay_lines(1, 1, "n", "ov9281", "t", "1280x800 MJPG", 240.0, 49.2, [],
                                    controls_error="COMError: 0x80070005")
check("no COM readout is one line saying so, in red, pointing at m",
      len(no_com) == 3 and no_com[2][1] and "unavailable (COMError: 0x80070005)" in no_com[2][0]
      and "press m" in no_com[2][0], str(no_com[2]))

failed = camera_check.failed_lines(1, 3, "Could not open camera at OpenCV index 2.\n  hint\n\n"
                                         "  more\n  and more\n  and yet more")
check("a camera that would not open still gets a pane, and it is all red",
      failed[0][0] == "cam 1  device 3   NOT OPEN" and all(a for _, a in failed), str(failed))
check("only the first few lines of the message; the rest belongs on the console",
      len(failed) == 5, str(len(failed)))

# --- the raw snapshot: full resolution, nothing drawn on it ---------------------------------
frame = np.full((800, 1280), 40, np.uint8)
frame[0, 0] = 200                       # a corner marker, so a rotation is visible in the data

tiled = camera_check.raw_tile([frame, None])
check("a missing camera contributes a black pane the size of the one that answered",
      tiled.shape == (800, 2560), str(tiled.shape))
check("the working camera's frame is copied through at full resolution, untouched",
      np.array_equal(tiled[:, :1280], frame))
check("and the missing one is black", not tiled[:, 1280:].any())
check("nothing is drawn on it: no third channel and no colour to draw in", tiled.ndim == 2)

rotated = camera_check.raw_tile([frame, None], rotate180=True)
check("rotation is applied per frame, exactly as the preview rotates each pane",
      np.array_equal(rotated[:, :1280], cv2.rotate(frame, cv2.ROTATE_180))
      and rotated[799, 1279] == 200, str(rotated[799, 1279]))

swapped = camera_check.raw_tile([None, frame])
check("the swap order is honoured -- the frame lands on the right",
      np.array_equal(swapped[:, 1280:], frame) and not swapped[:, :1280].any())
check("no frames at all falls back to the requested geometry rather than writing nothing",
      camera_check.raw_tile([None, None], fallback_shape=(800, 1280)).shape == (800, 2560))
check("and gives up only when there is no size to use either",
      camera_check.raw_tile([None, None]) is None)

# --- snapshot names -------------------------------------------------------------------------
over, raw = camera_check.snapshot_paths("C:/Temp/test/camera_check_20260907-101500.png")
check("the dialog's extension is stripped before the two names are built",
      over.name == "camera_check_20260907-101500_overlay.png"
      and raw.name == "camera_check_20260907-101500_raw.png", f"{over.name}, {raw.name}")
check("both land in the directory that was chosen", str(over.parent).endswith("test"),
      str(over.parent))
over2, raw2 = camera_check.snapshot_paths("C:/Temp/test/check")
check("a path with no extension works the same",
      (over2.name, raw2.name) == ("check_overlay.png", "check_raw.png"),
      f"{over2.name}, {raw2.name}")
check("the default name is stamped to the second",
      camera_check.default_snapshot_name(datetime(2026, 9, 7, 10, 15, 0))
      == "camera_check_20260907-101500",
      camera_check.default_snapshot_name(datetime(2026, 9, 7, 10, 15, 0)))

# --- session_video: the trigger clock, and which clock a session gets -----------------------
print("session_video")

# Guarded, and the only import in this file that is. session_video pulls in imageio_ffmpeg and,
# through strobe_timing, sonpy and matplotlib -- three packages nothing else offline needs. A
# machine without them should lose these checks and keep the rest of the suite.
try:
    import session_video
except ImportError as exc:
    session_video = None
    print(f"  skip   session_video is not importable here: {exc}")

if session_video is not None:
    pulses = np.array([10.00, 10.01, 10.02, 10.03])

    t = session_video.frame_times_from_trigger(pulses, 4)
    check("frame k is pulse k, as seconds since pulse 0",
          np.allclose(t, [0.00, 0.01, 0.02, 0.03]), str(t))
    t = session_video.frame_times_from_trigger(pulses, 6)
    check("frames past the last pulse are NaN, not clamped or wrapped",
          np.allclose(t[:4], [0.00, 0.01, 0.02, 0.03]) and np.isnan(t[4:]).all(), str(t))
    check("fewer frames than pulses takes the pulses from the START of the train",
          np.allclose(session_video.frame_times_from_trigger(pulses, 2), [0.00, 0.01]))
    check("no frames gives an empty array rather than failing",
          session_video.frame_times_from_trigger(pulses, 0).shape == (0,))

    # Stand-ins for the three clocks, so the SELECTOR can be tested without a .smrx, a ts.csv
    # or sonpy. Each records that it was called, which is what the checks below are about:
    # a trigger session must not touch the strobe anchors, and a free-run one must not touch
    # TTL2. What each clock returns is tested elsewhere.
    TRIG = {1: np.zeros(3), 2: np.zeros(3)}
    STROBE = {1: np.ones(3), 2: np.ones(3)}
    TSCSV = {1: np.full(3, 2.0), 2: np.full(3, 2.0)}
    calls = []
    mismatch_now = {}
    trigger_fails = [False]

    def fake_trigger(session_dir, sess, counts):
        calls.append("trigger")
        if trigger_fails[0]:
            raise session_video.StrobeUnusable("Ch12 (TTL2) holds no pulses")
        return TRIG, "fake TTL2", dict(mismatch_now)

    def fake_strobe(session_dir, sess, counts):
        calls.append("strobe")
        return STROBE, "fake strobe"

    def fake_tscsv(session_dir, counts):
        calls.append("ts.csv")
        return TSCSV, "fake ts.csv"

    real_clocks = (session_video.trigger_times, session_video.strobe_times,
                   session_video.ts_csv_times)
    try:
        session_video.trigger_times = fake_trigger
        session_video.strobe_times = fake_strobe
        session_video.ts_csv_times = fake_tscsv
        counts = {1: 3, 2: 3}
        TRIGGERED = {"trigger": {"mode": "trigger"}}

        del calls[:]
        t, label = session_video.camera_times("d", TRIGGERED, counts)
        check("a trigger session is timed from TTL2 and says so on screen",
              calls == ["trigger"] and t is TRIG and label == session_video.TRIGGER_LABEL,
              f"{calls} {label!r}")

        del calls[:]
        mismatch_now = {1: (3, 4)}
        t, label = session_video.camera_times("d", TRIGGERED, counts)
        check("a pulse count that does not match still renders, with the label marked",
              calls == ["trigger"] and t is TRIG
              and label == session_video.TRIGGER_MISMATCH_LABEL
              and label.startswith(session_video.TRIGGER_LABEL), f"{calls} {label!r}")
        mismatch_now = {}

        del calls[:]
        trigger_fails[0] = True
        t, label = session_video.camera_times("d", TRIGGERED, counts)
        check("an unusable TTL2 falls back to ts.csv and never tries the strobe anchor",
              calls == ["trigger", "ts.csv"] and t is TSCSV
              and label == session_video.TSCSV_LABEL, f"{calls} {label!r}")
        trigger_fails[0] = False

        del calls[:]
        t, label = session_video.camera_times("d", {"trigger": {"mode": "free-run",
                                                                "requested": True}}, counts)
        check("a run that fell back to free-run keeps the strobe path",
              calls == ["strobe"] and label == session_video.STROBE_LABEL,
              f"{calls} {label!r}")

        del calls[:]
        t, label = session_video.camera_times("d", {"perCamera": []}, counts)
        check("and so does an older session with no trigger block at all",
              calls == ["strobe"] and label == session_video.STROBE_LABEL,
              f"{calls} {label!r}")
    finally:
        (session_video.trigger_times, session_video.strobe_times,
         session_video.ts_csv_times) = real_clocks

# --- eyecal/dshow: the enumeration, on whatever this machine has -----------------------------
print("dshow")
devices = dshow.list_devices()
check("list_devices returns a list without raising", isinstance(devices, list),
      f"{len(devices)} video input device(s): "
      + ", ".join(str(d["friendlyName"]) for d in devices))
check("every entry has an index, a name and a device path",
      all(set(d) == {"index", "friendlyName", "devicePath"} for d in devices))
check("indices are 0-based and in enumeration order, so device N is index N-1",
      [d["index"] for d in devices] == list(range(len(devices))))
check("both interfaces use the same two mode flags",
      (dshow.FLAG_AUTO, dshow.FLAG_MANUAL) == (1, 2))
check("the property ids are the DirectShow ones",
      dict(dshow.CAMERA_CONTROL_PROPS)["Exposure"] == 4
      and dict(dshow.VIDEO_PROC_AMP_PROPS)["WhiteBalance"] == 7)


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


# --- the FSIN trigger, with no cameras anywhere near it ---------------------------------------
print("trigger")
check("the trigger settings are in the defaults",
      all(k in config.DEFAULTS for k in ("trigger", "trigger_wait_s", "trigger_end_margin_s",
                                         "pulse_hz")),
      str([k for k in config.DEFAULTS if "trigger" in k or "pulse" in k]))
check("and the shipped config.json still carries no unknown key",
      set(real) <= set(config.DEFAULTS))

# The driver's timeout frame is ALL ZERO; a real one carries the sensor's black-level pedestal in
# every pixel (near 28 counts, see cameras.PROBE_DARK), which is why the test needs no threshold
# and why a stride-64 subsample of it is enough.
for label, blank, lit in (("mono", np.zeros((800, 1280), np.uint8),
                           np.full((800, 1280), 28, np.uint8)),
                          ("colour", np.zeros((800, 1280, 3), np.uint8),
                           np.full((800, 1280, 3), 28, np.uint8))):
    one_pixel = blank.copy()
    one_pixel[0, 0] = 28
    check(f"{label}: an all-zero frame is the driver's timeout", trigger.is_timeout_frame(blank))
    check(f"{label}: a frame with the pedestal everywhere is real",
          not trigger.is_timeout_frame(lit))
    check(f"{label}: one lit pixel at [0,0] is enough to call it real",
          not trigger.is_timeout_frame(one_pixel))

WROTE_1 = {"requested": 1, "value": 1, "flags": 2, "attempts": 1}
check("all_agree only when EVERY device read the value back",
      trigger.all_agree({"1": WROTE_1, "2": dict(WROTE_1)}, 1))
check("a device that could not be written breaks it",
      not trigger.all_agree({"1": WROTE_1, "2": {"requested": 1, "error": "COMError"}}, 1))
check("so does a driver that answered with something else",
      not trigger.all_agree({"1": WROTE_1, "2": dict(WROTE_1, value=0)}, 1))
check("and no devices at all is not agreement", not trigger.all_agree({}, 1))

TEARDOWN = comtypes.COMError(-2147024865, "device not functioning", None)
NO_CONTROL = comtypes.COMError(-2147467259, "no such property", None)
check("the teardown HRESULT is recognised signed or unsigned",
      trigger._is_teardown_error(TEARDOWN) and 0x8007001F in trigger.TEARDOWN_HRESULTS)
check("a camera without the control is NOT a teardown, so it is not retried",
      not trigger._is_teardown_error(NO_CONTROL))
check("neither is anything that is not a COMError",
      not trigger._is_teardown_error(RuntimeError("no device at index 1")))


class ControlsFake:
    """dshow.Controls as trigger.write_ae_priority uses it: bind, set, read back, close.

    `raises` is popped once per attempt -- an exception, or None to let the write through -- so
    "it failed once and then worked" is written as a list. `binds` counts the attempts, which is
    the thing actually under test: a failure that is not the teardown window must cost one.
    """

    raises = []
    binds = 0
    closed = 0
    value = 0

    def __init__(self, index):
        ControlsFake.binds += 1
        self.index = index

    def set(self, interface, prop_id, value, flags=2):
        exc = ControlsFake.raises.pop(0) if ControlsFake.raises else None
        if exc is not None:
            raise exc
        ControlsFake.value = value

    def get(self, interface, prop_id):
        return ControlsFake.value, 2

    def close(self):
        ControlsFake.closed += 1


def with_fake_controls(raises, fn, retry_wait=None):
    """Run fn with dshow.Controls replaced, and the retry wait shortened only if asked.

    The wait is left at its real one second by default, deliberately: "this failure costs no
    sleep at all" cannot be asserted against a wait that has been patched down to nothing.
    """
    ControlsFake.raises, ControlsFake.binds, ControlsFake.closed = list(raises), 0, 0
    ControlsFake.value = 0
    real_controls, real_wait = dshow.Controls, trigger.WRITE_RETRY_WAIT_S
    dshow.Controls = ControlsFake
    if retry_wait is not None:
        trigger.WRITE_RETRY_WAIT_S = retry_wait
    try:
        return fn()
    finally:
        dshow.Controls, trigger.WRITE_RETRY_WAIT_S = real_controls, real_wait


entry = with_fake_controls([], lambda: trigger.write_ae_priority(1, 1))
check("a write that took is reported with the READBACK, not with the request",
      entry == {"requested": 1, "value": 1, "flags": 2, "attempts": 1}, str(entry))
check("the filter is closed again either way", ControlsFake.closed == 1)

t_begin = time.perf_counter()
entry = with_fake_controls([NO_CONTROL, NO_CONTROL, NO_CONTROL],
                           lambda: trigger.write_ae_priority(2, 1))
elapsed = time.perf_counter() - t_begin
# The ELP production cameras do not expose property 19 at all. Retried, that device would cost
# 2 s at each of the three points every run writes this control, three stderr lines each time.
check("a camera without the control fails on the FIRST attempt", ControlsFake.binds == 1,
      f"{ControlsFake.binds} attempts")
check("and costs no sleep at all", elapsed < 0.1, f"{elapsed * 1e3:.0f} ms")
check("the reason is carried, and no value is invented",
      "error" in entry and "value" not in entry, str(entry))

t_begin = time.perf_counter()
entry = with_fake_controls([TEARDOWN], lambda: trigger.write_ae_priority(1, 0), retry_wait=0.01)
elapsed = time.perf_counter() - t_begin
check("the teardown window IS waited out, and the second try works",
      entry.get("value") == 0 and entry.get("attempts") == 2, str(entry))
check("the earlier failure is dropped, not kept beside the answer", "error" not in entry,
      str(entry))
check("it really did bind twice, and waited between them",
      ControlsFake.binds == 2 and elapsed >= 0.01, f"{ControlsFake.binds} binds, {elapsed:.3f} s")


# --- waiting for the first pulse, and recording one frame per pulse ---------------------------

class FakeReader:
    """Just enough reader for the triggered loops: a queue, no error, and a thread that is not.

    capture.check_readers asks for exactly these -- error, ident and is_alive -- and ident None
    is its "never started", which is the one case it does not read as a thread that has died.
    """

    def __init__(self, name="cam"):
        self.q = queue.Queue()
        self.error = None
        self.ident = None
        self.name = name

    def is_alive(self):
        return True


def pulse_item(mark, h=8, w=8):
    """One queued REAL frame, marked at [0, 0] so the stored order can be read back off disk."""
    t = time.perf_counter()
    return (t, time.time(), np.full((h, w), mark, np.uint8))


def timeout_item(h=8, w=8):
    """One queued frame that is the driver's ~1000 ms timeout: all zero, and never stored."""
    t = time.perf_counter()
    return (t, time.time(), np.zeros((h, w), np.uint8))


readers = [FakeReader("cam1"), FakeReader("cam2")]
for item in (timeout_item(), pulse_item(10), pulse_item(11)):
    readers[0].q.put(item)
for item in (timeout_item(), pulse_item(20)):
    readers[1].q.put(item)
ok, pending, seen, timeouts = record.wait_for_pulses(readers, [0, 1], 2.0)
check("the wait succeeds once EVERY camera has delivered a real frame", ok)
check("every real frame is kept, in arrival order -- they are pulses 0 onwards",
      [int(i[2][0, 0]) for i in pending[0]] == [10, 11]
      and [int(i[2][0, 0]) for i in pending[1]] == [20], str(seen))
check("timeout frames are counted and dropped", seen == [2, 1] and timeouts == [1, 1],
      f"seen {seen}, timeouts {timeouts}")

readers = [FakeReader("cam1"), FakeReader("cam2")]
readers[0].q.put(pulse_item(10))            # one camera answers, the other never does
t_begin = time.perf_counter()
ok, pending, seen, timeouts = record.wait_for_pulses(readers, [0, 1], 0.3)
elapsed = time.perf_counter() - t_begin
check("one camera alone is not a pulse train: the wait gives up at its deadline",
      not ok and 0.3 <= elapsed < 1.5, f"{elapsed:.2f} s")
check("and says what each camera did see", seen == [1, 0], str(seen))

out_trig = tempfile.mkdtemp(prefix="eyecal_trigger_")
readers = [FakeReader("cam1"), FakeReader("cam2")]
pending = [[pulse_item(10)], [pulse_item(20)]]
for item in (pulse_item(11), pulse_item(12), timeout_item()):
    readers[0].q.put(item)
for item in (pulse_item(21), timeout_item()):
    readers[1].q.put(item)
sess_trig = session.Session(out_trig, "trigger_unit", 2)
sess_trig.mark_start()              # the caller does this BEFORE the wait; see record.py
stop_trig, report_trig = record.run_triggered(readers, [1, 2], [0, 1], sess_trig, seconds=5.0,
                                              pending=pending, preview_hz=0)
check("the pulse train ending is what stops it, not the clock",
      stop_trig.startswith("pulse train ended"), stop_trig)
check("the frames collected during the wait are stored FIRST, as pulse 0 onwards",
      report_trig["framesDuringWait"] == {"cam1": 1, "cam2": 1}, str(report_trig))
check("one timeout frame each, and not one of them stored",
      report_trig["timeouts"] == {"cam1": 1, "cam2": 1} and sess_trig.recorded == [3, 2],
      f"{report_trig['timeouts']}, recorded {sess_trig.recorded}")
stored = np.fromfile(Path(sess_trig.dir) / "frames" / "c1.bin", dtype=np.uint8)
check("c1.bin holds exactly the real frames, in pulse order",
      list(stored.reshape(-1, 8, 8)[:, 0, 0]) == [10, 11, 12],
      str(list(stored.reshape(-1, 8, 8)[:, 0, 0])))

# --- the warnings that only mean anything under the trigger -----------------------------------
TRIG_OK = {"requested": True, "mode": "trigger", "pulseHz": 100.0, "timeouts": {"cam1": 1},
           "endReason": "pulse train ended (1 s without a pulse on every camera)"}
PER_CAM = [{"cam": 1, "nFrames": 100, "fpsEffective": 100.0},
           {"cam": 2, "nFrames": 100, "fpsEffective": 100.0}]
check("nothing to say when the rate, the counts and the train all agree",
      session._trigger_warnings(PER_CAM, TRIG_OK) == [],
      str(session._trigger_warnings(PER_CAM, TRIG_OK)))

slow = [dict(PER_CAM[0], fpsEffective=80.0), PER_CAM[1]]
w = session._trigger_warnings(slow, TRIG_OK)
check("a delivered rate more than 10 percent off the pulse rate is a warning",
      len(w) == 1 and "10 percent apart" in w[0], str(w))

uneven = [PER_CAM[0], dict(PER_CAM[1], nFrames=97)]
w = session._trigger_warnings(uneven, TRIG_OK)
check("unequal frame counts are a warning: both cameras see the SAME pulses",
      len(w) == 1 and "different frame counts" in w[0], str(w))

w = session._trigger_warnings(PER_CAM, dict(TRIG_OK, timeouts={"cam1": 3}))
# One timeout is how a normal run ENDS -- the second of silence after the last pulse -- so three
# of them is two gaps, not three.
check("every timeout beyond the last one is a gap in the train",
      len(w) == 1 and "2 one-second gaps" in w[0], str(w))

w = session._trigger_warnings(PER_CAM, {"requested": True, "mode": "free-run",
                                        "reason": "no pulses within 3 s"})
check("a fallback is a warning of its own: the session is read a different way",
      len(w) == 1 and "fell back to free-run" in w[0] and "strobe anchors" in w[0], str(w))
check("but a run that never asked for the trigger says nothing",
      session._trigger_warnings(PER_CAM, {"requested": False, "mode": "free-run"}) == [])

# --- the pane that goes black when the pulses stop --------------------------------------------
blank = display.trigger_pane(None, 1, 7, 0, 2, float("nan"), True, (240, 320))
check("a pane with no pulse is drawn black, at the geometry it would have used",
      blank.shape == (240, 320, 3) and not blank[150:230, 8:300].any(), str(blank.shape))
live = display.trigger_pane(np.full((240, 320), 40, np.uint8), 1, 7, 123, 0, 99.5, False,
                            (240, 320))
check("a live pane is three channels with a red border",
      live.shape == (240, 320, 3) and tuple(live[0, 100]) == (0, 0, 255), str(live[0, 100]))


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
