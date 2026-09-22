"""Opening a camera that actually honours what was asked, and working out which one it is.

Everything here talks to one cv2.VideoCapture and knows nothing about threads -- see capture.py
for that. The rule that matters is stated in open_negotiated: the ORDER of the property writes
is load-bearing, and was paid for in bench time.
"""

import re
import sys
import time

import cv2

# ---------------------------------------------------------------------------------------
# Camera presets. Copied from the bench project's cameras/*.json; keep them in step if a profile
# is ever corrected. The format listed is the sensor's native full frame.
# ---------------------------------------------------------------------------------------
CAMERA_PRESETS = {
    "elp": {
        "title": "ELP-USBFHD01M-L21 (production eye camera)",
        "format": "MJPG_1920x1080",
        "grayscale": True,
        "source": {"ExposureMode": "manual", "Exposure": -10, "WhiteBalanceMode": "manual",
                   "BacklightCompensation": "off", "Contrast": 32, "FrameRate": "max"},
    },
    "ov2311": {
        "title": "Arducam B0322 / OV2311",
        "format": "MJPG_1600x1200",
        "grayscale": True,
        "source": {"ExposureMode": "manual", "Exposure": -11, "WhiteBalanceMode": "manual",
                   "BacklightCompensation": "off", "Gain": 0, "FrameRate": "max"},
    },
    "ov9281": {
        "title": "Arducam OV9281",
        "format": "MJPG_1280x800",
        "grayscale": True,
        "source": {"ExposureMode": "manual", "Exposure": -11, "FrameRate": "max"},
    },
}

AUTO = "auto"               # camera setting meaning "ask the device"; see identify_camera
MAX_RATE_SENTINEL = 240.0   # rate requested for FrameRate:"max"; above any camera in scope

# --- proving manual exposure; see secure_manual_exposure -------------------------------------
# Ordered by what actually works on this rig, not by the UVC convention. Measured on both
# OV2311s: 0.75 and 1.0 hand over manual control, 0.25 and 0.0 leave the camera in auto -- the
# opposite of the usual reading. Every candidate is still confirmed by measuring the image, so a
# camera that does follow the convention is handled without a special case.
MANUAL_AE_CANDIDATES = (0.75, 1.0, 0.25, 0.0)

# Tried only after every single candidate above has failed. WHY PAIRS RATHER THAN MORE VALUES:
# OpenCV's DirectShow backend maps AUTO_EXPOSURE 0.25 and 0.0 to IAMCameraControl::Set with the
# MANUAL flag, and 0.75 and 1.0 to Set(the driver's default value, AUTO flag) -- four values, two
# calls. So there is nothing left to try except the ORDER they are written in. The bench probe
# that found 0.75 escaping had written a manual value immediately before it, which makes a
# Manual-then-Auto transition, and the app's own candidate order never contains that transition
# because it starts at 0.75. Each pair is written first value then second, and measured once.
MANUAL_AE_PAIRS = ((0.25, 0.75), (0.0, 1.0))

# Waited out between closing a camera that is stuck in auto-exposure and opening it again; see
# open_camera. Measured on this rig 2026-09-07: the stuck state lives in the camera and survives
# the application that caused it, so a reopen is worth one try and is not certain to help.
REOPEN_WAIT_S = 2.0

# The response probe spans the driver's full usable exposure range rather than stepping a few
# stops around the operating point. Measured on this rig, a step down from the preset decides
# nothing: -10 -> -13 moves the mean only 1.20x, because a black-level pedestal of about 28
# counts does not scale with integration time and swamps the signal at the dark end. That is
# indistinguishable from an auto-exposure camera ignoring the write. Across the full range the
# same two cameras move 3.6x and 4.6x, against auto's 1.00x -- no ambiguity left. Above -6 the
# driver clamps, exposure having reached the frame period, so -6 is the top.
PROBE_DARK = -13.0
PROBE_BRIGHT = -6.0

# Manual control moves the mean by 3.6-4.6x across that range on these cameras; auto pins it to
# within 1 percent. 1.5 sits in the empty space between, well clear of both.
MIN_RESPONSE_RATIO = 1.5

# The probe ends on the BRIGHT endpoint on purpose, so a restore that silently failed leaves the
# camera conspicuously bright rather than subtly wrong: 3x here, against the 1.2x that separates
# the dark endpoint from a typical operating point. Only checked when the preset sits well below
# the bright endpoint, since otherwise "restored" and "still bright" are the same picture.
MIN_RESTORE_DROP = 1.2
RESTORE_CHECK_MARGIN = 2.0

_PROBE_SETTLE_S = 0.3       # discarded after each write, so the measurement is of the new value
_PROBE_FRAMES = 6

_FORMAT_RE = re.compile(r"^([A-Za-z0-9]{1,4})_(\d+)x(\d+)$")
_PROP_MAP = {
    "Exposure": cv2.CAP_PROP_EXPOSURE, "Gain": cv2.CAP_PROP_GAIN,
    "Brightness": cv2.CAP_PROP_BRIGHTNESS, "Contrast": cv2.CAP_PROP_CONTRAST,
    "Saturation": cv2.CAP_PROP_SATURATION, "Sharpness": cv2.CAP_PROP_SHARPNESS,
    "Gamma": cv2.CAP_PROP_GAMMA, "BacklightCompensation": cv2.CAP_PROP_BACKLIGHT,
    "Focus": cv2.CAP_PROP_FOCUS,
}
_MODE_MAP = {
    "ExposureMode": cv2.CAP_PROP_AUTO_EXPOSURE, "WhiteBalanceMode": cv2.CAP_PROP_AUTO_WB,
    "FocusMode": cv2.CAP_PROP_AUTOFOCUS,
}
_ONOFF = {"on": 1, "off": 0, "enable": 1, "disable": 0}

# Appended to every "the device is not really open" error. A geometry readback of -1 or 0 is not
# a media type the driver chose; it is the property query itself failing, which on DirectShow
# nearly always means something else still holds the device.
#
# Ordered by escalation, cheapest first, because the expensive remedies are the ones that get
# reached for. There is no "wait and retry" step: the entry point has already done that, three
# times, before any of this is ever printed.
_HELD_HINT = (
    "\n  No camera reports a non-positive frame size, so this is not a format problem -- the "
    "device is\n  almost certainly still claimed by something. In order of escalation:\n"
    "    - another application has it open: the Windows Camera app, Teams, OBS, a browser tab.\n"
    "      This is the most common cause and the easiest to miss, because nothing warns you\n"
    "    - a lingering interpreter from an earlier run (PowerShell: Get-Process python)\n"
    "    - MATLAB, if it has opened this camera in this Windows session -- run imaqreset\n"
    "    - disable and re-enable the camera in Device Manager, under Cameras. This resets the\n"
    "      driver without touching the cable, and clears a device Windows has not reclaimed\n"
    "    - unplug and replug only if none of the above clears it")

# How to read the per-attempt means printed by _stuck_message. The ratio alone cannot separate
# the two ways this check fails -- a camera pinned in auto and a camera looking at a lens cap
# both report about 1.00x -- so the operator is given the numbers the verdict was made on.
_MEANS_GUIDE = (
    "\n  Reading the means: the same mid value on every attempt, dark and bright alike, "
    "is the camera\n  running its own auto-exposure and holding its target against everything "
    "written here. Every mean\n  near the black level -- below about 35 counts -- is a dark "
    "scene instead, a lens cap or no light,\n  which this check cannot tell apart from auto. "
    "Every mean near 255 is a saturated scene.")

# Appended to the auto-exposure refusal, the way _HELD_HINT is appended to the busy one. Measured
# on this rig 2026-09-07, on two OV9281s, after the Windows Camera app had been opened on them.
_AUTO_STUCK_HINT = (
    "\n\n  Close anything else using this camera and retry. If nothing else has it, the cause "
    "is most likely\n  Windows itself: the Camera app -- and anything else that goes through "
    "Windows' camera pipeline,\n  the Frame Server -- can leave the camera in an auto-exposure "
    "mode that these writes cannot undo.\n  The state stays in the camera after that "
    "application is closed, and has appeared a minute or two\n  AFTER closing it, so the "
    "cause can look long finished. In order of escalation:\n"
    "    - do not use the Windows Camera app on the rig cameras at all\n"
    "    - if it was used, wait two minutes and retry: the state has cleared by itself\n"
    "    - disable and re-enable the camera in Device Manager, under Cameras\n"
    "    - unplug and replug it")


class CameraBusyError(RuntimeError):
    """The device never really came up -- it is almost certainly still claimed by a previous run.

    Kept distinct from every other opening failure because it is the only one that can clear
    ITSELF: Windows usually reclaims a leaked device within a few seconds. The entry point
    retries on this exception and on nothing else, so asking for the wrong preset, or for a
    geometry the camera cannot deliver, still fails on the first attempt with the message that
    says so -- retrying those would only bury the diagnosis under two more identical failures.
    """


class AutoExposureStuckError(RuntimeError):
    """Nothing this app can write took the camera out of its OWN auto-exposure.

    Carries `attempts`: the list of measurement dicts, exactly as it goes into session.json. That
    is what lets open_camera catch this, reopen the device, and then report BOTH rounds instead
    of only the round that raised last.

    Deliberately NOT busy-class. The device opened, streamed, and answered every property query,
    so the entry point's few-second retry would change nothing -- the state lives inside the
    camera and clears on the order of minutes, if at all (see _AUTO_STUCK_HINT).
    """

    def __init__(self, message, attempts=()):
        super().__init__(message)
        self.attempts = list(attempts)


def parse_format(fmt):
    """MJPG_1920x1080 -> (MJPG, 1920, 1080)."""
    m = _FORMAT_RE.match(fmt or "")
    if not m:
        raise ValueError(f"Cannot parse format {fmt!r}; expected <FOURCC>_<W>x<H>, "
                         f"e.g. MJPG_1920x1080.")
    return m.group(1).upper(), int(m.group(2)), int(m.group(3))


def fourcc_str(value):
    value = int(value)
    if value <= 0:
        return ""
    s = "".join(chr((value >> (8 * i)) & 0xFF) for i in range(4))
    return "".join(ch for ch in s if ch.isprintable()).strip()


def negotiated(cap):
    """(width, height, fourcc) as the driver currently reports them."""
    return (int(cap.get(cv2.CAP_PROP_FRAME_WIDTH)),
            int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT)),
            fourcc_str(cap.get(cv2.CAP_PROP_FOURCC)))


def _open_live(device_id):
    """Open the device and prove it really came up. Returns (cap, OpenCV index).

    isOpened() is not a liveness test on DSHOW: it returns True for a capture whose device was
    never really acquired, and every get() afterwards answers -1. Caught here so the failure is
    named, instead of surfacing later as a bogus "driver would not deliver 1920x1080".
    """
    index = int(device_id) - 1        # MATLAB winvideo is 1-based; OpenCV is 0-based
    cap = cv2.VideoCapture(index, cv2.CAP_DSHOW)
    if not cap.isOpened():
        raise CameraBusyError(f"Could not open camera at OpenCV index {index} "
                              f"(device {device_id}).{_HELD_HINT}")
    if int(cap.get(cv2.CAP_PROP_FRAME_WIDTH)) <= 0:
        cap.release()
        raise CameraBusyError(f"Camera at OpenCV index {index} (device {device_id}) reports "
                              f"itself open but answers a non-positive frame "
                              f"width.{_HELD_HINT}")
    return cap, index


def open_negotiated(device_id, fmt, source):
    """Open the device and put it into the requested media type. Returns cap.

    ORDER IS LOAD-BEARING, and not in the way you would guess. Two DirectShow traps, both
    measured on this rig with OpenCV 5.0.0:

      1. FOURCC set BEFORE the geometry is ignored. A full-frame MJPG request silently yields
         uncompressed YUY2 instead -- over the USB 2.0 budget that is 5.19 fps against 49.4.
         The geometry reads back correctly and nothing raises, so the run looks legitimate.
      2. FPS set AFTER the FOURCC flips the media type straight back to YUY2. Every value does
         it, including the value that was just set.

    So the sequence is width -> height -> FPS -> FOURCC and the rate is never touched again.
    Source properties are applied afterwards and the format re-verified, in case one of them has
    the same side effect.

    Leaving the rate unset is NOT the same as asking for the maximum: at MJPG_1600x1200 the
    OV2311 default-negotiates 31.8 fps against 49.4 available, so max asks for a sentinel above
    anything in scope and lets the driver clamp it.

    Split out of open_camera only so this whole sequence can be run a SECOND time, identically,
    when a camera has to be closed and reopened to get it out of auto-exposure. Nothing here may
    be reordered or skipped on that second run; that is the point of it being one function.

    Public because tools/camera_check.py opens the cameras THE SAME WAY and deliberately stops
    here, without proving manual exposure -- it exists to show what state a camera is already in,
    so it hands over a source dict holding only FrameRate and lets apply_source find nothing to
    write. Any other caller wanting a camera to record from wants open_camera, not this.
    """
    fourcc, width, height = parse_format(fmt)
    cap, _ = _open_live(device_id)

    rate = source.get("FrameRate", "max")
    fps_hint = MAX_RATE_SENTINEL if str(rate).lower() == "max" else float(rate)
    cap.set(cv2.CAP_PROP_FRAME_WIDTH, width)
    cap.set(cv2.CAP_PROP_FRAME_HEIGHT, height)
    cap.set(cv2.CAP_PROP_FPS, fps_hint)
    cap.set(cv2.CAP_PROP_FOURCC, cv2.VideoWriter_fourcc(*fourcc.ljust(4)[:4]))
    verify_format(cap, width, height, fourcc, "negotiation")

    apply_source(cap, source)
    verify_format(cap, width, height, fourcc, "applying source properties")
    return cap


def open_camera(device_id, fmt, source, anchor_exposure=None):
    """Open a camera, force the format, prove manual exposure. Returns (cap, text, evidence).

    The property order that makes the format stick is in open_negotiated, and is the part of
    this that was paid for in bench time.

    Opening is deliberately sequential across cameras even though it costs ~7 s each. It all
    happens before Spike2 starts sampling, so it is free in the only currency that matters here,
    and OpenCV's DirectShow backend runs through a shared global that is not documented as
    thread-safe.

    `anchor_exposure`, when given, is also checked for re-negotiating the media type.

    ONE REOPEN, AND ONLY ONE. A camera that will not leave its own auto-exposure is closed,
    waited out for REOPEN_WAIT_S and opened again from scratch, then measured again. Measured on
    this rig 2026-09-07: this state arrives from outside the app -- the Windows Frame Server
    leaves it behind -- and a fresh graph on a fresh device sometimes starts clean where property
    writes on the running one never do. It is one try because the state also persists across a
    reopen often enough that a loop would only spend the operator's time; the failure message
    then carries BOTH rounds' measurements, which is the evidence for what to do next.
    """
    fourcc, width, height = parse_format(fmt)
    cap = open_negotiated(device_id, fmt, source)

    # Always, not only when the anchor is on: a camera left in auto-exposure invalidates the
    # exposure setting itself, not just the anchor.
    reopened = False
    try:
        evidence = secure_manual_exposure(cap, source, width, height, fourcc, device_id,
                                          anchor_exposure)
    except AutoExposureStuckError as first:
        print(f"  Device {device_id}: exposure stayed in auto after "
              f"{_write_count(first.attempts)} writes; closing and reopening once before giving "
              f"up.", file=sys.stderr)
        # secure_manual_exposure already released it on that path. Released again here anyway,
        # because this function must not depend on where inside the call the failure happened,
        # and a device still claimed by this process is one the reopen cannot get back.
        cap.release()
        time.sleep(REOPEN_WAIT_S)
        # A CameraBusyError raised by the reopen propagates untouched: THAT one clears itself,
        # and the entry point already retries it.
        cap = open_negotiated(device_id, fmt, source)
        reopened = True
        try:
            evidence = secure_manual_exposure(cap, source, width, height, fourcc, device_id,
                                              anchor_exposure)
        except AutoExposureStuckError as second:
            cap.release()
            # from None: the chained message would be the first round printed twice, and the
            # message raised here already contains both rounds, side by side.
            raise AutoExposureStuckError(
                _stuck_message(device_id, (("first open", first.attempts),
                                           ("after close and reopen", second.attempts))),
                list(first.attempts) + list(second.attempts)) from None

    if evidence:
        evidence["reopened"] = reopened

    got_w, got_h, got_fourcc = negotiated(cap)
    description = (f"{got_w}x{got_h} {got_fourcc or '(codec not reported)'}, "
                   f"driver fps {cap.get(cv2.CAP_PROP_FPS):.4g}")
    if evidence:
        description += (f", manual exposure via AUTO_EXPOSURE={evidence['manualVia']} "
                        f"({evidence['responseRatio']:.2f}x)")
        if reopened:
            description += ", after one close and reopen"
    return cap, description, evidence


def identify_camera(device_id):
    """Work out WHICH preset this camera is, by asking it for each preset's native geometry.

    Returns (preset name, probes), where probes is the ask/got pair for every preset -- kept so
    session.json carries the evidence the identification was made on.

    WHY THE GEOMETRY AND NOT THE DEVICE NAME. Nothing Windows says about these devices can be
    trusted. The two Arducam models share a USB VID/PID, neither carries a serial number, and
    Windows caches the friendly name against the port, so moving a cable renames a camera. No PnP
    data is read here at all. What IS distinctive is the set of formats each sensor offers: every
    preset's native full frame is offered by that family and by no other -- 1920x1080 only on the
    ELP, 1600x1200 only on the OV2311, 1280x800 only on the OV9281. The bench project's
    findDevices.m identifies devices the same way.

    A geometry the camera does not have is NOT synthesised. DirectShow snaps the request to the
    nearest mode it does have and reports that back instead, so on an ELP asking for 1280x800
    reads back 1280x720 and asking for 1600x1200 reads back 1920x1080. Only the camera that
    really has the mode answers with exactly what was asked, which is what makes the readback a
    capability test rather than a guess.

    ONLY the geometry is set here -- no FOURCC, no FPS, no source properties. The ordering traps
    described in open_camera are exactly what this must not walk into; the camera is released
    again as soon as it is identified, and open_camera reopens it and does the full sequence.

    Costs about 5 s per camera, all of it in phase 1, before Spike2 is sampling.
    """
    cap, index = _open_live(device_id)
    try:
        probes, matches = [], []
        for name in sorted(CAMERA_PRESETS):
            _, width, height = parse_format(CAMERA_PRESETS[name]["format"])
            cap.set(cv2.CAP_PROP_FRAME_WIDTH, width)
            cap.set(cv2.CAP_PROP_FRAME_HEIGHT, height)
            got_w, got_h, _ = negotiated(cap)
            probes.append({"preset": name, "asked": [width, height], "got": [got_w, got_h]})
            if (got_w, got_h) == (width, height):
                matches.append(name)
    finally:
        # On every path: an identification that left the device claimed would make the open that
        # follows it fail for a reason this function invented.
        cap.release()

    # Built either way; three short strings, and both failure messages want it.
    asked = "\n".join(f"    {p['preset']:<8} asked {p['asked'][0]}x{p['asked'][1]}, "
                      f"got {p['got'][0]}x{p['got'][1]}" for p in probes)
    if not matches:
        raise RuntimeError(
            f"Could not identify the camera at OpenCV index {index} (device {device_id}): no "
            f"known preset's native geometry was delivered.\n{asked}\n"
            f"  Each preset's full frame is offered by one camera family and no other, so a "
            f"camera that delivers\n  none of them is not one this app has a preset for. Add "
            f"it to CAMERA_PRESETS, or pass --camera <preset>\n  to force a format on it "
            f"anyway.")
    if len(matches) > 1:
        raise RuntimeError(
            f"Camera at OpenCV index {index} (device {device_id}) answered to more than one "
            f"preset: {', '.join(matches)}.\n{asked}\n"
            f"  This cannot happen with the presets as shipped -- no two of them share a native "
            f"geometry -- so\n  either CAMERA_PRESETS has been edited or this camera really "
            f"offers both. Pass --camera <preset>\n  to say which one to use.")
    return matches[0], probes


def detect_preset(devices):
    """Identify every camera and insist they agree. Returns (preset name, evidence).

    The evidence goes to session.json as requested.cameraDetection, so a session can be checked
    afterwards for having identified the camera it says it recorded.

    This identifies the camera FAMILY, not the unit. Two cameras of the same model cannot be told
    apart by anything software can see here -- no serial number, and identical geometries by
    definition. Which physical camera is device 1 is still the operator's business, and is what
    the alignment stage is for.

    CameraBusyError propagates untouched: a device the previous run has not let go of is the one
    failure that clears itself, and the entry point retries identifying exactly as it retries
    opening.
    """
    evidence, found = [], []
    for device in devices:
        name, probes = identify_camera(device)
        evidence.append({"deviceId": device, "detected": name, "probes": probes})
        found.append(name)
        print(f"  device {device}: {name} ({CAMERA_PRESETS[name]['title']})")

    if len(set(found)) > 1:
        listing = "\n".join(f"    device {e['deviceId']}: {e['detected']}" for e in evidence)
        raise RuntimeError(
            f"The connected cameras are not the same model:\n{listing}\n"
            f"  One preset is applied to both, so a mixed pair cannot be recorded as it stands. "
            f"Either fix the\n  cabling so both devices are the intended pair, or pass "
            f"--camera <preset> to force one format\n  on both and accept what the other "
            f"camera makes of it.")
    return found[0], evidence


def verify_format(cap, width, height, fourcc, stage):
    """Fail loudly rather than quietly recording something else. Releases cap before raising."""
    got_w, got_h, got_fourcc = negotiated(cap)
    if got_w <= 0 or got_h <= 0:
        cap.release()
        raise CameraBusyError(f"Camera stopped answering property queries during {stage}: "
                              f"geometry reads back {got_w}x{got_h}.{_HELD_HINT}")
    if (got_w, got_h) != (width, height):
        cap.release()
        # Asking for the wrong preset is the easy mistake, and "it is at 1600x1200" is a lot more
        # useful when it also says which camera that is.
        other = next((n for n, p in CAMERA_PRESETS.items()
                      if p["format"].endswith(f"_{got_w}x{got_h}")), None)
        hint = ""
        if other:
            hint = (f"\n\n  {got_w}x{got_h} is the native format of the {other} preset; if that "
                    f"is the camera you\n  have connected, re-run with --camera {other}.")
        raise RuntimeError(f"Driver would not deliver {width}x{height}; it is at {got_w}x{got_h} "
                           f"after {stage}. Refusing to run: a cropped or rescaled frame is not "
                           f"comparable to a full-frame one, and the difference does not show up "
                           f"in the frame rate.{hint}")
    # An empty codec readback is not a mismatch -- some backends report nothing here even when
    # the media type is correct. Only a positively different codec is a failure.
    if fourcc and got_fourcc and got_fourcc.upper() != fourcc.upper():
        cap.release()
        raise RuntimeError(f"Asked for codec {fourcc} but the driver is on {got_fourcc} after "
                           f"{stage}. These have very different bandwidth costs -- uncompressed "
                           f"at full frame measures 5.19 fps against MJPG 49.4 over USB 2.0.")


def apply_source(cap, source):
    """Apply the preset's source properties.

    Auto modes go off BEFORE the manual values are written, or the driver overwrites them again
    on the next frame. Auto-exposure is the one that really matters: it lengthens integration as
    the scene darkens and the driver drops the frame rate to suit.

    FrameRate is absent on purpose -- it is set once during negotiation, in open_negotiated,
    and must never be re-applied. Its readback is worthless anyway; DirectShow echoes whatever it
    was handed, answering 1000.0 after a 1000 fps request on a 50 fps camera.
    """
    for name, prop in _MODE_MAP.items():
        # ExposureMode is NOT set here. It is the one mode whose effect can be measured, and
        # measuring is the only thing that works: see secure_manual_exposure, which owns it.
        if name in source and name != "ExposureMode":
            manual = str(source[name]).lower() == "manual"
            # 0.25/0.75 is the widely-honoured UVC convention and 0/1 is the other; both are
            # tried because neither is reliably right across backends. cap.set returning True
            # is NOT evidence the driver honoured it -- for white balance and focus there is
            # nothing cheap to measure, so this stays a best guess.
            for value in ((0.25, 0.0) if manual else (0.75, 1.0)):
                if cap.set(prop, value):
                    break

    for name, prop in _PROP_MAP.items():
        if name in source:
            raw = source[name]
            value = _ONOFF.get(str(raw).lower()) if isinstance(raw, str) else raw
            if value is None:
                print(f"  ignoring {name}={raw!r}: not a number or on/off", file=sys.stderr)
            else:
                cap.set(prop, float(value))


def set_exposure(cap, value):
    """Write one exposure value. Used for the recording anchor -- see record.py."""
    return bool(cap.set(cv2.CAP_PROP_EXPOSURE, float(value)))


def _mean_after(cap, settle_s=_PROBE_SETTLE_S, frames=_PROBE_FRAMES):
    """Mean intensity once the last control write has had time to take effect.

    The settle is not optional: a UVC control change lands several frames later, so measuring
    immediately would report the PREVIOUS setting and every verdict built on it would be wrong.
    """
    deadline = time.perf_counter() + settle_s
    while time.perf_counter() < deadline:
        cap.read()
    total, got = 0.0, 0
    for _ in range(frames):
        ok, frame = cap.read()
        if ok and frame is not None:
            plane = frame[:, :, 0] if frame.ndim == 3 else frame
            total += float(plane[::4, ::4].mean())      # 1/16 the pixels, same statistics
            got += 1
    return total / got if got else float("nan")


def _exposure_response(cap):
    """Drive the exposure across its full range and report how far the image actually moved.

    Dark first, bright second, so the sequence always ends on the bright endpoint -- see
    MIN_RESTORE_DROP for why that ordering is load-bearing.
    """
    set_exposure(cap, PROBE_DARK)
    at_dark = _mean_after(cap)
    set_exposure(cap, PROBE_BRIGHT)
    at_bright = _mean_after(cap)
    ratio = at_bright / at_dark if at_dark > 0 else float("nan")
    return ratio, at_dark, at_bright


def _write_count(attempts):
    """How many AUTO_EXPOSURE writes those attempts actually made.

    The first attempt measures what apply_source left behind and writes nothing; a pair attempt
    writes twice. Only used for the one-line note before a reopen, but a count that matches what
    the code did is what makes that note worth printing.
    """
    return sum(2 if " then " in str(a["autoExposure"]) else 1
               for a in attempts if a["autoExposure"] != "as applied")


def _attempt_lines(attempts):
    """One line per attempt: what was written, and what the image did about it."""
    return "\n".join(
        f"    AUTO_EXPOSURE={a['autoExposure']}: mean {a['meanDark']:.1f} at exposure "
        f"{PROBE_DARK:g} -> {a['meanBright']:.1f} at {PROBE_BRIGHT:g}  ({a['ratio']:.2f}x)"
        for a in attempts)


def _stuck_message(device_id, rounds):
    """The refusal: every attempt's measurements, how to read them, and what to do about it.

    `rounds` is (label, attempts) pairs -- one round for a plain failure, two once the camera has
    also been closed and reopened. EVERY attempt is printed, not just the last ratio, because the
    ratio on its own does not say which failure this is: a camera pinned in its own auto-exposure
    and a camera looking at a lens cap both measure about 1.00x, and only the means separate them.
    """
    blocks = []
    for label, attempts in rounds:
        head = f"  {label}:\n" if label else ""
        blocks.append(head + _attempt_lines(attempts))
    return (
        f"Device {device_id}: exposure is not under manual control. Driving it across its whole "
        f"range, {PROBE_DARK:g} to {PROBE_BRIGHT:g}, never moved the image by the "
        f"{MIN_RESPONSE_RATIO:g}x that manual control gives -- on any setting tried:\n"
        + "\n".join(blocks) + "\n"
        + _MEANS_GUIDE
        + "\n\n  If this is auto-exposure, the preset Exposure value does nothing and the "
          "recording anchor cannot\n  mark anything. Refusing to run either way: the frames "
          "would look plausible and be untrustworthy."
        + _AUTO_STUCK_HINT)


def secure_manual_exposure(cap, source, width, height, fourcc, device_id, anchor_exposure=None):
    """Prove this camera is in MANUAL exposure, and that the anchor cannot break the format.

    Returns evidence for session.json, or None if the preset pins no exposure.

    WHY THIS IS MEASURED RATHER THAN SET. CAP_PROP_AUTO_EXPOSURE reads back -1 on these cameras,
    so the mode cannot be queried, and cap.set() returns True for values the driver ignores.
    Measured on both OV2311s (2026-08-18, tools/probe_exposure_mode.py in the bench project):
    writing 0.25 -- the value the UVC convention calls manual, and the value this app used to
    write -- leaves the camera in AUTO. So does 0.0. Only 0.75 and 1.0 hand over control, which
    is inverted from the usual reading of the standard.

    The consequence of getting it wrong is not obvious in the data: auto-exposure holds the image
    near its own target, so frames look reasonable while the exposure setting does nothing, the
    two cameras drift apart as each chases its own scene, and the recording anchor writes an
    exposure change that never happens. A whole session looks fine and is not.

    Nothing here trusts a convention. The exposure is driven across its full usable range and the
    image is checked for having followed -- the full range, not a step around the operating point,
    because near the dark end a black-level pedestal compresses the response to 1.2x, which is
    indistinguishable from a camera ignoring the write (see PROBE_DARK). Candidates are ordered by
    what works on this rig, but every one is confirmed by measurement, so a camera that does
    follow the standard is handled without a special case.

    Raises AutoExposureStuckError, carrying every attempt's measurements, when nothing written
    here moves the image. open_camera catches that one and tries a closed-and-reopened device
    before the refusal reaches the operator.

    All of this happens in phase 1, before the ready flag exists and before Spike2 is sampling,
    so the second or two it costs is free.
    """
    normal = source.get("Exposure")
    if normal is None:
        return None
    normal = float(normal)

    attempts = []
    ratio, at_dark, at_bright = _exposure_response(cap)
    attempts.append({"autoExposure": "as applied", "ratio": ratio,
                     "meanDark": at_dark, "meanBright": at_bright})
    winner = "as applied" if ratio >= MIN_RESPONSE_RATIO else None

    if winner is None:
        for candidate in MANUAL_AE_CANDIDATES:
            cap.set(cv2.CAP_PROP_AUTO_EXPOSURE, candidate)
            ratio, at_dark, at_bright = _exposure_response(cap)
            attempts.append({"autoExposure": candidate, "ratio": ratio,
                             "meanDark": at_dark, "meanBright": at_bright})
            if ratio >= MIN_RESPONSE_RATIO:
                winner = candidate
                break

    # The last thing left to try before giving up: two values in sequence rather than one. See
    # MANUAL_AE_PAIRS -- the backend has only a manual flag and an auto flag to hand the driver,
    # so the TRANSITION between them is the one thing the single candidates cannot say.
    if winner is None:
        for first, second in MANUAL_AE_PAIRS:
            cap.set(cv2.CAP_PROP_AUTO_EXPOSURE, first)
            cap.set(cv2.CAP_PROP_AUTO_EXPOSURE, second)
            ratio, at_dark, at_bright = _exposure_response(cap)
            label = f"{first} then {second}"
            attempts.append({"autoExposure": label, "ratio": ratio,
                             "meanDark": at_dark, "meanBright": at_bright})
            if ratio >= MIN_RESPONSE_RATIO:
                winner = label
                break

    if winner is None:
        cap.release()
        raise AutoExposureStuckError(_stuck_message(device_id, ((None, attempts),)), attempts)

    # Restore, then CHECK the restore -- an unverified restore is how a camera ends up running a
    # whole session at the probe value. The probe ended bright, so a failed restore is a big,
    # obvious difference rather than a subtle one.
    set_exposure(cap, normal)
    restored = _mean_after(cap)
    if normal <= PROBE_BRIGHT - RESTORE_CHECK_MARGIN and restored > 0:
        if at_bright / restored < MIN_RESTORE_DROP:
            cap.release()
            raise RuntimeError(
                f"Device {device_id}: the exposure did not go back to {normal:g} after the "
                f"manual-control check -- the image is still as bright as it was at "
                f"{PROBE_BRIGHT:g} (mean {restored:.1f} against {at_bright:.1f}). Refusing to "
                f"run rather than record at the wrong exposure.")

    # The anchor's own value is the one that must not re-negotiate the media type: several
    # property writes on this backend revert MJPG to YUY2, which would drop the recording from
    # 49 fps to 5 without raising anything.
    anchor_ok = None
    if anchor_exposure is not None:
        set_exposure(cap, anchor_exposure)
        got = negotiated(cap)
        set_exposure(cap, normal)
        reverted = (got[0], got[1]) != (width, height)
        if fourcc and got[2] and got[2].upper() != fourcc.upper():
            reverted = True
        if reverted:
            cap.release()
            raise RuntimeError(
                f"Device {device_id}: changing the exposure re-negotiated the media type "
                f"({width}x{height} {fourcc} -> {got[0]}x{got[1]} {got[2] or '?'}). The recording "
                f"anchor cannot be used on this camera. Set anchor to false in config.json and "
                f"fall back to timestamp matching for the strobe alignment.")
        verify_format(cap, width, height, fourcc, "the exposure-anchor probe")
        anchor_ok = True

    return {"manualVia": winner, "responseRatio": ratio,
            "probeDark": PROBE_DARK, "probeBright": PROBE_BRIGHT,
            "normalExposure": normal, "meanAtDark": at_dark, "meanAtBright": at_bright,
            "meanAfterRestore": restored, "anchorFormatSafe": anchor_ok,
            "minRatioRequired": MIN_RESPONSE_RATIO, "attempts": attempts}
