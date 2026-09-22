r"""Camera check: are both rig cameras there, what state are they in, and take a picture.

    python tools\camera_check.py
    python tools\camera_check.py --devices 1 2 --camera ov9281
    python tools\camera_check.py --quit-after 15 --snapshot "C:/Temp/test/check"
    python "eye_sensor_calibration/recording_app/tools/camera_check.py"

USE THIS INSTEAD OF THE WINDOWS CAMERA APP. That is the whole reason it exists. The Camera app
-- and anything else that goes through Windows' own camera pipeline, the Frame Server -- can
leave a camera in an auto-exposure mode that the recording app's property writes cannot undo,
and the state stays in the camera after the app is closed. See cameras._AUTO_STUCK_HINT: it has
appeared a minute or two AFTER the Camera app was shut, so the cause looks long finished by the
time a session refuses to run. Nothing here goes through the Frame Server, and nothing here
writes a single camera control unless you ask for it with `m`.

What it shows, per camera, live:

    which device it is, and the friendly name Windows has for it (not to be trusted for
    identification -- see cameras.identify_camera -- but useful for saying which one you mean)
    the negotiated media type, the driver's frame rate, and the rate frames really arrive at
    every control the driver exposes: value, range, step, default, and whether it is on AUTO

The AUTO readout is the part the app cannot give you. It comes from DirectShow directly, through
eyecal/dshow.py, and it is instant -- against the couple of seconds cameras.secure_manual_exposure
spends driving the exposure and measuring the image. It is a HINT, not the proof: `m` runs the
app's real measurement on demand, and its verdict is the one that decides whether a session runs.

A tool, not a test, for the same reason session_video.py is one: it looks at the rig rather than
at the code, it needs a camera and an operator, and a test_ prefix would make pytest collect it.
Its pure parts -- the overlay lines, the raw tiling, the snapshot names -- are covered in
tests/test_offline.py under "camera_check".

WHAT IT DOES NOT DO. It does not force the exposure and it does not prove manual control at
start-up, both of which the app does on every open. The cameras are opened by the app's own
cameras.open_negotiated, so the media type is negotiated in exactly the load-bearing order the
app uses, but the source dict handed to it holds ONLY FrameRate. apply_source therefore finds
nothing to write, and each camera keeps whatever state it was already in -- which is the state
this tool exists to show. Press `m` to change that deliberately.

Keys:  Esc quit | p snapshot | m prove manual exposure | o reopen failed | c colour | r rotate
       | s swap sides | h hide overlay

A camera that will not open keeps its pane and shows the failure in red; the other camera goes
on running. `o` retries them. Camera detection (`camera: auto`) identifies EACH device on its
own, so a MIXED PAIR is fine here -- an OV9281 on device 1 and an OV2311 on device 2 each get
their own preset, their own format and their own pane size. The app cannot do that: it applies
one preset to both cameras, so cameras.detect_preset makes both devices agree and refuses a mixed
pair. A device that cannot be identified at all keeps its pane and says so, but there is no
preset to open it with and `o` cannot retry it; pass `--camera <preset>` for that.
"""

import argparse
import sys
import time
from datetime import datetime
from pathlib import Path

import cv2
import numpy as np

# eyecal/ is the recording app itself, one level up. Everything about opening a camera and
# drawing on a frame is used from there, never copied, so this tool and the app can never
# disagree about what "opened correctly" means.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from eyecal import cameras, capture, config, display, dshow          # noqa: E402

KEYS_HELP = ("Esc quit | p snapshot | m prove manual exposure | o reopen failed | c colour | "
             "r rotate | s swap sides | h hide overlay")

# Same retry as the app's entry point, for the same narrow reason: a device a previous run has
# not been released from is the one opening failure that clears itself, given a few seconds.
OPEN_ATTEMPTS = 3
OPEN_RETRY_WAIT_S = 5.0

# The driver properties are re-read on this period rather than every frame. They change only
# when something writes them, a COM round trip per property is not free, and a readout that
# flickered at 25 Hz would be harder to read, not fresher.
PROPERTY_PERIOD_S = 1.0

OVERLAY_SCALE = 0.4                 # measured: the longest line is 475 px in a 640 px pane
OVERLAY_TOP = 16                    # y of the first baseline
OVERLAY_STEP = 15                   # px between baselines; the glyphs are 11 px tall
# The image is dimmed to this fraction behind the text. Measured on the first snapshot: green
# text over a white test pattern was unreadable, and this tool's readout is the thing being
# looked at. The picture stays visible through it, and `h` removes both.
OVERLAY_DIM = 0.45

# The measured frame rate is the reader's frame count differenced over this window. One second
# at 120 fps is 120 frames, which is plenty for a two-digit readout that updates visibly.
RATE_WINDOW_S = 1.0

# --snapshot waits this long after every camera's first frame before writing. The measured rate
# on the overlay reads nan until the first window has elapsed, and a snapshot taken for the
# record should carry the real number.
SNAPSHOT_SETTLE_S = 3.0


# The format a device that could not be identified falls back to, and only when no other device
# was identified either. Nothing is ever opened with it -- a pane with no preset cannot be opened
# at all -- so it only has to give the black pane a plausible shape.
FALLBACK_FORMAT = "MJPG_1280x800"


# --- the pure parts: overlay text, raw tiling, snapshot names -------------------------------

def property_line(row):
    """One driver property as one line: what it is, what it is set to, and which mode it is in.

    `row` is one dict from dshow.Controls.read(). A value of None means GetRange succeeded but
    Get did not, which is rare and worth showing as "?" rather than as a number that was never
    read.
    """
    value = "?" if row["value"] is None else row["value"]
    mode = {True: "AUTO", False: "MANUAL"}.get(row["auto"], "mode unknown")
    return (f"{row['name']:<22}{value!s:>6}   "
            f"[{row['min']}..{row['max']} step {row['step']}, default {row['default']}]   {mode}")


def overlay_lines(pane, device, friendly, preset_name, preset_title, media_type, driver_fps,
                  measured_fps, rows, controls_error=None, note=""):
    """The whole readout for one working camera, as (line, is_auto) pairs.

    is_auto is what the caller colours by, so the decision is made once, here, where the property
    row is still in hand -- the drawing code never has to parse a string back apart.

    Properties the driver does not expose are already absent from `rows`; dshow.Controls.read
    drops anything whose GetRange fails, because a line reading "Pan 0" on a camera with no pan
    would be a lie with a number on it.

    BOTH RATES ARE SHOWN because they are not the same number and neither replaces the other.
    `driver_fps` is what the driver ANSWERS, and on these cameras that is the sentinel the
    negotiation asked for -- DirectShow echoes whatever it was handed, so a "max" request reads
    back as 240 (see cameras.apply_source). `measured_fps` is the rate frames really arrive at,
    and is the one to believe.
    """
    lines = [(f"cam {pane}  device {device}   {friendly or '(no name)'}", False),
             (f"preset {preset_name} ({preset_title})   {media_type}   "
              f"driver {driver_fps:.4g} fps   measured {measured_fps:.1f} fps", False)]
    if controls_error:
        lines.append((f"driver properties: unavailable ({controls_error}) -- press m to prove "
                      f"exposure by measurement", True))
    else:
        for row in rows:
            # The mode word is already in the line; the flag is what makes it RED, which is what
            # carries across a room. Repeating "<-- AUTO" after a column that says AUTO would
            # only look like a bug.
            lines.append((property_line(row), bool(row["auto"])))
    if note:
        lines.append((note, False))
    return lines


def failed_lines(pane, device, message, max_lines=4):
    """What a pane shows when its camera would not open: the device, and why.

    The first few lines only. These messages carry a whole escalation list (see
    cameras._HELD_HINT) which belongs on the console, where it can be read and scrolled; the
    pane says which device failed and what the first line of the reason was.
    """
    body = [ln.rstrip() for ln in str(message).splitlines() if ln.strip()][:max_lines]
    return [(f"cam {pane}  device {device}   NOT OPEN", True)] + [(ln, True) for ln in body]


def raw_tile(frames, rotate180=False, fallback_shape=None):
    """The full-resolution frames, side by side, with NOTHING drawn on them.

    `frames` is in display order, left to right, and a None is a camera that is not delivering.
    A missing camera contributes a BLACK pane the size of a frame that did arrive, so the two
    halves of the file stay where the operator saw them on screen.

    Each frame is rotated on its own, exactly as the preview rotates each pane on its own.
    Rotating the finished canvas instead would also swap the panes left to right, which is not
    what the screen showed.

    Returns None when there is nothing at all to write and no fallback shape was given.
    """
    shapes = [f.shape for f in frames if f is not None] or ([fallback_shape]
                                                            if fallback_shape else [])
    if not shapes:
        return None
    height = max(s[0] for s in shapes)
    colour = any(len(s) == 3 for s in shapes)

    panes = []
    for f in frames:
        if f is None:
            panes.append(np.zeros(shapes[0][:2] + ((3,) if colour else ()), np.uint8))
            continue
        f = cv2.rotate(f, cv2.ROTATE_180) if rotate180 else f
        if colour and f.ndim == 2:
            f = cv2.cvtColor(f, cv2.COLOR_GRAY2BGR)
        panes.append(f)

    width = sum(p.shape[1] for p in panes)
    canvas = np.zeros((height, width) + ((3,) if colour else ()), np.uint8)
    x = 0
    for p in panes:
        canvas[:p.shape[0], x:x + p.shape[1]] = p
        x += p.shape[1]
    return canvas


def snapshot_paths(chosen):
    """The two files one snapshot writes, from whatever path the dialog handed back.

    Any extension is stripped first, so accepting the dialog's own ".png" does not produce
    "check.png_raw.png". The two files are the canvas as shown and the frames as captured, and
    they are deliberately named so they sort next to each other.
    """
    base = Path(chosen)
    if base.suffix:
        base = base.with_suffix("")
    return (base.with_name(base.name + "_overlay.png"),
            base.with_name(base.name + "_raw.png"))


def default_snapshot_name(now=None):
    """camera_check_<YYYYmmdd-HHMMSS>, the name the save dialog opens with."""
    return (now or datetime.now()).strftime("camera_check_%Y%m%d-%H%M%S")


# --- drawing --------------------------------------------------------------------------------

def _bgr(img):
    """A drawable 3-channel copy. These sensors are monochrome, so this is usually a convert."""
    return cv2.cvtColor(img, cv2.COLOR_GRAY2BGR) if img.ndim == 2 else img.copy()


def dim_band(canvas, top, bottom):
    """Darken the rows top..bottom in place, so text drawn there reads against any picture."""
    top, bottom = max(int(top), 0), min(int(bottom), canvas.shape[0])
    if bottom > top:
        band = canvas[top:bottom]
        cv2.addWeighted(band, OVERLAY_DIM, np.zeros_like(band), 0.0, 0, dst=band)


def draw_pane(img, lines, colour, show_overlay):
    """One pane: the display copy, the readout, and a border. No crosshair.

    This is not the alignment view. A crosshair here would invite someone to aim a camera by a
    reticle that this tool makes no promise about, and aiming is what align.py is for.

    A line marked auto is drawn in RED whatever the current overlay colour is: the operator is
    looking for exactly that, and a green "AUTO" among green everything is a line to be missed.
    """
    canvas = _bgr(img)
    if show_overlay:
        dim_band(canvas, 0, OVERLAY_TOP + OVERLAY_STEP * len(lines) - 8)
        y = OVERLAY_TOP
        for text, is_auto in lines:
            display.text(canvas, text, (8, y), display.RED if is_auto else colour,
                         OVERLAY_SCALE)
            y += OVERLAY_STEP
    # The border stays regardless of the overlay toggle: tiled edge to edge, it is the only
    # thing showing where one camera's frame ends and the next begins.
    cv2.rectangle(canvas, (0, 0), (canvas.shape[1] - 1, canvas.shape[0] - 1), colour, 1)
    return canvas


# --- one camera -----------------------------------------------------------------------------

class Pane:
    """Everything about ONE device: its capture, its reader, its readout and its last failure.

    A pane exists for every requested device whether or not the camera opened, because a missing
    camera is the single most useful thing this tool can show. `error` set means cap and reader
    are both None and the pane draws red.

    The preset is PER PANE, not per run: this tool identifies each device on its own and is happy
    with a mixed pair, so the two panes may be different cameras with different formats and
    different frame sizes. Everything downstream -- the open, the overlay, `m`, `o` and the size
    of the black pane a silent camera leaves behind -- reads it from here.
    """

    def __init__(self, device):
        self.device = device
        self.preset_name = ""               # which preset this device is; see assign_preset
        self.preset = None                  # the CAMERA_PRESETS entry, or None if unidentified
        self.fmt = ""                       # the format this pane is actually opened with
        self.cap = None
        self.reader = None
        self.error = None                   # str: why this camera is not open
        self.latest = None                  # newest FULL-RESOLUTION frame, for the raw snapshot
        self.fps = float("nan")             # frames the reader really read per second
        self._rate_mark = None              # (preview time, reader.frames_read) at the last rate
        self.friendly = ""
        self.media_type = ""
        self.driver_fps = float("nan")
        self.controls = None                # dshow.Controls, or None
        self.controls_error = ""
        self.rows = []                      # last property readout
        self.next_property_read = 0.0
        self.note = ""                      # what the last `m` proved, until the next one

    @property
    def open(self):
        return self.cap is not None

    def reset_rate(self):
        self.fps, self._rate_mark = float("nan"), None

    def update_rate(self, t, window_s=RATE_WINDOW_S):
        """The rate from the reader's OWN frame counter, not from what the preview drained.

        The preview drains to the newest frame once per loop, so counting drains -- which is
        what a RateMeter ticked from drain_newest measures -- reports the loop rate and reads the
        same number for every camera. capture.CameraReader counts every successful read in
        frames_read; the difference over a window is the rate frames really arrive at.
        """
        if self.reader is None:
            return
        n = self.reader.frames_read
        if self._rate_mark is None:
            self._rate_mark = (t, n)
        elif t - self._rate_mark[0] >= window_s:
            t0, n0 = self._rate_mark
            self.fps = (n - n0) / (t - t0)
            self._rate_mark = (t, n)


def assign_preset(pane, preset_name, fmt_override=""):
    """Give one pane the preset it will be opened with, and the format that follows from it.

    Per pane rather than once for the whole run, because the two devices need not be the same
    camera in this tool. A --format override still applies to EVERY pane, as it always has: it is
    a deliberate "open them all like this", and the only way to ask a camera for a geometry that
    is not its own.
    """
    pane.preset_name = preset_name
    pane.preset = cameras.CAMERA_PRESETS[preset_name]
    pane.fmt = fmt_override or pane.preset["format"]


def identify_pane(pane, fmt_override=""):
    """Identify ONE camera from the geometries it offers, and give the pane that preset.

    Each device on its own, and a MIXED PAIR IS ALLOWED here. cameras.detect_preset -- what the
    app uses -- identifies both devices and then insists they agree, because the app applies one
    preset to both cameras and cannot record a pair that is not one. This tool records nothing,
    so two different cameras are a thing to report rather than a thing to refuse.

    A device that cannot be identified is marked failed and keeps its pane, for the same reason a
    device that will not open does: which camera is missing is half of what this tool is run for.
    It cannot be opened afterwards, since there is no preset to open it with; --camera <preset>
    is what forces one on it. A busy device propagates after the retries, exactly as it does for
    the app's entry point -- the retry is the whole remedy for that one.

    Sets pane.error and returns False on failure.
    """
    try:
        preset_name, _ = _retry_busy(lambda: cameras.identify_camera(pane.device))
    except cameras.CameraBusyError:
        raise
    except Exception as exc:
        pane.error = str(exc)
        print(f"\n  Device {pane.device} could not be identified:\n{exc}\n\n  Carrying on "
              f"without it; pass --camera <preset> to open it anyway.", file=sys.stderr)
        return False
    assign_preset(pane, preset_name, fmt_override)
    return True


def fill_unidentified(panes, fallback=FALLBACK_FORMAT):
    """Give every unidentified pane a format, purely so its black pane is a sensible size.

    Nothing is opened with it -- a pane with no preset is never opened -- but run() sizes the
    black pane from pane.fmt, and that needs a number. The first identified pane's format is the
    best guess available, and is right whenever the pair really is a pair.
    """
    known = next((p.fmt for p in panes if p.preset is not None), fallback)
    for pane in panes:
        if pane.preset is None:
            pane.fmt = known


def open_pane(pane, fmt, source, want_gray):
    """Open one camera the app's way, WITHOUT writing a single control, and start reading it.

    The source dict is cut down to FrameRate alone on purpose. cameras.open_negotiated applies
    whatever it is given through cameras.apply_source, which for a real preset writes Exposure,
    Gain, the auto modes and the rest -- that would change the very state this tool exists to
    report. FrameRate is kept because it is set during negotiation, not by apply_source, and
    leaving it out would negotiate a different media type from the one the app gets (see
    cameras.open_negotiated: unset is NOT the same as maximum).

    Retries a busy device exactly as the app's entry point does, and for the same reason: that
    is the one failure that clears itself within a few seconds. Every other failure is a
    statement about the configuration and is reported at once.

    Sets pane.error and returns False on failure, so one dead camera does not stop the tool.
    """
    rate = {"FrameRate": source.get("FrameRate", "max")}
    for attempt in range(1, OPEN_ATTEMPTS + 1):
        try:
            pane.cap = cameras.open_negotiated(pane.device, fmt, rate)
            got_w, got_h, got_fourcc = cameras.negotiated(pane.cap)
            pane.media_type = f"{got_w}x{got_h} {got_fourcc or '(codec not reported)'}"
            pane.driver_fps = pane.cap.get(cv2.CAP_PROP_FPS)
            pane.error = None
            pane.reader = capture.CameraReader(pane.cap, pane.device, want_gray)
            pane.reader.start()
            print(f"  device {pane.device}: {pane.media_type}, driver fps "
                  f"{pane.driver_fps:.4g}  (nothing written to this camera)")
            return True
        except cameras.CameraBusyError as exc:
            if attempt == OPEN_ATTEMPTS:
                pane.error = str(exc)
                print(f"\n  Device {pane.device} still claimed after {OPEN_ATTEMPTS} attempts. "
                      f"Giving up on it; press o to retry.\n{exc}", file=sys.stderr)
                return False
            print(f"\n  Device {pane.device}, attempt {attempt} of {OPEN_ATTEMPTS}:\n{exc}",
                  file=sys.stderr)
            print(f"\n  Waiting {OPEN_RETRY_WAIT_S:g} s for Windows to reclaim the device, then "
                  f"retrying.", file=sys.stderr)
            time.sleep(OPEN_RETRY_WAIT_S)
        except Exception as exc:
            pane.error = str(exc)
            print(f"\n  Device {pane.device} would not open:\n{exc}", file=sys.stderr)
            return False


def open_controls(pane):
    """Bind this device's DirectShow filter for the property readout.

    Measured on this rig: this coexists with an OpenCV capture on the same device in either
    order, so it is done after the open and no filter has to be held across it (see
    eyecal/dshow.py). A failure here costs the readout and nothing else -- the preview, the
    snapshot and `m` all still work -- so it is reported on the pane and never raised.
    """
    try:
        pane.controls = dshow.Controls(pane.device - 1)
        pane.friendly = pane.controls.friendly_name or ""
        pane.controls_error = ""
    except Exception as exc:
        pane.controls, pane.controls_error = None, f"{type(exc).__name__}: {exc}"
        print(f"  device {pane.device}: driver properties unavailable ({pane.controls_error})",
              file=sys.stderr)


def stop_reader(pane, timeout=5.0):
    """Stop and join this pane's reader so the MAIN thread can have the capture to itself.

    Nothing may call cap.read() on a capture another thread is also reading -- OpenCV makes no
    such promise, and capture.close_all documents what tearing a graph down under an active read
    does to the device. So `m` joins first and only then measures.

    Returns True once the thread is really gone. False means it is still inside a read that has
    not returned, and the caller must leave that capture alone.
    """
    if pane.reader is None:
        return True
    pane.reader.stop()
    pane.reader.join(timeout=timeout)
    if pane.reader.is_alive():
        print(f"  device {pane.device}: its reader is still inside a read after {timeout:g} s; "
              f"leaving the camera alone.", file=sys.stderr)
        return False
    pane.reader = None
    return True


def prove_manual(pane, source, want_gray):
    """Apply the preset and PROVE manual exposure by measurement, on one already-open camera.

    This is the app's own check, cameras.secure_manual_exposure, run on demand -- the one thing
    in this tool that writes to a camera. It exists because the AUTO flag from DirectShow says
    what the driver was last told, and only the measurement says what the sensor actually does.

    The reader is stopped first and started again afterwards, so the capture has exactly one
    thread reading it throughout.

    On failure the capture has ALREADY been released by secure_manual_exposure, on every path it
    raises from, so the pane is marked failed and its cap dropped; `o` reopens it. Returns the
    evidence dict, or None.
    """
    if not pane.open:
        return None
    if not stop_reader(pane):
        return None

    width, height, fourcc = cameras.negotiated(pane.cap)
    try:
        cameras.apply_source(pane.cap, source)
        evidence = cameras.secure_manual_exposure(pane.cap, source, width, height, fourcc,
                                                  pane.device)
    except RuntimeError as exc:
        # AutoExposureStuckError and the failed-restore refusal both land here, and both have
        # released the capture themselves. Everything the operator needs is in the message.
        print(f"\n{exc}\n", file=sys.stderr)
        pane.cap, pane.error = None, str(exc)
        pane.note = ""
        return None

    if evidence is None:
        pane.note = "preset pins no exposure -- nothing to prove"
        print(f"  device {pane.device}: the preset pins no exposure; nothing measured.")
    else:
        pane.note = (f"manual proven via AUTO_EXPOSURE={evidence['manualVia']} "
                     f"({evidence['responseRatio']:.2f}x)")
        print(f"  device {pane.device}: manual exposure via AUTO_EXPOSURE="
              f"{evidence['manualVia']}, response {evidence['responseRatio']:.2f}x "
              f"(mean {evidence['meanAtDark']:.1f} at {evidence['probeDark']:g} -> "
              f"{evidence['meanAtBright']:.1f} at {evidence['probeBright']:g}), restored to "
              f"{evidence['normalExposure']:g} at mean {evidence['meanAfterRestore']:.1f}")

    pane.reader = capture.CameraReader(pane.cap, pane.device, want_gray)
    pane.reader.start()
    pane.reset_rate()                       # a new reader counts from zero
    return evidence


def close_pane(pane):
    """Stop the reader, THEN release the capture. The join is what makes the release safe.

    Released even when the join timed out, exactly as capture.close_all does and for the same
    reason: a capture that is never released leaves the device claimed after this process exits,
    which is worse than the risk of tearing the graph down under a read. stop_reader has already
    said so on stderr. `m` makes the other choice, because it is not exiting.
    """
    stop_reader(pane)
    if pane.cap is not None:
        try:
            pane.cap.release()
        except Exception:
            pass
        pane.cap = None
    if pane.controls is not None:
        pane.controls.close()
        pane.controls = None


# --- the snapshot ---------------------------------------------------------------------------

def ask_snapshot_path(out_dir):
    """Where to save. A withdrawn Tk root, so no empty window is left behind.

    The preview freezes while this is up, which is fine -- the capture threads keep running and
    fill their queues, and capture.CameraReader BLOCKS on a full queue rather than discarding, so
    nothing is lost and nothing has to be told to wait.
    """
    import tkinter
    from tkinter import filedialog
    out_dir = Path(out_dir).expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)
    root = tkinter.Tk()
    root.withdraw()
    picked = filedialog.asksaveasfilename(
        title="Save camera check snapshot", initialdir=str(out_dir),
        initialfile=default_snapshot_name(), defaultextension=".png",
        filetypes=[("PNG image", "*.png")])
    root.destroy()
    return picked or None


def write_snapshot(chosen, canvas, frames, rotate180, fallback_shape):
    """Write the two PNGs: the canvas as shown, and the frames as captured.

    Two files because they answer different questions. The overlay one is evidence of what the
    operator was looking at, readout and all. The raw one is the pixels the sensor delivered, at
    full resolution with nothing drawn on them, so it can be measured rather than merely looked
    at -- which a green overlay burnt into it would ruin.

    Returns (written paths, failures), both lists, so the caller can say so on screen as well as
    on the console. cv2.imwrite returns False rather than raising for an unwritable path.
    """
    overlay_path, raw_path = snapshot_paths(chosen)
    raw = raw_tile(frames, rotate180, fallback_shape)
    written, failed = [], []
    for path, image in ((overlay_path, canvas), (raw_path, raw)):
        if image is None:
            failed.append(f"{path.name}: no frames to write")
            continue
        path.parent.mkdir(parents=True, exist_ok=True)
        if cv2.imwrite(str(path), image):
            written.append(path)
        else:
            failed.append(f"{path.name}: cv2.imwrite refused the path")
    for path in written:
        print(f"  wrote {path}")
    for message in failed:
        print(f"  SNAPSHOT FAILED -- {message}", file=sys.stderr)
    return written, failed


# --- the loop -------------------------------------------------------------------------------

def run(panes, cfg, snapshot_path=None, quit_after=None):
    """The live view. Mirrors align.run_alignment, minus the crosshair and the accept key.

    Every pane brings its OWN preset and format, because the two devices need not be the same
    camera here (see identify_pane). So the overlay, `m`, `o` and the size of the black pane a
    camera that is not delivering leaves behind are each taken from that pane, never from one
    preset shared between them. Frames of two different sizes are fine on screen: display.tile
    sizes the canvas to the tallest pane and leaves the rest black.

    Returns nothing: this tool has no verdict to hand back. Its output is what the operator saw
    and whatever snapshots were written.
    """
    downsample = max(int(cfg["downsample"]), 1)
    preview_hz = float(cfg["align_preview_hz"])
    # Once, not per redraw: a pane's format cannot change while the tool is running.
    full_shapes, pane_shapes = [], []
    for pane in panes:
        _, full_w, full_h = cameras.parse_format(pane.fmt)
        full_shapes.append((full_h, full_w))
        pane_shapes.append((full_h // downsample, full_w // downsample))

    colour_idx, show_overlay, rot = 0, True, bool(cfg["rotate180"])
    order = list(range(len(panes)))
    period = 1.0 / preview_hz if preview_hz > 0 else 0.0
    t0, next_at = time.perf_counter(), 0.0
    snapshot_pending = snapshot_path is not None
    frames_since = None                 # preview time at which every open camera had delivered
    # A refused write is shown on the panes as well as on the console, because the console is
    # behind the preview window and an operator who pressed p is looking at the preview.
    snapshot_error = ""

    window = display.Window("camera check -- records nothing", cfg["fullscreen"])
    print(f"\nKeys: {KEYS_HELP}")
    print("Nothing is written to a camera unless you press m.\n")

    try:
        while True:
            t = time.perf_counter() - t0
            for pane in panes:
                if pane.reader is None:
                    continue
                if pane.reader.error is not None:
                    # A capture thread that died takes its camera down with it; the pane says so
                    # and `o` is what tries again.
                    pane.error = f"capture thread failed: {pane.reader.error}"
                    print(f"\n  device {pane.device}: {pane.error}", file=sys.stderr)
                    close_pane(pane)
                    continue
                item = capture.drain_newest(pane.reader)
                if item is not None:
                    pane.latest = item[2]
                pane.update_rate(t)

            if t >= next_at:
                for pane in panes:
                    if pane.controls is not None and t >= pane.next_property_read:
                        pane.rows = _read_rows(pane)
                        pane.next_property_read = t + PROPERTY_PERIOD_S

                drawn = []
                for position, k in enumerate(order, start=1):
                    pane = panes[k]
                    img = (display.display_copy(pane.latest, downsample, rot)
                           if pane.latest is not None else np.zeros(pane_shapes[k], np.uint8))
                    if pane.open and pane.latest is not None:
                        lines = overlay_lines(
                            position, pane.device, pane.friendly, pane.preset_name,
                            pane.preset["title"], pane.media_type, pane.driver_fps, pane.fps,
                            pane.rows, pane.controls_error, pane.note)
                    elif pane.open:
                        lines = [(f"cam {position}  device {pane.device}   waiting for the "
                                  f"first frame...", False)]
                    else:
                        lines = failed_lines(position, pane.device, pane.error or "not open")
                    if snapshot_error:
                        lines.append((f"SNAPSHOT FAILED -- {snapshot_error}", True))
                    drawn.append(draw_pane(img, lines, display.COLOURS[colour_idx][1],
                                           show_overlay))
                canvas = display.tile(drawn)
                if show_overlay:
                    dim_band(canvas, canvas.shape[0] - 24, canvas.shape[0])
                    display.text(canvas, KEYS_HELP, (8, canvas.shape[0] - 12),
                                 display.COLOURS[colour_idx][1], 0.45)
                window.show(canvas)
                next_at = t + period

                # Once every open camera has delivered AND the rate windows have had time to
                # fill, so --snapshot never writes a black pane for a camera that was merely
                # still starting up, nor a "measured nan fps" for one that had just started.
                if snapshot_pending and all(p.latest is not None for p in panes if p.open):
                    if frames_since is None:
                        frames_since = t
                    elif t - frames_since >= SNAPSHOT_SETTLE_S:
                        snapshot_pending = False
                        snapshot_error = _snapshot(snapshot_path, panes, order, canvas, rot,
                                                   full_shapes[order[0]])

            key = window.pump()
            if key == 27:                                       # Esc
                print("\nQuit (Esc).")
                return
            elif key == ord("p"):
                chosen = ask_snapshot_path(cfg["out"])
                if chosen is None:
                    print("  snapshot cancelled")
                else:
                    snapshot_error = _snapshot(chosen, panes, order, canvas, rot,
                                               full_shapes[order[0]])
                t0, next_at = time.perf_counter() - t, 0.0      # the dialog is not preview time
            elif key == ord("m"):
                for pane in panes:
                    # A pane with no preset was never identified, so there is nothing to apply
                    # and nothing to prove; it is not open either.
                    if pane.preset is not None:
                        prove_manual(pane, pane.preset["source"], pane.preset["grayscale"])
                t0, next_at = time.perf_counter() - t, 0.0
            elif key == ord("o"):
                for pane in panes:
                    if not pane.open and pane.preset is not None:
                        pane.latest, pane.note = None, ""
                        pane.reset_rate()
                        if open_pane(pane, pane.fmt, pane.preset["source"],
                                     pane.preset["grayscale"]):
                            open_controls(pane)
                t0, next_at = time.perf_counter() - t, 0.0
            elif key == ord("c"):
                colour_idx = (colour_idx + 1) % len(display.COLOURS)
                print(f"  overlay colour: {display.COLOURS[colour_idx][0]}")
            elif key == ord("r"):
                rot = not rot
                print(f"  rotation: {'180 deg' if rot else 'none'}")
            elif key == ord("s"):
                order.reverse()
                print(f"  sides swapped -- left to right is now devices "
                      f"{[panes[i].device for i in order]}")
            elif key == ord("h"):
                show_overlay = not show_overlay
                print(f"  overlay: {'shown' if show_overlay else 'hidden'}")

            # AFTER the waitKey that delivers the close message and BEFORE the next redraw --
            # see display.Window.closed for why both orderings are load-bearing.
            if window.closed():
                print("\nQuit (window closed).")
                return
            if quit_after is not None and t >= quit_after:
                print(f"\nQuit ({quit_after:g} s elapsed, --quit-after).")
                return
    finally:
        window.close()


def _read_rows(pane):
    """The property readout, or the reason there is not one. Never raises.

    A camera unplugged while the tool is running fails here first, and losing the whole view to
    that would be the opposite of useful.
    """
    try:
        return pane.controls.read()
    except Exception as exc:
        pane.controls_error = f"{type(exc).__name__}: {exc}"
        pane.controls = None
        return []


def _snapshot(chosen, panes, order, canvas, rot, fallback_shape):
    """Write both files for the current view. Returns "" or what went wrong, for the panes."""
    frames = [panes[k].latest for k in order]
    _, failed = write_snapshot(chosen, canvas, frames, rot, fallback_shape)
    return "; ".join(failed)


def main(argv=None):
    """Identify each camera, open what can be opened, and show what state it is in.

    EACH DEVICE IS IDENTIFIED ON ITS OWN and a mixed pair is accepted: device 1 an OV9281 and
    device 2 an OV2311 get a preset, a format and a pane size each. The app refuses that pair --
    cameras.detect_preset makes both devices agree, because one preset is applied to both cameras
    there and two different geometries are not a pair to record from. Nothing is recorded here,
    so what is plugged in is a thing to report.

    --camera <preset> skips identification and applies that one preset to every device, which is
    also how a camera too broken to be identified is looked at. --format overrides the format on
    every pane, as it always has.
    """
    args = parse_args(argv)
    cfg = config.apply_cli(config.load(args.config, required=args.config is not None), args)

    if cfg["camera"] not in cameras.CAMERA_PRESETS and cfg["camera"] != cameras.AUTO:
        raise KeyError(f"Unknown camera {cfg['camera']!r}. "
                       f"Known: {', '.join(sorted(cameras.CAMERA_PRESETS))}, "
                       f"or {cameras.AUTO} to identify it")
    devices = [int(d) for d in cfg["devices"]]

    print("\n" + "=" * 88)
    print("CAMERA CHECK -- shows what state the cameras are in. Use this, not the Camera app.")
    print("=" * 88)
    for entry in _enumerate_devices():
        print(f"  index {entry['index']} (device {entry['index'] + 1}): "
              f"{entry['friendlyName']}")

    panes = [Pane(d) for d in devices]
    if cfg["camera"] == cameras.AUTO:
        print("Identifying each camera from the geometries it offers (auto)...")
        for pane in panes:
            identify_pane(pane, cfg["format"])
        fill_unidentified(panes)
    else:
        for pane in panes:
            assign_preset(pane, cfg["camera"], cfg["format"])
    for pane in panes:
        if pane.preset is None:
            print(f"  device {pane.device}: NOT IDENTIFIED (the reason is above)")
        else:
            print(f"  device {pane.device}: {pane.preset_name} ({pane.preset['title']})  "
                  f"{pane.fmt}")
    print(f"Devices     : {devices}")
    print("Opening cameras. NOTHING is written to any camera control.")

    try:
        for pane in panes:
            # Sequential, as the app opens them: OpenCV's DirectShow backend runs through a
            # shared global that is not documented as thread-safe. A pane with no preset is
            # skipped: it already carries the reason it was not identified, and there is no
            # format to ask that camera for.
            if pane.preset is None:
                continue
            if open_pane(pane, pane.fmt, pane.preset["source"], pane.preset["grayscale"]):
                open_controls(pane)
        run(panes, cfg, args.snapshot, args.quit_after)
    finally:
        for pane in panes:
            close_pane(pane)
    return 0


def _enumerate_devices():
    """Everything Windows lists as a video input, or nothing if the enumeration fails.

    Printed before anything is opened, because "device 2 does not exist" is the answer to half
    the questions this tool gets run for, and it is available without touching a camera.
    """
    try:
        return dshow.list_devices()
    except Exception as exc:
        print(f"  could not enumerate video devices ({exc})", file=sys.stderr)
        return []


def _retry_busy(fn):
    """Retry only a busy device, exactly as the app's entry point does and for its reasons."""
    for attempt in range(1, OPEN_ATTEMPTS + 1):
        try:
            return fn()
        except cameras.CameraBusyError as exc:
            if attempt == OPEN_ATTEMPTS:
                raise
            print(f"\n  Attempt {attempt} of {OPEN_ATTEMPTS} failed:\n{exc}", file=sys.stderr)
            print(f"\n  Waiting {OPEN_RETRY_WAIT_S:g} s for Windows to reclaim the device, then "
                  f"retrying.", file=sys.stderr)
            time.sleep(OPEN_RETRY_WAIT_S)


def parse_args(argv):
    """The app's own options for the keys this tool uses, plus two for running it unattended.

    Everything defaults to None so config.json wins unless the flag was actually given, which is
    the same rule config.apply_cli applies for the app.
    """
    p = argparse.ArgumentParser(
        description="Show what state the rig cameras are in, and save snapshots. Use this "
                    "instead of the Windows Camera app, which can leave a camera stuck in "
                    "auto-exposure.")
    p.add_argument("--config", default=None,
                   help="path to config.json (default: the app's, beside recording_app/)")
    p.add_argument("--camera", default=None,
                   choices=sorted(cameras.CAMERA_PRESETS) + [cameras.AUTO],
                   help="preset name, or auto to identify it from the geometries offered")
    p.add_argument("--devices", type=int, nargs="+", default=None,
                   help="1-based device ids; the ORDER is the initial left-to-right arrangement")
    p.add_argument("--format", default=None, help="e.g. MJPG_1280x720")
    p.add_argument("--downsample", type=int, default=None, help="preview pixel stride")
    p.add_argument("--align-preview-hz", type=float, default=None, help="preview redraw ceiling")
    p.add_argument("--rotate180", action="store_true", default=None)
    p.add_argument("--no-rotate180", dest="rotate180", action="store_false", default=None)
    p.add_argument("--fullscreen", action="store_true", default=None)

    p.add_argument("--snapshot", default=None,
                   help="write the two snapshot files to this base path, with no dialog, a few "
                        "seconds after every open camera has delivered a frame, then keep "
                        "running")
    p.add_argument("--quit-after", type=float, default=None,
                   help="exit after this many seconds (unattended runs and tests)")

    args = p.parse_args(argv)
    if args.config is None:
        beside = Path(__file__).resolve().parent.parent / "config.json"
        args.config = str(beside) if beside.exists() else None
    return args


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\nInterrupted.")
        sys.exit(1)
    except Exception as exc:
        print(f"\nFAILED: {exc}", file=sys.stderr)
        sys.exit(1)
