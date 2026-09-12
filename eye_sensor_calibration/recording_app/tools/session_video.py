r"""Session video: watch a recorded session play back with its own timestamps on it.

    python tools\session_video.py "C:\Temp\test\20260821-205753Z_ov2311_2cam_10s"
    python tools\session_video.py                     # asks for the folder with a dialog
    python eye_sensor_calibration\recording_app\tools\session_video.py

Renders one recorded session to calibration_video.mp4 in the session folder: both cameras
side by side at the top, and underneath them three growing histograms -- each camera's own
frame interval, and the moment-to-moment simultaneity between the two.

A tool, not a test, for the same reason strobe_timing.py is one: it looks at the app's
OUTPUT rather than at the app's code, and a test_ prefix would make pytest collect it.

What it is for. strobe_timing.py answers "is the timing good?" with numbers and two static
figures. This answers "does the recording LOOK right?" -- whether the two cameras really
show the same moment, whether either one stutters, and where in the run any oddity sits.
Some faults are obvious in motion and invisible in a summary statistic.

Three clocks, and session.json picks between them. A session recorded under the FSIN
trigger -- session.json trigger.mode == "trigger" -- is timed from the Spike2 TTL2 channel,
the pulse train that drove the cameras: stored frame k of BOTH c1.bin and c2.bin is pulse k,
so there is nothing to align and nothing to verify (see eyecal/trigger.py). Every other
session, free-run or one that asked for the trigger and fell back, keeps the older path: the
strobe train when its anchor verifies, and the host arrival times in ts.csv when it does not.
A hardware clock beats host arrivals by a wide margin, so ts.csv is only ever the fallback,
and the console always says which source was used and why. All three are converted to seconds
since that camera's own first frame before anything is drawn, so they look the same on screen
and none is quietly flattering.

Two things about the trigger path are stated rather than worked around. The WHOLE TTL2 channel
is used and its first spike is frame 0, so a test train recorded earlier in the same Spike2
file would read every frame against an earlier pulse, with nothing on screen to show it --
record the trigger train once, for the run. And both cameras take their times from the SAME
spike array, so the simultaneity panel reads a flat zero by construction: that panel documents
the shared clock there, it does not measure anything.

A pulse count that does not match the stored frame count never refuses the render. Frame k is
still shown against spike k from the START of the train, frames past the last pulse draw '--',
the difference is printed per camera, and the on-screen label says COUNT MISMATCH so a still
from the video cannot be mistaken for a clean run.

Three decisions worth stating up front, because each one is a deliberate refusal:

No OFF-gap correction on the strobe indices. strobe_timing.anchor_map measures, but does
not assume, what a camera does at the anchor's OFF transition: device 1 on this rig stores
the frame and emits no strobe, device 2 appears to drop both. The correction is therefore
camera- and run-dependent, and applying a guessed one would silently shift every frame
after the anchor. Frame k is shown against spike s0+k, raw. The visible consequence is
honest and small: the one doubled interval at the start anchor shows up as a single count
out at about 2P in the interval histograms. Better a datapoint you can see and explain
than a correction you cannot check.

Fixed histogram axes, set from the BULK of the data in a full pre-pass before encoding
starts. Fixed, because autoscaling per frame would make the axes crawl and the bars rescale
continuously, so a bar that grew and one that was merely re-normalised would look identical
-- the one thing a growing histogram exists to distinguish. From the bulk rather than the
extremes, because the first version set them from min and max, and the two ~40 ms anchor
intervals then stretched the interval axes twentyfold and put 99.7% of the data in one bar.
Outliers are not dropped for this: they are counted at whichever end of the axis they fell
off, and drawn there. See make_panel.

ffmpeg rather than cv2.VideoWriter. OpenCV's H.264 writer is not usable on this machine:
the "avc1" fourcc needs an OpenH264 DLL that is not installed, and rather than failing it
writes a valid-looking file with no frames in it. imageio-ffmpeg ships a static ffmpeg with
libx264 built in, so raw bgr24 frames are piped to that instead. It cannot fail silently:
the return code is checked and the encoder's log is kept.
"""

import argparse
import csv
import json
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import cv2
import imageio_ffmpeg

# strobe_timing.py sits beside this file and is a script rather than a package, so its
# directory goes on the path before it is imported. Everything to do with reading the
# .smrx and finding the anchor lives there and is used from there, not copied.
sys.path.insert(0, str(Path(__file__).resolve().parent))
import strobe_timing                                                    # noqa: E402

# eyecal/ is the recording app itself, one level up. read_frames is its own reader for the
# raw .bin files, so the video and the app can never disagree about the frame layout.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from eyecal.session import read_frames                                  # noqa: E402

DEFAULT_ROOT = Path("C:/Temp/test")     # where sessions land on this rig; only a dialog hint
HIST_BINS = 50
# The x-range comes from the BULK of the data, not its extremes: median +/- IQR_MULT x IQR,
# never narrower than +/- MIN_HALF_MS. Anything outside is counted in an overflow marker at
# that end of the axis instead of stretching the axis to reach it. See make_panel.
HIST_IQR_MULT = 4.0
HIST_MIN_HALF_MS = 0.25                 # 50 bins over 0.5 ms = one 10 us Spike2 tick per bin
HIST_OVERFLOW_ZONE = 34                 # px reserved at each axis end for the marker + count
CROSSHAIR = (0, 255, 0)                 # BGR green, blended to ~35% over the image

FONT = cv2.FONT_HERSHEY_SIMPLEX

# Every string drawn on the canvas is ASCII on purpose: the Hershey fonts OpenCV ships hold
# no glyphs above 127, so a plus-minus sign or a middle dot comes out as a hollow box.
STROBE_LABEL = "strobe (Spike2), +/-2 frames"
TSCSV_LABEL = "ts.csv arrival (no valid strobe)"
TRIGGER_LABEL = "trigger (Spike2 TTL2), frame k = pulse k"
TRIGGER_MISMATCH_LABEL = TRIGGER_LABEL + " -- COUNT MISMATCH, see console"

# Where the trigger train is recorded: Spike2's TTL2 input. The channel number comes from
# spike2_experimental_protocols/utils/sampling_window_config.s2s, which sets TTL2_ch% := 12
# with the title "TTL2" and the comment "TTL2: Camera Trigger".
TRIGGER_CHAN = 12


class StrobeUnusable(Exception):
    """A Spike2 clock -- strobe or trigger -- cannot time this session.

    Carries the reason, which gets printed. One exception covers both clocks because the
    handling is identical: say why on the console and fall back to ts.csv. The name is from
    the strobe, which was the only Spike2 clock when this was written.
    """


# --- small drawing helpers --------------------------------------------------------------

def bgr(hex_colour):
    """'#3987e5' -> (229, 135, 57), because OpenCV is blue-first."""
    h = hex_colour.lstrip("#")
    r, g, b = (int(h[i:i + 2], 16) for i in (0, 2, 4))
    return (b, g, r)


def text(img, s, org, scale=0.5, colour=(255, 255, 255), weight=1, anchor="left",
         outline=False):
    """Draw a string and return its pixel width, so the caller can place things beside it.

    outline=True strokes the glyphs in black first. That is only needed over the camera
    images, where the background is whatever the camera happened to be looking at and white
    on white is a real possibility; over the flat panel surface it is just mud.
    """
    (w, _h), _ = cv2.getTextSize(s, FONT, scale, weight)
    x, y = org
    if anchor == "center":
        x -= w // 2
    elif anchor == "right":
        x -= w
    if outline:
        cv2.putText(img, s, (x, y), FONT, scale, (0, 0, 0), weight + 2, cv2.LINE_AA)
    cv2.putText(img, s, (x, y), FONT, scale, colour, weight, cv2.LINE_AA)
    return w


def fmt_seconds(t):
    return f"{t:.5f} s" if np.isfinite(t) else "  --   s"


# --- timestamps: the trigger train, the strobe train, or the host arrival times ---------

def fit(values, n):
    """Force an array to exactly n entries, padding short ones with NaN.

    A camera can store more or fewer frames than the other, and ts.csv can hold fewer rows
    than frames if the run was cut short. Both are real and neither should stop the render;
    a missing timestamp draws as '--' and contributes nothing to the histograms.
    """
    out = np.full(n, np.nan)
    m = min(n, len(values))
    out[:m] = values[:m]
    return out


def relative(t):
    """Seconds since this camera's own first frame.

    Done for BOTH timestamp sources, before anything is drawn, so a strobe-timed render and
    a ts.csv-timed one are read the same way and neither source gets a flattering origin.
    """
    finite = np.where(np.isfinite(t))[0]
    return t - t[finite[0]] if len(finite) else t


def frame_times_from_anchor(times, anchor, n):
    """Frame k -> spike s0+k, raw: no correction for the strobes lost at the OFF gaps.

    See the module docstring. The correction differs per camera on this rig and is not
    knowable from the strobe train alone, so it is left off and its one visible artefact --
    the doubled interval at the start anchor -- is allowed to show in the histogram.

    An index off either end of the channel gives NaN rather than a wrapped or clamped time.
    """
    idx = anchor["s0"] + np.arange(n)
    t = np.full(n, np.nan)
    inside = (idx >= 0) & (idx < len(times))
    t[inside] = np.asarray(times)[idx[inside]]
    return relative(t)


def find_smrx(session_dir):
    """The session's single .smrx, or raise StrobeUnusable saying how many there were.

    Shared by both Spike2 clocks: neither can be read without the file, and both treat none
    and more-than-one the same way -- say so, and let the caller fall back to ts.csv.
    """
    found = sorted(Path(session_dir).glob("*.smrx"))
    if len(found) != 1:
        raise StrobeUnusable(f"the session holds {len(found)} .smrx files, not 1"
                             + (f" ({', '.join(p.name for p in found)})" if found else ""))
    return found[0]


def strobe_times(session_dir, sess, counts):
    """Per-camera timestamps from the Spike2 strobe channels, or raise StrobeUnusable.

    Only the strobe-INTRINSIC checks in anchor_map gate this, exactly as they gate the
    analysis window in strobe_timing.choose_bracket. They are properties of the TTL train by
    itself, so failing one means the anchor is not really an anchor. a["measured"] is
    deliberately NOT gated on: it records a benign per-camera hardware difference in how the
    OFF-transition frame is handled, and gating on it rejects a healthy camera.
    """
    smrx = find_smrx(session_dir)
    f, open_note = strobe_timing.open_sonfile(smrx)     # copies it if Spike2 has it locked
    channels = strobe_timing.read_all_channels(f)
    ttl5 = strobe_timing.find_channel(channels, strobe_timing.TTL5_CHAN, "TTL5")
    ttl6 = strobe_timing.find_channel(channels, strobe_timing.TTL6_CHAN, "TTL6")
    # The Keyboard channel is deliberately not looked up. It plays no part in timing frames
    # -- this render covers the whole file, not the bracket between its marks -- and
    # requiring it would reject a perfectly good strobe train over an unrelated channel.
    tick_ms = f.GetTimeBase() * 1e3

    ttl5_cam, why, _ = strobe_timing.infer_assignment(ttl5, ttl6, sess)
    if ttl5_cam is None:
        raise StrobeUnusable(f"which .bin each TTL channel feeds is undecided: {why}")

    out = {}
    for cam in (1, 2):
        ch = ttl5 if cam == ttl5_cam else ttl6
        a = strobe_timing.anchor_for(ch, sess, cam, tick_ms)
        if a is None:
            raise StrobeUnusable(f"no strobe anchor on Ch{ch['chan']}, the channel feeding "
                                 f"c{cam}.bin")
        failed = [f"Ch{ch['chan']} (c{cam}.bin) {name}: {detail}"
                  for name, (passed, detail) in a["checks"].items() if not passed]
        if failed:
            raise StrobeUnusable("anchor check failed -- " + "; ".join(failed))
        if a["s0"] is None:
            raise StrobeUnusable(f"the c{cam}.bin anchor has no frame index "
                                 f"(no bright frames recorded in session.json)")
        out[cam] = frame_times_from_anchor(ch["times"], a, counts[cam])

    return out, (f"{smrx.name} {open_note}; TTL5 = c{ttl5_cam}.bin because {why}; "
                 f"every anchor check passed on both channels")


def frame_times_from_trigger(times, n):
    """Frame k -> pulse k, the same pulses for both cameras.

    Under the trigger there is nothing to align. The camera exposes only when a pulse tells it
    to, so stored frame k IS pulse k of the train on every camera (see eyecal/trigger.py): no
    anchor, no per-camera offset, and no OFF-gap question to refuse to guess at.

    Frames past the last pulse get NaN, like any other missing timestamp, so a count that does
    not match is reported rather than patched. See trigger_times.
    """
    return relative(fit(np.asarray(times, dtype=float), n))


def trigger_times(session_dir, sess, counts):
    """Per-camera timestamps from the Spike2 trigger train, or raise StrobeUnusable.

    Returns (times_per_cam, why, mismatch). `mismatch` carries (frames, pulses) for each camera
    whose stored frame count differs from the number of pulses, and is empty when they agree;
    camera_times prints it and marks the label, and the render goes ahead either way.

    THE WHOLE TTL2 CHANNEL IS USED, and its FIRST spike is frame 0. That is exact for a file
    holding the run's trigger train and nothing else, and wrong for one that also holds a test
    train recorded earlier -- every frame would then be read against an earlier pulse. Record
    the trigger train once, for the run.

    `sess` is not read. It is in the signature so camera_times can call this and strobe_times
    the same way; under the trigger, session.json has nothing left to contribute to the timing.
    """
    smrx = find_smrx(session_dir)
    f, open_note = strobe_timing.open_sonfile(smrx)     # copies it if Spike2 has it locked
    channels = strobe_timing.read_all_channels(f)
    try:
        ttl2 = strobe_timing.find_channel(channels, TRIGGER_CHAN, "TTL2")
    except KeyError as exc:
        raise StrobeUnusable(f"no trigger train: {exc.args[0]}") from exc
    times = ttl2["times"]
    if times is None or not len(times):
        raise StrobeUnusable(f"Ch{TRIGGER_CHAN} ({ttl2['title']}) holds no pulses, so there is "
                             f"no trigger train to read the frames against")

    # Both cameras from the SAME array, deliberately: that is what the trigger means. It also
    # makes the simultaneity panel a flat zero, which documents the shared clock.
    out = {cam: frame_times_from_trigger(times, counts[cam]) for cam in sorted(counts)}
    mismatch = {cam: (counts[cam], len(times)) for cam in sorted(counts)
                if counts[cam] != len(times)}
    why = (f"{smrx.name} {open_note}; TTL2 (Ch{TRIGGER_CHAN}) {len(times)} spikes; "
           + ", ".join(f"c{cam}.bin {counts[cam]} frames" for cam in sorted(counts)))
    return out, why, mismatch


def ts_csv_times(session_dir, counts):
    """Per-camera host arrival times from ts.csv, relative to each camera's first row.

    Two schemas exist on disk. Current sessions write eyecal.session.TS_COLUMNS, whose
    arrival column is t_arrive_s and is already relative to the run start; older ones wrote
    t_host_s, a raw perf_counter reading. Either is fine here, because relative() removes
    the origin anyway -- only the column name has to be found.
    """
    path = Path(session_dir) / "ts.csv"
    if not path.exists():
        raise FileNotFoundError(f"no such file: {path}")

    with path.open(newline="", encoding="utf-8") as fh:
        rows = csv.reader(fh)
        head = next(rows)
        col = next((c for c in ("t_arrive_s", "t_host_s") if c in head), None)
        if col is None:
            raise ValueError(f"{path} has no t_arrive_s or t_host_s column; found {head}")
        i_cam, i_idx, i_t = head.index("cam"), head.index("frame_idx"), head.index(col)
        per = {cam: [] for cam in counts}
        for r in rows:
            cam = int(r[i_cam])
            if cam in per:
                per[cam].append((int(r[i_idx]), float(r[i_t])))

    out = {}
    for cam, pairs in per.items():
        pairs.sort()                                   # written per thread, so not in order
        out[cam] = relative(fit(np.array([t for _, t in pairs], dtype=float), counts[cam]))
    return out, f"{path.name} column {col}"


def camera_times(session_dir, sess, counts):
    """The timestamps to draw, plus the on-screen label saying where they came from.

    session.json chooses the Spike2 clock: trigger.mode == "trigger" means the TTL2 pulse
    train, anything else -- free-run, a run that fell back to free-run, or an older session
    with no trigger block at all -- means the strobe anchors. A trigger session is never tried
    against the strobe anchors: its anchor is switched off, so there is nothing there to verify.

    Any failure on either Spike2 path -- a locked or absent .smrx, a channel that is not there,
    an anchor that does not verify, sonpy throwing -- falls back to ts.csv and prints the
    reason. There is always a render; there is never a silent substitution.
    """
    triggered = bool(sess) and (sess.get("trigger") or {}).get("mode") == "trigger"
    try:
        if triggered:
            t, why, mismatch = trigger_times(session_dir, sess, counts)
            label = TRIGGER_MISMATCH_LABEL if mismatch else TRIGGER_LABEL
            print(f"  timestamps    {label}")
            print(f"                {why}")
            for cam in sorted(mismatch):
                n_frames, n_pulses = mismatch[cam]
                print(f"                COUNT MISMATCH cam {cam}: {n_frames} frames stored, "
                      f"{n_pulses} TTL2 spikes -- frame k is still shown against spike k "
                      f"from the start")
            return t, label
        t, why = strobe_times(session_dir, sess, counts)
        print(f"  timestamps    {STROBE_LABEL}")
        print(f"                {why}")
        return t, STROBE_LABEL
    except StrobeUnusable as exc:
        reason = str(exc)
    except Exception as exc:
        reason = f"{type(exc).__name__}: {exc}"

    t, why = ts_csv_times(session_dir, counts)
    print(f"  timestamps    {TSCSV_LABEL}")
    print(f"                from {why}")
    print(f"                fell back because {reason}")
    return t, TSCSV_LABEL


def real_fps(sess, times):
    """The rate the session was actually recorded at, and where that number came from."""
    if sess:
        rates = [c.get("fpsEffective") for c in (sess.get("perCamera") or [])
                 if c.get("fpsEffective")]
        if rates:
            return float(np.mean(rates)), "mean of perCamera fpsEffective in session.json"
    # No session.json: fall back to the frames themselves -- frames over elapsed time, the
    # same arithmetic fpsEffective uses. NOT the median interval: host arrivals bunch under
    # the OS scheduler, so their median gap (measured 16.3 ms on this rig) sits well below
    # the real 20.27 ms period and would play the video a quarter too fast. The mean is
    # right to a few tens of ppm because bunching delays frames without losing any.
    t = times[1][np.isfinite(times[1])]
    if len(t) < 2 or t[-1] <= t[0]:
        raise ValueError("cannot work out the recorded frame rate: no session.json and no "
                         "usable cam 1 timestamps")
    return ((len(t) - 1) / float(t[-1] - t[0]),
            "cam 1 frames over elapsed time (no session.json)")


# --- the growing histograms --------------------------------------------------------------

def make_panel(title, values, colour, as_rate=False):
    """One bottom panel, with its axes fixed from the whole series before rendering starts.

    values[k] is the datapoint that arrives at frame k, NaN where there is none -- so the
    interval panels are empty at k=0 and every panel gains at most one count per frame.

    as_rate=True also shows the median as a frequency, which is meaningful for an interval
    (1000 / 20.27 ms = 49.33 Hz) and nonsense for a skew, so the simultaneity panel leaves
    it off.

    The range is ROBUST: centred on the median, half-width the larger of IQR_MULT x IQR and
    MIN_HALF_MS. Taking it from the min and max instead -- the obvious thing, and what this
    first did -- let the two ~40 ms anchor intervals set a 20 ms axis, so each of 50 bins
    was 0.4 ms wide and 99.7% of the data sat in one bar. The panel then said nothing about
    consistency, which is its whole job. With the bulk setting the range, a bin is one
    Spike2 tick wide and the quantised structure at 20.26 / 20.27 / 20.28 ms is visible
    from the first few frames. The outliers are not dropped: each is counted at whichever
    end it fell off, and draw_panel puts a marker with the count there.

    The edges are laid so the median sits at the CENTRE of a bin. Values quantised to a tick
    then land mid-bin rather than on an edge, where floating-point rounding would scatter
    them between two bins and comb the histogram.

    The bin index of every datapoint is worked out here too. Rendering then only has to
    increment one counter per frame instead of re-binning the whole history 1500 times.
    """
    ok = np.isfinite(values)
    if ok.any():
        v = values[ok]
        med = float(np.median(v))
        q1, q3 = np.percentile(v, [25, 75])
        half = max(HIST_IQR_MULT * float(q3 - q1), HIST_MIN_HALF_MS)
    else:
        v, med, half = np.empty(0), 0.0, HIST_MIN_HALF_MS
    binw = 2.0 * half / HIST_BINS
    edges = med + (np.arange(HIST_BINS + 1) - HIST_BINS / 2 - 0.5) * binw
    lo, hi = float(edges[0]), float(edges[-1])

    # -1 no datapoint; -2 fell off the low end; -3 off the high end; otherwise the bin.
    which = np.full(len(values), -1, dtype=np.int64)
    idx = np.where(ok)[0]
    which[idx] = np.clip(np.digitize(v, edges) - 1, 0, HIST_BINS - 1)
    which[idx[v < lo]] = -2
    which[idx[v > hi]] = -3
    inside = which >= 0
    # The tallest bar the finished histogram will ever have, so the y scale never moves.
    tallest = (int(np.max(np.bincount(which[inside], minlength=HIST_BINS)))
               if inside.any() else 1)

    return {"title": title, "colour": colour, "values": values, "which": which,
            "edges": edges, "ymax": max(1, tallest), "as_rate": as_rate,
            "counts": np.zeros(HIST_BINS, dtype=np.int64), "below": 0, "above": 0}


def advance(panel, k):
    """Take in frame k's datapoint, if it has one: into a bin, or off one end of the axis."""
    b = panel["which"][k]
    if b >= 0:
        panel["counts"][b] += 1
    elif b == -2:
        panel["below"] += 1
    elif b == -3:
        panel["above"] += 1


def draw_panel(img, box, panel, k, theme):
    """Draw one panel's histogram as it stands at frame k, into box = (x, y, w, h)."""
    x0, y0, w, h = box
    left, right = x0 + 12, x0 + w - 12
    top, bottom = y0 + 44, y0 + h - 34
    # A strip at each end of the axis is kept clear of bars, so a count that fell off the
    # axis has somewhere to stand without sitting on top of the edge bin.
    bl, br = left + HIST_OVERFLOW_ZONE, right - HIST_OVERFLOW_ZONE

    text(img, panel["title"], (left, y0 + 18), 0.45, theme["secondary"])

    seen = panel["values"][:k + 1]
    seen = seen[np.isfinite(seen)]
    if len(seen):
        med = float(np.median(seen))
        median = f"median {med:.3f}"
        if panel["as_rate"] and med > 0:
            median += f" ({1000.0 / med:.2f} Hz)"
        sd = f"{np.std(seen, ddof=1):.3f}" if len(seen) > 1 else "--"
        stat = f"n {len(seen)}   {median}   sd {sd} ms"
    else:
        stat = "n 0"
    outside = panel["below"] + panel["above"]
    if outside:
        stat += f"   off axis {outside}"
    text(img, stat, (left, y0 + 34), 0.40, theme["muted"])

    cv2.line(img, (bl, bottom), (br, bottom), theme["grid"], 1, cv2.LINE_AA)

    counts, ymax, span = panel["counts"], panel["ymax"], br - bl
    for i, c in enumerate(counts):
        if c <= 0:
            continue
        bx0 = bl + int(i * span / HIST_BINS)
        # One pixel of air between neighbouring bars: without it a full histogram reads as
        # a single filled blob and the shape of the distribution disappears.
        bx1 = max(bx0, bl + int((i + 1) * span / HIST_BINS) - 1)
        bh = max(1, int(round(c / ymax * (bottom - top))))
        cv2.rectangle(img, (bx0, bottom - bh), (bx1, bottom), panel["colour"], -1)

    # Overflow: an arrowhead pointing off the axis at the end the values fell off, with
    # the count beside it. Drawn only when there is something to report.
    for n, x, d in ((panel["below"], bl - 8, -1), (panel["above"], br + 8, 1)):
        if n:
            tri = np.array([[x, bottom - 9], [x, bottom - 1], [x + d * 8, bottom - 5]])
            cv2.fillPoly(img, [tri], panel["colour"], cv2.LINE_AA)
            text(img, str(n), (x + d * 12, bottom - 2), 0.40, panel["colour"],
                 anchor="left" if d > 0 else "right")

    e = panel["edges"]
    ty = y0 + h - 14
    text(img, f"{e[0]:.2f}", (bl, ty), 0.38, theme["muted"])
    text(img, f"{0.5 * (e[0] + e[-1]):.2f}", ((bl + br) // 2, ty), 0.38,
         theme["muted"], anchor="center")
    text(img, f"{e[-1]:.2f}", (br, ty), 0.38, theme["muted"], anchor="right")


# --- the render --------------------------------------------------------------------------

def draw_camera(img, box, frame, geom, cam, t_now, t_prev, label, theme):
    """One camera panel: the frame, letterboxed, with its timestamps written over it."""
    x0, y0, w, h = box
    dw, dh, ox, oy = geom
    small = cv2.resize(frame, (dw, dh), interpolation=cv2.INTER_AREA)
    roi = img[y0 + oy:y0 + oy + dh, x0 + ox:x0 + ox + dw]
    roi[:] = cv2.cvtColor(small, cv2.COLOR_GRAY2BGR)

    # Crosshair through the image centre, for judging where each camera is pointed.
    # Blended rather than drawn solid: it stays visible over a bright frame and never hides
    # the pixel it is sitting on. `roi` is a view into img, so the blend lands in place.
    lines = roi.copy()
    cx, cy = dw // 2, dh // 2
    cv2.line(lines, (cx, 0), (cx, dh - 1), CROSSHAIR, 1, cv2.LINE_AA)
    cv2.line(lines, (0, cy), (dw - 1, cy), CROSSHAIR, 1, cv2.LINE_AA)
    cv2.addWeighted(lines, 0.35, roi, 0.65, 0, dst=roi)

    text(img, f"cam {cam} - c{cam}.bin", (x0 + 16, y0 + 62), 0.55, theme["primary"],
         outline=True)
    text(img, label, (x0 + w - 16, y0 + h - 16), 0.45, theme["secondary"], anchor="right",
         outline=True)

    # 15% up from the panel's bottom edge: clear of the letterbox bar, and low enough in the
    # image that it sits under the eye rather than across it.
    ty = y0 + h - int(0.15 * h)
    now = fmt_seconds(t_now)
    width = cv2.getTextSize(now, FONT, 1.0, 2)[0][0]
    text(img, now, (x0 + w // 2, ty), 1.0, theme["primary"], weight=2, anchor="center",
         outline=True)
    if t_prev is not None:
        # The frame just gone, dim and small, to its left: it reads as already passed, and
        # the gap between the two numbers is the interval the left histogram is counting.
        text(img, fmt_seconds(t_prev), (x0 + w // 2 - width // 2 - 18, ty), 0.7,
             (110, 110, 110), anchor="right", outline=True)


def draw_progress(img, k, n, t, seam_y, width, theme):
    """A bar along the seam between the two rows, filled to how far through the run we are.

    Filled by frame fraction, not time fraction, so it is exact and monotonic even across a
    frame whose timestamp is missing. The elapsed / total seconds beside it come from cam 1
    and fall back to a frame count if either end has no timestamp.
    """
    y0, y1 = seam_y - 2, seam_y + 2
    cv2.rectangle(img, (0, y0), (width - 1, y1), theme["grid"], -1)
    fill = int(round((k + 1) / n * width))
    cv2.rectangle(img, (0, y0), (max(0, fill - 1), y1), theme["secondary"], -1)

    finite = np.isfinite(t)
    total = t[finite][-1] if finite.any() else np.nan
    label = (f"{t[k]:.2f} / {total:.2f} s" if np.isfinite(t[k]) and np.isfinite(total)
             else f"frame {k + 1} / {n}")
    text(img, label, (width - 16, seam_y + 18), 0.42, theme["secondary"], anchor="right")


def render(session_dir, out_path, size, fps, mm, times, label, panels, theme):
    """Draw every frame and pipe it to ffmpeg. Returns (frames written, seconds taken)."""
    W, H = size
    top_h = int(round(0.70 * H))                       # 756 of 1080: the cameras dominate
    cam_w, bot_h, bot_w = W // 2, H - top_h, W // 3
    n = len(panels[0]["values"])

    fh, fw = mm[1].shape[1], mm[1].shape[2]
    scale = min(cam_w / fw, top_h / fh)                # 1600x1200 into 960x756 -> 960x720
    dw, dh = int(fw * scale), int(fh * scale)
    geom = (dw, dh, (cam_w - dw) // 2, (top_h - dh) // 2)

    exe = imageio_ffmpeg.get_ffmpeg_exe()
    # stderr goes to a file rather than a pipe: an undrained pipe deadlocks a long encode,
    # and DEVNULL would throw away the only explanation of a failure.
    log = tempfile.TemporaryFile()
    proc = subprocess.Popen(
        [exe, "-y", "-f", "rawvideo", "-vcodec", "rawvideo", "-pix_fmt", "bgr24",
         "-s", f"{W}x{H}", "-r", str(fps), "-i", "-", "-an", "-c:v", "libx264",
         "-pix_fmt", "yuv420p", "-crf", "20", "-preset", "medium",
         "-movflags", "+faststart", str(out_path)],
        stdin=subprocess.PIPE, stdout=subprocess.DEVNULL, stderr=log)

    canvas = np.empty((H, W, 3), dtype=np.uint8)
    name = Path(session_dir).name
    started = time.perf_counter()
    written = 0
    try:
        for k in range(n):
            canvas[:] = theme["surface"]
            for cam in (1, 2):
                draw_camera(canvas, ((cam - 1) * cam_w, 0, cam_w, top_h),
                            np.asarray(mm[cam][k]), geom, cam, times[cam][k],
                            times[cam][k - 1] if k else None, label, theme)
            text(canvas, name, (W // 2, 28), 0.7, theme["primary"], weight=2,
                 anchor="center", outline=True)

            for i, panel in enumerate(panels):
                advance(panel, k)
                draw_panel(canvas, (i * bot_w, top_h, bot_w, bot_h), panel, k, theme)
            draw_progress(canvas, k, n, times[1], top_h, W, theme)

            proc.stdin.write(canvas.tobytes())
            written += 1
            if written % 100 == 0 or written == n:
                print(f"  {written}/{n} frames  ({time.perf_counter() - started:.1f} s)")
    except BrokenPipeError:
        pass                                           # ffmpeg died; the log below says why
    finally:
        if proc.stdin and not proc.stdin.closed:
            proc.stdin.close()
        proc.wait()

    if proc.returncode != 0:
        log.seek(0)
        tail = log.read().decode("utf-8", "replace").strip().splitlines()[-12:]
        log.close()
        raise RuntimeError(f"ffmpeg exited {proc.returncode} after {written} frames:\n  "
                           + "\n  ".join(tail))
    log.close()
    return written, time.perf_counter() - started


# --- driver -------------------------------------------------------------------------------

def choose_session_dir(given):
    """The session to render. None means the operator cancelled the dialog."""
    if given:
        return Path(given).expanduser().resolve()

    import tkinter
    from tkinter import filedialog
    root = tkinter.Tk()
    root.withdraw()
    opts = {"title": "Choose a recorded session folder"}
    if DEFAULT_ROOT.exists():
        opts["initialdir"] = str(DEFAULT_ROOT)
    picked = filedialog.askdirectory(**opts)
    root.destroy()
    return Path(picked).resolve() if picked else None


def main(argv=None):
    args = parse_args(argv)
    session_dir = choose_session_dir(args.session_dir)
    if session_dir is None:
        print("No folder chosen; nothing to render.")
        return 0
    if not session_dir.is_dir():
        raise NotADirectoryError(f"no such session directory: {session_dir}")

    print()
    print("=" * 88)
    print("SESSION VIDEO")
    print("=" * 88)
    print(f"  session       {session_dir}")

    mm = {}
    for cam in (1, 2):
        try:
            mm[cam] = read_frames(session_dir, cam)
        except (OSError, KeyError, ValueError, json.JSONDecodeError) as exc:
            raise FileNotFoundError(f"cannot read c{cam}.bin from {session_dir}: "
                                    f"{exc}") from exc
    counts = {cam: int(mm[cam].shape[0]) for cam in mm}
    n = min(counts.values())
    if n < 1:
        raise ValueError(f"nothing to render: c1.bin has {counts[1]} frames, c2.bin "
                         f"{counts[2]}")
    print(f"  frames        c1.bin {counts[1]}, c2.bin {counts[2]}  ->  rendering {n}")

    sess = strobe_timing.read_session(session_dir)
    times, label = camera_times(session_dir, sess, counts)
    times = {cam: fit(v, n) for cam, v in times.items()}

    fps_real, fps_why = real_fps(sess, times)
    fps = args.speed * fps_real
    if fps <= 0:
        raise ValueError(f"--speed {args.speed} gives a video frame rate of {fps}")
    print(f"  recorded at   {fps_real:.4f} fps  ({fps_why})")
    print(f"  playing at    {fps:.4f} fps  ({args.speed:g} x real time, "
          f"{n / fps:.1f} s long)")

    theme = {k: (bgr(v) if isinstance(v, str) else v)
             for k, v in strobe_timing.THEMES["dark"].items()}
    c1, c2 = (bgr(c) for c in strobe_timing.THEMES["dark"]["series"])
    pair = bgr(strobe_timing.THEMES["dark"]["pair"])

    # values[k] is the datapoint that appears at frame k. The interval panels have none at
    # k=0 by definition; simultaneity is defined from the first frame, so it starts with one.
    d1 = np.concatenate([[np.nan], np.diff(times[1]) * 1e3])
    d2 = np.concatenate([[np.nan], np.diff(times[2]) * 1e3])
    panels = [
        make_panel("cam 1 interval   t1[k] - t1[k-1]   (ms)", d1, c1, as_rate=True),
        make_panel("simultaneity   t1[k] - t2[k]   (ms)", (times[1] - times[2]) * 1e3, pair),
        make_panel("cam 2 interval   t2[k] - t2[k-1]   (ms)", d2, c2, as_rate=True),
    ]

    out_path = Path(args.out).expanduser() if args.out \
        else session_dir / "calibration_video.mp4"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"  writing       {out_path}")

    written, elapsed = render(session_dir, out_path, (args.width, args.height), fps, mm,
                              times, label, panels, theme)
    mb = out_path.stat().st_size / 1e6
    print(f"\nwrote {out_path}\n  {written} frames at {fps:.4f} fps, {mb:.1f} MB, "
          f"{elapsed:.1f} s to render")
    print("=" * 88)
    return 0


def parse_args(argv):
    p = argparse.ArgumentParser(
        description="Render a recorded camera session to an H.264 mp4, with each camera's "
                    "timestamps and growing interval histograms drawn on it.")
    p.add_argument("session_dir", nargs="?", default=None,
                   help="the session folder holding frames/, ts.csv and session.json. "
                        "Omit it to pick one with a folder dialog")
    p.add_argument("--speed", type=float, default=0.3,
                   help="playback speed as a fraction of real time; 0.3 means the video "
                        "runs at 30%% of the recorded rate (default: %(default)s)")
    p.add_argument("--out", default=None,
                   help="output file (default: <session_dir>/calibration_video.mp4)")
    p.add_argument("--width", type=int, default=1920, help="default: %(default)s")
    p.add_argument("--height", type=int, default=1080, help="default: %(default)s")
    return p.parse_args(argv)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"\nFAILED: {exc}", file=sys.stderr)
        sys.exit(1)
