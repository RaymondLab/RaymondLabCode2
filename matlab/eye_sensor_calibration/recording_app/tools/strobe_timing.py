r"""Strobe timing: do the camera strobes in Spike2 match the frames on disk?

    python tools\strobe_timing.py "C:\Temp\test\20260819-011856Z_ov2311_2cam_35s\calibration.smrx"

Name the .smrx; frames/ and session.json are read from beside it. A session directory is
accepted too when it holds exactly one .smrx.

A tool, not a test. It lives in tools/ rather than tests/ because it evaluates the app's
OUTPUT -- a recorded session, against the Spike2 file -- rather than the app's code, and
because a test_ prefix is a discovery convention: pytest would collect a 1200-line module
that asserts nothing.

Reads the Spike2 .smrx with sonpy, extracts every channel, and analyses the two camera
strobe trains -- TTL5 (Ch15) and TTL6 (Ch16) -- over the segment bracketed by the two
Keyboard (Ch31) markers, which mark the start and end of image acquisition. Two figures
come out:

    fig1_strobe_intervals.png   interval between consecutive spikes per channel over
                                time, plus the marginal histogram of both
                                -> how steady each camera's own sampling clock is
    fig2_pair_skew.png          TTL5-TTL6 skew per matched pair over time, plus its
                                marginal histogram
                                -> how simultaneous the two cameras are with each other

This is a measurement script, not a pass/fail unit test: it prints numbers and draws
pictures. tests/test_offline.py is the one that asserts.

Two things this file works around, both real and both easy to trip over:

Spike2 keeps an exclusive lock on a .smrx it has open, and sonpy then fails to open it
with No_File (-1) no matter which OpenFlags you pass -- Shared included. Rather than
requiring the operator to close Spike2 first, the reader copies the file and reads the
copy. See open_sonfile.

Which .bin belongs to which TTL channel is NOT knowable from the channel titles alone.
The Spike2 channel comments name a wiring position ("Camera 1 Strobe"), while c1.bin and
c2.bin are named by the operator's ACCEPTED left-to-right order, which the alignment step
may have swapped. So the assignment is inferred from evidence and printed with its
reasoning; --assign overrides it. See infer_assignment.
"""

import argparse
import json
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

import sonpy

# Channel numbers as Spike2 shows them (1-based); sonpy indexes from 0.
TTL5_CHAN, TTL6_CHAN, KEYBOARD_CHAN = 15, 16, 31

MAX_EVENTS = 5_000_000        # per event channel; far above any real session
MAX_WAVE_POINTS = 2_000_000   # waveform channels are summarised, not plotted

# Categorical slots 1-3 of the reference palette, in their fixed order, validated against
# both surfaces: lightness band and chroma floor pass; adjacent-pair CVD dE 9.2 (light) /
# 9.4 (dark) against a target of 8; normal-vision dE 27.6 / 26.5 against a floor of 15.
#
# Slot 1 is always TTL5 and slot 2 always TTL6, in both figures -- colour follows the
# channel, never the position. Slot 3 (aqua) is the pair-skew series in figure 2, which is
# a relation BETWEEN the two channels rather than either one of them, so it must not wear
# either channel's hue. Aqua sits at 2.74:1 on the light surface, a contrast WARN: it is
# carried by a direct median label on the plot and the full stats table on the console,
# never by colour alone.
THEMES = {
    "light": {"series": ("#2a78d6", "#eb6834"), "pair": "#1baf7a", "surface": "#fcfcfb",
              "primary": "#0b0b0b", "secondary": "#52514e", "muted": "#8a8880",
              "grid": "#e3e2dd"},
    "dark":  {"series": ("#3987e5", "#d95926"), "pair": "#199e70", "surface": "#1a1a19",
              "primary": "#ffffff", "secondary": "#c3c2b7", "muted": "#8a8880",
              "grid": "#333331"},
}


# --- reading the Spike2 file ----------------------------------------------------------

def error_name(code):
    """Name a sonpy open error.

    GetErrorString is not a lookup for these codes -- handed -1 it returns "This object
    does not own a file handle or any resources", which says nothing about the file. The
    module-level constants are the real mapping, so read them instead.
    """
    for name in ("Son_OK", "No_File", "No_Access", "Read_Only", "Bad_Read", "Bad_Write",
                 "Wrong_File", "Corrupt_File", "No_Channel", "No_Memory", "Bad_Param",
                 "Past_EOF", "Past_SOF", "No_Extra", "No_Block"):
        value = getattr(sonpy, name, None)
        if value is not None and int(value) == int(code):
            return f"{code} {name}"
    return str(code)


def open_sonfile(path):
    """Open a .smrx read-only, copying it first if Spike2 has it locked.

    A .smrx open in Spike2 cannot be opened by sonpy at all: every OpenFlags value returns
    No_File (-1), including Shared. Copying sidesteps it, because the lock blocks opening
    rather than reading. Sessions are small -- a ten second one is 320 KB -- so the copy
    costs nothing and the operator never has to close Spike2 to run this.

    Returns (SonFile, note), the note describing what happened for the console header.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"no such file: {path}")

    f = sonpy.SonFile(str(path), True)
    if f.GetOpenError() == 0:
        return f, "opened directly"
    locked_err = f.GetOpenError()

    tmp = Path(tempfile.mkdtemp(prefix="spiketiming_")) / path.name
    try:
        shutil.copy2(path, tmp)
    except OSError as exc:
        # Spike2 normally leaves the file readable, so the copy succeeds where the open
        # failed. If even the copy is refused, something is holding it far more tightly
        # and there is nothing this script can do about it from here.
        raise OSError(
            f"sonpy could not open {path} ({error_name(locked_err)}) and copying "
            f"it failed too: {exc}. "
            f"Close the file in Spike2 (or any other program holding it) and re-run."
        ) from exc

    f = sonpy.SonFile(str(tmp), True)
    err = f.GetOpenError()
    if err != 0:
        raise OSError(f"sonpy could not open {path} ({error_name(locked_err)}) nor a "
                      f"copy of it at {tmp} ({error_name(err)})")
    return f, (f"original was locked ({error_name(locked_err)}, probably open in "
               f"Spike2); read a copy at {tmp}")


def read_all_channels(f):
    """Every channel that is not Off, as a list of dicts, dispatched by data type.

    Waveform channels are read but only summarised -- nothing here plots them, and a long
    session's ADC data is large enough that carrying it around costs more than it is worth.
    """
    tb = f.GetTimeBase()
    events = {sonpy.DataType.EventRise, sonpy.DataType.EventFall, sonpy.DataType.EventBoth}
    out = []
    for idx in range(f.MaxChannels()):
        kind = f.ChannelType(idx)
        if kind == sonpy.DataType.Off:
            continue

        ch = {
            "chan": idx + 1, "index": idx, "type": str(kind).replace("DataType.", ""),
            "title": f.GetChannelTitle(idx), "units": f.GetChannelUnits(idx),
            "comment": f.GetChannelComment(idx), "rate": f.GetIdealRate(idx),
            "max_time_s": f.ChannelMaxTime(idx) * tb, "times": None, "values": None,
        }

        if kind in events:
            ch["times"] = np.asarray(f.ReadEvents(idx, MAX_EVENTS, 0), dtype=np.int64) * tb
        elif kind == sonpy.DataType.Marker:
            marks = f.ReadMarkers(idx, MAX_EVENTS, 0)
            ch["times"] = np.array([m.Tick for m in marks], dtype=np.int64) * tb
            ch["values"] = [(m.Code1, m.Code2, m.Code3, m.Code4) for m in marks]
        elif kind == sonpy.DataType.TextMark:
            marks = f.ReadTextMarks(idx, MAX_EVENTS, 0)
            ch["times"] = np.array([m.Tick for m in marks], dtype=np.int64) * tb
            ch["values"] = [getattr(m, "Text", "") for m in marks]
        elif kind in (sonpy.DataType.Adc, sonpy.DataType.RealWave):
            reader = f.ReadInts if kind == sonpy.DataType.Adc else f.ReadFloats
            ch["values"] = np.asarray(reader(idx, MAX_WAVE_POINTS, 0))
            ch["times"] = None      # regularly sampled; rate and divide describe it
            ch["divide"] = f.ChannelDivide(idx)
        else:
            # AdcMark / RealMark: recorded so the channel is not silently dropped.
            ch["values"] = "not read (extended marker channel)"

        ch["n"] = (len(ch["times"]) if ch["times"] is not None
                   else (len(ch["values"]) if hasattr(ch["values"], "__len__") else 0))
        out.append(ch)
    return out


def find_channel(channels, chan_no, want_title=None):
    for ch in channels:
        if ch["chan"] == chan_no:
            if want_title and want_title.lower() not in ch["title"].lower():
                print(f"  note: Ch{chan_no} is titled {ch['title']!r}, expected something "
                      f"like {want_title!r} -- using it anyway.", file=sys.stderr)
            return ch
    raise KeyError(f"Ch{chan_no} is not present (or is Off) in this file. "
                   f"Present: {sorted(c['chan'] for c in channels)}")


def keyboard_window(kb):
    """The two Keyboard marks bracketing image acquisition, as (t0, t1, label).

    More than two marks is plausible if the operator typed anything else while sampling,
    so take the first and last rather than insisting on exactly two.
    """
    t = kb["times"]
    if t is None or len(t) < 2:
        raise ValueError(f"Ch{kb['chan']} ({kb['title']}) has {0 if t is None else len(t)} "
                         f"marker(s); need 2 to bracket the acquisition.")

    def code(i):
        c = kb["values"][i][0] if kb["values"] else 0
        return chr(c) if 32 <= c < 127 else f"0x{c:02x}"

    if len(t) > 2:
        print(f"  note: Ch{kb['chan']} has {len(t)} markers "
              f"({', '.join(code(i) for i in range(len(t)))}); using the first and last.",
              file=sys.stderr)
    return float(t[0]), float(t[-1]), f"{code(0)!r} at {t[0]:.4f} s -> {code(len(t) - 1)!r}"


# --- frames on disk -------------------------------------------------------------------

def read_frame_count(frames_dir, cam):
    """Frame count for c{cam}.bin, taken from the file's real size.

    The sidecar count is what the recorder believed it wrote; the size is what is actually
    on disk. They should agree, and a disagreement means a truncated write, so both are
    returned and the report flags any mismatch.
    """
    frames_dir = Path(frames_dir)
    binary, sidecar = frames_dir / f"c{cam}.bin", frames_dir / f"c{cam}.json"
    if not binary.exists():
        raise FileNotFoundError(f"no such file: {binary}")
    if not sidecar.exists():
        raise FileNotFoundError(f"no such file: {sidecar} (needed for the frame geometry)")

    meta = json.loads(sidecar.read_text(encoding="utf-8"))
    frame_bytes = meta["width"] * meta["height"] * np.dtype(meta["dtype"]).itemsize
    size = binary.stat().st_size
    return {
        "cam": cam, "path": binary, "declared": int(meta["count"]),
        "actual": size // frame_bytes, "partial": size % frame_bytes,
        "width": meta["width"], "height": meta["height"], "dtype": meta["dtype"],
    }


def read_session(session_dir):
    """session.json if present. It carries the accepted device order, which is the key
    piece of evidence for matching TTL channels to .bin files. Absent is not fatal."""
    p = Path(session_dir) / "session.json"
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        print(f"  note: could not read {p}: {exc}", file=sys.stderr)
        return None


# --- the strobe anchor ------------------------------------------------------------------
#
# The recorder kicks the exposure up four stops for holdS seconds at each end of the run
# (see eyecal/record.py). On this rig that does NOT widen the strobe -- the pulse width is
# unchanged either side of it. What it leaves in the TTL train is a pair of timing marks
# per anchor, and their arithmetic is exact:
#
#     ON  gap = P + d           the frame period stretches once as exposure lengthens
#     OFF gap = 2*P - d         a DOUBLED period: the camera stores a frame here but
#                               emits no strobe pulse for it
#
# with P the frame period and d the same constant in both (about 1.63 ms at 49.3 Hz). The
# EventBoth recording confirms the OFF gap is a true absence rather than a runt pulse: the
# line simply stays low across it, with no anomalous width anywhere in the file.
#
# So a recording loses exactly TWO strobes per channel, one at each OFF transition, and
# frames = strobes + 2 over the whole run. Nothing else goes missing. Every "off by one"
# in the counts traces back to how many OFF gaps a given interval happens to enclose.
#
# The four marks are found by their SPACING, not their shape: the ON/OFF gap signature on
# its own occurs dozens of times per run while the cameras are being set up (97 ISI
# anomalies in one measured file), but only one quadruple anywhere matches the schedule.

ANCHOR_HOLD_TOL_S = 0.12    # holdS is nominal; measured 0.2246-0.2651 across runs
ANCHOR_SPAN_TOL_S = 0.25
ANCHOR_ISI_TOL_MS = 1.0     # an ISI this far off the median is a candidate mark


def find_anchor(times, hold_s, seconds):
    """The four anchor spike indices, or None.

    Returns indices into `times` of the spike BEFORE each gap, in schedule order:
    (start_on, start_off, end_on, end_off). The gap itself is times[i+1] - times[i].
    """
    if len(times) < 8 or not hold_s or not seconds:
        return None
    isi = np.diff(times) * 1e3
    med = float(np.median(isi))
    cand = np.where(np.abs(isi - med) > ANCHOR_ISI_TOL_MS)[0]
    span = seconds - 2 * hold_s
    for i in range(len(cand) - 3):
        k = cand[i:i + 4]
        t = times[k]
        if (abs((t[1] - t[0]) - hold_s) < ANCHOR_HOLD_TOL_S
                and abs((t[3] - t[2]) - hold_s) < ANCHOR_HOLD_TOL_S
                and abs((t[2] - t[0]) - span) < ANCHOR_SPAN_TOL_S):
            return tuple(int(x) for x in k)
    return None


def anchor_map(times, idx, bright, tick_ms):
    """Frame-to-spike mapping from an anchor, plus the checks that say whether to trust it.

    The mapping is

        frame k  <->  spike (s0 + k) - (OFF gaps passed before frame k)

    where s0 comes from the start anchor: the first bright frame is exposed by the first
    strobe after the ON gap, so s0 = (start_on + 1) - first_bright_frame. The end anchor
    gives s0 a second time, independently, most of a recording downstream.

    Only the STROBE-INTRINSIC checks gate. They are properties of the Spike2 data by
    itself -- the gap arithmetic and the absence of stray anomalies -- so failing one means
    the train really is not what an anchor should look like, and the window cannot be
    trusted.

    s0 agreement does NOT gate, and the reason is worth recording. It was written as a
    gating check on the strength of a single run, on the assumption that at the OFF
    transition a camera stores a frame but emits no strobe, so spikes = frames - 1. Five
    subsequent 30 s runs showed that is one camera's behaviour and not a law: on this rig
    device 1 gives exactly that in 5 of 5 runs, while device 2 gives spikes = frames, as
    though it drops the frame along with the strobe. Gating on it rejected a healthy camera
    for a benign hardware difference. The offset is stable per camera, so it is measured
    and reported (`off_offset`) rather than assumed -- calibrate it per camera if you need
    a frame-exact index.

    d is measured from the ON gap rather than assumed, so a change of frame rate does not
    silently invalidate the OFF test.

    Returns a dict; `ok` reflects the gating checks only.
    """
    isi = np.diff(times) * 1e3
    P = float(np.median(isi))
    on_s, off_s, on_e, off_e = idx
    g_on_s, g_off_s = float(isi[on_s]), float(isi[off_s])
    g_on_e, g_off_e = float(isi[on_e]), float(isi[off_e])

    d = g_on_s - P                                   # the stretch, measured not assumed
    tol = max(4 * tick_ms, 0.25)
    checks = {
        "ON gaps are P+d": (abs(g_on_e - (P + d)) < tol,
                            f"{g_on_s:.2f} / {g_on_e:.2f} ms vs P+d = {P + d:.2f}"),
        "OFF gaps are 2P-d": (abs(g_off_s - (2 * P - d)) < tol
                              and abs(g_off_e - (2 * P - d)) < tol,
                              f"{g_off_s:.2f} / {g_off_e:.2f} ms vs 2P-d = {2 * P - d:.2f}"),
    }

    seg = isi[off_s:off_e]                           # between the two OFF gaps
    expected = {off_s, on_e}                         # the enclosed anchor events
    unexpected = [off_s + int(j) for j in np.where(np.abs(seg - P) > tol)[0]
                  if off_s + int(j) not in expected]
    checks["no stray anomalies in the bracket"] = (
        not unexpected, f"{len(unexpected)} unexpected"
        + (f" at spikes {unexpected[:5]}" if unexpected else ""))

    # Measured, not gated. off_offset is 0 for a camera that stores the OFF-transition
    # frame without a strobe, and +1 for one that drops the frame too.
    s0 = s0_end = off_offset = None
    if bright and bright.get("start") and bright.get("end"):
        s0 = (on_s + 1) - bright["start"][0]
        s0_end = (on_e + 1) - bright["end"][0] + 1   # +1: one OFF gap lies between them
        off_offset = s0_end - s0
        note = (f"start anchor {s0}, end anchor {s0_end}, offset {off_offset:+d}"
                + ("  (stores the OFF frame, no strobe)" if off_offset == 0 else
                   "  (drops the OFF frame as well as the strobe)" if off_offset == 1 else
                   "  -- neither known pattern; treat the frame index as uncertain"))
    else:
        note = "no bright frames in session.json"
    measured = {"s0 from each anchor": note}

    return {
        "idx": idx, "P_ms": P, "d_ms": d, "s0": s0, "s0_end": s0_end,
        "off_offset": off_offset,
        "t": [float(times[i]) for i in idx],
        "gaps_ms": [g_on_s, g_off_s, g_on_e, g_off_e],
        "checks": checks, "measured": measured, "ok": all(v[0] for v in checks.values()),
        "off_bracket": (float(times[off_s]), float(times[off_e])),
        # the steady-state span: past the start-OFF gap, short of the end-ON gap
        "inner": (float(times[off_s + 1]), float(times[on_e])),
        "spikes_between_off": off_e - off_s,
    }


def anchor_for(ch, sess, cam, tick_ms):
    """find_anchor + anchor_map for one channel, or None if the schedule is unknown."""
    if not sess:
        return None
    anc = (sess.get("strobeAnchors") or {})
    if not anc.get("enabled"):
        return None
    hold = anc.get("holdS")
    seconds = (sess.get("requested") or {}).get("seconds")
    idx = find_anchor(ch["times"], hold, seconds)
    if idx is None:
        return None
    bright = (anc.get("brightFrames") or {}).get(f"cam{cam}")
    return anchor_map(ch["times"], idx, bright, tick_ms)


# --- which TTL channel is which camera ------------------------------------------------

def infer_assignment(ttl5, ttl6, sess):
    """Decide whether TTL5 is c1.bin or c2.bin, and say why.

    The problem: c1.bin/c2.bin are numbered by the operator's ACCEPTED left-to-right order
    (session.json acceptedDeviceOrder), which alignment may have swapped relative to the
    device ids, while the TTL channels are wired to fixed physical inputs. So the titles
    and comments in the .smrx describe devices, not files.

    The evidence used:

    1. Open order, which is decisive here. open_all() opens the requested devices in order,
       sequentially, and each camera strobes from the moment its graph starts -- so the TTL
       channel whose FIRST event is earlier belongs to requested.deviceIds[0]. Map that
       device through acceptedDeviceOrder to get its cam number, hence its .bin. The margin
       is enormous: the cameras open seconds apart while their strobes are ~20 ms apart, so
       this is not a close call. It is still guarded by a 10-interval threshold, and gives
       up rather than guessing if the two trains start too close together.

    2. Channel comments, as corroboration only. "Camera 1 Strobe" names a wiring position,
       and whether that means device 1 or c1.bin is exactly the ambiguity in question -- so
       it is reported for the operator to weigh, never trusted.

    Returns (ttl5_cam, reason, corroboration), ttl5_cam in (1, 2) or None if undecidable.
    """
    comments = f"Ch15 comment {ttl5['comment']!r}; Ch16 comment {ttl6['comment']!r}"

    if ttl5["times"] is None or ttl6["times"] is None \
            or not len(ttl5["times"]) or not len(ttl6["times"]):
        return None, "one of the TTL channels is empty", comments
    if sess is None:
        return None, "no session.json, so the accepted device order is unknown", comments

    requested = list(sess.get("requested", {}).get("deviceIds") or [])
    accepted = list(sess.get("acceptedDeviceOrder") or [])
    if len(requested) < 2 or len(accepted) < 2:
        return None, "session.json has no usable device order", comments

    first5, first6 = float(ttl5["times"][0]), float(ttl6["times"][0])
    gap = abs(first5 - first6)
    isi = float(np.median(np.diff(ttl5["times"]))) if len(ttl5["times"]) > 1 else 0.0
    if gap < 10 * max(isi, 1e-3):
        return None, (f"the two strobe trains start only {gap * 1e3:.1f} ms apart, too "
                      f"close to read the open order from"), comments

    early_device = requested[0] if first5 < first6 else requested[1]
    if early_device not in accepted:
        return None, (f"device {early_device} is not in acceptedDeviceOrder "
                      f"{accepted}"), comments

    # c{k}.bin holds whichever camera ended up in accepted position k-1.
    early_cam = accepted.index(early_device) + 1
    ttl5_cam = early_cam if first5 < first6 else (3 - early_cam)

    earlier, later = ("TTL5", "TTL6") if first5 < first6 else ("TTL6", "TTL5")
    reason = (f"{earlier} starts {gap:.2f} s before {later}, so {earlier} is device "
              f"{early_device} (opened first, of requested {requested}); accepted order "
              f"{accepted} puts device {early_device} in c{early_cam}.bin")
    return ttl5_cam, reason, comments


# --- the two analyses -----------------------------------------------------------------

def in_window(times, t0, t1):
    return times[(times >= t0) & (times <= t1)]


def intervals(times):
    """Consecutive-spike intervals in ms, with the time each is plotted at -- the later
    spike of the pair, so an interval is drawn where it finished."""
    if len(times) < 2:
        return np.empty(0), np.empty(0)
    return times[1:], np.diff(times) * 1e3


def pair_spikes(a, b, tol_s):
    """Match TTL5 spikes to TTL6 spikes one-to-one, in time order.

    Why this rather than pairing by index. The two trains almost never hold the same count:
    one camera opens earlier, one keeps strobing a few frames longer as the app tears down,
    and either can drop a pulse mid-run. Zipping by index makes every pair after the first
    discrepancy meaningless -- it silently compares spike n of one camera against spike n+k
    of the other, and the plot then shows a staircase of multiples of the frame interval
    instead of the real sub-millisecond skew.

    So: a monotonic two-pointer merge. Walk both trains together; if the heads are within
    tol_s that is a pair and both advance, otherwise the earlier head is unmatched and only
    it advances. This is the order-preserving matching for two nearly simultaneous point
    processes, and with tol_s at half the frame interval it cannot pair across a frame
    boundary -- a spike is matched either to its true partner or to nothing. Order-preserving
    is the property that matters: pairs can never cross, so the unmatched spikes stay
    interpretable as real drops or as tail overhang rather than as bookkeeping artefacts.

    Returns (ta, dt_ms, unmatched_a, unmatched_b): the pair time (taken as the TTL5 spike),
    the signed skew a-b in ms, and the times of whatever did not pair on each side.
    """
    i = j = 0
    ta, dt, ua, ub = [], [], [], []
    while i < len(a) and j < len(b):
        delta = a[i] - b[j]
        if abs(delta) <= tol_s:
            ta.append(a[i])
            dt.append(delta * 1e3)
            i += 1
            j += 1
        elif delta < 0:
            ua.append(a[i])
            i += 1
        else:
            ub.append(b[j])
            j += 1
    ua.extend(a[i:])
    ub.extend(b[j:])
    return np.array(ta), np.array(dt), np.array(ua), np.array(ub)


def describe(x, unit="ms"):
    if len(x) == 0:
        return {"n": 0, "unit": unit}
    return {
        "n": len(x), "unit": unit, "mean": float(np.mean(x)), "median": float(np.median(x)),
        "sd": float(np.std(x, ddof=1)) if len(x) > 1 else float("nan"),
        "min": float(np.min(x)), "max": float(np.max(x)),
        "p1": float(np.percentile(x, 1)), "p99": float(np.percentile(x, 99)),
    }


def anomalies(times, isi, tol_ms):
    """Intervals more than tol_ms off the median -- the dropped or doubled pulses.

    Everything else sits within a tick or two of the median, so this is what "how
    consistent is it" actually reduces to once the bulk is known to be flat.
    """
    if not len(isi):
        return np.empty(0), np.empty(0)
    off = np.abs(isi - np.median(isi)) > tol_ms
    return times[off], isi[off]


def drift_fit(t, skew_ms):
    """Least-squares slope of the pair skew, as ms/s and as ppm.

    A non-zero slope means the two cameras' pixel clocks run at slightly different rates,
    so the pair separation walks steadily instead of sitting still. That is a different
    defect from jitter and the standard deviation alone will not show it -- a skew that
    ramps smoothly across the window has a large sd but perfect short-term stability.
    """
    if len(t) < 3:
        return None
    slope, intercept = np.polyfit(t - t[0], skew_ms, 1)
    return {"slope_ms_per_s": float(slope), "ppm": float(slope) * 1e3,
            "intercept_ms": float(intercept), "span_ms": float(slope) * (t[-1] - t[0])}


def analyse(ttl5, ttl6, t0, t1, window_label, tick_ms):
    w5 = in_window(ttl5["times"], t0, t1)
    w6 = in_window(ttl6["times"], t0, t1)
    if len(w5) < 2 or len(w6) < 2:
        raise ValueError(f"only {len(w5)} TTL5 and {len(w6)} TTL6 spikes fall between the "
                         f"Keyboard markers ({t0:.4f} s .. {t1:.4f} s); nothing to analyse.")

    isi5_t, isi5 = intervals(w5)
    isi6_t, isi6 = intervals(w6)

    # Half the median interval: wide enough for any real skew, narrow enough that a spike
    # can never be paired with its partner's neighbour.
    tol_s = 0.5 * float(np.median(np.concatenate([np.diff(w5), np.diff(w6)])))
    skew_t, skew, un5, un6 = pair_spikes(w5, w6, tol_s)

    # One tick is the finest difference the file can express, so "within a tick" is the
    # honest definition of "identical" here; anything past two ticks is a real deviation.
    near_tol = 2 * tick_ms
    an5_t, an5 = anomalies(isi5_t, isi5, near_tol)
    an6_t, an6 = anomalies(isi6_t, isi6, near_tol)

    return {
        "t0": t0, "t1": t1, "window_label": window_label, "w5": w5, "w6": w6,
        "isi5_t": isi5_t, "isi5": isi5, "isi6_t": isi6_t, "isi6": isi6,
        "skew_t": skew_t, "skew": skew, "unmatched5": un5, "unmatched6": un6,
        "tol_ms": tol_s * 1e3, "tick_ms": tick_ms, "near_tol_ms": near_tol,
        "anom5_t": an5_t, "anom5": an5, "anom6_t": an6_t, "anom6": an6,
        "drift": drift_fit(skew_t, skew),
        "s_isi5": describe(isi5), "s_isi6": describe(isi6), "s_skew": describe(skew),
    }


# --- figures ---------------------------------------------------------------------------

def style(theme):
    t = THEMES[theme]
    matplotlib.rcParams.update({
        "figure.facecolor": t["surface"], "axes.facecolor": t["surface"],
        "savefig.facecolor": t["surface"],
        "text.color": t["primary"], "axes.labelcolor": t["secondary"],
        "xtick.color": t["secondary"], "ytick.color": t["secondary"],
        "axes.edgecolor": t["grid"], "axes.linewidth": 0.8,
        "grid.color": t["grid"], "grid.linewidth": 0.8,
        "xtick.direction": "out", "ytick.direction": "out",
        "font.size": 9, "axes.titlesize": 10, "figure.dpi": 110,
        "legend.frameon": False,
    })
    return t


def tidy(ax, t, xlabel=None, ylabel=None, title=None):
    """Recessive frame: only the spines that carry meaning, grid behind the marks."""
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.set_axisbelow(True)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title, color=t["secondary"], loc="left", pad=8)


def full_limits(*arrays, pad=0.08, floor=None):
    """y-range covering ALL the data.

    Percentile clipping was the wrong instinct here and is worth spelling out. These
    intervals are quantised to the file's 10 us tick and over 99% of them land on a single
    value, so a 0.5-99.5 percentile range collapses the axis onto one tick -- turning
    meaningless quantisation dither into a picket fence of full-height spikes while hiding
    the handful of genuinely short intervals, which are the only interesting events in the
    series. Showing everything puts the anomalies on screen and lets the flat band read as
    what it is: flat. The marginal histogram carries the shape of the bulk, on a log count
    axis so a bin of one is still visible beside a bin of five hundred.
    """
    present = [a for a in arrays if len(a)]
    pooled = np.concatenate(present) if present else np.array([0.0, 1.0])
    lo, hi = float(np.min(pooled)), float(np.max(pooled))
    if floor is not None and hi - lo < floor:
        mid = 0.5 * (lo + hi)
        lo, hi = mid - floor / 2, mid + floor / 2
    m = (hi - lo) * pad
    return lo - m, hi + m


def tick_edges(lims, tick_ms, max_bins=70):
    """Bin edges straddling the tick grid, so each bar is a whole number of ticks.

    Data quantised to a tick and binned on an arbitrary linear grid produces a comb: some
    bins swallow two tick values and their neighbours none. Snapping the edges to
    half-tick offsets makes every bar a real count.
    """
    span = lims[1] - lims[0]
    if tick_ms <= 0 or span / tick_ms > max_bins:
        return np.linspace(lims[0], lims[1], max_bins + 1)
    step = tick_ms * max(1, int(np.ceil(span / tick_ms / max_bins)))
    first = np.floor(lims[0] / step) * step - step / 2
    return np.arange(first, lims[1] + step, step)


def marginal_hist(ax, series, lims, t, edges, log=False):
    """Horizontal histogram sharing the time-series y axis: literally its marginal."""
    for label, values, colour in series:
        if not len(values):
            continue
        ax.hist(values, bins=edges, orientation="horizontal", color=colour,
                alpha=0.6 if len(series) > 1 else 0.85, label=label,
                edgecolor=t["surface"], linewidth=0.4, log=log)
    ax.set_ylim(*lims)
    if log:
        ax.set_xlim(left=0.7)
    tidy(ax, t, xlabel="count (log scale)" if log else "count", title="distribution")
    ax.tick_params(labelleft=False)
    ax.grid(axis="x", alpha=0.6)


def panels(theme, suptitle, subtitle):
    """The shared two-panel frame: wide time series, narrow marginal histogram.

    Title, context line and legend stack above the axes in that order and nothing is
    written inside the plot area except the marks and their direct labels -- the earlier
    arrangement put an axes title and a corner note in the frame, where both landed on top
    of the data as soon as the y-scale changed.
    """
    t = style(theme)
    fig = plt.figure(figsize=(12.5, 5.1))
    gs = GridSpec(1, 2, width_ratios=[3.3, 1.0], wspace=0.06, figure=fig,
                  left=0.07, right=0.985, top=0.76, bottom=0.16)
    ax = fig.add_subplot(gs[0])
    axh = fig.add_subplot(gs[1], sharey=ax)
    fig.suptitle(suptitle, x=0.07, ha="left", fontsize=13, color=t["primary"], y=0.985)
    fig.text(0.07, 0.905, subtitle, ha="left", fontsize=8.5, color=t["muted"])
    return t, fig, ax, axh


def footnote(fig, t, text):
    fig.text(0.07, 0.015, text, ha="left", fontsize=8, color=t["muted"])


def figure_intervals(res, out_path, theme, subtitle):
    """Figure 1 -- how steady each camera's own strobe interval is."""
    t, fig, ax, axh = panels(theme, "1. Sampling consistency within each camera", subtitle)
    c5, c6 = t["series"]
    t0, tick = res["t0"], res["tick_ms"]

    lims = full_limits(res["isi5"], res["isi6"], floor=8 * tick)
    # Markers, not a line. The values are quantised to one tick, so a connecting line
    # draws a vertical stroke between two indistinguishable values and the plot fills with
    # spurious full-height spikes. Unconnected marks show the band and the outliers.
    for name, times, isi, colour in (
            ("TTL5 (Ch15)", res["isi5_t"], res["isi5"], c5),
            ("TTL6 (Ch16)", res["isi6_t"], res["isi6"], c6)):
        ax.plot(times - t0, isi, linestyle="none", marker="o", markersize=3.2,
                markeredgewidth=0, color=colour, alpha=0.65, label=name)

    for isi, colour in ((res["isi5"], c5), (res["isi6"], c6)):
        if len(isi):
            ax.axhline(np.median(isi), color=colour, lw=1.0, ls=(0, (4, 3)), alpha=0.5)

    # Ring the anomalies: at this y-scale they are single marks among a thousand.
    for times, isi, colour in ((res["anom5_t"], res["anom5"], c5),
                               (res["anom6_t"], res["anom6"], c6)):
        for x, y in zip(times - t0, isi):
            ax.plot([x], [y], marker="o", markersize=9, markerfacecolor="none",
                    markeredgecolor=colour, markeredgewidth=1.4)
            ax.annotate(f"{y:.2f} ms at {x:.2f} s", xy=(x, y), xytext=(-8, 10),
                        textcoords="offset points", ha="right", fontsize=8, color=colour)

    ax.set_ylim(*lims)
    ax.set_xlim(0, res["t1"] - t0)
    ax.grid(axis="y", alpha=0.6)
    tidy(ax, t, xlabel=f"time from the start of the analysis window at {t0:.3f} s  (s)",
         ylabel="interval between consecutive spikes (ms)")
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.01), ncol=2, labelcolor=t["secondary"])

    steady = sum(int(np.sum(np.abs(isi - np.median(isi)) <= res["near_tol_ms"]))
                 for isi in (res["isi5"], res["isi6"]))
    total = len(res["isi5"]) + len(res["isi6"])
    rest = total - steady
    footnote(fig, t,
             f"Strobe interval across the analysis window.  {steady} of {total} intervals "
             f"lie within {res['near_tol_ms']:.2f} ms of their channel median "
             f"({100.0 * steady / total:.2f}%)"
             + (f"; the {rest} ringed mark{'s' if rest > 1 else ''} "
                f"{'are' if rest > 1 else 'is'} the rest.  " if rest else ".  ")
             + f"File tick {tick:.2f} ms, so the flat band is as flat as this file "
               f"can show.")

    marginal_hist(axh, [("TTL5", res["isi5"], c5), ("TTL6", res["isi6"], c6)],
                  lims, t, tick_edges(lims, tick), log=True)
    m5 = np.median(res["isi5"]) if len(res["isi5"]) else None
    m6 = np.median(res["isi6"]) if len(res["isi6"]) else None
    if m5 is not None and m6 is not None and abs(m5 - m6) < tick / 2:
        axh.annotate(f"both medians {m5:.2f} ms", xy=(0.5, m5), xytext=(0, 9),
                     xycoords=("axes fraction", "data"), textcoords="offset points",
                     ha="center", va="bottom", fontsize=8, color=t["secondary"])
    else:
        for m, colour, dy, va in ((m5, c5, 9, "bottom"), (m6, c6, -9, "top")):
            if m is not None:
                axh.annotate(f"median {m:.2f} ms", xy=(0.5, m), xytext=(0, dy),
                             xycoords=("axes fraction", "data"),
                             textcoords="offset points", ha="center", va=va,
                             fontsize=8, color=colour)

    fig.savefig(out_path, bbox_inches="tight")
    return fig


def figure_skew(res, out_path, theme, subtitle):
    """Figure 2 -- how simultaneous the two cameras are with each other."""
    t, fig, ax, axh = panels(theme, "2. Simultaneity between the two cameras", subtitle)
    c5, c6, pair = t["series"][0], t["series"][1], t["pair"]
    t0, dt, tick = res["t0"], res["skew"], res["tick_ms"]

    lims = full_limits(dt, floor=8 * tick)
    ax.plot(res["skew_t"] - t0, dt, lw=1.3, color=pair, alpha=0.9, label="pair skew")
    if len(dt):
        med = float(np.median(dt))
        ax.axhline(med, color=pair, lw=1.0, ls=(0, (4, 3)), alpha=0.55)
        ax.annotate(f"median {med:+.3f} ms", xy=(0.004, med),
                    xycoords=("axes fraction", "data"), va="bottom", fontsize=8,
                    color=pair)

    # The trend line is the point of this panel as much as the scatter: a steady ramp is a
    # clock-rate difference, which the standard deviation alone would misreport as jitter.
    d = res["drift"]
    if d:
        x = res["skew_t"] - res["skew_t"][0]          # the origin drift_fit regressed on
        ax.plot(res["skew_t"] - t0, d["intercept_ms"] + d["slope_ms_per_s"] * x,
                lw=1.4, color=t["primary"],
                alpha=0.55, ls=(0, (6, 3)),
                label=f"drift {d['slope_ms_per_s'] * 1e3:+.1f} us/s  ({d['ppm']:+.1f} ppm)")

    # Unmatched spikes as a rug along the foot, in their own channel's colour: they are the
    # whole reason the two counts differ, and where they sit in time says whether they are
    # teardown overhang or a pulse dropped mid-run.
    y = lims[0] + (lims[1] - lims[0]) * 0.04
    for times, colour, label in ((res["unmatched5"], c5, "TTL5 unmatched"),
                                 (res["unmatched6"], c6, "TTL6 unmatched")):
        if len(times):
            ax.plot(times - t0, np.full(len(times), y), linestyle="none", marker="|",
                    markersize=10, markeredgewidth=1.3, color=colour,
                    label=f"{label} ({len(times)})")

    ax.set_ylim(*lims)
    ax.set_xlim(0, res["t1"] - t0)
    ax.grid(axis="y", alpha=0.6)
    tidy(ax, t, xlabel=f"time from the start of the analysis window at {t0:.3f} s  (s)",
         ylabel="TTL5 - TTL6 skew within a matched pair (ms)")
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.01), ncol=4, labelcolor=t["secondary"])

    frame_ms = res["s_isi5"]["median"] or 1.0
    note = (f"Skew between paired strobes, in time order.  {len(dt)} pairs matched by "
            f"monotonic nearest-neighbour within {res['tol_ms']:.2f} ms; "
            f"|skew| never exceeds {np.max(np.abs(dt)):.3f} ms "
            f"({np.max(np.abs(dt)) / frame_ms * 100:.1f}% of a {frame_ms:.2f} ms frame "
            f"interval).")
    if d:
        resid = dt - (d["intercept_ms"] + d["slope_ms_per_s"] * x)
        note += (f"  The ramp is a clock-rate difference, not jitter: about the trend the "
                 f"sd is {np.std(resid, ddof=1):.4f} ms against {np.std(dt, ddof=1):.4f} "
                 f"ms raw.")
    footnote(fig, t, note)

    marginal_hist(axh, [("skew", dt, pair)], lims, t, tick_edges(lims, tick))
    fig.savefig(out_path, bbox_inches="tight")
    return fig


# --- console report ---------------------------------------------------------------------

def row(label, s):
    if not s["n"]:
        return f"  {label:<24} (no data)"
    return (f"  {label:<24} {s['n']:>7} {s['mean']:>10.4f} {s['median']:>10.4f} "
            f"{s['sd']:>10.4f} {s['min']:>10.4f} {s['max']:>10.4f}")


def header(label):
    return f"\n-- {label} " + "-" * max(3, 84 - len(label))


def report_anchor(anchors, bracket, frames, assignment, res):
    """The anchor section: what was found, whether it verified, and the frame accounting."""
    print(header("strobe anchor"))
    if not anchors or not any(anchors.values()):
        print(f"  none usable    {bracket['reason'] or 'no anchor detected'}")
        print(f"  bracket        {bracket['label']}")
        return

    ttl5_cam = assignment[0]
    print(f"  {'chan':<6}{'start ON':>11}{'start OFF':>11}{'end ON':>11}{'end OFF':>11}"
          f"{'P':>8}{'d':>7}{'s0':>7}")
    for nm in ("TTL5", "TTL6"):
        a = anchors.get(nm)
        if not a:
            print(f"  {nm:<6}   not found")
            continue
        t, g = a["t"], a["gaps_ms"]
        print(f"  {nm:<6}{t[0]:>11.5f}{t[1]:>11.5f}{t[2]:>11.5f}{t[3]:>11.5f}"
              f"{a['P_ms']:>8.2f}{a['d_ms']:>7.2f}"
              f"{(a['s0'] if a['s0'] is not None else 0):>7}")
        print(f"  {'':<6}{'gaps ms':>11}" + "".join(f"{x:>11.2f}" for x in g))
    print("  ON = P+d (the period stretches once); OFF = 2P-d (a frame whose strobe never "
          "fired)")

    print("  checks (strobe-intrinsic; these gate the bracket):")
    for nm in ("TTL5", "TTL6"):
        a = anchors.get(nm)
        if not a:
            continue
        for name, (passed, detail) in a["checks"].items():
            print(f"    {'ok  ' if passed else 'FAIL'}  {nm}  {name:<34} {detail}")

    print("  measured (reported, does NOT gate -- see anchor_map on why):")
    for nm in ("TTL5", "TTL6"):
        a = anchors.get(nm)
        if not a:
            continue
        for name, detail in a["measured"].items():
            print(f"          {nm}  {name:<34} {detail}")

    print(f"  bracket        {bracket['label']}")
    if bracket["reason"] and not bracket["anchor_ok"]:
        print(f"  fell back      {bracket['reason']}")

    # frames = strobes + 2 (one strobe lost at each OFF transition). Reported, not asserted.
    print("  frame accounting (frames = strobes + 2, one lost at each OFF transition):")
    for nm in ("TTL5", "TTL6"):
        a = anchors.get(nm)
        if not a or ttl5_cam is None:
            continue
        cam = ttl5_cam if nm == "TTL5" else 3 - ttl5_cam
        fr = next((x for x in frames if x["cam"] == cam), None)
        w = res["w5"] if nm == "TTL5" else res["w6"]
        n_off = a["spikes_between_off"]
        print(f"    {nm} -> c{cam}.bin: {n_off} spikes between the OFF gaps, "
              f"+1 for the enclosed OFF gap = {n_off + 1} frames spanned")
        if fr is not None and a["s0"] is not None:
            print(f"      whole run: {fr['actual']} frames stored, "
                  f"{len(w)} spikes in the analysis window; "
                  f"frame k <-> spike {a['s0']} + k, less 1 per OFF gap passed")


def report(res, channels, frames, assignment, open_note, smrx, anchors=None, bracket=None):
    t0, t1 = res["t0"], res["t1"]
    stat_head = (f"  {'':<24} {'n':>7} {'mean':>10} {'median':>10} {'sd':>10} "
                 f"{'min':>10} {'max':>10}")

    print()
    print("=" * 88)
    print("SPIKE TIMING ACCURACY")
    print("=" * 88)
    print(f"  file          {smrx}")
    print(f"  access        {open_note}")

    print(header("channels found in the .smrx"))
    print(f"  {'chan':>5}  {'title':<12} {'type':<11} {'n':>8}  {'max t (s)':>10}  comment")
    for c in channels:
        print(f"  Ch{c['chan']:<3}  {c['title']:<12} {c['type']:<11} {c['n']:>8}  "
              f"{c['max_time_s']:>10.4f}  {c['comment']}")

    print(header("analysis window"))
    print(f"  source        {bracket['label'] if bracket else 'Keyboard marks'}")
    print(f"  markers       {res['window_label']}")
    print(f"  window        {t0:.6f} s  ->  {t1:.6f} s")
    print(f"  duration      {t1 - t0:.6f} s")

    print(header("frames on disk"))
    print(f"  {'file':<9} {'frames':>8} {'declared':>9}  {'geometry':<18} path")
    for fr in frames:
        flag = "" if fr["actual"] == fr["declared"] else "   <-- MISMATCH"
        if fr["partial"]:
            flag += f"   <-- {fr['partial']} trailing bytes"
        geom = f"{fr['width']}x{fr['height']} {fr['dtype']}"
        print(f"  c{fr['cam']}.bin{'':<3} {fr['actual']:>8} {fr['declared']:>9}  "
              f"{geom:<18} {fr['path']}{flag}")

    print(header("strobe counts inside the window"))
    print(f"  TTL5 (Ch15)   {len(res['w5']):>6} spikes   "
          f"{res['w5'][0]:.4f} s -> {res['w5'][-1]:.4f} s   "
          f"(spans {res['w5'][-1] - res['w5'][0]:.4f} s)")
    print(f"  TTL6 (Ch16)   {len(res['w6']):>6} spikes   "
          f"{res['w6'][0]:.4f} s -> {res['w6'][-1]:.4f} s   "
          f"(spans {res['w6'][-1] - res['w6'][0]:.4f} s)")
    print(f"  difference    {len(res['w5']) - len(res['w6']):>+6} (TTL5 - TTL6)")

    ttl5_cam, reason, corroboration = assignment
    print(header("which .bin belongs to which channel"))
    if ttl5_cam is None:
        print(f"  UNDECIDED     {reason}")
        print( "                pass --assign ttl5-c1 or --assign ttl5-c2 to force it")
    else:
        print(f"  TTL5 (Ch15) <-> c{ttl5_cam}.bin       "
              f"TTL6 (Ch16) <-> c{3 - ttl5_cam}.bin")
        print(f"  because       {reason}")
    print(f"  corroborate   {corroboration}")
    if ttl5_cam is not None:
        for fr in frames:
            chan = "TTL5" if fr["cam"] == ttl5_cam else "TTL6"
            spikes = len(res["w5"]) if fr["cam"] == ttl5_cam else len(res["w6"])
            # Deliberately NOT called a frame/strobe mismatch. The analysis window is a
            # sub-span of the recording -- always for the anchor bracket, which stops short
            # of both anchors, and by a different amount for the Keyboard marks, which
            # overrun it. The real reconciliation is in the anchor section, which compares
            # like with like; this line is only here to show the scale of each.
            print(f"  c{fr['cam']}.bin        {fr['actual']:>5} frames stored in total, "
                  f"{spikes:>5} {chan} spikes in the analysis window "
                  f"({spikes - fr['actual']:+d}; the window is not the whole recording)")

    if bracket is not None:
        report_anchor(anchors, bracket, frames, assignment, res)

    print(header("1. interval between consecutive spikes  (ms)"))
    print(stat_head)
    print(row("TTL5 (Ch15)", res["s_isi5"]))
    print(row("TTL6 (Ch16)", res["s_isi6"]))
    for name, s in (("TTL5", res["s_isi5"]), ("TTL6", res["s_isi6"])):
        if s["n"]:
            print(f"  {name} rate {1000.0 / s['median']:.4f} Hz from the median interval;"
                  f"   1-99 pct {s['p1']:.4f} - {s['p99']:.4f} ms")
    print(f"  file tick     {res['tick_ms']:.4f} ms -- the finest interval difference "
          f"this file can express")
    for name, at, vals in (("TTL5", res["anom5_t"], res["anom5"]),
                           ("TTL6", res["anom6_t"], res["anom6"])):
        n_ok = (len(res["isi5"]) if name == "TTL5" else len(res["isi6"])) - len(vals)
        total = n_ok + len(vals)
        print(f"  {name} steady    {n_ok} of {total} intervals within "
              f"{res['near_tol_ms']:.2f} ms of the median "
              f"({100.0 * n_ok / total:.2f}%); {len(vals)} anomal"
              f"{'y' if len(vals) == 1 else 'ies'}"
              + (":  " + ", ".join(f"{v:.2f} ms at {x - t0:.3f} s"
                                   for x, v in zip(at[:4], vals[:4]))
                 + (" ..." if len(vals) > 4 else "") if len(vals) else ""))

    print(header("2. skew between paired spikes, TTL5 - TTL6  (ms)"))
    print(stat_head)
    print(row("matched pairs", res["s_skew"]))
    s = res["s_skew"]
    if s["n"]:
        print(f"  1-99 pct {s['p1']:.4f} - {s['p99']:.4f} ms;   "
              f"largest |skew| {max(abs(s['min']), abs(s['max'])):.4f} ms   "
              f"({max(abs(s['min']), abs(s['max'])) / res['s_isi5']['median'] * 100:.2f}% "
              f"of a frame interval)")
    print(f"  pairing       {s['n']} pairs matched within {res['tol_ms']:.2f} ms "
          f"(half the median interval)")
    d = res["drift"]
    if d:
        print(f"  drift         {d['slope_ms_per_s'] * 1e3:+.2f} us/s "
              f"({d['ppm']:+.1f} ppm), i.e. {d['span_ms']:+.4f} ms across the window -- "
              f"a clock-rate difference, not jitter")
        resid = res["skew"] - (d["intercept_ms"] + d["slope_ms_per_s"]
                               * (res["skew_t"] - res["skew_t"][0]))
        print(f"  jitter        sd about that trend is {np.std(resid, ddof=1):.4f} ms, "
              f"against {s['sd']:.4f} ms raw")
    for name, un in (("TTL5", res["unmatched5"]), ("TTL6", res["unmatched6"])):
        line = f"  {name} unmatched  {len(un):>5}"
        if len(un):
            shown = ", ".join(f"{x - t0:.3f}" for x in un[:6])
            line += (f"   at {shown}{' ...' if len(un) > 6 else ''} s into the window")
        print(line)
    print("=" * 88)


# --- driver ------------------------------------------------------------------------------

def choose_bracket(kb, anchors, mode):
    """Pick the analysis window: the anchor OFF..OFF span, or the Keyboard marks.

    The anchor is the better window when it is trustworthy -- it is cut from the strobe
    train itself, so it starts and ends on a real camera event rather than on a software
    mark that brackets rather more than the recording. But a mis-detected anchor would
    silently analyse the wrong span, which is worse than a loose window, so `auto` demands
    that BOTH channels found an anchor and that every check passed. Anything less falls
    back and says why.

    Returns (t0, t1, window_label, bracket) where bracket carries the reason for the log.
    """
    kb_t0, kb_t1, kb_label = keyboard_window(kb)
    keyboard = {"kind": "keyboard", "label": "window between the Keyboard marks",
                "reason": "", "anchor_ok": False}
    if mode == "keyboard":
        keyboard["reason"] = "forced with --bracket keyboard"
        return kb_t0, kb_t1, kb_label, keyboard

    have = {nm: a for nm, a in anchors.items() if a}
    failures = []
    if len(have) < 2:
        failures.append(f"anchor found on {len(have)} of 2 channels")
    for nm, a in sorted(have.items()):
        for name, (passed, detail) in a["checks"].items():
            if not passed:
                failures.append(f"{nm}: {name} -- {detail}")

    if failures and mode == "auto":
        keyboard["reason"] = "; ".join(failures)
        return kb_t0, kb_t1, kb_label, keyboard
    if failures and mode == "anchor":
        print("  WARNING: --bracket anchor was forced, but the anchor did not verify:",
              file=sys.stderr)
        for why in failures:
            print(f"    {why}", file=sys.stderr)
        if len(have) < 2:
            print("    falling back to the Keyboard marks anyway -- there is no anchor "
                  "to use.", file=sys.stderr)
            keyboard["reason"] = "; ".join(failures)
            return kb_t0, kb_t1, kb_label, keyboard

    # Strictly INSIDE the anchors, and inside them on BOTH channels: start at the first
    # spike after the later start-OFF gap, end at the last spike before the earlier end-ON
    # gap. Bracketing on the OFF marks themselves would leave each channel's own 38.91 ms
    # gap sitting in the window as its first interval, which is a real event but not part
    # of the steady-state timing being measured -- it alone lifts the ISI sd from ~0.02 ms
    # to ~0.70 ms and drags the mean off the median. The frame accounting is unaffected:
    # that is done per channel from the anchor indices, in report_anchor.
    t0 = max(a["inner"][0] for a in have.values())
    t1 = min(a["inner"][1] for a in have.values())
    label = ("strobe-anchor bracket (inside both anchors)"
             + ("" if not failures else " (FORCED -- checks failed)"))
    return t0, t1, f"anchor-bracketed ({t1 - t0:.4f} s)", {
        "kind": "anchor", "label": label, "anchor_ok": not failures,
        "reason": "; ".join(failures) if failures else "all checks passed",
    }


def resolve_paths(args):
    """Work out (smrx, session_dir) from whatever the operator pointed at.

    The .smrx is the primary thing to name: it is the file being analysed, and it is the
    one whose name the operator controls in Spike2. Everything else -- frames/,
    session.json -- sits beside it in the session directory, so that can be derived rather
    than typed. Naming the DIRECTORY and deriving the file cannot work: it needs a rule for
    the filename, and any such rule breaks the moment the Spike2 file is renamed.

    A directory is still accepted, in which case the single .smrx inside it is used. Two or
    more is ambiguous and says so rather than picking one.
    """
    if args.smrx:                                   # explicit flag always wins
        smrx = Path(args.smrx).expanduser().resolve()
    elif args.path:
        target = Path(args.path).expanduser().resolve()
        if target.is_dir():
            found = sorted(target.glob("*.smrx"))
            if not found:
                raise FileNotFoundError(f"no .smrx in {target}. Name the file directly.")
            if len(found) > 1:
                raise ValueError(f"{len(found)} .smrx files in {target}: "
                                 f"{', '.join(p.name for p in found)}. Name one directly.")
            smrx = found[0]
        else:
            smrx = target
    else:
        raise ValueError("nothing to analyse: give the path to a .smrx file "
                         "(or a session directory containing one)")

    session_dir = (Path(args.session_dir).expanduser().resolve() if args.session_dir
                   else smrx.parent)
    return smrx, session_dir


def main(argv=None):
    args = parse_args(argv)
    smrx, session_dir = resolve_paths(args)
    frames_dir = Path(args.frames_dir) if args.frames_dir else session_dir / "frames"
    out_dir = Path(args.out) if args.out else session_dir / "analysis"
    out_dir.mkdir(parents=True, exist_ok=True)

    f, open_note = open_sonfile(smrx)
    channels = read_all_channels(f)

    ttl5 = find_channel(channels, args.ttl5_chan, "TTL5")
    ttl6 = find_channel(channels, args.ttl6_chan, "TTL6")
    kb = find_channel(channels, args.keyboard_chan, "Keyboard")

    tick_ms = f.GetTimeBase() * 1e3
    frames = [read_frame_count(frames_dir, 1), read_frame_count(frames_dir, 2)]
    sess = read_session(session_dir)
    if args.assign == "auto":
        assignment = infer_assignment(ttl5, ttl6, sess)
    else:
        cam = 1 if args.assign == "ttl5-c1" else 2
        assignment = (cam, f"forced with --assign {args.assign}",
                      f"Ch15 comment {ttl5['comment']!r}; "
                      f"Ch16 comment {ttl6['comment']!r}")

    # The anchor needs to know which .bin a channel feeds, because the bright frames it
    # checks itself against are recorded per stored camera, not per TTL channel.
    ttl5_cam = assignment[0]
    anchors = {}
    if ttl5_cam is not None:
        for ch, nm, cam in ((ttl5, "TTL5", ttl5_cam), (ttl6, "TTL6", 3 - ttl5_cam)):
            anchors[nm] = anchor_for(ch, sess, cam, tick_ms)

    t0, t1, window_label, bracket = choose_bracket(kb, anchors, args.bracket)
    res = analyse(ttl5, ttl6, t0, t1, window_label, tick_ms)

    report(res, channels, frames, assignment, open_note, smrx, anchors, bracket)

    if ttl5_cam is None:
        binmap = "TTL5/TTL6 to c1.bin/c2.bin assignment undetermined"
    else:
        binmap = f"TTL5 = c{ttl5_cam}.bin, TTL6 = c{3 - ttl5_cam}.bin"
    subtitle = (f"{smrx.name}   |   {bracket['label']} {t0:.3f}-{t1:.3f} s "
                f"({t1 - t0:.3f} s)   |   {binmap}")

    p1 = out_dir / "fig1_strobe_intervals.png"
    p2 = out_dir / "fig2_pair_skew.png"
    figure_intervals(res, p1, args.theme, subtitle)
    figure_skew(res, p2, args.theme, subtitle)
    print(f"\nwrote {p1}\nwrote {p2}")

    if args.show:
        plt.show()
    else:
        plt.close("all")
    return 0


def parse_args(argv):
    p = argparse.ArgumentParser(
        description="Characterise camera strobe timing in a Spike2 .smrx against the "
                    "frames the recording app stored.")
    p.add_argument("path", nargs="?", default=None,
                   help="the .smrx to analyse. A session directory is also accepted, in "
                        "which case the single .smrx inside it is used. frames/ and "
                        "session.json are read from beside the .smrx unless "
                        "--session-dir says otherwise")
    p.add_argument("--session-dir", default=None,
                   help="session directory holding frames/ and session.json "
                        "(default: the directory containing the .smrx)")
    p.add_argument("--smrx", default=None,
                   help="path to the .smrx, if it is not the positional argument")
    p.add_argument("--frames-dir", default=None,
                   help="directory holding c1.bin/c2.bin (default: <session-dir>/frames)")
    p.add_argument("--out", default=None,
                   help="where the figures go (default: <session-dir>/analysis)")

    p.add_argument("--ttl5-chan", type=int, default=TTL5_CHAN, help="default: %(default)s")
    p.add_argument("--ttl6-chan", type=int, default=TTL6_CHAN, help="default: %(default)s")
    p.add_argument("--keyboard-chan", type=int, default=KEYBOARD_CHAN,
                   help="default: %(default)s")

    p.add_argument("--bracket", choices=("auto", "anchor", "keyboard"), default="auto",
                   help="which window to analyse: 'anchor' is the OFF..OFF strobe-anchor "
                        "span, 'keyboard' the S/s marks. 'auto' uses the anchor when both "
                        "channels find one and every check passes, else falls back to the "
                        "Keyboard marks and says why (default: %(default)s)")
    p.add_argument("--assign", choices=("auto", "ttl5-c1", "ttl5-c2"), default="auto",
                   help="which .bin TTL5 belongs to; auto infers it from the camera open "
                        "order and session.json (default: %(default)s)")
    p.add_argument("--theme", choices=tuple(THEMES), default="light",
                   help="figure colour scheme (default: %(default)s)")
    p.add_argument("--show", action="store_true", help="also open the figures in a window")
    return p.parse_args(argv)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"\nFAILED: {exc}", file=sys.stderr)
        sys.exit(1)
