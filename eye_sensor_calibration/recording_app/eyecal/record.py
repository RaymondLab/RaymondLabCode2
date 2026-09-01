"""The recording phase: store every frame, and bracket the run with a strobe anchor.

Runs on cameras that are ALREADY open and streaming -- alignment handed them over without a
release, so there is no second start-up and no second burst of settling pulses in the Spike2
file.

THE ANCHOR. Each camera's strobe pin is wired to a Spike2 digital input, so Spike2 timestamps
every exposure on its own timebase. What Spike2 cannot see is which pulse is stored frame 0: the
camera has been free-running since long before recording began, and the pulse train carries no
edge at the moment storage starts.

So the recording marks itself. At the start and again at the end, the exposure is briefly raised
for a fraction of a second. Two things happen, and either one is enough:

    in the Spike2 file   the strobe pulses go conspicuously WIDE (the strobe tracks exposure)
    in the stored data   those same frames come out saturated, and their indices are recorded

The bright frames are STORED, not discarded, so the anchor is a frame index on one side and a
pulse index on the other, with nothing counted or assumed in between. Two anchors also bracket
the run: the number of pulses between them must equal the number of stored frames between them,
which is a complete drop check that no queue-depth heuristic can give you.

Whether writing CAP_PROP_EXPOSURE also restarts the graph (a GAP in the strobe instead of a wide
pulse) is a property of the driver. Either way it is a marker. cameras.secure_manual_exposure proves
at start-up that it at least does not re-negotiate the media type.
"""

import collections
import queue
import time

import cv2

from . import capture, display

Anchor = collections.namedtuple("Anchor", "enabled normal bright hold_s")

# Frames are sampled for brightness over a window twice the hold, so the frames on either side of
# the transition are measured too and the threshold has something dark to sit against.
_SAMPLE_FACTOR = 2


def make_anchor(enabled, normal_exposure, bright_exposure=None, hold_s=0.25):
    """Default the bright exposure to four stops up from normal.

    DirectShow exposure is log2 seconds, so +4 is 16x the integration time: at -10 (about 1 ms)
    that gives about 16 ms, long enough to saturate and to widen the strobe obviously, while
    still fitting inside a 20 ms frame period so the FRAME RATE does not change. Dropping the
    rate would work as a marker too, but it would perturb the very timing being anchored.
    """
    if not enabled or normal_exposure is None:
        return Anchor(False, normal_exposure, None, hold_s)
    bright = normal_exposure + 4 if bright_exposure is None else bright_exposure
    return Anchor(True, normal_exposure, bright, hold_s)


def run_recording(readers, devices, order, sess, seconds, rotate180=False, anchor=None,
                  preview_hz=10.0, downsample=4, fullscreen=False, window_scale=0.5):
    """Record `seconds` of frames from every camera. Returns (stop_reason, anchor_report).

    `order[cam_index]` is the reader that becomes cam 1, cam 2 ... and therefore c1.bin, c2.bin.
    """
    n = len(order)
    anchor = anchor or Anchor(False, None, None, 0.25)
    if anchor.enabled and seconds < 6 * anchor.hold_s:
        print(f"  recording is too short for both anchors ({seconds:g} s); disabling them")
        anchor = anchor._replace(enabled=False)

    # Everything queued during alignment is stale by definition. Dropping it here means stored
    # frame 0 is the first frame that ARRIVES after storage begins, not one that was already
    # waiting -- which is what makes the anchor land where it claims to.
    for r in readers:
        capture.drain_newest(r)

    sess.open()
    schedule = _anchor_schedule(anchor, seconds)
    samples = {cam: [] for cam in range(n)}     # cam index -> [(frame_idx, mean, "start"/"end")]
    latest_disp = [None] * n
    depths = [0] * n
    period = 1.0 / preview_hz if preview_hz > 0 else 0.0
    # Smaller than the alignment window, and downsampled to match. The two go together: a smaller
    # window with the alignment phase's downsample would leave the canvas LARGER than the window,
    # which flips the letterbox onto its expensive interpolation and costs more than the full-size
    # window did. Measured 1.15 ms per refresh at this pairing against 3.53 ms at the alignment
    # settings -- ergonomics first, the couple of milliseconds is incidental.
    window = display.Window(f"RECORDING -- {sess.session_id}", fullscreen, window_scale) \
        if preview_hz > 0 else None
    stop_reason = "duration reached"
    t0 = sess.mark_start()          # stamps perf_counter AND the wall clock, and records both
    next_preview_at = 0.0

    print(f"\nRecording {seconds:g} s... (q, Esc, closing the preview, or Ctrl-C stops early)")
    try:
        while True:
            capture.check_readers(readers)
            t = time.perf_counter() - t0
            if t >= seconds:
                break

            while schedule and schedule[0][0] <= t:
                _, value, label = schedule.pop(0)
                for r in readers:
                    r.request_exposure(value)
                print(f"  [{t:6.2f}s] anchor {label}: exposure -> {value:g}")

            got_any = False
            for cam, k in enumerate(order):
                r = readers[k]
                while True:
                    try:
                        t_arrive, t_abs, img = r.q.get_nowait()
                    except queue.Empty:
                        break
                    t_dequeue = time.perf_counter()     # closes the queue wait, opens the write
                    got_any = True
                    # Measured on the CONSUMER side, and it has to be: sampled where the producer
                    # puts, it would report the depth before its own frame was added and read 0
                    # at steady state, asserting headroom rather than measuring it.
                    depths[cam] = r.q.qsize() + 1
                    if rotate180:
                        img = cv2.rotate(img, cv2.ROTATE_180)
                    idx = sess.write(cam, img, t_arrive, t_abs, depths[cam], t_dequeue)
                    which = _sample_window(anchor, seconds, t_arrive - t0)
                    if which:
                        samples[cam].append((idx, float(img.mean()), which))
                    latest_disp[cam] = img

            if window is not None and t >= next_preview_at and \
                    all(f is not None for f in latest_disp):
                panes = [display.recording_pane(
                    display.display_copy(latest_disp[cam], downsample), cam + 1,
                    devices[order[cam]], sess.recorded[cam], depths[cam]) for cam in range(n)]
                canvas = display.tile(panes)
                display.text(canvas, f"REC  {t:5.1f}s / {seconds:g}s   "
                                     f"{max(seconds - t, 0.0):5.1f}s left   q or close to stop",
                             (8, canvas.shape[0] - 12), display.RED, 0.5)
                window.show(canvas)
                next_preview_at = t + period

            if window is not None:
                key = window.pump()
                if key in (ord("q"), 27):
                    stop_reason = "stopped by operator (q/Esc)"
                    break
                # AFTER the waitKey that delivers the close message -- see Window.closed.
                if window.closed():
                    stop_reason = "stopped by operator (preview window closed)"
                    break
            elif not got_any:
                # Without a window there is no waitKey to pace the loop, and a spin at full tilt
                # would steal a core from the capture threads.
                time.sleep(0.001)
    except KeyboardInterrupt:
        stop_reason = "interrupted by operator (Ctrl-C)"
        print("\nInterrupted -- closing files cleanly.")
    finally:
        # Before anything slow: the closing clock pair should bracket the frames, not the
        # teardown that follows them.
        sess.mark_finish()
        if anchor.enabled:
            for r in readers:                   # never leave a camera at the bright value
                r.request_exposure(anchor.normal)
        sess.close_files()
        if window is not None:
            window.close()

    return stop_reason, _anchor_report(anchor, samples)


def _anchor_schedule(anchor, seconds):
    """(elapsed_s, exposure, label) events, in order. Empty when the anchor is off."""
    if not anchor.enabled:
        return []
    h = anchor.hold_s
    return [(0.0, anchor.bright, "start on"),
            (h, anchor.normal, "start off"),
            (seconds - 2 * h, anchor.bright, "end on"),
            (seconds - h, anchor.normal, "end off")]


def _sample_window(anchor, seconds, elapsed):
    """"start", "end", or None. Only frames near an anchor get their mean intensity measured."""
    if not anchor.enabled:
        return None
    w = _SAMPLE_FACTOR * anchor.hold_s
    if elapsed <= w:
        return "start"
    if elapsed >= seconds - w:
        return "end"
    return None


def _anchor_report(anchor, samples):
    """Turn the sampled intensities into the frame indices that came out bright.

    Thresholded at the midpoint between the darkest and brightest sample within each window: the
    anchor is a 16x exposure change, so the two groups are nowhere near each other and nothing
    subtler is warranted. The raw samples are kept too, so the split can be redone offline
    without re-reading any frames.
    """
    if not anchor.enabled:
        return {"enabled": False}
    bright, raw = {}, {}
    for cam, rows in samples.items():
        name = f"cam{cam + 1}"
        raw[name] = [[int(i), round(m, 2), w] for i, m, w in rows]
        marks, thresholds = {}, {}
        for which in ("start", "end"):
            means = [m for _, m, w in rows if w == which]
            if not means:
                marks[which] = []
                continue
            cut = (min(means) + max(means)) / 2.0
            thresholds[which] = round(cut, 2)
            marks[which] = [int(i) for i, m, w in rows if w == which and m > cut]
        bright[name] = {"start": marks.get("start", []), "end": marks.get("end", []),
                        "thresholds": thresholds}
    return {"enabled": True, "normalExposure": anchor.normal, "brightExposure": anchor.bright,
            "holdS": anchor.hold_s, "sampleWindowS": _SAMPLE_FACTOR * anchor.hold_s,
            "brightFrames": bright, "samples": raw,
            "note": "Bright stored frames correspond to the WIDE strobe pulses on that camera's "
                    "Spike2 digital input. The pulse count between the start and end anchors "
                    "must equal the stored frame count between them."}
