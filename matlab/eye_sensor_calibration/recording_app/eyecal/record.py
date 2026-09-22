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

THE TRIGGERED PATH, at the bottom of this file, needs no anchor at all. With the cameras taking
their frame timing from Spike2's FSIN pulse train, stored frame k IS pulse k on every camera, so
there is nothing to mark and the exposure is left alone -- it has to stay well inside the pulse
period anyway. run_recording above is still the fallback, and runs exactly as it always has when
no pulses arrive. See eyecal/trigger.py for the control, the ghost frames and the measurements.
"""

import collections
import queue
import time

import cv2

from . import capture, display, trigger
from .trigger import HOLD_MS

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


# --- the triggered path: one stored frame per FSIN pulse --------------------------------------

# The stored frame rate on the preview is counted over this window, exactly as
# tools/camera_check.py counts the reader's own rate. One second at 100 Hz is 100 frames, plenty
# for a readout that updates visibly, and it costs two comparisons per redraw rather than any
# per-frame work.
RATE_WINDOW_S = 1.0

IDLE_SLEEP_S = 0.001        # when nothing was queued; a spin would steal a core from the readers


def wait_for_pulses(readers, order, wait_s):
    """Wait up to `wait_s` for the first REAL frame from EVERY camera.

    Returns (ok, pending, seen, timeouts); all but `ok` are per camera index, in `order`.

    This is the test for "is the pulse train actually running?", and it is the only test there
    is: in trigger mode a camera with no pulses does not fail, it hands back one ALL-ZERO frame
    about once a second after the driver's ~1000 ms timeout (see eyecal/trigger.py). So every
    item is classified, timeout frames are counted and dropped, and only real frames are an answer.

    EVERY REAL FRAME IS KEPT, not just the first. Between the first pulse and the moment this
    returns -- the other camera's first frame may be a pulse period later -- more pulses arrive,
    and each one is a frame the recording must store. They are handed back in `pending`, in
    arrival order, for run_triggered to write as pulse 0, 1, 2... Dropping them would silently
    shift every pulse index in the session by however many were lost.

    ON A FAILED WAIT THE REAL FRAMES COLLECTED SO FAR ARE DISCARDED, deliberately. The caller
    falls back to run_recording, which drains the queues itself and stamps its own origin, so
    frames from before that origin have no place in its ledger -- and a handful of frames from a
    train that never got going is not a recording of anything.
    """
    n = len(order)
    pending = [[] for _ in range(n)]
    seen, timeouts = [0] * n, [0] * n
    deadline = time.perf_counter() + wait_s

    while True:
        capture.check_readers(readers)
        got_any = False
        for cam, k in enumerate(order):
            r = readers[k]
            while True:
                try:
                    item = r.q.get_nowait()
                except queue.Empty:
                    break
                got_any = True
                if trigger.is_timeout_frame(item[2]):
                    timeouts[cam] += 1
                    continue
                pending[cam].append(item)
                seen[cam] += 1
        if all(seen):
            return True, pending, seen, timeouts
        if time.perf_counter() >= deadline:
            return False, pending, seen, timeouts
        if not got_any:
            time.sleep(IDLE_SLEEP_S)


def run_triggered(readers, devices, order, sess, seconds, pending, rotate180=False,
                  preview_hz=10.0, downsample=4, fullscreen=False, window_scale=0.5,
                  end_margin_s=5.0):
    """Record one stored frame per FSIN pulse. Returns (stop_reason, report).

    `pending` is what wait_for_pulses already took off the queues, per camera, in arrival order.
    Those frames are pulse 0 onwards and are written FIRST, before anything else is drained,
    through the same sess.write path and with the same rotation as run_recording.

    sess.mark_start() must already have been called -- by the caller, at the START of the wait --
    so the pending frames and everything after them share one origin. That is also why this calls
    sess.open() itself and the caller does not: nothing may be created on disk before the wait,
    because the wait is what decides whether this path runs at all.

    NO ANCHOR HERE. The anchor exists to say which pulse is stored frame 0, and under the trigger
    that is not a question: frame k is pulse k. Raising the exposure would also fight the pulse
    period, which the exposure has to stay well inside.

    THE PULSE TRAIN ENDING IS WHAT STOPS THIS, not the clock. A camera that has delivered real
    frames and then produces a timeout frame has gone a second without a pulse; once every such
    camera has, the train is over and the recording is complete. `seconds` is only the ceiling:
    `seconds + end_margin_s` stops a run whose train never ended, so a Spike2 left pulsing cannot
    record until the disk fills.
    """
    n = len(order)
    if sess.t0 is None:
        raise RuntimeError("sess.mark_start() must be called before run_triggered(); it is the "
                           "origin the pending frames were already stamped against")
    sess.open()
    t0 = sess.t0
    pending = pending or [[] for _ in range(n)]
    timeouts = [0] * n
    real_seen = [False] * n
    ended = [False] * n                 # a timeout frame since this camera's last real one
    last_real_t = [None] * n            # newest real frame, relative to t0
    latest_disp = [None] * n
    pane_shapes = [None] * n            # geometry of the display copy, for the black pane
    rate_fps = [float("nan")] * n
    rate_mark = [None] * n
    frames_during_wait = [len(items) for items in pending]

    for cam, items in enumerate(pending):
        for item in items:
            img = _store_real(sess, cam, readers[order[cam]], item, rotate180)
            real_seen[cam] = True
            last_real_t[cam] = item[0] - t0
            latest_disp[cam] = img
            pane_shapes[cam] = _pane_shape(img, downsample)

    period = 1.0 / preview_hz if preview_hz > 0 else 0.0
    window = display.Window(f"RECORDING (FSIN trigger) -- {sess.session_id}", fullscreen,
                            window_scale) if preview_hz > 0 else None
    ceiling = seconds + end_margin_s
    stop_reason = "pulse train ended (1 s without a pulse on every camera)"
    next_preview_at = 0.0

    print(f"\nRecording one frame per FSIN pulse. It stops when the train does, or at "
          f"{ceiling:g} s. (q, Esc, closing the preview, or Ctrl-C stops early)")
    print(f"  {sum(frames_during_wait)} frame(s) arrived while waiting for the first pulse and "
          f"are stored as pulse 0 onwards: {per_cam_counts(frames_during_wait)}")
    try:
        while True:
            capture.check_readers(readers)
            t = time.perf_counter() - t0

            for cam, k in enumerate(order):
                r = readers[k]
                while True:
                    try:
                        item = r.q.get_nowait()
                    except queue.Empty:
                        break
                    # The ONLY per-frame work this adds to run_recording's loop: one subsampled
                    # test for the driver's all-zero timeout frame. Storing one would put a frame
                    # in the ledger that the sensor never took and shift every pulse after it.
                    if trigger.is_timeout_frame(item[2]):
                        timeouts[cam] += 1
                        if real_seen[cam]:
                            ended[cam] = True
                        continue
                    img = _store_real(sess, cam, r, item, rotate180)
                    real_seen[cam] = True
                    ended[cam] = False
                    last_real_t[cam] = item[0] - t0
                    latest_disp[cam] = img
                    pane_shapes[cam] = _pane_shape(img, downsample)

            live = [cam for cam in range(n) if real_seen[cam]]
            if live and all(ended[cam] for cam in live):
                break
            if t >= ceiling:
                stop_reason = "ceiling reached (seconds + margin) before the pulse train ended"
                break

            if window is not None and t >= next_preview_at:
                panes = []
                for cam in range(n):
                    _update_rate(rate_mark, rate_fps, cam, t, sess.recorded[cam])
                    waiting = (last_real_t[cam] is None
                               or (t - last_real_t[cam]) * 1e3 > HOLD_MS)
                    img = (None if waiting or latest_disp[cam] is None
                           else display.display_copy(latest_disp[cam], downsample))
                    panes.append(display.trigger_pane(
                        img, cam + 1, devices[order[cam]], sess.recorded[cam], timeouts[cam],
                        rate_fps[cam], waiting, pane_shapes[cam] or _BLANK_SHAPE))
                canvas = display.tile(panes)
                display.text(canvas, f"REC  {t:5.1f}s   ceiling {ceiling:g}s   stops when the "
                                     f"pulse train does   q or close to stop",
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
            else:
                # Without a window there is no waitKey to pace the loop. Slept whether or not
                # anything was drained, unlike run_recording: under the trigger most iterations
                # drain nothing even while the train runs, and at 100 Hz a millisecond is a tenth
                # of a pulse period.
                time.sleep(IDLE_SLEEP_S)
    except KeyboardInterrupt:
        stop_reason = "interrupted by operator (Ctrl-C)"
        print("\nInterrupted -- closing files cleanly.")
    finally:
        # Before anything slow: the closing clock pair should bracket the frames, not the
        # teardown that follows them. There is no exposure to put back -- the anchor never ran.
        sess.mark_finish()
        sess.close_files()
        if window is not None:
            window.close()

    return stop_reason, {
        "timeouts": {f"cam{cam + 1}": int(timeouts[cam]) for cam in range(n)},
        "framesDuringWait": {f"cam{cam + 1}": int(frames_during_wait[cam]) for cam in range(n)},
        "lastRealFrameS": {f"cam{cam + 1}": (None if last_real_t[cam] is None
                                             else float(last_real_t[cam])) for cam in range(n)},
        "endReason": stop_reason,
    }


# Only ever drawn for a camera that has not delivered one real frame, which under the trigger
# means a pane that never showed anything. Small on purpose: it is a placeholder, and tile()
# pads the rest.
_BLANK_SHAPE = (240, 320)


def _store_real(sess, cam, reader, item, rotate180):
    """Write one real frame exactly as run_recording's inner loop does. Returns the stored image.

    The queue depth is measured HERE, on the consumer side, for the reason run_recording gives:
    sampled where the producer puts, it would report the depth before its own frame was added.
    """
    t_arrive, t_abs, img = item
    t_dequeue = time.perf_counter()     # closes the queue wait, opens the write
    if rotate180:
        img = cv2.rotate(img, cv2.ROTATE_180)
    sess.write(cam, img, t_arrive, t_abs, reader.q.qsize() + 1, t_dequeue)
    return img


def _pane_shape(img, downsample):
    """The geometry display_copy would produce, without producing it.

    Needed because a pane with no pulses is drawn BLACK at that geometry, and at exactly that
    moment there is no image to take it from. A strided slice keeps the CEILING of the division,
    so this is that and not a floor.
    """
    h, w = img.shape[:2]
    step = max(int(downsample), 1)
    return (-(-h // step), -(-w // step))


def _update_rate(mark, fps, cam, t, count):
    """Stored frames per second over the last RATE_WINDOW_S, from the ledger's own count.

    Differenced over a window rather than measured per frame, exactly as
    tools/camera_check.Pane.update_rate does: the recording loop must not grow per-frame work for
    a readout, and under the trigger this number IS the delivered pulse rate.
    """
    if mark[cam] is None:
        mark[cam] = (t, count)
    elif t - mark[cam][0] >= RATE_WINDOW_S:
        fps[cam] = (count - mark[cam][1]) / (t - mark[cam][0])
        mark[cam] = (t, count)


def per_cam_counts(counts):
    """Per-camera counts on one line -- "cam1 3, cam2 1" -- for the console and for the reason
    recorded in session.json when the trigger falls back."""
    return ", ".join(f"cam{k + 1} {n}" for k, n in enumerate(counts))
