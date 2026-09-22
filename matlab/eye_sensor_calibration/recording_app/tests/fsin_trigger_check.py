r"""FSIN trigger check: do the OV9281s really take their frame timing from an external pulse?

    python tests\fsin_trigger_check.py --seconds 5
    python tests\fsin_trigger_check.py --seconds 5 --mode free
    python tests\fsin_trigger_check.py --seconds 5 --devices 1 2 --camera ov9281 --out C:/Temp/test
    python "eye_sensor_calibration/recording_app/tests/fsin_trigger_check.py" --seconds 5

A bench proof of concept for the FSIN external-trigger input on two Arducam B0332 (OV9281) USB
cameras. THE PULSE TRAIN IS ALREADY RUNNING WHEN THIS STARTS. The operator brings up the ~50 Hz
train -- 20 ms period, 1 ms high, from a Power1401 digital output, divided down to about 3.3 V --
into both cameras' FSIN pins BEFORE running this. The script then opens both cameras, puts them
into external trigger mode, records for N seconds into the app's usual flat .bin files, and
reports how many frames each camera delivered, plus a PNG of the first and the last frame of each
camera so the picture itself can be looked at.

HOW TRIGGER MODE IS ENTERED ON THIS CAMERA. MEASURED ON THIS RIG, 2026-09-11, on BOTH devices,
MJPG 1280x800, with FSIN disconnected:

  1. BACKLIGHT COMPENSATION IS NOT THE SWITCH, whatever Arducam's application note's name for it
     ("low-brightness compensation") suggests. Its range here is 0..2 with default 1, and 0, 1
     and 2 all leave the camera free-running at about 121 fps; only the brightness changed. The
     first version of this script wrote that control and proved nothing: "on" is 1, which is the
     DEFAULT, so the camera never left free-run and the run looked like a trigger that failed.
  2. THE SWITCH IS IAMCameraControl PROPERTY 19, AUTO_EXPOSURE_PRIORITY -- the UVC
     CT_AE_PRIORITY_CONTROL, the "Low Light Compensation" checkbox on Windows' Camera Control
     tab, `exposure_dynamic_framerate` on Linux. See eyecal/dshow.py, AE_PRIORITY. Written 1 the
     camera delivers one frame per rising edge on FSIN; written 0 it free-runs again at once.
  3. THE WRITE GOES THROUGH A FILTER OF ITS OWN, dshow.Controls, WHILE OpenCV's capture graph is
     streaming and the reader thread is still running. Nothing is stopped, no source is
     re-applied, and the media type does not change.
  4. THE SETTING PERSISTS IN THE CAMERA across a release and a reopen. So this script writes 0
     BEFORE it opens anything -- a camera left in trigger mode by an earlier run would otherwise
     run capture.open_all's settle wait out -- and writes 0 again on the way out.
  5. IN TRIGGER MODE WITH NO PULSES, cap.read() DOES NOT FAIL. It answers ok=True about once a
     second, after the driver's ~1000 ms timeout, with an ALL-BLACK frame. Those are counted as
     `timeouts` and never stored: a real triggered frame is never all zero, because the sensor's
     black-level pedestal sits near 28 counts (see cameras.PROBE_DARK).

The mechanism itself -- the control write, its one retry, the all-zero test and the settle before
the drain -- lives in eyecal/trigger.py, which the recording app uses too. This script only drives
it and reports what the cameras then did.

THE TWO DECISIVE OBSERVATIONS:

  1. with the pulse train running, the delivered frame rate should equal the PULSE rate (~50 Hz)
     instead of the camera's own free-running rate;
  2. if the operator stops the pulse train mid-run (Spike2 key `t`), frames should stop arriving.

--mode free is the control condition: the same run with the camera left free-running, for the
rate to be compared against. It is no longer needed as a repair -- every run now writes free-run
back on the way out, and writes it again before the next open -- but it stays, because the
control condition is half the proof.

Exit codes: 0 every camera delivered frames, 1 it failed, 3 a camera delivered nothing.

A bench check, not a test. It needs two real cameras and an operator with a pulse train, so the
name carries no test_ prefix and pytest never collects it.
"""

import argparse
import json
import queue
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import cv2
import numpy as np

# eyecal/ is the recording app itself, one level up. The cameras are opened, read and stored by
# the app's own code, never by a copy of it, so this check and the app can never disagree about
# what "opened correctly" or "one stored frame" means.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from eyecal import cameras, capture, session, trigger        # noqa: E402

EXIT_OK, EXIT_FAILED, EXIT_NO_FRAMES = 0, 1, 3

# Same retry as the app's entry point, for the same narrow reason: a device a previous run has not
# been released from is the one opening failure that clears itself, given a few seconds.
OPEN_ATTEMPTS = 3
OPEN_RETRY_WAIT_S = 5.0

# How far the delivered rate may sit from the pulse rate and still be called the pulse rate. The
# two candidate answers here are far apart -- the pulse train against the camera's free-running
# rate -- so this only has to exclude a coincidence, not resolve anything fine.
RATE_TOLERANCE = 0.10

IDLE_SLEEP_S = 0.002        # when nothing was queued; a spin would steal a core from the readers
PROGRESS_PERIOD_S = 1.0


def main(argv=None):
    args = parse_args(argv)
    seconds = float(args.seconds)
    if seconds <= 0:
        raise ValueError(f"--seconds must be positive, got {seconds:g}")
    devices = [int(d) for d in args.devices]
    pulse_hz = float(args.pulse_hz)

    camera, detection = args.camera, None
    if camera == cameras.AUTO:
        print("Identifying cameras from the geometries they offer (auto)...")
        camera, detection = _retry_busy(lambda: cameras.detect_preset(devices))
    preset = cameras.CAMERA_PRESETS[camera]
    camera_note = (preset["title"] if detection is None
                   else f"{cameras.AUTO} -> {camera} ({preset['title']})")

    # The preset is used exactly as it ships. The trigger is NOT one of its source properties --
    # OpenCV cannot reach it at all -- so it is written separately, through dshow.
    want = 1 if args.mode == "trigger" else 0
    fmt = args.format or preset["format"]

    started = datetime.now(timezone.utc)
    session_id = (started.strftime("%Y%m%d-%H%M%SZ")
                  + f"_fsin_{args.mode}_{camera}_{len(devices)}cam_{seconds:g}s")
    session_dir = Path(args.out).expanduser().resolve() / session_id

    mode_note = ("trigger -- AE priority 1, one frame per rising edge on FSIN" if want
                 else "free running -- AE priority 0, the control condition")
    print(f"Mode        : {mode_note}")
    print(f"Camera      : {camera_note}")
    print(f"Format      : {fmt}")
    print(f"Devices     : {devices}")
    print(f"Duration    : {seconds:g} s")
    print(f"Pulse train : {pulse_hz:g} Hz expected on FSIN (started by the operator already)")
    print(f"Session dir : {session_dir}")

    # BEFORE the open, not after it: the setting persists in the camera, so a camera left in
    # trigger mode by an earlier run would deliver one black frame a second and nothing else, and
    # capture.open_all would sit out its whole settle timeout waiting for it to stream.
    print("Putting both cameras into free-run before opening them...")
    trigger.free_run_all(devices)
    print("Opening cameras (nothing is recorded yet)...")

    readers = open_cameras(preset, devices, args)
    try:
        print(f"Writing the trigger control (AE priority {want})...")
        ae = trigger.set_all(devices, want)
        # Before the drain, so the drain really empties the queues of everything the cameras
        # produced BEFORE the write. See trigger.MODE_SETTLE_S for the measurement.
        time.sleep(trigger.MODE_SETTLE_S)
        dropped = drop_stale(readers)
        print(f"  dropped {dropped} frame(s) per camera, queued before the mode was written")

        sess = session.Session(args.out, session_id, len(devices))
        sess.open()
        last_frames, timeouts, stopped_early, elapsed = record_fsin(readers, sess, seconds)

        # Released here rather than only in the finally, so the cameras are free before the
        # sidecars, the PNGs and the report are written -- and so close_all's verdict can go into
        # fsin_check.json, which is the only record that outlives this console.
        stuck = capture.close_all(readers)
        sess.mark_finish()
        sess.close_files()
        sess.write_sidecars()

        pngs = write_pngs(sess, last_frames)
        # The rate is frames over the window they were actually collected in. Normally that is
        # exactly --seconds; after a Ctrl-C it is the shorter run, and dividing by --seconds then
        # would understate the rate and invent a trigger failure that did not happen.
        window = elapsed if stopped_early else seconds
        stats = per_camera_stats(sess, window, timeouts)

        report(sess, stats, args.mode, pulse_hz, pngs, stuck, stopped_early)
        write_json(sess, args, camera, fmt, seconds, window, pulse_hz, want, ae, stats,
                   stuck, pngs, stopped_early, started)
        return EXIT_OK if all(s["frames"] for s in stats) else EXIT_NO_FRAMES
    finally:
        capture.close_all(readers)
        # The setting PERSISTS IN THE CAMERA, so a camera left in trigger mode delivers nothing
        # but one black frame a second to the next ordinary recording, and nothing in the
        # recording app would explain why. The restore is part of exiting, on every path.
        if not args.keep_trigger:
            print("\nrestoring free-run")
            trigger.free_run_all(devices)


def open_cameras(preset, devices, args):
    """Open every camera through the app's own path, retrying only a busy device.

    This is the ONLY place a cv2.VideoCapture is touched from the main thread. Everything after it
    happens on the reader threads, which is capture.py's rule: nothing may call into a capture
    another thread is reading.

    Every camera is put into free-run before this runs, so the open never happens in trigger mode
    and a settle timeout out of capture.open_all means exactly what it says there.
    """
    return _retry_busy(lambda: capture.open_all(preset, devices, args.format, args.exposure,
                                                None, float(args.settle_timeout)))


def drop_stale(readers):
    """Empty every reader queue and say how many frames went. See trigger.MODE_SETTLE_S.

    The frames are dropped, not stored somewhere else: they are ordinary free-running video from
    before the run, and the only thing this check does with a frame is count it.
    """
    dropped = []
    for r in readers:
        n = 0
        while True:
            try:
                r.q.get_nowait()
            except queue.Empty:
                break
            n += 1
        dropped.append(n)
    return dropped


def record_fsin(readers, sess, seconds):
    """Store what both cameras deliver for `seconds`. Returns (last, timeouts, early, span).

    Drains the reader queues IN ORDER and keeps every real frame, exactly as record.run_recording
    does -- the same storage path, without the anchor and without a preview. Neither says anything
    about the trigger, and the anchor would move the exposure while the camera is being asked to
    hold it under the trigger period.

    `timeouts` is the per-camera count of the driver's one-second black frames; see _store.
    """
    n = len(readers)
    last_frames = [None] * n
    timeouts = [0] * n
    stopped_early = False
    t0 = sess.mark_start()
    t_end = t0 + seconds
    next_progress = t0 + PROGRESS_PERIOD_S

    print(f"\nRecording {seconds:g} s... (Ctrl-C stops early and still writes everything)")
    try:
        while True:
            capture.check_readers(readers)
            got_any = False
            for k, r in enumerate(readers):
                while True:
                    try:
                        item = r.q.get_nowait()
                    except queue.Empty:
                        break
                    got_any = True
                    _store(sess, k, r, item, last_frames, timeouts)

            now = time.perf_counter()
            if now >= next_progress:
                counts = "  ".join(f"cam{k + 1} {sess.recorded[k]} ({timeouts[k]} timeouts)"
                                   for k in range(n))
                print(f"  [{now - t0:6.2f}s] {counts}")
                next_progress += PROGRESS_PERIOD_S
            if now >= t_end:
                break
            if not got_any:
                time.sleep(IDLE_SLEEP_S)
    except KeyboardInterrupt:
        stopped_early = True
        print("\nStopped early by the operator -- writing what was recorded.")

    elapsed = time.perf_counter() - t0
    # One last drain, because a frame can arrive between the final get_nowait and the deadline.
    # Anything stamped AFTER the deadline is not part of this run and would stretch the window the
    # delivered rate is computed over.
    for k, r in enumerate(readers):
        while True:
            try:
                item = r.q.get_nowait()
            except queue.Empty:
                break
            if item[0] <= t_end:
                _store(sess, k, r, item, last_frames, timeouts)
    return last_frames, timeouts, stopped_early, elapsed


def _store(sess, k, reader, item, last_frames, timeouts):
    """Store one queued frame, or count it as a timeout when it is the driver's black one.

    An ALL-ZERO frame is the ~1000 ms timeout a camera in trigger mode with no pulses hands back
    (fact 5 in the module docstring), not an exposure, and storing it would put a frame into the
    ledger that the sensor never took. trigger.is_timeout_frame is the test, subsampled and exact,
    and the reason it needs no threshold at all is written out there.
    """
    t_host, t_abs, frame = item
    if trigger.is_timeout_frame(frame):
        timeouts[k] += 1
        return
    sess.write(k, frame, t_host, t_abs, reader.q.qsize(), time.perf_counter())
    last_frames[k] = frame


def write_pngs(sess, last_frames):
    """First and last stored frame of each camera, as plain PNGs. Returns the names written.

    The point is visual confirmation: a camera delivering under the trigger should still deliver a
    normal picture, not a torn or half-exposed one. The frames come straight off the queue and are
    never reused by anything, so no copy is needed.
    """
    written = {}
    for k in range(sess.n_cams):
        first = sess.first_frames[k]
        if not sess.recorded[k] or first is None:
            print(f"  cam {k + 1}: no frames recorded, so no PNGs were written")
            continue
        names = []
        for label, frame in (("first", first), ("last", last_frames[k])):
            if frame is None:
                continue
            name = f"{label}_frame_cam{k + 1}.png"
            if not cv2.imwrite(str(sess.dir / name), frame):
                raise RuntimeError(f"cv2.imwrite failed to write {name}")
            names.append(name)
        written[f"cam{k + 1}"] = names
    return written


def per_camera_stats(sess, window, timeouts):
    """Frames, timeouts, delivered rate and arrival-interval spread per camera, from the ledger.

    The interval is the second number that separates a triggered camera from a free-running one: a
    camera taking its timing from the pulse train inherits the train's own steadiness.
    """
    stats = []
    for k in range(sess.n_cams):
        t = np.array([r[2] for r in sess.rows if r[0] == k + 1])    # t_arrive_s
        ifi = np.diff(t) * 1e3 if t.size > 1 else np.array([])
        frames = int(sess.recorded[k])
        stats.append({
            "cam": k + 1,
            "frames": frames,
            "timeouts": int(timeouts[k]),
            "rateHz": float(frames / window) if window > 0 else float("nan"),
            # Guarded rather than handed to nanmedian: a camera that delivered nothing is the
            # very result this check is looking for, and an empty-slice warning printed over its
            # own report helps nobody.
            "ifiMedianMs": float(np.median(ifi)) if ifi.size else float("nan"),
            "ifiSdMs": float(np.std(ifi)) if ifi.size else float("nan"),
            "maxQueueDepth": int(sess.max_queue[k]),
        })
    return stats


def verdict(stat, mode, pulse_hz):
    """One plain line per camera: is this the pulse rate, or is the camera running its own?"""
    if not stat["frames"]:
        if stat["timeouts"]:
            return (f"cam {stat['cam']}: no triggered frames; {stat['timeouts']} one-second "
                    f"timeouts -- no FSIN pulses reached this camera")
        return f"cam {stat['cam']}: no frames at all -- this camera delivered nothing"
    rate = stat["rateHz"]
    if mode != "trigger":
        return f"cam {stat['cam']}: rate {rate:.1f} fps -- free running, no trigger was requested"
    if pulse_hz > 0 and abs(rate - pulse_hz) <= RATE_TOLERANCE * pulse_hz:
        return (f"cam {stat['cam']}: rate {rate:.1f} fps -- consistent with the {pulse_hz:g} Hz "
                f"FSIN train")
    return (f"cam {stat['cam']}: rate {rate:.1f} fps -- NOT the pulse rate; the camera is "
            f"probably free-running")


def report(sess, stats, mode, pulse_hz, pngs, stuck, stopped_early):
    print("\n" + "=" * 78)
    print(f"FSIN check: {sess.session_id}")
    print(f"  {sess.dir}")
    if stopped_early:
        print("  stopped early by the operator")
    print("-" * 78)
    print(f"{'cam':>3} {'frames':>8} {'timeouts':>9} {'fps':>8} {'IFI med':>9} {'IFI SD':>8} "
          f"{'peak q':>8}")
    for s in stats:
        print(f"{s['cam']:>3} {s['frames']:>8} {s['timeouts']:>9} {s['rateHz']:>8.2f} "
              f"{s['ifiMedianMs']:>8.2f}m {s['ifiSdMs']:>7.2f}m {s['maxQueueDepth']:>8}")
    print("-" * 78)
    for s in stats:
        print(f"  {verdict(s, mode, pulse_hz)}")
    for cam in sorted(pngs):
        for name in pngs[cam]:
            print(f"  {cam}: {sess.dir / name}")
    if stuck:
        print(f"  teardown: {', '.join(stuck)} was still reading when its capture was released. "
              f"THIS RUN'S FRAMES ARE FINE; the NEXT open may fail, and waiting a few seconds "
              f"clears it.")
    print("=" * 78)


def write_json(sess, args, camera, fmt, seconds, window, pulse_hz, want, ae, stats, stuck,
               pngs, stopped_early, started):
    """Everything the report printed, in the session directory, because the console will go."""
    info = {
        "sessionId": sess.session_id,
        "sessionDir": str(sess.dir),
        "startedUtc": started.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "mode": args.mode,
        "aePriorityRequested": want,
        "devices": [int(d) for d in args.devices],
        "camera": camera,
        "format": fmt,
        "seconds": float(seconds),
        "rateWindowS": float(window),
        "pulseHz": float(pulse_hz),
        "stoppedEarly": bool(stopped_early),
        "aePriorityReadback": ae,
        "perCamera": [dict(s, verdict=verdict(s, args.mode, pulse_hz)) for s in stats],
        "stuckReaders": list(stuck),
        "png": pngs,
    }
    with open(sess.dir / "fsin_check.json", "w", encoding="utf-8") as fh:
        json.dump(info, fh, indent=2)
    return info


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
    p = argparse.ArgumentParser(
        description="Bench check of the FSIN external trigger on the Arducam OV9281 cameras. "
                    "Start the pulse train BEFORE running this.")
    p.add_argument("--seconds", type=float, default=5.0, help="default: %(default)s")
    p.add_argument("--devices", type=int, nargs="+", default=[1, 2],
                   help="1-based device ids (default: 1 2)")
    p.add_argument("--camera", default="ov9281",
                   choices=sorted(cameras.CAMERA_PRESETS) + [cameras.AUTO],
                   help="preset name, or auto to identify it from the geometries offered "
                        "(default: %(default)s)")
    p.add_argument("--mode", choices=("trigger", "free"), default="trigger",
                   help="trigger writes AE priority 1, which is this camera's external-trigger "
                        "switch; free writes 0, which is the control condition "
                        "(default: %(default)s)")
    p.add_argument("--exposure", type=float, default=None,
                   help="override the preset exposure; it must stay shorter than the trigger "
                        "period")
    p.add_argument("--format", default=None, help="e.g. MJPG_1280x800")
    p.add_argument("--out", default="C:/Temp/test",
                   help="parent directory for the session (default: %(default)s)")
    p.add_argument("--pulse-hz", type=float, default=50.0,
                   help="the FSIN pulse rate the operator is running, for the verdict "
                        "(default: %(default)s)")
    p.add_argument("--settle-timeout", type=float, default=30.0,
                   help="how long a camera may deliver nothing before the open gives up "
                        "(default: %(default)s)")
    p.add_argument("--keep-trigger", action="store_true",
                   help="do NOT restore free-run on exit; a camera left in trigger mode gives "
                        "the next ordinary recording nothing but a black frame a second")
    return p.parse_args(argv)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\nInterrupted.")
        sys.exit(EXIT_FAILED)
    except Exception as exc:
        print(f"\nFAILED: {exc}", file=sys.stderr)
        sys.exit(EXIT_FAILED)
