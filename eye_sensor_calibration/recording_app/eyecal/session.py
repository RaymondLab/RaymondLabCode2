"""The session on disk: flat frame binaries, timestamps, sidecars and session.json.

    <out>/<sessionId>/
        frames/c1.bin, c1.json     flat frames + sidecar giving shape, count and dtype
        frames/c2.bin, c2.json
        first_frame_cam1.tif       first frame of each camera, one plain TIFF each,
        first_frame_cam2.tif       double-clickable from Windows Explorer
        ts.csv                     one row per frame: arrival timestamps and queue depth
        session.json               what was asked, what happened, and the strobe anchors

session.json is also the SUCCESS SIGNAL for Spike2: it is written last, only on a run that
completed, so its existence means the recording finished. See README.md.

Storage is a flat binary of uncompressed frames, deliberately. Measured on this rig, PNG
encoding costs ~1164 ms per frame against ~1 ms for a raw write; the capture queue saturates and
the storage path falls seconds behind. Convert to images offline if you need them.
"""

import json
import time
from datetime import datetime, timezone
from pathlib import Path

import cv2
import numpy as np

# One row per stored frame. Timestamps are seconds, durations milliseconds, and the suffix
# says which -- the previous schema had two columns sharing a `t_*_s` pattern while carrying
# different origins (one raw perf_counter, one relative to the start of recording), so their
# difference was meaningless and nothing in the file said so.
#
#   t_arrive_s     when the frame reached Python, from cap.read(), relative to the start of
#                  recording. session.json's clock block pins that origin to the wall clock.
#   t_abs_posix    the SAME instant on the wall clock. Kept per frame, not just as an
#                  endpoint pair: two endpoints cannot tell a clock STEP from a rate change,
#                  and only the per-frame series can. summarise fits it and stores the
#                  residual, which is the one detector for the OS adjusting the clock mid-run.
#   queue_depth    frames waiting at the moment this one was taken off the queue.
#   queue_wait_ms  arrival -> taken off the queue. The consumer falling behind.
#   write_ms       taken off the queue -> bytes written. Everything the consumer did with the
#                  frame: the rotation if one was accepted, and the file write.
#
# Durations rather than a second timestamp, deliberately. t_arrive_s + queue_wait_ms +
# write_ms reconstructs when the write finished, so nothing is lost -- and a duration has no
# epoch, so it cannot be misread the way the old t_sink_s was. Splitting the two also
# separates their causes: a busy consumer loop and a slow disk look identical when merged.
TS_COLUMNS = ["cam", "frame_idx", "t_arrive_s", "t_abs_posix", "queue_depth",
              "queue_wait_ms", "write_ms"]


class Session:
    """Holds the open output files and the per-camera ledger while a recording runs.

    The directory is NOT created until open() is called, which the recorder does only once the
    cameras are known good. A failed run should not leave an empty session directory behind for
    someone to find later and wonder about.
    """

    def __init__(self, out, session_id, n_cams):
        self.session_id = session_id
        self.dir = Path(out).expanduser().resolve() / session_id
        self.frames_dir = self.dir / "frames"
        self.n_cams = n_cams
        # (cam, idx, t_arrive_s, t_abs, depth, queue_wait_ms, write_ms)
        self.rows = []
        self.t0 = None                          # perf_counter at the start of recording
        self.clock = None                       # the perf_counter <-> wall clock bridge
        self.recorded = [0] * n_cams
        self.max_queue = [0] * n_cams
        self.bytes_out = [0] * n_cams
        self.dims = [None] * n_cams
        self.dtypes = [None] * n_cams
        self.first_frames = [None] * n_cams
        self._files = []

    def open(self):
        self.frames_dir.mkdir(parents=True, exist_ok=True)
        self._files = [open(self.frames_dir / f"c{k + 1}.bin", "wb", buffering=1 << 20)
                       for k in range(self.n_cams)]

    def mark_start(self):
        """Stamp the recording origin on both clocks and return it. Call once, at t0.

        Reading perf_counter and time.time back to back pins the arbitrary perf_counter epoch
        (typically boot) to a real moment. Without it, t_arrive_s is a number with no meaning
        outside this process -- and t0 was previously recorded nowhere at all, which is what
        made the old t_sink_s irrecoverable rather than merely wrong.
        """
        self.t0 = time.perf_counter()
        self.clock = {"t0Perf": self.t0, "t0Posix": time.time()}
        return self.t0

    def mark_finish(self):
        """The same pair again at the end, so the two clocks' rate difference is on record."""
        if self.clock is None:
            return
        self.clock["t1Perf"] = time.perf_counter()
        self.clock["t1Posix"] = time.time()
        span = self.clock["t1Perf"] - self.t0
        if span > 0:
            drift = (self.clock["t1Posix"] - self.clock["t0Posix"]) - span
            self.clock["driftPpm"] = float(drift / span * 1e6)

    def write(self, k, img, t_arrive, t_abs, depth, t_dequeue):
        """Append one frame for camera index k. Returns its 0-based frame index.

        Takes the closing timestamp ITSELF, after the bytes are down. Passing it in as an
        argument -- as this used to -- evaluates it before the call runs, so the write cost
        the caller thought it was measuring was never in the number.

        Writes the array's buffer directly rather than tobytes(): at 50 fps a full frame is a
        2 MB copy per camera per frame, and there is nothing to gain by making it.
        """
        if self.t0 is None:
            raise RuntimeError("Session.mark_start() must be called before write()")
        if self.dims[k] is None:
            self.dims[k], self.dtypes[k] = img.shape[:2], str(img.dtype)
        elif img.shape[:2] != self.dims[k]:
            raise RuntimeError(f"Frame size changed mid-run on camera {k + 1} "
                               f"({self.dims[k]} -> {img.shape[:2]}); the flat binary assumes "
                               f"fixed geometry.")

        buf = np.ascontiguousarray(img)         # no-op unless the frame was rotated
        self._files[k].write(buf.data)
        self.bytes_out[k] += buf.nbytes
        if self.first_frames[k] is None:
            self.first_frames[k] = buf.copy()   # same bytes as frame 0 in the .bin

        idx = self.recorded[k]
        self.rows.append((k + 1, idx, t_arrive - self.t0, t_abs, depth,
                          (t_dequeue - t_arrive) * 1e3,
                          (time.perf_counter() - t_dequeue) * 1e3))
        self.recorded[k] += 1
        self.max_queue[k] = max(self.max_queue[k], depth)
        return idx

    def close_files(self):
        for fh in self._files:
            try:
                fh.close()
            except Exception:
                pass
        self._files = []

    # -- finalising ----------------------------------------------------------------------

    def write_sidecars(self):
        self.rows.sort(key=lambda r: (r[0], r[1]))
        with open(self.dir / "ts.csv", "w", encoding="utf-8", newline="") as fh:
            fh.write(",".join(TS_COLUMNS) + "\n")
            for cam, idx, t_arrive, t_abs, depth, q_ms, w_ms in self.rows:
                fh.write(f"{cam},{idx},{t_arrive:.9f},{t_abs:.6f},{depth:.0f},"
                         f"{q_ms:.4f},{w_ms:.4f}\n")

        for k in range(self.n_cams):
            if self.dims[k] is None:
                continue
            h, w = self.dims[k]
            with open(self.frames_dir / f"c{k + 1}.json", "w", encoding="utf-8") as fh:
                json.dump({
                    "camera": k + 1, "count": int(self.recorded[k]),
                    "height": int(h), "width": int(w), "dtype": self.dtypes[k],
                    # "mode" is not needed by read_frames() below. It is written because the
                    # lab's existing readers -- pycamrig.frames.RunFrames and MATLAB's
                    # camrig.io.FrameReader -- key off it, and being able to open a recording
                    # with the tools already in the repo costs one line here.
                    "mode": "raw",
                    # numpy's native order. Stated rather than assumed so a reader never has to
                    # guess -- reading row-major data as column-major yields a transposed image
                    # with exactly the right byte count, which passes every other sanity check.
                    "order": "row-major",
                    "bytes": int(self.bytes_out[k]),
                    "layout": "count consecutive frames of height*width elements, no padding, "
                              "no header",
                }, fh, indent=2)

    def write_first_frame_tiffs(self):
        """ONE single-page TIFF per camera, so they preview in Explorer and open with a
        double-click. Uncompressed and lossless, so each is byte-identical to the first frame in
        the matching frames/c<K>.bin -- a faithful preview, not a re-rendering."""
        written = []
        for k, f in enumerate(self.first_frames):
            if f is None:
                continue
            name = f"first_frame_cam{k + 1}.tif"
            if not cv2.imwrite(str(self.dir / name), f):
                raise RuntimeError(f"cv2.imwrite failed to write {name}")
            written.append(name)
        return written


def clock_bridge(session):
    """Tie t_arrive_s to the wall clock, and say whether the wall clock behaved.

    Two things live here. The endpoint pair from mark_start/mark_finish converts t_arrive_s
    to a real moment, which perf_counter alone cannot do -- its epoch is arbitrary.

    The fit is the other half, and it needs the per-frame t_abs_posix rather than just those
    endpoints: two endpoints cannot distinguish the OS STEPPING the clock from the two clocks
    simply running at slightly different rates, because both look like a changed slope.
    Fitting every frame and keeping the residual does distinguish them -- a step leaves a
    residual far above the few microseconds a smooth slew leaves. Measured on this rig the
    two clocks diverge 64-177 us over a run, cleanly linear, residual under 4 us and no steps
    seen, so this is a cheap alarm rather than a correction anyone needs to apply.
    """
    info = dict(session.clock or {})
    t = np.array([r[2] for r in session.rows])          # t_arrive_s
    a = np.array([r[3] for r in session.rows])          # wall clock, same instants
    if t.size > 2:
        a0 = a - a[0]                                   # centre: 1.79e9 s is ill-conditioned
        slope, intercept = np.polyfit(t, a0, 1)
        resid = a0 - (intercept + slope * t)
        info.update({
            "absVsArriveSlope": float(slope),
            "absVsArrivePpm": float((slope - 1.0) * 1e6),
            "absVsArriveResidualUs": float(np.std(resid, ddof=1) * 1e6),
            "absVsArriveMaxStepUs": float(np.max(np.abs(np.diff(a0 - t))) * 1e6),
        })
    return info


def summarise(session, started, request, readers, devices, order, anchors, stop_reason,
              tiffs, rotate180, stuck_readers=(), trigger=None):
    """Build session.json and write it. Written LAST, because it is the success signal.

    `stuck_readers` is whatever capture.close_all reported -- see warnings_for.

    `trigger` is the FSIN trigger block, or None when the feature is off entirely. It is the one
    record of WHICH path a session took: under the trigger, stored frame k is pulse k, and a
    session that fell back to free-run is read the old way, through the strobe anchors. Nothing
    downstream can tell the two apart from the frames themselves.
    """
    per_cam = []
    for k in range(session.n_cams):
        rows = [r for r in session.rows if r[0] == k + 1]
        t = np.array([r[2] for r in rows])
        ifi = np.diff(t) * 1000.0 if t.size > 1 else np.array([np.nan])
        duration = (t[-1] - t[0]) if t.size > 1 else float("nan")
        reader = readers[order[k]]

        # Sink cost: arrival -> bytes down. p95 rather than the peak, because every run has
        # one 30-60 ms outlier and the peak therefore separates a healthy session from a
        # struggling one by under 2x, while p95 separates them by nearly 8x. slowFrames is
        # sharper still, and is what warnings_for actually gates on.
        q_ms = np.array([r[5] for r in rows]) if rows else np.array([np.nan])
        w_ms = np.array([r[6] for r in rows]) if rows else np.array([np.nan])
        sink = q_ms + w_ms
        slow_over = 0.5 * np.nanmedian(ifi) if np.isfinite(np.nanmedian(ifi)) else np.inf
        per_cam.append({
            "cam": k + 1,
            "deviceId": devices[order[k]],
            "nFrames": int(session.recorded[k]),
            "framesReadByThread": int(reader.frames_read),
            "durationS": float(duration),
            "fpsEffective": float(session.recorded[k] / duration)
                            if duration and duration > 0 else float("nan"),
            "ifiMedianMs": float(np.nanmedian(ifi)),
            "ifiSdMs": float(np.nanstd(ifi)),
            "ifiP99Ms": float(np.nanpercentile(ifi, 99)) if ifi.size > 1 else float("nan"),
            "ifiMaxMs": float(np.nanmax(ifi)),
            "maxQueueDepth": int(session.max_queue[k]),
            "queueCap": request.get("queueCap"),
            "queueWaitMedianMs": float(np.nanmedian(q_ms)),
            "writeMedianMs": float(np.nanmedian(w_ms)),
            "sinkMedianMs": float(np.nanmedian(sink)),
            "sinkP95Ms": float(np.nanpercentile(sink, 95)) if sink.size > 1 else float("nan"),
            "sinkMaxMs": float(np.nanmax(sink)),
            "sinkSlowThresholdMs": float(slow_over),
            "slowFrames": int(np.sum(sink > slow_over)),
            "bytes": int(session.bytes_out[k]),
            "height": int(session.dims[k][0]) if session.dims[k] else None,
            "width": int(session.dims[k][1]) if session.dims[k] else None,
            "exposureChanges": [[float(t_h), float(v)] for t_h, v in reader.exposure_changes],
            # How manual exposure was proven at open, measured rather than assumed. Without this
            # there is no way to tell afterwards whether a session ran under manual control or
            # under the camera's own auto-exposure, which looks entirely plausible and is not.
            "exposureControl": getattr(reader, "exposure_evidence", None),
        })

    info = {
        "sessionId": session.session_id,
        "sessionDir": str(session.dir),
        "clock": clock_bridge(session),
        "method": "single process: alignment then recording, cameras never released between",
        "backend": "dshow",
        "startedUtc": started.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "finishedUtc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "requested": request,
        "acceptedDeviceOrder": [devices[i] for i in order],
        "rotate180Stored": bool(rotate180),
        "stopReason": stop_reason,
        "perCamera": per_cam,
        "strobeAnchors": anchors,
        "trigger": trigger,
        "firstFrameTiffs": list(tiffs or []),
        "warnings": warnings_for(per_cam, session.recorded, anchors, stuck_readers, trigger),
    }
    with open(session.dir / "session.json", "w", encoding="utf-8") as fh:
        json.dump(info, fh, indent=2)
    return info


def warnings_for(per_cam, recorded, anchors, stuck_readers=(), trigger=None):
    """The guard rails. This strategy has no hardware drop counter -- the DirectShow backend
    exposes none -- so the queue depth, the frame ledger and the anchors are the available
    evidence that a session is trustworthy.

    `stuck_readers` is the one entry here that says nothing about THIS session's data: those
    frames are already stored and are as good as any other. It is a warning about the NEXT run,
    which is the one that will fail to open the camera, and it is recorded because by then the
    console that carried capture.close_all's warning has scrolled away or been closed."""
    w = []
    for c in per_cam:
        cap = c.get("queueCap") or 0
        # A queue at its cap means the storage path could not keep up. Sustained, backpressure
        # blocks the capture thread and frames go missing below Python, where DirectShow cannot
        # count them.
        if cap and c["maxQueueDepth"] >= cap:
            w.append(f"cam {c['cam']}: queue reached its {cap}-frame cap. The storage path did "
                     f"not keep up; frames may have been lost below Python, where the dshow "
                     f"backend cannot count them. Treat this session as suspect.")
        elif cap and c["maxQueueDepth"] > cap // 2:
            w.append(f"cam {c['cam']}: peak queue {c['maxQueueDepth']} of {cap} -- over half the "
                     f"buffer. Headroom is thinner than it should be.")

        # Queue depth says the buffer is filling; sink cost says why. Gated on slowFrames and
        # p95 rather than the maximum: every run, healthy or not, contains one 30-60 ms
        # outlier, so a peak-based gate fires on everything and separates nothing.
        thr = c.get("sinkSlowThresholdMs")
        slow, p95 = c.get("slowFrames"), c.get("sinkP95Ms")
        if thr and np.isfinite(thr) and slow and c["nFrames"]:
            frac = slow / c["nFrames"]
            if frac > 0.05:
                w.append(f"cam {c['cam']}: {slow} of {c['nFrames']} frames ({frac:.1%}) took "
                         f"longer than {thr:.1f} ms to reach disk, p95 {p95:.1f} ms. The "
                         f"storage path is the bottleneck, not the camera.")
            elif np.isfinite(p95 or np.nan) and p95 > thr:
                w.append(f"cam {c['cam']}: sink cost p95 {p95:.1f} ms against a "
                         f"{thr:.1f} ms half-frame budget. Storage is keeping up, but "
                         f"without much room.")
        if np.isfinite(c["ifiMaxMs"]) and np.isfinite(c["ifiMedianMs"]) and \
                c["ifiMedianMs"] > 0 and c["ifiMaxMs"] > 3 * c["ifiMedianMs"]:
            w.append(f"cam {c['cam']}: worst frame gap {c['ifiMaxMs']:.1f} ms against a "
                     f"{c['ifiMedianMs']:.1f} ms median -- roughly "
                     f"{c['ifiMaxMs'] / c['ifiMedianMs']:.1f} frame periods missing.")
    w.extend(_trigger_warnings(per_cam, trigger))
    if anchors and anchors.get("enabled"):
        for cam, marks in sorted(anchors.get("brightFrames", {}).items()):
            if not marks.get("start") or not marks.get("end"):
                w.append(f"{cam}: the strobe anchor did not produce detectable bright frames at "
                         f"both ends. Frame-to-pulse mapping must fall back to timestamp "
                         f"matching for this session.")
    if stuck_readers:
        w.append(f"teardown: {', '.join(stuck_readers)} was still reading when its capture was "
                 f"released, so this run may have left the camera claimed. THIS SESSION'S FRAMES "
                 f"ARE FINE. The risk is to the NEXT run: if it cannot open the camera, this is "
                 f"the run that caused it, and waiting a few seconds before retrying clears it.")
    return w


def _trigger_warnings(per_cam, trigger):
    """The guard rails that only mean anything under the FSIN trigger.

    Under the trigger the pulse train is the ground truth, and it says three things no free-run
    recording can say: both cameras see the SAME pulses, so their frame counts must match; the
    delivered rate must be the pulse rate; and a gap in the train shows up directly, as the
    driver's one-second timeout frame (see eyecal/trigger.py). Each is a drop check that needs no
    anchor and no queue-depth reasoning.

    The fallback warning is the other half. A run that asked for the trigger and recorded
    free-running instead is a perfectly good recording read an entirely different way, and the
    operator must not find that out by counting pulses in Spike2 afterwards.
    """
    if not trigger:
        return []
    if trigger.get("mode") != "trigger":
        if trigger.get("requested"):
            return [f"FSIN trigger fell back to free-run: {trigger.get('reason')}. This session "
                    f"has no pulse-to-frame identity; use the strobe anchors as usual."]
        return []

    w = []
    pulse_hz = trigger.get("pulseHz")
    if pulse_hz:
        for c in per_cam:
            fps = c["fpsEffective"]
            if np.isfinite(fps) and abs(fps - pulse_hz) > 0.10 * float(pulse_hz):
                w.append(f"cam {c['cam']}: stored {fps:.1f} frames per second against a "
                         f"{float(pulse_hz):g} Hz pulse train -- more than 10 percent apart. "
                         f"Either pulses were missed or the train was not the rate given.")
    counts = {c["cam"]: c["nFrames"] for c in per_cam}
    if len(set(counts.values())) > 1:
        w.append("the cameras stored different frame counts "
                 + ", ".join(f"cam {k} {n}" for k, n in sorted(counts.items()))
                 + ". Every camera sees the SAME FSIN pulses, so under the trigger these must "
                   "match; the difference is frames one camera lost.")
    # One timeout frame is how a normally ended run ENDS -- it is the second of silence after the
    # last pulse -- so it is not a gap. Every other one is a second in which the train stopped and
    # started again, which the frame indices cannot show on their own.
    ended_normally = str(trigger.get("endReason") or "").startswith("pulse train ended")
    for cam, n in sorted((trigger.get("timeouts") or {}).items()):
        gaps = n - 1 if ended_normally else n
        if gaps > 0:
            w.append(f"{cam}: {gaps} one-second gaps in the pulse train while recording. The "
                     f"frames on either side of a gap are consecutive in c<K>.bin but are NOT "
                     f"consecutive pulses.")
    return w


def report(info):
    print("\n" + "=" * 78)
    print(f"Session: {info['sessionId']}")
    print(f"  {info['sessionDir']}")
    print(f"  stop reason: {info['stopReason']}")
    print(f"  device order (left to right / cam1, cam2): {info['acceptedDeviceOrder']}")
    print("-" * 78)
    print(f"{'cam':>3} {'device':>7} {'frames':>8} {'span s':>8} {'fps':>8} {'IFI med':>9} "
          f"{'IFI SD':>8} {'worst gap':>10} {'peak q':>8} {'GB':>7}")
    for c in info["perCamera"]:
        print(f"{c['cam']:>3} {c['deviceId']:>7} {c['nFrames']:>8} {c['durationS']:>8.2f} "
              f"{c['fpsEffective']:>8.2f} {c['ifiMedianMs']:>8.2f}m {c['ifiSdMs']:>7.2f}m "
              f"{c['ifiMaxMs']:>9.1f}m {c['maxQueueDepth']:>8} {c['bytes'] / 1e9:>7.2f}")

    for c in info["perCamera"]:
        ec = c.get("exposureControl")
        if ec:
            print(f"  cam {c['cam']}: manual exposure via AUTO_EXPOSURE={ec['manualVia']}, "
                  f"image moved {ec['responseRatio']:.2f}x over "
                  f"{ec['probeDark']:g} -> {ec['probeBright']:g}")

    anchors = info.get("strobeAnchors") or {}
    if anchors.get("enabled"):
        print("\n  strobe anchors (bright frames -- match these to the wide strobe pulses):")
        for cam, marks in sorted(anchors.get("brightFrames", {}).items()):
            print(f"    {cam}: start {marks.get('start')}  end {marks.get('end')}")
    trig = info.get("trigger")
    if trig:
        note = f" -- {trig['reason']}" if trig.get("reason") else ""
        print(f"\n  FSIN trigger: {trig.get('mode')}{note}")
    if trig and trig.get("mode") == "trigger":
        rate = (f" against {float(trig['pulseHz']):g} Hz asked for" if trig.get("pulseHz")
                else "")
        for c in info["perCamera"]:
            name = f"cam{c['cam']}"
            print(f"    {name}: {c['nFrames']} frames at {c['fpsEffective']:.1f} fps{rate}, "
                  f"{(trig.get('timeouts') or {}).get(name, 0)} one-second timeouts")
        print(f"    stored frame k is pulse k. Ended: {trig.get('endReason')}")
    if info["warnings"]:
        print("\n  WARNINGS:")
        for w in info["warnings"]:
            print(f"    ! {w}")
    else:
        print("\n  No warnings: queue stayed well inside its cap and no stalls were seen.")
    print("=" * 78)
    print("Read the frames back with:")
    print("    import json, numpy as np")
    print(f"    d = r'{info['sessionDir']}'")
    print("    s = json.load(open(d + r'/frames/c1.json'))")
    print("    f = np.memmap(d + r'/frames/c1.bin', dtype=s['dtype'], mode='r',")
    print("                  shape=(s['count'], s['height'], s['width']))")


def read_frames(session_dir, cam=1):
    """Memory-mapped (count, height, width) view of one camera's frames.

    Included so a recording is never trapped in a format only this app understands. memmap
    rather than a read: a two-camera session is tens of gigabytes, and nothing that reads a
    session should have to hold all of it at once.
    """
    d = Path(session_dir)
    s = json.loads((d / "frames" / f"c{cam}.json").read_text(encoding="utf-8"))
    return np.memmap(d / "frames" / f"c{cam}.bin", dtype=s["dtype"], mode="r",
                     shape=(s["count"], s["height"], s["width"]))
