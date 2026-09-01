"""One capture thread per camera, feeding one bounded queue each.

The thread has exactly ONE behaviour: read a frame, timestamp it, convert it, queue it. It has
no idea whether the app is currently aligning or recording. That distinction lives entirely in
whoever drains the queue -- align.py throws everything away but the newest frame, record.py
keeps every one. Keeping the phase out of the thread is what lets the cameras stay open and
streaming across the whole session without a mode switch that could go wrong mid-run.

Threads are worth it because cv2.VideoCapture.read releases the GIL while it waits on the driver
and decodes. In a single-threaded loop the next read cannot start until the current frame has
been converted, written and drawn, so per-frame work does not merely add latency -- it delays
the next acquisition.
"""

import collections
import queue
import sys
import threading
import time

import cv2

from . import cameras

# Generous enough that normal operation never touches it, small enough that a consumer which has
# stopped keeping up applies backpressure instead of eating RAM. On overflow the capture thread
# BLOCKS rather than discarding: a frame this code threw away would be indistinguishable from
# one the camera never delivered, and DirectShow has no counter to tell them apart.
QUEUE_MAX = 256

# Discarded after the stream starts, by time rather than by frame count. The old 15-frame rule
# existed to let AUTO-exposure settle; exposure here is manual, and by the time anything reads
# this camera it has usually been free-running at its final settings for seconds already. This
# is a small safety margin, not a warm-up.
SETTLE_S = 0.3


class CameraReader(threading.Thread):
    """Reads one camera as fast as it will deliver and timestamps at the point of arrival.

    The timestamp is taken immediately after read() returns and before any per-frame work, so it
    reflects when the frame arrived rather than when the consumer got to it.

    Exposure changes are requested from other threads but APPLIED HERE, between two reads. All
    access to the cv2.VideoCapture therefore happens on the thread that owns it -- OpenCV makes
    no thread-safety promise about set() racing a read() on the same capture, and the recording
    anchor would otherwise do exactly that.
    """

    def __init__(self, cap, device, want_gray, settle_s=SETTLE_S, maxsize=QUEUE_MAX):
        super().__init__(daemon=True, name=f"cam{device}")
        self.cap, self.device, self.want_gray = cap, device, want_gray
        self.settle_s = settle_s
        self.q = queue.Queue(maxsize=maxsize)
        self.warmed = threading.Event()
        self.error = None
        self.frames_read = 0
        self.exposure_changes = []          # (t_host, value) applied, for session.json
        self._pending_exposure = None
        self._lock = threading.Lock()
        self._stopping = threading.Event()
        self.released = False               # set by close_all, so a second call is a no-op

    def run(self):
        try:
            settle_until = time.perf_counter() + self.settle_s
            while time.perf_counter() < settle_until:
                if self._stopping.is_set():
                    return
                self.cap.read()
            self.warmed.set()

            failures = 0
            while not self._stopping.is_set():
                ok, frame = self.cap.read()
                t_host = time.perf_counter()
                t_abs = time.time()
                if not ok or frame is None:
                    failures += 1
                    if failures > 50:
                        # Busy-class: a device that opened but answers nothing is wedged in the
                        # same way one that never opened is, and clears the same way.
                        raise cameras.CameraBusyError(
                            f"camera {self.device}: 50 consecutive failed reads")
                    continue
                failures = 0
                self.frames_read += 1
                self._apply_pending_exposure(t_host)

                # OpenCV hands back BGR even for a monochrome sensor, because the MJPG payload
                # decodes to three channels. Converting here keeps the array the size the sensor
                # actually produced, stores one byte per pixel instead of three, and keeps the
                # work off the consumer thread.
                if self.want_gray and frame.ndim == 3:
                    frame = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)

                while not self._stopping.is_set():
                    try:
                        self.q.put((t_host, t_abs, frame), timeout=0.1)
                        break
                    except queue.Full:
                        continue
        except BaseException as exc:        # surfaced by the main thread, never swallowed
            self.error = exc

    def _apply_pending_exposure(self, t_host):
        with self._lock:
            value, self._pending_exposure = self._pending_exposure, None
        if value is not None:
            cameras.set_exposure(self.cap, value)
            self.exposure_changes.append((t_host, float(value)))

    def request_exposure(self, value):
        """Ask the capture thread to write this exposure after its next read."""
        with self._lock:
            self._pending_exposure = value

    def exposure_settled(self):
        """True once every requested exposure change has actually been written."""
        with self._lock:
            return self._pending_exposure is None

    def stop(self):
        self._stopping.set()


def drain_newest(reader):
    """Empty the queue and return only the newest item, or None.

    What a live preview wants: a view showing a frame from two seconds ago is worse than useless.
    The recording path deliberately does the opposite and consumes in order.
    """
    item = None
    while True:
        try:
            item = reader.q.get_nowait()
        except queue.Empty:
            return item


def check_readers(readers):
    """Surface a capture-thread failure on the main thread.

    The busy/not-busy distinction is preserved across the thread boundary, because it is what
    decides whether the entry point retries the open or gives up (see cameras.CameraBusyError).
    """
    for r in readers:
        if r.error is not None:
            kind = (cameras.CameraBusyError
                    if isinstance(r.error, cameras.CameraBusyError) else RuntimeError)
            raise kind(f"capture thread {r.name} failed: {r.error}") from r.error
        if r.ident is not None and not r.is_alive():
            raise RuntimeError(f"capture thread {r.name} exited unexpectedly")


def open_all(preset, devices, fmt=None, exposure=None, anchor_exposure=None,
             settle_timeout=30.0):
    """Open every camera, start a reader each, and block until all of them are streaming.

    A camera that never streams is the most common bench failure and must not be able to trap
    the operator forever, hence the timeout. Any failure closes whatever was already opened, so
    a half-built set of cameras is never left claimed.
    """
    source = dict(preset["source"])
    if exposure is not None:
        source["Exposure"] = exposure

    readers = []
    try:
        for device in devices:
            cap, description, evidence = cameras.open_camera(
                device, fmt or preset["format"], source, anchor_exposure)
            reader = CameraReader(cap, device, preset["grayscale"])
            # Carried on the reader purely so session.json can record HOW manual exposure was
            # secured on this device. Nothing in the capture loop reads it.
            reader.exposure_evidence = evidence
            readers.append(reader)
            print(f"  device {device}: {description}")
        for r in readers:
            r.start()

        deadline = time.perf_counter() + settle_timeout
        for r in readers:
            while not r.warmed.wait(timeout=0.05):
                check_readers(readers)
                if time.perf_counter() > deadline:
                    raise cameras.CameraBusyError(
                        f"Device {r.device} delivered nothing in {settle_timeout:g} s. It is not "
                        f"streaming -- unplug and replug it, or run imaqreset if MATLAB is still "
                        f"holding it.")
        print(f"  all {len(readers)} cameras streaming")
    except BaseException:
        close_all(readers)
        raise
    return readers


def close_all(readers, timeout=10.0):
    """Stop every reader, join it, THEN release. The join is what makes the release safe.

    A reader spends almost all its life blocked inside cap.read() and only sees the stop flag
    once that returns -- normally within a frame period, far longer if the device has stalled.
    Releasing a VideoCapture while another thread is still inside read() on it tears the
    DirectShow graph down underneath an active read, and the camera can stay claimed after this
    process exits. That is what makes the NEXT run open a capture whose properties all answer -1.

    Returns the names of the readers that were STILL RUNNING when their capture was released --
    that is, the ones this run may have leaked. The entry point puts that list into session.json,
    because the run it predicts a failure for is the NEXT one, and by then this console is gone.

    Called both from the entry point's normal path (before session.json is written, so the list
    can go into it) and from its one finally, so every exit -- accept, cancel, exception, Ctrl-C
    -- cleans up by exactly this path. The second call is a no-op: a reader is only closed once.
    """
    readers = [r for r in readers if not r.released]
    if not readers:
        return []

    for r in readers:
        r.stop()

    deadline = time.perf_counter() + timeout
    for r in readers:
        if r.ident is not None:                     # never started; nothing to join
            r.join(timeout=max(deadline - time.perf_counter(), 0.0))
    stuck = [r.name for r in readers if r.is_alive()]

    for r in readers:
        try:
            r.cap.release()
        except Exception:
            pass
        r.released = True
    if stuck:
        print(f"  WARNING: released the capture while {', '.join(stuck)} was still reading. The "
              f"camera may stay claimed after this process exits -- if the next run cannot open "
              f"it, wait a few seconds and retry.", file=sys.stderr)
    return stuck


class RateMeter:
    """Rolling frame rate over the last `window` timestamps, for on-screen readout only."""

    def __init__(self, window=60):
        self._t = collections.deque(maxlen=window)

    def tick(self, t):
        self._t.append(t)

    @property
    def fps(self):
        if len(self._t) < 2:
            return float("nan")
        span = self._t[-1] - self._t[0]
        return (len(self._t) - 1) / span if span > 0 else float("nan")
