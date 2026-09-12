"""External trigger mode on the Arducam OV9281s: one frame per rising edge on FSIN.

The recording app normally lets the cameras free-run and marks the stored frames with the strobe
anchor (see record.py). Under the trigger the mapping needs no marking at all: stored frame k of
c<K>.bin IS pulse k, on every camera, because the camera exposes only when Spike2's pulse train
tells it to.

THE SWITCH IS ONE CONTROL WRITE, AND OPENCV CANNOT REACH IT. It is IAMCameraControl property 19,
AUTO_EXPOSURE_PRIORITY -- the UVC CT_AE_PRIORITY_CONTROL, Windows' "Low Light Compensation"
checkbox, Linux's exposure_dynamic_framerate. There is no CAP_PROP for it, so it is written
through a dshow.Controls filter of this module's own, beside whatever OpenCV is doing. See
eyecal/dshow.py, AE_PRIORITY, for the property itself.

MEASURED ON THIS RIG, 2026-09-11, two Arducam B0332 (OV9281), MJPG 1280x800, on BOTH devices:

  - AE priority 1 = one frame per rising edge on FSIN. AE priority 0 = free-running at ~121 fps
    again, immediately, with the media type unchanged.
  - THE WRITE WORKS WHILE OPENCV IS STREAMING, through a separately bound dshow.Controls, with
    capture.CameraReader still running. No reader stop, no re-apply of the source and no format
    check are needed.
  - IT PERSISTS IN THE CAMERA across cap.release() and a reopen. A camera left at 1 comes back
    at 1, which is why free_run_all is called BEFORE anything is identified or opened and again
    on EVERY exit path of the app. Without that, the next ordinary recording gets one black
    frame a second out of capture.open_all and sits out its whole settle timeout looking dead,
    and nothing in the app would explain why.
  - A WRITE FAILS WITH COMError 0x8007001F ("a device attached to the system is not
    functioning") for about half a second after cap.release(), while the driver tears the
    capture graph down. The exit restore runs exactly there, so THAT ONE failure is waited out --
    WRITE_ATTEMPTS tries, WRITE_RETRY_WAIT_S apart -- rather than reported. Every other failure
    is returned immediately; see write_ae_priority on the camera that has no such control.
  - IN TRIGGER MODE WITH NO PULSES THE READ DOES NOT FAIL. The driver answers ok=True about once
    a second, after its ~1000 ms timeout, with an ALL-ZERO frame. Those ghost frames are counted
    and dropped, never stored: one in the ledger would be a frame the sensor never took, and it
    would shift every pulse index after it.

WHY THE ALL-ZERO TEST IS EXACT, AND WHY A STRIDE OF 64 IS SAFE. A really triggered frame cannot
be all zero: the sensor's black-level pedestal sits near 28 counts even at the darkest exposure
the driver allows (see the comment on cameras.PROBE_DARK), so a real frame has that pedestal in
every pixel, not merely somewhere. The test therefore needs no threshold, and it does not need
every pixel either -- a stride-64 subsample of a 1280x800 frame is 20x13 pixels, and a real frame
carries the pedestal in all of them. The subsample is what makes the test cheap enough to run on
every frame in the recording loop.

The bench tools that established all of this, and the place to go back to when the hardware
changes, are tests/fsin_trigger_check.py (records and reports rates) and tests/fsin_trigger_view.py
(watch the train drive both cameras live).

COM lives on the main thread, like everything in dshow.py: nothing here may be called from a
capture thread.
"""

import sys
import time

import comtypes

from . import dshow

# A control write in the teardown window after cap.release() answers 0x8007001F and succeeds
# again about half a second later. Three tries a second apart covers it with room to spare, and
# the restore on exit is the write that must not quietly fail.
WRITE_ATTEMPTS = 3
WRITE_RETRY_WAIT_S = 1.0

# THE ONE FAILURE WORTH RETRYING, and the reason _is_teardown_error exists: 0x8007001F, "a device
# attached to the system is not functioning", is what the driver answers for about half a second
# after cap.release() while it tears the capture graph down. comtypes reports an HRESULT signed,
# so the same number is compared both ways round.
TEARDOWN_HRESULTS = (0x8007001F, -2147024865)

# Waited out after trigger mode is written, before the queues are drained, so the drain really
# empties them of everything the cameras produced BEFORE the switch. MEASURED BY
# tests/fsin_trigger_check.py, 2026-09-11: 8 to 16 free-running frames per camera still land after
# the write, at about 9 ms spacing, which is the free-running rate and not the trigger's. Drained
# too early they would be stored as pulse 0, 1, 2 and shift every pulse index after them.
MODE_SETTLE_S = 0.25

# Every 64th pixel in each direction: 20x13 of a 1280x800 frame. See the module docstring on why
# that is enough to tell the driver's all-zero timeout frame from a real one.
TIMEOUT_STRIDE = 64

# A preview pane blanks when no real frame has arrived within this long. At a 100 Hz pulse train
# that is ten pulse periods: a continuous picture while the train runs, and black within a tenth
# of a second of it stopping.
HOLD_MS = 100.0


def is_timeout_frame(img):
    """True when this frame is the driver's ~1000 ms timeout rather than an exposure.

    Subsampled rather than tested whole, because this runs on every frame of every camera in the
    recording loop. Exact, not a threshold: see the module docstring on the black-level pedestal.
    """
    return not img[::TIMEOUT_STRIDE, ::TIMEOUT_STRIDE].any()


def _is_teardown_error(exc):
    """True only for a write that failed because the capture graph is still coming down.

    That failure clears itself within about half a second; every other one is a statement about
    the device that will be exactly as true on the third try as on the first. See
    TEARDOWN_HRESULTS and write_ae_priority.
    """
    return (isinstance(exc, comtypes.COMError)
            and getattr(exc, "hresult", None) in TEARDOWN_HRESULTS)


def write_ae_priority(device, value):
    """Write the trigger control on ONE device and read it straight back. Never raises.

    Returns {requested, value, flags, attempts} on success, or {requested, error} once it has
    given up. THE READBACK IS THE EVIDENCE: a write that returned without an error is not, which
    is the same rule cameras.secure_manual_exposure follows for exposure.

    ONLY THE TEARDOWN FAILURE IS RETRIED. Anything else returns at once, with no sleep, because
    it is not going to change: a camera that does not expose camera-control property 19 at all --
    the ELP production cameras are the expected case -- makes Controls.set raise immediately, and
    retrying that would cost 2 s per device at each of the three points every run writes this
    control, with three stderr lines each time, for a device that will never answer.

    The filter is bound fresh for each attempt and closed again in a finally. It is not a claim
    on the camera -- it works while OpenCV is streaming from that same device -- and a filter
    bound before a failure is no use afterwards anyway.
    """
    entry = {"requested": value}
    for attempt in range(1, WRITE_ATTEMPTS + 1):
        try:
            controls = dshow.Controls(device - 1)
            try:
                controls.set("CameraControl", dshow.AE_PRIORITY, value)
                readback, flags = controls.get("CameraControl", dshow.AE_PRIORITY)
            finally:
                controls.close()
        except Exception as exc:
            entry["error"] = f"{type(exc).__name__}: {exc}"
            if attempt < WRITE_ATTEMPTS and _is_teardown_error(exc):
                time.sleep(WRITE_RETRY_WAIT_S)
                continue
            print(f"  device {device}: AE priority unavailable after {attempt} attempt(s) "
                  f"({entry['error']})", file=sys.stderr)
            return entry

        # An earlier attempt's error is dropped, not kept beside the answer: the write that
        # matters is the one that worked, and all_agree reads "error" as "this device was not
        # written".
        entry.pop("error", None)
        entry.update({"value": readback, "flags": flags, "attempts": attempt})
        print(f"  device {device}: AE priority -> {readback} (requested {value})")
        return entry


def set_all(devices, value):
    """Write the control on every device. {str(device): entry}, as write_ae_priority returns.

    Every device is attempted even after one has failed: the caller decides what a partial write
    means, and a device left unwritten is a device left in whatever state the last run put it in.
    """
    return {str(device): write_ae_priority(device, value) for device in devices}


def all_agree(readback, value):
    """True when EVERY device in a set_all result was written and read that value back."""
    entries = list(readback.values())
    return bool(entries) and all("error" not in e and e.get("value") == value for e in entries)


def free_run_all(devices):
    """Put every device back to free-running. The restore, and the only safe starting state.

    Called from the entry point in three places, all for the same reason -- the setting PERSISTS
    IN THE CAMERA between runs (module docstring): once before the cameras are identified and
    opened, so an earlier run that died in trigger mode cannot hand this one a black frame a
    second; once when a trigger attempt is abandoned, so the run falls back to an ordinary
    free-running recording; and once in the outer finally, after capture.close_all, so no exit
    path leaves a camera triggered.
    """
    return set_all(devices, 0)
