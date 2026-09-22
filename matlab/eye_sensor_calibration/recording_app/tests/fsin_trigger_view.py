r"""FSIN trigger view: watch the pulse train drive both cameras, live, and see it stop.

    python tests\fsin_trigger_view.py
    python tests\fsin_trigger_view.py --mode free
    python tests\fsin_trigger_view.py --devices 1 2 --camera ov9281 --hold-ms 150
    python "eye_sensor_calibration/recording_app/tests/fsin_trigger_view.py" --camera ov9281

A live two-camera window, like tools/camera_check.py, with both cameras put into EXTERNAL
TRIGGER MODE. A pane shows a picture only while frames keep arriving, so the screen follows the
FSIN pulse train directly: the panes are the pulses. Nothing is recorded. The measuring version
of this check is tests/fsin_trigger_check.py, which records N seconds and reports rates; this one
is for looking at the train while somebody changes it.

HOW TRIGGER MODE IS ENTERED ON THIS CAMERA. MEASURED ON THIS RIG, 2026-09-11, on BOTH devices:

  1. BACKLIGHT COMPENSATION IS NOT THE SWITCH, whatever Arducam's name for it ("low-brightness
     compensation") suggests. Its range is 0..2 with default 1, and 0, 1 and 2 all leave the
     camera free-running at about 121 fps. An earlier version of this view wrote it and proved
     nothing.
  2. THE SWITCH IS IAMCameraControl PROPERTY 19, AUTO_EXPOSURE_PRIORITY -- the "Low Light
     Compensation" checkbox on Windows, exposure_dynamic_framerate on Linux. See
     eyecal/dshow.py, AE_PRIORITY. Written 1 the sensor delivers ONE FRAME PER RISING EDGE on
     FSIN; written 0 it free-runs again immediately.
  3. THE WRITE GOES THROUGH THE PANE'S OWN dshow.Controls FILTER while its reader thread keeps
     streaming. Nothing is stopped, no source is re-applied and the media type does not change,
     so `t` is one COM call and the picture does not even blink.
  4. THE SETTING PERSISTS IN THE CAMERA across a release and a reopen, so this view writes
     free-run on every device BEFORE it opens anything, and again on the way out.
  5. IN TRIGGER MODE WITH NO PULSES THE READ DOES NOT FAIL. The driver hands OpenCV one
     ALL-BLACK frame about once a second, after its ~1000 ms timeout. This view counts those as
     `timeouts` and never displays one, so a pane with no pulses still goes black. A really
     triggered frame is never all zero -- the sensor's black-level pedestal sits near 28 counts
     (see the comment on cameras.PROBE_DARK). The full account is in tests/fsin_trigger_check.py.

THE OPERATOR DRIVES THE TRAIN, not this script. The ~50 Hz train -- 20 ms period, 1 ms high,
from a Power1401 digital output -- is started and stopped from Spike2 with ITS keys T and t,
in the Spike2 window. The `t` key IN THIS WINDOW is a different thing: it toggles trigger mode
on the cameras.

WHAT A WORKING TRIGGER LOOKS LIKE. While the train runs, both panes show a continuous picture
and the frame counters climb at about the pulse rate. When the operator stops the train, both
panes go BLACK within --hold-ms (100 ms by default, five pulse periods), the status line turns
red and says WAITING FOR FSIN PULSE, and the counters freeze. Starting the train again brings
the picture straight back.

WHAT A NON-WORKING ONE LOOKS LIKE. The picture never blanks and the measured rate stays at the
camera's own free-running rate whatever the operator does with the train. That is the camera
ignoring the control, i.e. still free-running, and it is the same picture --mode free gives on
purpose as the control condition.

Keys:  Esc quit | t trigger mode on/off | m prove manual exposure | p snapshot | c colour
       | r rotate | s swap sides | h hide overlay

ON EXIT BOTH CAMERAS ARE PUT BACK TO FREE-RUN, because the control PERSISTS IN THE CAMERA
between runs: a camera left in trigger mode gives the next ordinary recording nothing but one
black frame a second, and there is nothing in the recording app that would explain why. The
restore goes through a filter of its own, so it works even for a pane that never opened.
--keep-trigger skips it, for the one case where the next thing to run is another trigger check.

`m` and the opening exposure measurement both need frames, so they only do anything while the
train is running. The READER survives a camera with no pulses: the driver's timeout frames come
back as successful reads, so capture.CameraReader's give-up-after-50-failed-reads rule does not
fire in trigger mode. There is no reopen key here; restart the script.

The mechanism itself -- the control write, the all-zero timeout test and the blanking period --
lives in eyecal/trigger.py, which the recording app uses too; this view only drives it live.

A bench check, not a test. It needs two real cameras and an operator with a pulse train, so the
name carries no test_ prefix and pytest never collects it. Everything that can be checked
without a camera -- opening, the panes, the property readout, the snapshot -- lives in
tools/camera_check.py and is imported from there rather than copied.
"""

import argparse
import queue
import sys
import time
from pathlib import Path

import comtypes
import numpy as np

# tools/ is where camera_check lives, and camera_check is most of this script: the Pane, the
# open, the property readout, the drawing and the snapshot are all used from there so this view
# and the check can never disagree about what a camera in a pane is.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "tools"))
import camera_check                                             # noqa: E402

# eyecal/ is the recording app itself, one level up.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from eyecal import cameras, config, display, dshow, trigger     # noqa: E402

KEYS_HELP = ("Esc quit | t trigger mode on/off | m prove manual exposure | p snapshot | "
             "c colour | r rotate | s swap sides | h hide overlay")

# A pane blanks when no frame has arrived within this long. The app's own trigger preview blanks
# on the same number, for the same reason, so it is kept in one place: see trigger.HOLD_MS. At
# 50 Hz the frames are 20 ms apart, so it is five pulse periods.
HOLD_MS = trigger.HOLD_MS


# --- the pure parts: the readout ------------------------------------------------------------

def trigger_lines(position, pane, t, hold_ms=HOLD_MS):
    """The whole readout for one camera under the trigger, as (line, is_red) pairs.

    Pure, so it can be checked without a camera: everything it needs is already on the pane.
    is_red is what the caller colours by, decided here where the numbers are still in hand.

    The driver frame rate camera_check shows is left out. It is the sentinel the negotiation
    asked for (see cameras.apply_source), and under an external trigger the driver's opinion of
    the rate means even less than usual -- the pulse train sets it. The MEASURED rate is the one
    that answers the question this view exists for.

    THAT MEASURED RATE IS THE READER'S RAW RATE, timeout frames and all. In trigger mode with no
    pulses it settles near 1 fps, which is the driver handing back one black frame a second, not
    the camera exposing anything. `frames` and `timeouts` on the next line are the split.
    """
    blanked = pane.last_frame_t is None or t - pane.last_frame_t > hold_ms / 1000.0
    want = 1 if pane.trigger_on else 0
    value = pane.ae_readback
    # Unknown is not a disagreement: a readback that never arrived says nothing about the camera,
    # and colouring it red would cry wolf on every camera whose COM bind failed.
    disagrees = value is not None and value != want
    age = "never" if pane.last_frame_t is None else f"{(t - pane.last_frame_t) * 1e3:.0f} ms ago"

    lines = [(f"cam {position}  device {pane.device}   {pane.friendly or '(no name)'}", False),
             (f"preset {pane.preset_name} ({pane.preset['title']})   {pane.media_type}   "
              f"measured {pane.fps:.1f} fps", False),
             (f"trigger mode: {'ON' if pane.trigger_on else 'OFF'} requested   AE priority "
              f"reads {'?' if value is None else value}", disagrees),
             (f"frames {pane.frames_seen}   timeouts {pane.timeouts}   last frame {age}   "
              f"[{'*' if pane.blink else '.'}]", blanked)]
    if not pane.trigger_on:
        lines.append(("free running -- the camera sets its own frame timing", False))
    elif blanked:
        lines.append(("WAITING FOR FSIN PULSE", True))
    else:
        lines.append(("FRAMES ARRIVING", False))
    if pane.note:
        lines.append((pane.note, False))
    return lines


# --- one camera under the trigger -------------------------------------------------------------

def new_pane(device):
    """A camera_check.Pane with this view's own trigger bookkeeping bolted on.

    Pane is a plain class, so the four extra attributes simply live on the instance. They are
    set here rather than in set_trigger so that a pane which never opened, or whose first write
    failed, still draws instead of raising on a missing attribute.
    """
    pane = camera_check.Pane(device)
    pane.trigger_on = False             # what the last set_trigger REQUESTED, not a readback
    pane.ae_readback = None             # what the driver said AE priority is, or None
    reset_trigger_state(pane)
    return pane


def reset_trigger_state(pane):
    """Forget what this pane has seen. Called whenever it gets a new reader, which counts from
    zero -- keeping the old count would show a frame arriving that never did, and worse, would
    hold the pane's picture up for one more hold window."""
    pane.last_frame_t = None            # preview time of the newest REAL frame, or None
    pane.frames_seen = 0                # real frames drained since the last reset
    pane.timeouts = 0                   # all-black frames drained; see drain_pane
    pane.blink = False                  # flipped per frame, so the overlay ticks visibly


def drain_pane(pane, t):
    """Take everything this pane's reader has queued and split real frames from timeouts.

    Everything, not capture.drain_newest: a timeout frame and a triggered frame come off the same
    queue and only the pixels tell them apart, so every item has to be looked at. The preview
    still ends up showing the newest real frame, which is all a live view wants.

    AN ALL-ZERO FRAME IS THE DRIVER'S ~1000 ms TIMEOUT, not an exposure (fact 5 in the module
    docstring). It is counted and dropped: latest and last_frame_t are left alone, so a pane with
    no pulses still goes black on time instead of showing a black frame and calling it live.
    trigger.is_timeout_frame is the test, and why it needs no threshold is written out there.

    frames_seen counts REAL frames here, not reader.frames_read, because the reader counts the
    timeouts too. The rate on the overlay is still the reader's own; trigger_lines says so.
    """
    if pane.reader is None:
        return
    while True:
        try:
            _t_host, _t_abs, frame = pane.reader.q.get_nowait()
        except queue.Empty:
            return
        if trigger.is_timeout_frame(frame):
            pane.timeouts += 1
            continue
        pane.latest = frame
        pane.last_frame_t = t
        pane.frames_seen += 1
        pane.blink = not pane.blink


def refresh_ae(pane):
    """Re-read AE priority onto the pane. Never raises; unavailable reads as None.

    Read on the property readout's own period, so the overlay follows a change made from OUTSIDE
    this script -- somebody ticking Low Light Compensation in AMCap, say -- instead of only ever
    repeating what set_trigger last asked for.
    """
    if pane.controls is None:
        return
    try:
        pane.ae_readback = pane.controls.get("CameraControl", dshow.AE_PRIORITY)[0]
    except Exception:
        pane.ae_readback = None


def set_trigger(pane, on):
    """Put ONE camera into external trigger mode, or take it back out. THE ONLY CONTROL WRITE.

    One COM call on this pane's own filter, with the reader thread left running. MEASURED ON THIS
    RIG, 2026-09-11: the write takes effect on a camera OpenCV is streaming from, and the media
    type is unchanged afterwards. So there is no reader stop, no re-apply of the preset source and
    no verify_format here, and the picture does not even blink (see eyecal/dshow.py, Controls.set).

    The exposure is left exactly as the open set it. It has to stay shorter than the trigger
    period, or the sensor cannot finish a frame before the next edge arrives; the ov9281 preset's
    -11 is about 0.5 ms against the train's 20 ms, with room to spare.

    The frame bookkeeping is reset because counts from either side of a mode change are not
    comparable. Returns False, having written nothing, when there are no driver controls.
    """
    if pane.controls is None:
        print(f"  device {pane.device}: no driver controls, cannot switch trigger mode",
              file=sys.stderr)
        return False
    try:
        pane.controls.set("CameraControl", dshow.AE_PRIORITY, 1 if on else 0)
        pane.ae_readback = pane.controls.get("CameraControl", dshow.AE_PRIORITY)[0]
    except comtypes.COMError as exc:
        print(f"  device {pane.device}: the trigger write failed ({exc})", file=sys.stderr)
        return False

    pane.trigger_on = on
    reset_trigger_state(pane)
    print(f"  device {pane.device}: trigger mode {'ON' if on else 'OFF'} requested "
          f"(AE priority reads {pane.ae_readback})")
    return True


def free_run_device(device):
    """Put one device back to free-running, through a filter bound just for it. Never raises.

    No pane and no capture needed, which is the point. It runs BEFORE anything is opened, because
    the setting persists in the camera and one left in trigger mode by a crashed run would hand
    open_pane's reader a black frame a second and look dead; and it runs again on the way out,
    where it has to work for a pane that never opened or whose reader has gone.

    trigger.write_ae_priority does the write, the readback and the reporting -- the same code the
    recording app restores with, including waiting out the teardown window after a release.
    """
    return "error" not in trigger.write_ae_priority(device, 0)


# --- the loop -------------------------------------------------------------------------------

def run(panes, cfg, hold_ms=HOLD_MS, quit_after=None):
    """The live view. camera_check.run's loop, with the picture gated on the pulse train.

    THE GATE IS THE WHOLE POINT. camera_check keeps showing the last frame it got, which is the
    right thing for a tool that reports camera state and the wrong thing entirely here: a frozen
    picture and a live one look identical, and "did the frames stop?" is the question. So a pane
    shows its latest frame only while one has arrived within hold_ms, and black otherwise.

    Returns nothing. Its output is what the operator saw.
    """
    downsample = max(int(cfg["downsample"]), 1)
    preview_hz = float(cfg["align_preview_hz"])
    hold_s = hold_ms / 1000.0
    # Once, not per redraw: a pane's format cannot change while the view is running.
    full_shapes, pane_shapes = [], []
    for pane in panes:
        _, full_w, full_h = cameras.parse_format(pane.fmt)
        full_shapes.append((full_h, full_w))
        pane_shapes.append((full_h // downsample, full_w // downsample))

    colour_idx, show_overlay, rot = 0, True, bool(cfg["rotate180"])
    order = list(range(len(panes)))
    period = 1.0 / preview_hz if preview_hz > 0 else 0.0
    t0, next_at = time.perf_counter(), 0.0
    snapshot_error = ""

    window = display.Window("FSIN trigger view -- records nothing", cfg["fullscreen"])
    print(f"\nKeys: {KEYS_HELP}")
    print(f"A pane goes black when no frame has arrived for {hold_ms:g} ms.")
    print("Start and stop the pulse train from Spike2 (its own T and t).\n")

    try:
        while True:
            t = time.perf_counter() - t0
            for pane in panes:
                if pane.reader is None:
                    continue
                if pane.reader.error is not None:
                    # A capture thread that died takes its camera down with it. The trigger is
                    # NOT the cause: the driver's timeout reads come back successful, so a camera
                    # with no pulses never trips capture.CameraReader's failed-read rule. The
                    # pane says what failed and the script has to be restarted for that camera.
                    pane.error = f"capture thread failed: {pane.reader.error}"
                    print(f"\n  device {pane.device}: {pane.error}", file=sys.stderr)
                    camera_check.close_pane(pane)
                    continue
                drain_pane(pane, t)
                pane.update_rate(t)

            if t >= next_at:
                for pane in panes:
                    if pane.controls is not None and t >= pane.next_property_read:
                        pane.rows = camera_check._read_rows(pane)
                        # Alongside the rest of the readout, and on its own period for the same
                        # reason: one COM round trip per property is not free.
                        refresh_ae(pane)
                        pane.next_property_read = t + camera_check.PROPERTY_PERIOD_S

                drawn = []
                for position, k in enumerate(order, start=1):
                    pane = panes[k]
                    live = (pane.latest is not None and pane.last_frame_t is not None
                            and t - pane.last_frame_t <= hold_s)
                    img = (display.display_copy(pane.latest, downsample, rot) if live
                           else np.zeros(pane_shapes[k], np.uint8))
                    if pane.open:
                        lines = trigger_lines(position, pane, t, hold_ms)
                    else:
                        lines = camera_check.failed_lines(position, pane.device,
                                                          pane.error or "not open")
                    if snapshot_error:
                        lines.append((f"SNAPSHOT FAILED -- {snapshot_error}", True))
                    drawn.append(camera_check.draw_pane(img, lines,
                                                        display.COLOURS[colour_idx][1],
                                                        show_overlay))
                canvas = display.tile(drawn)
                if show_overlay:
                    camera_check.dim_band(canvas, canvas.shape[0] - 24, canvas.shape[0])
                    display.text(canvas, KEYS_HELP, (8, canvas.shape[0] - 12),
                                 display.COLOURS[colour_idx][1], 0.45)
                window.show(canvas)
                next_at = t + period

            key = window.pump()
            if key == 27:                                       # Esc
                print("\nQuit (Esc).")
                return
            elif key == ord("t"):
                for pane in panes:
                    if pane.open:
                        set_trigger(pane, not pane.trigger_on)
                t0, next_at = time.perf_counter() - t, 0.0      # the writes are not preview time
            elif key == ord("m"):
                for pane in panes:
                    if pane.preset is None:
                        continue
                    camera_check.prove_manual(pane, pane.preset["source"],
                                              pane.preset["grayscale"])
                    # Not because of the control -- prove_manual cannot reach it, and the camera
                    # holds it anyway -- but because prove_manual gives the pane a NEW reader.
                    # This puts the frame bookkeeping back in step with it.
                    set_trigger(pane, pane.trigger_on)
                t0, next_at = time.perf_counter() - t, 0.0
            elif key == ord("p"):
                chosen = camera_check.ask_snapshot_path(cfg["out"])
                if chosen is None:
                    print("  snapshot cancelled")
                else:
                    # The raw file holds the last frame each camera DELIVERED, even for a pane
                    # that has gone black: the blanking is a statement about timing, and the
                    # pixels are still the pixels.
                    _, failed = camera_check.write_snapshot(
                        chosen, canvas, [panes[k].latest for k in order], rot,
                        full_shapes[order[0]])
                    snapshot_error = "; ".join(failed)
                t0, next_at = time.perf_counter() - t, 0.0      # the dialog is not preview time
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


def main(argv=None):
    """Open both cameras the app's way, put them under the trigger, and show what arrives.

    Each device is identified on its own and a mixed pair is accepted, exactly as in
    camera_check.main -- nothing is recorded here either. The open itself writes no control
    (camera_check.open_pane), so the only camera control this view ever writes is AE priority:
    free-run before the open, the requested mode after it, free-run again on the way out.
    """
    args = parse_args(argv)
    cfg = config.apply_cli(config.load(args.config, required=args.config is not None), args)

    if cfg["camera"] not in cameras.CAMERA_PRESETS and cfg["camera"] != cameras.AUTO:
        raise KeyError(f"Unknown camera {cfg['camera']!r}. "
                       f"Known: {', '.join(sorted(cameras.CAMERA_PRESETS))}, "
                       f"or {cameras.AUTO} to identify it")
    devices = [int(d) for d in cfg["devices"]]
    want_trigger = args.mode == "trigger"

    print("\n" + "=" * 88)
    print("FSIN TRIGGER VIEW -- the panes follow the pulse train. Nothing is recorded.")
    print("=" * 88)
    for entry in camera_check._enumerate_devices():
        print(f"  index {entry['index']} (device {entry['index'] + 1}): "
              f"{entry['friendlyName']}")

    panes = [new_pane(d) for d in devices]
    if cfg["camera"] == cameras.AUTO:
        print("Identifying each camera from the geometries it offers (auto)...")
        for pane in panes:
            camera_check.identify_pane(pane, cfg["format"])
        camera_check.fill_unidentified(panes)
    else:
        for pane in panes:
            camera_check.assign_preset(pane, cfg["camera"], cfg["format"])
    for pane in panes:
        if pane.preset is None:
            print(f"  device {pane.device}: NOT IDENTIFIED (the reason is above)")
        else:
            print(f"  device {pane.device}: {pane.preset_name} ({pane.preset['title']})  "
                  f"{pane.fmt}")
    print(f"Devices     : {devices}")
    mode_note = ("trigger -- AE priority 1, one frame per rising edge on FSIN" if want_trigger
                 else "free running -- AE priority 0, the control condition")
    print(f"Mode        : {mode_note}")
    print(f"Blank after : {args.hold_ms:g} ms with no frame")

    # BEFORE the open, on every device, whether or not its pane could be identified: the setting
    # persists in the camera, so one left in trigger mode by a crashed run would deliver nothing
    # but a black frame a second and its pane would simply look dead.
    print("Putting every device into free-run before opening it...")
    for device in devices:
        free_run_device(device)
    print("Opening cameras (the open writes nothing; the trigger write comes after it).")

    try:
        for pane in panes:
            # Sequential, as the app opens them: OpenCV's DirectShow backend runs through a
            # shared global that is not documented as thread-safe.
            if pane.preset is None:
                continue
            if camera_check.open_pane(pane, pane.fmt, pane.preset["source"],
                                      pane.preset["grayscale"]):
                camera_check.open_controls(pane)
        for pane in panes:
            if pane.open:
                set_trigger(pane, want_trigger)
        run(panes, cfg, float(args.hold_ms), args.quit_after)
    finally:
        # The control PERSISTS IN THE CAMERA, so a camera left in trigger mode gives the next
        # ordinary recording nothing but a black frame a second. Restoring free-run is therefore
        # part of exiting, on every path, and --keep-trigger is the only way out of it. Through a
        # filter of its own, not through the pane, so it still happens for a pane that never
        # opened or whose reader died.
        if not args.keep_trigger:
            for device in devices:
                print(f"  restoring free-run on device {device}")
                free_run_device(device)
        for pane in panes:
            camera_check.close_pane(pane)
    return 0


def parse_args(argv):
    """camera_check's options, plus the three this view adds.

    Everything defaults to None so config.json wins unless the flag was actually given, which is
    the same rule config.apply_cli applies for the app. --camera included: the app's own
    "auto" then decides, per device.
    """
    p = argparse.ArgumentParser(
        description="Live view of two cameras under the FSIN external trigger. The panes go "
                    "black when the pulse train stops. Records nothing.")
    p.add_argument("--config", default=None,
                   help="path to config.json (default: the app's, beside recording_app/)")
    p.add_argument("--camera", default=None,
                   choices=sorted(cameras.CAMERA_PRESETS) + [cameras.AUTO],
                   help="preset name, or auto to identify it from the geometries offered")
    p.add_argument("--devices", type=int, nargs="+", default=None,
                   help="1-based device ids; the ORDER is the initial left-to-right arrangement")
    p.add_argument("--format", default=None, help="e.g. MJPG_1280x800")
    p.add_argument("--downsample", type=int, default=None, help="preview pixel stride")
    p.add_argument("--align-preview-hz", type=float, default=None, help="preview redraw ceiling")
    p.add_argument("--rotate180", action="store_true", default=None)
    p.add_argument("--no-rotate180", dest="rotate180", action="store_false", default=None)
    p.add_argument("--fullscreen", action="store_true", default=None)

    p.add_argument("--mode", choices=("trigger", "free"), default="trigger",
                   help="the mode written right after the open. trigger writes AE priority 1, "
                        "which is this camera's external-trigger switch; free writes 0, which "
                        "is the control condition (default: %(default)s)")
    p.add_argument("--hold-ms", type=float, default=HOLD_MS,
                   help="a pane goes black when no frame has arrived for this long; at 50 Hz "
                        "the frames are 20 ms apart (default: %(default)s)")
    p.add_argument("--keep-trigger", action="store_true",
                   help="do NOT restore free-run on exit; a camera left in trigger mode gives "
                        "the next ordinary recording nothing but a black frame a second")
    p.add_argument("--quit-after", type=float, default=None,
                   help="exit after this many seconds (unattended runs)")

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
