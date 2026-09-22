"""Align a camera pair, then record them. One process, cameras never released in between.

    python "eye_sensor_calibration/recording_app/eye-calibration-recording.py" --seconds 10 --out "C:/Temp/test"

Run by Spike2 as a single step. The phases are:

    1. open both cameras and negotiate the format          Spike2 NOT sampling
    2. live alignment view; the operator presses 'a'       Spike2 NOT sampling
    3. put both cameras into FSIN trigger mode             Spike2 NOT sampling
    4. create the ready flag                               Spike2 sees it and marks the file
    5. record one frame per FSIN pulse                     Spike2 sampling and pulsing
    6. write session.json                                  Spike2 sees it and continues

Phase 5 falls back to a free-running recording, bracketed by the strobe anchor, if no pulses
arrive within --trigger-wait-s. Under the trigger the mapping is free -- stored frame k IS pulse
k -- and the anchor is not used; free-running, the anchor is how a frame is matched to a pulse.
Both paths are recorded in session.json, and a fallback is a warning (exit code 3), because the
two are read in completely different ways. See eyecal/trigger.py.

Why one process. The cameras strobe once per exposure from the moment their graph starts, and
those pulses are wired to Spike2 digital inputs. Anything that happens before SampleStart() is
therefore invisible -- which is why all of the opening, negotiating and aligning is done first,
and why the cameras are NOT closed and reopened afterwards. Keeping them open also means the
accepted device order and rotation never have to be handed between processes: they stay in
memory.

Alignment keys:  a accept | Esc cancel | c colour | r rotate 180 | s swap sides | h hide overlay
Closing the alignment window cancels. Closing the recording preview stops early, which still
finalises the session.

Process exit codes -- useful from a terminal, but NOT what Spike2 reads (ProgStatus cannot carry
them; see eyecal/spike2.py):

    0  recorded and finalised
    1  failed (camera would not open or could not be identified, wrong format, not streaming)
    2  bad command line (argparse)
    3  recorded, but a guard tripped -- read the warnings in session.json
    4  cancelled by the operator at the alignment stage
"""

import argparse
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from eyecal import align, cameras, capture, config, record, session, spike2, trigger

EXIT_OK, EXIT_FAILED, EXIT_USAGE, EXIT_WARNINGS, EXIT_CANCELLED = 0, 1, 2, 3, 4

# A device left claimed by a previous run is the one opening failure that clears itself, given a
# few seconds (see cameras.CameraBusyError and capture.close_all). Retrying it here is free in the
# only currency this app spends: every attempt happens in phase 1, before the ready flag exists,
# so Spike2 is not sampling and has not been told anything yet. The alternative is what the bench
# project measured -- one run in 394 lost to a camera that would have been ready shortly after.
OPEN_ATTEMPTS = 3
OPEN_RETRY_WAIT_S = 5.0


def main(argv=None):
    args = parse_args(argv)
    cfg = config.apply_cli(config.load(args.config, required=args.config is not None), args)

    if cfg["camera"] not in cameras.CAMERA_PRESETS and cfg["camera"] != cameras.AUTO:
        raise KeyError(f"Unknown camera {cfg['camera']!r}. "
                       f"Known: {', '.join(sorted(cameras.CAMERA_PRESETS))}, "
                       f"or {cameras.AUTO} to identify it")
    devices = [int(d) for d in cfg["devices"]]
    seconds = float(cfg["seconds"])
    if seconds <= 0:
        raise ValueError(f"seconds must be positive, got {seconds}")

    # BEFORE anything opens a camera, because the trigger control PERSISTS IN THE CAMERA between
    # runs: one left in trigger mode by a crashed run delivers a single black frame a second, and
    # capture.open_all would sit out its whole settle timeout on a camera that looks dead. ONE
    # call is enough for both stages that follow -- identification releases each camera and the
    # open reopens it, but nothing writes trigger mode in between, and the camera keeps the 0.
    if cfg["trigger"]:
        print("Putting every camera into free-run before opening it (the control persists)...")
        trigger.free_run_all(devices)

    # Identified FIRST, because everything downstream is built out of the preset: the recording
    # anchor takes its normal exposure from it, and the format cannot be negotiated until it is
    # known which camera is on the other end. It costs about 5 s a camera and is still phase 1,
    # so Spike2 is not sampling and has been told nothing yet. Retried on a busy device exactly
    # as the open is -- see _retry_busy.
    camera, detection = cfg["camera"], None
    if camera == cameras.AUTO:
        print("Identifying cameras from the geometries they offer (auto)...")
        camera, detection = _retry_busy(lambda: cameras.detect_preset(devices))
    preset = cameras.CAMERA_PRESETS[camera]
    camera_note = (preset["title"] if detection is None
                   else f"{cameras.AUTO} -> {camera} ({preset['title']})")

    started = datetime.now(timezone.utc)
    session_id = args.session_id or (
        started.strftime("%Y%m%d-%H%M%SZ") + f"_{camera}_{len(devices)}cam_{seconds:g}s"
        + (f"_{args.tag}" if args.tag else ""))
    session_dir = Path(cfg["out"]).expanduser().resolve() / session_id

    anchor = record.make_anchor(
        cfg["anchor"],
        cfg["exposure"] if cfg["exposure"] is not None else preset["source"].get("Exposure"),
        cfg["anchor_exposure"], float(cfg["anchor_hold_s"]))

    anchor_note = "off"
    if anchor.enabled:
        anchor_note = (f"exposure {anchor.normal:g} -> {anchor.bright:g} for "
                       f"{anchor.hold_s:g} s at each end")
    trigger_note = "off (free-running recording with the strobe anchor)"
    if cfg["trigger"]:
        trigger_note = (f"on (wait {float(cfg['trigger_wait_s']):g} s for the first pulse, then "
                        f"fall back to free-run)")
    pulse_note = ("not given (--pulse-hz; no rate check)" if cfg["pulse_hz"] is None
                  else f"{float(cfg['pulse_hz']):g} Hz (from --pulse-hz)")
    handshake_note = "off (recording starts straight after accept)"
    if cfg["handshake"]:
        handshake_note = cfg["ready_flag"]
        if float(cfg["handshake_wait_s"]) > 0:
            handshake_note += f"  (wait up to {float(cfg['handshake_wait_s']):g} s for delete)"

    print(f"Camera      : {camera_note}")
    print(f"Format      : {cfg['format'] or preset['format']}")
    print(f"Devices     : {devices}")
    print(f"Duration    : {seconds:g} s")
    print(f"Session dir : {session_dir}")
    print(f"Anchor      : {anchor_note}")
    print(f"Trigger     : {trigger_note}")
    print(f"Pulse rate  : {pulse_note}")
    print(f"Handshake   : {handshake_note}")
    print("Opening cameras (nothing is recorded yet)...")

    # Deliberately NOT cleared here. Spike2 must delete the flag before ProgRun -- this process
    # takes seconds to import cv2, so anything it cleared would be cleared long after a Spike2
    # poll loop had already seen the stale file. All this can usefully do is say so.
    if cfg["handshake"] and spike2.flag_exists(cfg["ready_flag"]):
        print(f"  WARNING: {cfg['ready_flag']} already exists. Spike2 should FileDelete() it "
              f"before ProgRun, or it will read this run as ready before it is.",
              file=sys.stderr)

    readers = open_cameras(preset, devices, cfg, anchor)
    try:
        accepted = align.run_alignment(
            readers, devices, float(cfg["align_preview_hz"]), int(cfg["downsample"]),
            bool(cfg["rotate180"]), bool(cfg["fullscreen"]), seconds)
        if accepted is None:
            print("Nothing was recorded. Spike2 should halt.")
            return EXIT_CANCELLED

        # BEFORE the ready flag, deliberately: the flag is what makes Spike2 start its pulse
        # train, and a camera that is not yet in trigger mode when the first pulses arrive would
        # miss them. The wait for the first pulse comes after the flag, for the same reason.
        trig = attempt_trigger(readers, devices, cfg)

        if cfg["handshake"]:
            signal_ready(readers, cfg, session_dir)

        sess = session.Session(cfg["out"], session_id, len(devices))
        pending = None
        if trig["mode"] == "trigger":
            # mark_start BEFORE the wait, not inside the recorder: the real frames the wait
            # collects ARE pulse 0 onwards and are stored, so they and everything after them must
            # be stamped against one origin.
            sess.mark_start()
            ok, pending, seen, timeouts = record.wait_for_pulses(readers, accepted.order,
                                                                 trig["waitS"])
            if not ok:
                trig["aePriorityRestored"] = trigger.free_run_all(devices)
                trig["mode"] = "free-run"
                trig["reason"] = (f"no pulses within {trig['waitS']:g} s (real frames: "
                                  f"{record.per_cam_counts(seen)}; timeout frames: "
                                  f"{record.per_cam_counts(timeouts)})")
                print(f"  {trig['reason']}\n  Recording free-running with the strobe anchor "
                      f"instead.", file=sys.stderr)

        if trig["mode"] == "trigger":
            # No anchor under the trigger: there is nothing to mark when frame k is pulse k, and
            # the exposure has to stay well inside the pulse period.
            anchors = {"enabled": False, "reason": "trigger mode: frame k is pulse k"}
            stop_reason, trig_report = record.run_triggered(
                readers, devices, accepted.order, sess, seconds, pending, accepted.rotate180,
                float(cfg["rec_preview_hz"]), int(cfg["rec_downsample"]),
                bool(cfg["fullscreen"]), float(cfg["rec_window_scale"]), trig["endMarginS"])
            trig.update(trig_report)
        else:
            stop_reason, anchors = record.run_recording(
                readers, devices, accepted.order, sess, seconds, accepted.rotate180, anchor,
                float(cfg["rec_preview_hz"]), int(cfg["rec_downsample"]), bool(cfg["fullscreen"]),
                float(cfg["rec_window_scale"]))

        # Release the cameras HERE rather than leaving it to the finally, for two reasons. It is
        # the only point at which close_all's verdict can still reach session.json, which is
        # written a few lines below and is the only record that outlives this console. And it
        # frees the devices before the sidecars and TIFFs are written, which on a long session is
        # seconds of work with the cameras needlessly still claimed. The finally still runs;
        # close_all is a no-op the second time.
        stuck = capture.close_all(readers)

        info = finalise(sess, started, cfg, camera, detection, preset, devices, accepted,
                        readers, anchors, stop_reason, seconds, args.tag, stuck, trig)
        session.report(info)
        return EXIT_WARNINGS if info["warnings"] else EXIT_OK
    finally:
        capture.close_all(readers)
        # EVERY exit path, including the operator cancelling at the alignment stage (which
        # returns from inside the try, so this still runs) and any exception. The control
        # persists in the camera, so a run that left one triggered would give the NEXT ordinary
        # recording one black frame a second and nothing to explain it. After close_all, not
        # before: the write fails for about half a second after a release, and trigger.py waits
        # that out.
        if cfg["trigger"]:
            trigger.free_run_all(devices)


def _retry_busy(fn):
    """Run something that claims the cameras, waiting out a device a previous run still holds.

    Only CameraBusyError is retried, and deliberately so. That is the narrow class of failure
    where nothing is wrong with the request and nothing is wrong with the hardware -- the device
    simply has not been released yet, and Windows reclaims it on its own within a few seconds.
    Every other failure -- the wrong preset, a geometry the camera cannot deliver, a camera that
    matches no preset at all, an exposure change that re-negotiates the media type -- is a
    statement about the configuration that will be exactly as true on the third attempt as it was
    on the first, so it is raised immediately and the operator reads the message that says what
    to fix.

    Used by both stages that touch the devices before recording: identifying them and opening
    them. Neither leaves anything to clean up between attempts -- cameras.identify_camera
    releases the capture on every path, and capture.open_all closes whatever it managed to open
    before it re-raises, so a half-built set of cameras is never carried into the next try.
    """
    for attempt in range(1, OPEN_ATTEMPTS + 1):
        try:
            return fn()
        except cameras.CameraBusyError as exc:
            if attempt == OPEN_ATTEMPTS:
                print(f"\n  Still claimed after {OPEN_ATTEMPTS} attempts. Giving up.",
                      file=sys.stderr)
                raise
            print(f"\n  Attempt {attempt} of {OPEN_ATTEMPTS} failed:\n{exc}", file=sys.stderr)
            print(f"\n  Waiting {OPEN_RETRY_WAIT_S:g} s for Windows to reclaim the device, then "
                  f"retrying. Nothing has been recorded and Spike2 has not been signalled.",
                  file=sys.stderr)
            time.sleep(OPEN_RETRY_WAIT_S)


def open_cameras(preset, devices, cfg, anchor):
    """Open every camera, waiting out a device the previous run has not finished letting go of."""
    return _retry_busy(lambda: capture.open_all(preset, devices, cfg["format"], cfg["exposure"],
                                                anchor.bright if anchor.enabled else None))


def signal_ready(readers, cfg, session_dir):
    """Tell Spike2 that alignment is accepted and recording is about to start.

    Creating the flag is the whole signal. Spike2 sees it with FileStatus(), drops its
    SampleKey("S") marker, and carries on -- it does not have to answer. Line 1 of the flag is
    `session_dir`, so Spike2 can save its own recording alongside the frames; see spike2.py.

    If handshake_wait_s is set, the app additionally waits that long for Spike2 to DELETE the
    flag before storing anything, which lets Spike2 gate the start (see README). It then records
    regardless: a Spike2 that never answers must not be able to hang a session, and the strobe
    anchor identifies frame 0 either way.

    The queues are drained while waiting. The capture threads never stop, so without this they
    would hit their cap in about five seconds and start blocking -- harmless, since these frames
    are discarded anyway, but it leaves the driver dropping frames at exactly the moment the
    recording is about to start.
    """
    flag, wait_s = cfg["ready_flag"], float(cfg["handshake_wait_s"])
    spike2.write_flag(flag, session_dir, note="alignment accepted; recording starts now")
    print(f"\nCameras hot. Told Spike2 via:\n    {flag}")
    if wait_s <= 0:
        return

    print(f"  waiting up to {wait_s:g} s for Spike2 to delete it...")
    deadline = time.perf_counter() + wait_s
    while spike2.flag_exists(flag) and time.perf_counter() < deadline:
        capture.check_readers(readers)
        for r in readers:
            capture.drain_newest(r)
        time.sleep(0.01)
    print("  Spike2 acknowledged; recording now." if not spike2.flag_exists(flag)
          else "  no acknowledgement, recording anyway.")


def attempt_trigger(readers, devices, cfg):
    """Put every camera into FSIN trigger mode. Returns the trigger block for session.json.

    `mode` is "trigger" only when EVERY device was written and read the value back. Anything else
    -- the feature switched off, a COM failure, a driver that answered with a different value --
    puts every device back to free-run and returns "free-run" with the reason, and the run
    records exactly as it always did, with the strobe anchor. A trigger that half worked is the
    one outcome that must not happen: one camera pulsing and one free-running produces two files
    whose frame indices mean different things.

    Called AFTER the operator accepts and BEFORE the ready flag, because the flag is what starts
    Spike2's pulse train. The queues are drained here too, for the reason run_recording drains
    its own: everything in them is free-running video from before the switch, and stored under
    the trigger it would claim to be a pulse.

    THE DRAIN WAITS trigger.MODE_SETTLE_S FIRST. Free-running frames keep landing for a moment
    after the write (see the constant), so a drain that ran straight away would leave them
    queued and the recording would store them as pulse 0, 1, 2. The quarter second is free:
    Spike2 is not sampling until the ready flag, which is written after this returns.

    Whether any pulses actually ARRIVE is a separate question, answered later by
    record.wait_for_pulses -- this only proves the camera was told.
    """
    trig = {"requested": bool(cfg["trigger"]), "mode": "free-run", "reason": "",
            "aePriority": None, "pulseHz": cfg["pulse_hz"],
            "waitS": float(cfg["trigger_wait_s"]),
            "endMarginS": float(cfg["trigger_end_margin_s"])}
    if not cfg["trigger"]:
        trig["reason"] = "trigger disabled in config"
        return trig

    print("\nWriting the FSIN trigger control on every camera...")
    readback = trigger.set_all(devices, 1)
    trig["aePriority"] = readback
    if not trigger.all_agree(readback, 1):
        trig["reason"] = _write_failure(devices, readback)
        trig["aePriorityRestored"] = trigger.free_run_all(devices)
        print(f"  {trig['reason']}\n  Recording free-running with the strobe anchor instead.",
              file=sys.stderr)
        return trig

    time.sleep(trigger.MODE_SETTLE_S)
    for r in readers:
        capture.drain_newest(r)
    trig["mode"] = "trigger"
    return trig


def _write_failure(devices, readback):
    """Name the FIRST device the trigger write did not take on, and what it answered.

    The first rather than all of them: the run falls back whichever it was, and one device and
    one reason is what fits on the console line the operator actually reads.
    """
    for device in devices:
        entry = readback.get(str(device), {})
        if "error" in entry:
            detail = entry["error"]
        elif entry.get("value") != 1:
            detail = f"the driver read it back as {entry.get('value')}, not 1"
        else:
            continue
        return f"trigger write failed on device {device}: {detail}"
    return "trigger write failed: no device was written"


def finalise(sess, started, cfg, camera, detection, preset, devices, accepted, readers,
             anchors, stop_reason, seconds, tag, stuck=(), trig=None):
    """Write everything except the frames, session.json last of all."""
    sess.write_sidecars()
    tiffs = sess.write_first_frame_tiffs() if cfg["first_frame_tiff"] else []
    request = {
        # "camera" is the preset that was USED, always a preset name, because that is what
        # every downstream reader looks up. What was asked for is beside it.
        "camera": camera, "cameraRequested": cfg["camera"], "cameraDetection": detection,
        "format": cfg["format"] or preset["format"],
        "deviceIds": devices, "seconds": seconds,
        "recPreviewHz": cfg["rec_preview_hz"], "tag": tag,
        "source": preset["source"], "queueCap": capture.QUEUE_MAX,
        "settleS": capture.SETTLE_S, "handshake": bool(cfg["handshake"]),
        "trigger": bool(cfg["trigger"]), "pulseHz": cfg["pulse_hz"],
    }
    return session.summarise(sess, started, request, readers, devices, accepted.order,
                             anchors, stop_reason, tiffs, accepted.rotate180, stuck, trig)


def parse_args(argv):
    """Everything defaults to None so config.json wins unless the flag was actually given."""
    p = argparse.ArgumentParser(
        description="Align a camera pair, then record it. One process, driven by Spike2.")
    p.add_argument("--config", default=None,
                   help="path to config.json (default: the one beside this script, if present)")
    p.add_argument("--session-id", default=None,
                   help="name of the session directory. PASS THIS FROM SPIKE2: it is how Spike2 "
                        "knows where to look for session.json afterwards")
    p.add_argument("--tag", default=None, help="free text appended to a generated session id")

    p.add_argument("--camera", default=None,
                   choices=sorted(cameras.CAMERA_PRESETS) + [cameras.AUTO],
                   help="preset name, or auto to identify the family from the geometries the "
                        "device offers (default from config.json)")
    p.add_argument("--devices", type=int, nargs="+", default=None,
                   help="1-based device ids; the ORDER is the initial left-to-right arrangement")
    p.add_argument("--format", default=None, help="e.g. MJPG_1280x720")
    p.add_argument("--exposure", type=float, default=None)
    p.add_argument("--seconds", type=float, default=None)
    p.add_argument("--out", default=None, help="parent directory for the session")

    p.add_argument("--rotate180", action="store_true", default=None,
                   help="start with the view rotated; if accepted, STORED frames are rotated too")
    p.add_argument("--no-rotate180", dest="rotate180", action="store_false", default=None)
    p.add_argument("--downsample", type=int, default=None,
                   help="alignment preview pixel stride only")
    p.add_argument("--rec-downsample", type=int, default=None,
                   help="recording preview pixel stride; coupled to --rec-window-scale")
    p.add_argument("--rec-window-scale", type=float, default=None,
                   help="recording window as a fraction of the screen (alignment uses 0.92)")
    p.add_argument("--align-preview-hz", type=float, default=None)
    p.add_argument("--rec-preview-hz", type=float, default=None)
    p.add_argument("--no-preview", dest="rec_preview_hz", action="store_const", const=0.0,
                   default=None, help="no window while recording")
    p.add_argument("--fullscreen", action="store_true", default=None)

    p.add_argument("--no-handshake", dest="handshake", action="store_false", default=None,
                   help="skip the Spike2 flag and record immediately after accept (testing)")
    p.add_argument("--ready-flag", default=None, help="path of the handshake flag file")
    p.add_argument("--handshake-wait-s", type=float, default=None,
                   help="wait this long for Spike2 to DELETE the flag before storing, then "
                        "record regardless; 0 records as soon as the flag is written")

    p.add_argument("--trigger", action="store_true", default=None,
                   help="after accept, put both cameras into FSIN external trigger mode: one "
                        "stored frame per pulse (default from config.json)")
    p.add_argument("--no-trigger", dest="trigger", action="store_false", default=None,
                   help="record free-running with the strobe anchor, as before")
    p.add_argument("--trigger-wait-s", type=float, default=None,
                   help="wait this long for the first FSIN pulse on every camera; if none "
                        "arrives, fall back to a free-running recording")
    p.add_argument("--trigger-end-margin-s", type=float, default=None,
                   help="how long past --seconds a triggered recording may run before it stops "
                        "without the pulse train having ended")
    p.add_argument("--pulse-hz", type=float, default=None,
                   help="the FSIN pulse rate Spike2 is driving. PASS THIS FROM SPIKE2: it is "
                        "recorded, and the delivered rate is warned about if it disagrees")

    p.add_argument("--no-anchor", dest="anchor", action="store_false", default=None,
                   help="do not bracket the recording with the exposure/strobe anchor")
    p.add_argument("--anchor-exposure", type=float, default=None)
    p.add_argument("--anchor-hold-s", type=float, default=None)

    args = p.parse_args(argv)
    if args.config is None:
        beside = Path(__file__).resolve().parent / "config.json"
        args.config = str(beside) if beside.exists() else None
    return args


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\nInterrupted.")
        sys.exit(EXIT_CANCELLED)
    except Exception as exc:
        print(f"\nFAILED: {exc}", file=sys.stderr)
        sys.exit(EXIT_FAILED)
