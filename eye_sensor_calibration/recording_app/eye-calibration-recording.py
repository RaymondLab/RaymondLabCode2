"""Align a camera pair, then record them. One process, cameras never released in between.

    python eye-calibration-recording.py --session-id m123_run01 --seconds 30

Run by Spike2 as a single step. The phases are:

    1. open both cameras and negotiate the format          Spike2 NOT sampling
    2. live alignment view; the operator presses 'a'       Spike2 NOT sampling
    3. create the ready flag                               Spike2 sees it and marks the file
    4. record, bracketed by the strobe anchor              Spike2 sampling
    5. write session.json                                  Spike2 sees it and continues

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
    1  failed (camera would not open, wrong format, device not streaming)
    2  bad command line (argparse)
    3  recorded, but a guard tripped -- read the warnings in session.json
    4  cancelled by the operator at the alignment stage
"""

import argparse
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from eyecal import align, cameras, capture, config, record, session, spike2

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

    if cfg["camera"] not in cameras.CAMERA_PRESETS:
        raise KeyError(f"Unknown camera {cfg['camera']!r}. "
                       f"Known: {', '.join(sorted(cameras.CAMERA_PRESETS))}")
    preset = cameras.CAMERA_PRESETS[cfg["camera"]]
    devices = [int(d) for d in cfg["devices"]]
    seconds = float(cfg["seconds"])
    if seconds <= 0:
        raise ValueError(f"seconds must be positive, got {seconds}")

    started = datetime.now(timezone.utc)
    session_id = args.session_id or (
        started.strftime("%Y%m%d-%H%M%SZ") + f"_{cfg['camera']}_{len(devices)}cam_{seconds:g}s"
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
    handshake_note = "off (recording starts straight after accept)"
    if cfg["handshake"]:
        handshake_note = cfg["ready_flag"]
        if float(cfg["handshake_wait_s"]) > 0:
            handshake_note += f"  (wait up to {float(cfg['handshake_wait_s']):g} s for delete)"

    print(f"Camera      : {preset['title']}")
    print(f"Format      : {cfg['format'] or preset['format']}")
    print(f"Devices     : {devices}")
    print(f"Duration    : {seconds:g} s")
    print(f"Session dir : {session_dir}")
    print(f"Anchor      : {anchor_note}")
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

        if cfg["handshake"]:
            signal_ready(readers, cfg, session_dir)

        sess = session.Session(cfg["out"], session_id, len(devices))
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

        info = finalise(sess, started, cfg, preset, devices, accepted, readers, anchors,
                        stop_reason, seconds, args.tag, stuck)
        session.report(info)
        return EXIT_WARNINGS if info["warnings"] else EXIT_OK
    finally:
        capture.close_all(readers)


def open_cameras(preset, devices, cfg, anchor):
    """Open every camera, waiting out a device the previous run has not finished letting go of.

    Only CameraBusyError is retried, and deliberately so. That is the narrow class of failure
    where nothing is wrong with the request and nothing is wrong with the hardware -- the device
    simply has not been released yet, and Windows reclaims it on its own within a few seconds.
    Every other opening failure -- the wrong preset, a geometry the camera cannot deliver, an
    exposure change that re-negotiates the media type -- is a statement about the configuration
    that will be exactly as true on the third attempt as it was on the first, so it is raised
    immediately and the operator reads the message that says what to fix.

    Nothing has to be cleaned up between attempts: capture.open_all closes whatever it managed to
    open before it re-raises, so a half-built set of cameras is never carried into the next try.
    """
    for attempt in range(1, OPEN_ATTEMPTS + 1):
        try:
            return capture.open_all(preset, devices, cfg["format"], cfg["exposure"],
                                    anchor.bright if anchor.enabled else None)
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


def finalise(sess, started, cfg, preset, devices, accepted, readers, anchors, stop_reason,
             seconds, tag, stuck=()):
    """Write everything except the frames, session.json last of all."""
    sess.write_sidecars()
    tiffs = sess.write_first_frame_tiffs() if cfg["first_frame_tiff"] else []
    request = {
        "camera": cfg["camera"], "format": cfg["format"] or preset["format"],
        "deviceIds": devices, "seconds": seconds,
        "recPreviewHz": cfg["rec_preview_hz"], "tag": tag,
        "source": preset["source"], "queueCap": capture.QUEUE_MAX,
        "settleS": capture.SETTLE_S, "handshake": bool(cfg["handshake"]),
    }
    return session.summarise(sess, started, request, readers, devices, accepted.order,
                             anchors, stop_reason, tiffs, accepted.rotate180, stuck)


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

    p.add_argument("--camera", default=None, choices=sorted(cameras.CAMERA_PRESETS))
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
