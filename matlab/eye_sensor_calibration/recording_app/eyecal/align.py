"""The alignment phase: position the cameras, then accept or cancel.

Records nothing. Draws a crosshair at the exact centre of each camera's OWN frame so the
operator can position each camera until the feature of interest sits on the intersection.

Ends one of two ways, and the caller distinguishes them by the return value:

    'a'                     accept   -> returns the arrangement, recording follows
    Esc, or closing window  cancel   -> returns None, the app exits and Spike2 halts

The cameras are NOT released on accept. They stay open and streaming straight into the recording
phase, which is the whole reason this app is one process: every strobe pulse from opening,
negotiating and aligning happens before Spike2 starts sampling, so none of it lands in the file.
"""

import collections
import time

from . import capture, display

Alignment = collections.namedtuple("Alignment", "order rotate180 seconds")

KEYS_HELP = "a accept | Esc cancel | c colour | r rotate | s swap sides | h hide overlay"


def run_alignment(readers, devices, preview_hz=25.0, downsample=2, rotate180=False,
                  fullscreen=False, seconds=None):
    """Show the live alignment view. Returns an Alignment on accept, or None on cancel.

    `order[j]` is the reader drawn in pane j, left to right. It is kept apart from the reader
    list so 's' rearranges the display without disturbing the threads, the meters or the device
    ids -- all of which stay indexed by reader throughout. That same order then decides which
    camera becomes cam 1 / c1.bin, which is the one thing this phase exists to decide.
    """
    downsample = max(int(downsample), 1)
    colour_idx, show_overlay, rot = 0, True, bool(rotate180)     # colour_idx 0 == green
    order = list(range(len(readers)))
    meters = [capture.RateMeter() for _ in readers]
    # Shorter window than the camera meters: at 25 Hz, 60 samples would smooth over 2.4 s and
    # hide exactly the brief stalls this readout exists to make visible.
    display_meter = capture.RateMeter(window=30)
    latest = [None] * len(readers)
    period = 1.0 / preview_hz if preview_hz > 0 else 0.0
    t0, next_at = time.perf_counter(), 0.0

    window = display.Window("camera alignment -- records nothing", fullscreen)
    print("\nAligning. Put the target on the crosshair intersection of BOTH cameras.")
    print(f"Keys: {KEYS_HELP}")
    print("Closing the window cancels the whole run; Spike2 will stop.")
    print(f"Left to right: devices {[devices[i] for i in order]}\n")

    try:
        while True:
            capture.check_readers(readers)
            # Drain and discard all but the newest. This is what keeps the capture threads from
            # ever hitting their queue cap during a long alignment session.
            for k, r in enumerate(readers):
                item = capture.drain_newest(r)
                if item is not None:
                    meters[k].tick(item[0])
                    latest[k] = item[2]

            t = time.perf_counter() - t0
            if t >= next_at and all(f is not None for f in latest):
                display_meter.tick(t)      # before drawing, so the rate shown includes this one
                panes = []
                for pane, k in enumerate(order, start=1):
                    img = latest[k]
                    panes.append(display.alignment_pane(
                        display.display_copy(img, downsample, rot), pane, devices[k],
                        display.COLOURS[colour_idx][1], show_overlay, meters[k].fps,
                        display_meter.fps, img.shape))
                canvas = display.tile(panes)
                if show_overlay:
                    label = "NOT RECORDING -- alignment only.  " + KEYS_HELP
                    if seconds:
                        label = f"NOT RECORDING.  'a' accepts and records {seconds:g} s.  " \
                                f"{KEYS_HELP}"
                    display.text(canvas, label, (8, canvas.shape[0] - 12),
                                 display.COLOURS[colour_idx][1], 0.45)
                window.show(canvas)
                next_at = t + period

            key = window.pump()
            if key == ord("a"):
                print(f"\nAccepted: left to right is devices {[devices[i] for i in order]}"
                      f"{', display rotated 180 deg' if rot else ''}.")
                return Alignment(order=list(order), rotate180=rot, seconds=seconds)
            elif key == 27:                                     # Esc
                print("\nCancelled (Esc).")
                return None
            elif key == ord("c"):
                colour_idx = (colour_idx + 1) % len(display.COLOURS)
                print(f"  overlay colour: {display.COLOURS[colour_idx][0]}")
            elif key == ord("r"):
                rot = not rot
                print(f"  rotation: {'180 deg' if rot else 'none'} "
                      f"(carries into the STORED frames)")
            elif key == ord("s"):
                order.reverse()
                print(f"  sides swapped -- left to right is now devices "
                      f"{[devices[i] for i in order]}")
            elif key == ord("h"):
                show_overlay = not show_overlay
                print(f"  overlay: {'shown' if show_overlay else 'hidden (crosshair too)'}")

            # AFTER the waitKey that delivers the close message, and BEFORE the next redraw --
            # see Window.closed for why both orderings are load-bearing.
            if window.closed():
                print("\nCancelled (window closed).")
                return None
    finally:
        window.close()
