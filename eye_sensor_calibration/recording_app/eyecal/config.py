"""Settings: built-in defaults, overridden by config.json, overridden by the command line.

Spike2 should pass only what varies between runs (--session-id, usually --seconds), and leave
the rig's setup in config.json where it can be changed without editing a command line inside a
.s2s script.
"""

import json
from pathlib import Path

DEFAULTS = {
    # what to open
    "camera": "elp",                    # preset name; see eyecal/cameras.py
    "devices": [1, 2],                  # 1-based, as MATLAB's winvideo numbers them
    "format": None,                     # override the preset's native format, e.g. MJPG_1280x720
    "exposure": None,                   # override the preset's Exposure value

    # what to record
    "seconds": 30.0,
    "out": "./recordings",
    "rotate180": False,                 # starting state; 'r' during alignment changes it, and
                                        # whatever is accepted rotates the STORED frames
    "first_frame_tiff": True,

    # what to show
    "align_preview_hz": 25.0,
    "rec_preview_hz": 10.0,             # 0 disables the recording preview entirely
    "downsample": 2,                    # alignment display pixel stride; never affects frames

    # The recording preview is deliberately smaller than the alignment one -- it is a progress
    # indicator, not something to aim a camera by, and a half-size window leaves room for Spike2.
    # THESE TWO ARE COUPLED. rec_window_scale shrinks the window; rec_downsample must shrink the
    # canvas at least as much, or the canvas ends up larger than the window and every refresh
    # takes the expensive downscale path (measured 7x per output pixel). Change one, check the
    # other -- the window says so on stderr if the pairing is wrong. At 0.5 on a 1920x1080
    # screen, a two-camera 1600x1200 pair needs rec_downsample 4.
    "rec_downsample": 4,
    "rec_window_scale": 0.5,            # fraction of the screen; alignment uses 0.92
    "fullscreen": False,

    # Spike2 handshake (see eyecal/spike2.py)
    "handshake": True,
    "ready_flag": "C:/Temp/eyecal_ready.flag",
    "handshake_wait_s": 0.0,            # >0 waits for Spike2 to DELETE the flag before
                                        # storing, then records anyway; 0 records immediately

    # strobe anchor (see eyecal/record.py)
    "anchor": True,
    "anchor_exposure": None,            # None = normal + 4 stops
    "anchor_hold_s": 0.25,
}


def load(path=None, required=False):
    """Defaults, with config.json merged over the top.

    An unknown key is an error rather than a shrug. A typo in a config file that is silently
    ignored is the classic way to spend an afternoon wondering why a setting did nothing.
    """
    cfg = dict(DEFAULTS)
    if path is None:
        return cfg
    p = Path(path)
    if not p.exists():
        if required:
            raise FileNotFoundError(f"config file not found: {p}")
        return cfg

    loaded = json.loads(p.read_text(encoding="utf-8"))
    unknown = sorted(set(loaded) - set(DEFAULTS) - {"_comment"})
    if unknown:
        raise KeyError(f"{p}: unknown setting(s) {', '.join(unknown)}. "
                       f"Known settings: {', '.join(sorted(DEFAULTS))}")
    cfg.update({k: v for k, v in loaded.items() if k != "_comment"})
    return cfg


def apply_cli(cfg, args):
    """Overlay the command line. Only arguments the user actually gave (not None) win."""
    for key in DEFAULTS:
        value = getattr(args, key, None)
        if value is not None:
            cfg[key] = value
    return cfg
