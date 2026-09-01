"""Everything drawn on screen: pane overlays, tiling, and one window wrapper.

Nothing in here is a measurement. The preview exists so an operator can aim a camera and see
that a recording is running -- the stored frames are untouched by any of it. It still has to be
geometrically honest, because a stretched alignment view moves the apparent centre of a feature
away from the crosshair that exists to mark it, which is the one error this app must not cause.
"""

import sys

import cv2
import numpy as np

# GREEN FIRST: it is the default overlay colour and the one the operator expects on startup.
# 'c' cycles to the others for backgrounds green disappears against.
COLOURS = [("green", (0, 255, 0)), ("magenta", (255, 0, 255)),
           ("white", (255, 255, 255)), ("amber", (0, 191, 255))]

RED = (0, 0, 255)


def display_copy(img, downsample=1, rotate180=False):
    """A contiguous, cheap-to-draw copy of one frame for display only.

    ascontiguousarray because a column-strided slice is not a layout cv2 accepts, and cv2.rotate
    rather than np.rot90 because the latter returns a negative-stride view that some OpenCV calls
    reject outright.
    """
    out = img if downsample <= 1 else np.ascontiguousarray(img[::downsample, ::downsample])
    return cv2.rotate(out, cv2.ROTATE_180) if rotate180 else out


def text(canvas, s, org, colour, scale=0.5):
    """LINE_AA is right for text: smoothing helps, and no single pixel here carries meaning."""
    cv2.putText(canvas, s, org, cv2.FONT_HERSHEY_SIMPLEX, scale, colour, 1, cv2.LINE_AA)


def line(canvas, p0, p1, colour):
    """LINE_8, NOT LINE_AA. Anti-aliasing spreads a 1-pixel line across three columns at roughly
    20/90/20 intensity, so the brightest pixel is no longer unambiguously the centre. For a
    reticle whose whole job is to mark one exact column and one exact row, that is a defect.

    There is no black underlay either: it bought contrast against a bright iris at the cost of a
    3-pixel-wide mark where a 1-pixel one was meant. Cycling the colour is the answer to a
    background the reticle disappears against.
    """
    cv2.line(canvas, p0, p1, colour, 1, cv2.LINE_8)


def _canvas_from(img):
    """A drawable 3-channel copy. cvtColor already returns a fresh array; .copy() covers the
    case where the frame arrives in colour and must not be drawn on in place."""
    return cv2.cvtColor(img, cv2.COLOR_GRAY2BGR) if img.ndim == 2 else img.copy()


def alignment_pane(img, pane, device, colour, show, cam_fps, display_hz, full_shape):
    """Crosshair and readout for ONE camera's frame, at THAT frame's own centre.

    Drawing one crosshair across the composited image would put both marks in the wrong place,
    which is exactly the mistake this view exists to prevent.

    `img` is the already-downsampled display copy and `full_shape` the true sensor geometry, so
    the readout reports real pixel coordinates. `pane` is the position (1 = leftmost) and
    `device` the physical camera filling it; both are drawn, because 's' pulls them apart and
    seeing which device you just moved where is the entire point of that key. `cam_fps` is the
    frame arrival rate and `display_hz` the redraw rate -- unrelated, routinely 2x apart, so both
    are labelled.
    """
    canvas = _canvas_from(img)
    h, w = canvas.shape[:2]
    # Integer centre of the displayed image. For an even dimension the true centre falls between
    # two pixels; w // 2 puts the line on the first pixel right of it, identically for every
    # camera, which is what matters for a relative alignment.
    cx, cy = w // 2, h // 2

    if show:
        full_h, full_w = full_shape[:2]
        # Readout first, so on a small display the text can never paint over the reticle: the
        # reticle is the one element that must be exactly where it claims to be.
        text(canvas, f"cam {pane}  (device {device})   {full_w}x{full_h}   "
                     f"center ({full_w // 2}, {full_h // 2}) px", (8, 20), colour)
        text(canvas, f"{cam_fps:5.1f} fps camera   {display_hz:5.1f} Hz display", (8, 40), colour)
        line(canvas, (0, cy), (w, cy), colour)
        line(canvas, (cx, 0), (cx, h), colour)

    # The border stays regardless of the overlay toggle: with panes tiled edge to edge it is what
    # shows where one camera's frame ends and the next begins, and it obscures nothing.
    cv2.rectangle(canvas, (0, 0), (w - 1, h - 1), colour, 1)
    return canvas


def recording_pane(img, cam, device, frames, queue_depth):
    """One camera's pane while recording. Red, and says where the frames are going.

    cam 1 is only cam 1 by virtue of the accepted device order, and that order is the thing that
    silently changes between sessions -- so the pane states both, and the file it feeds.
    """
    canvas = _canvas_from(img)
    text(canvas, f"cam {cam} (device {device}) -> c{cam}.bin   {frames} frames   "
                 f"queue {queue_depth}", (8, 20), RED)
    cv2.rectangle(canvas, (0, 0), (canvas.shape[1] - 1, canvas.shape[0] - 1), RED, 1)
    return canvas


def tile(panes):
    """Panes side by side in ONE window. n windows would be n blits, n sets of window chrome and
    n compositor surfaces."""
    if len(panes) == 1:
        return panes[0]
    height = max(p.shape[0] for p in panes)
    canvas = np.zeros((height, sum(p.shape[1] for p in panes), 3), np.uint8)
    x = 0
    for p in panes:
        canvas[:p.shape[0], x:x + p.shape[1]] = p
        x += p.shape[1]
    return canvas


class Window:
    """One highgui window: sized once, letterboxed on every frame, closable by the operator.

    A class rather than loose functions because the app has two of these -- alignment and
    recording -- and every call needs the window name. Bundling the name with the "have I been
    sized yet" flag keeps both out of the phase loops.
    """

    def __init__(self, name, fullscreen=False, screen_fraction=0.92):
        self.name = name
        self.fullscreen = fullscreen
        self.screen_fraction = screen_fraction
        self._sized = False
        self._warned_downscale = False
        # WINDOW_KEEPRATIO is requested but not relied on: measured on OpenCV 5.0.0 Win32 it has
        # no effect. Aspect is preserved by _letterbox instead, where it can be verified. The
        # flag costs nothing and is correct on backends that do honour it.
        cv2.namedWindow(name, cv2.WINDOW_NORMAL | cv2.WINDOW_KEEPRATIO)
        if fullscreen:
            cv2.setWindowProperty(name, cv2.WND_PROP_FULLSCREEN, cv2.WINDOW_FULLSCREEN)

    def show(self, canvas):
        """Size on the first canvas -- its shape depends on the camera count and downsample, so
        it is not known at construction -- then letterbox and blit. Sized once only, so a window
        the operator has resized stays resized."""
        if not self._sized:
            self._size_to(canvas.shape[1], canvas.shape[0])
            self._sized = True
        cv2.imshow(self.name, self._letterbox(canvas))

    def pump(self):
        """EXACTLY ONE waitKey per loop iteration, and the key it returns must be dispatched by
        the caller. waitKey POPS a key rather than peeking at it, so a second call that only
        recognises one key does not politely ignore the rest -- it swallows them. It must also
        run on iterations that drew nothing: waitKey is what pumps the window's event loop and
        blits, and without it imshow queues work that never happens."""
        return cv2.waitKey(1) & 0xFF

    def closed(self):
        """True once the operator has closed the window with the title-bar X.

        Measured on OpenCV 5.0.0 Win32: WND_PROP_VISIBLE reads 1.0 from the moment namedWindow
        returns, before any imshow, so there is no startup false positive. Must be called AFTER
        the waitKey that delivers the close message and BEFORE the next show(): imshow recreates
        a closed window on this build, so a check after the redraw would find a resurrected
        window and never report the close. Guarded because a destroyed window is a NULL window to
        highgui and other calls raise rather than answer.
        """
        try:
            return cv2.getWindowProperty(self.name, cv2.WND_PROP_VISIBLE) < 1
        except cv2.error:
            return True

    def close(self):
        try:
            cv2.destroyWindow(self.name)
        except cv2.error:
            pass
        cv2.waitKey(1)      # highgui needs one more pump to actually destroy the window

    def _letterbox(self, canvas):
        """Scale the canvas to fill the window without distorting it, padding with black.

        Sizing the output to the window's image area EXACTLY leaves highgui's own scaling a 1:1
        copy with nothing left to stretch. Upscaling uses INTER_NEAREST for the same reason line()
        uses LINE_8: interpolating a 1-pixel reticle spreads it over several rows and destroys the
        unambiguous centre. Only the downscale path, where aliasing is the greater risk, uses
        INTER_AREA.
        """
        try:
            _, _, w, h = cv2.getWindowImageRect(self.name)
        except cv2.error:
            return canvas                               # window gone; let imshow deal with it
        ch, cw = canvas.shape[:2]
        if w <= 0 or h <= 0:
            return canvas                               # minimised, or not yet mapped

        scale = min(w / cw, h / ch)
        tw, th = max(int(round(cw * scale)), 1), max(int(round(ch * scale)), 1)
        if scale < 1.0 and not self._warned_downscale:
            # Measured: INTER_AREA costs 5.75 ms per output megapixel against INTER_NEAREST's
            # 0.78, so landing on the downscale path is roughly 7x dearer PER PIXEL and easily
            # outweighs having fewer pixels to write. It happens when the canvas is bigger than
            # the window, which is what a shrunken window and an unchanged downsample produce
            # together. Said once, naming the fix, because it is a config mistake and not a
            # runtime fault.
            self._warned_downscale = True
            print(f"  note: {self.name} is showing a {cw}x{ch} canvas in a {w}x{h} window, so "
                  f"every refresh downscales (the expensive interpolation). Raise the downsample "
                  f"for this phase until the canvas fits inside the window.", file=sys.stderr)
        scaled = cv2.resize(canvas, (tw, th),
                            interpolation=cv2.INTER_AREA if scale < 1.0 else cv2.INTER_NEAREST)
        if (tw, th) == (w, h):
            return scaled                               # already at the canvas aspect; no bars
        padded = np.zeros((h, w, 3), canvas.dtype)
        y0, x0 = (h - th) // 2, (w - tw) // 2
        padded[y0:y0 + th, x0:x0 + tw] = scaled
        return padded

    def _size_to(self, canvas_w, canvas_h):
        """Largest box fitting `screen_fraction` of the screen, at the canvas's own aspect.

        Purely cosmetic, so nothing here may propagate -- losing a session to window geometry
        would be absurd. ctypes rather than a GUI toolkit because it is stdlib and this app is
        Windows-only anyway through the DirectShow backend. SetProcessDPIAware is deliberately
        NOT called: left DPI-unaware, GetSystemMetrics and highgui's window coordinates stay in
        the same virtual-pixel space.

        `screen_fraction` is not free to choose alone: shrinking the window without shrinking the
        canvas to match pushes _letterbox onto its expensive path. See the note there.
        """
        if self.fullscreen:
            return
        try:
            import ctypes
            user32 = ctypes.windll.user32
            sw, sh = int(user32.GetSystemMetrics(0)), int(user32.GetSystemMetrics(1))
            if sw <= 0 or sh <= 0:
                sw, sh = 1920, 1080
            frac = self.screen_fraction
            scale = min(sw * frac / canvas_w, sh * frac / canvas_h)
            cv2.resizeWindow(self.name, max(int(canvas_w * scale), 320),
                             max(int(canvas_h * scale), 240))
        except Exception as exc:
            print(f"  could not size the preview window ({exc}); using its default size",
                  file=sys.stderr)
