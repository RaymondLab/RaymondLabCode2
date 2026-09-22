"""Reading a camera's DirectShow controls DIRECTLY, including whether each one is on AUTO.

OpenCV cannot answer this. CAP_PROP_AUTO_EXPOSURE reads back -1 on these cameras and cap.set()
returns True for values the driver ignores, which is why cameras.secure_manual_exposure proves
manual exposure by measuring the image instead. That measurement costs a couple of seconds and
speaks about exposure only. The two DirectShow control interfaces underneath OpenCV --
IAMCameraControl and IAMVideoProcAmp -- carry a per-property Auto/Manual flag, plus the driver's
range and default, for every property the device exposes, and answer instantly. This reads them.

ONE WRITE, and it is Controls.set. Everything else here reads. The write exists for exactly one
control -- AE_PRIORITY, the external-trigger switch below -- which no other path in this code base
can reach: OpenCV has no CAP_PROP for it. The app's exposure path keeps its own mechanism, which
it has proven by measurement, and nothing here touches it.

DEVICE NUMBERING. Index k here is the k-th moniker of CLSID_VideoInputDeviceCategory, which is
the same enumeration, in the same order, that OpenCV's DirectShow backend counts its capture
index with. cameras._open_live turns the app's 1-based device id into an OpenCV index with
`device_id - 1`, so the app's device N is index N-1 here too. Nothing cross-checks that at run
time and nothing can: the two Arducams share a USB VID/PID and carry no serial number, so the
enumeration holds no field that would confirm it (cameras.identify_camera says why no PnP data
is trusted anywhere in this app).

MEASURED ON THIS RIG, 2026-09-07, two Arducam OV9281s, OpenCV 5.0.0 DSHOW, comtypes 1.4.16:

  - BINDING WORKS WHILE OPENCV IS STREAMING. The filter binds and every Get/GetRange answers
    before an OpenCV capture is opened, while it is open and delivering frames, and after it is
    released, with the same numbers in all three. Both orders were tried -- bound first and kept,
    and bound afterwards -- and neither disturbed the other. So nothing has to be held across the
    OpenCV open, and this is not a claim on the device in the sense that makes a CameraBusyError.
  - WHAT AN OV9281 EXPOSES. IAMCameraControl: Exposure only. IAMVideoProcAmp: Brightness,
    Contrast, Hue, Saturation, Sharpness, Gamma, WhiteBalance, BacklightCompensation, Gain. Pan,
    Tilt, Roll, Zoom, Iris, Focus and ColorEnable all fail GetRange and are therefore not
    reported. Exposure reads -13..-1 step 1 default -6, in log2 seconds exactly as OpenCV's
    CAP_PROP_EXPOSURE does. Out of the box both cameras report WhiteBalance on AUTO (flags=1)
    and everything else MANUAL (flags=2).
  - THE FLAG SAYS WHAT THE DRIVER WAS LAST TOLD, at the instant it is read. Writing OpenCV's
    CAP_PROP_AUTO_EXPOSURE 0.75 or 1.0 reads back AUTO here immediately afterwards; 0.25 or 0.0
    reads back MANUAL; and ANY CAP_PROP_EXPOSURE write puts it back to MANUAL, because the
    backend sets a value and the Manual flag together. So a readout taken after an exposure
    write always says MANUAL and says nothing about what the sensor is doing.
  - AGAINST THE MEASUREMENT. On both cameras as they stand, the image responded 2.6x and 3.0x
    across the -13..-6 probe on every AUTO_EXPOSURE setting -- manual control by
    cameras.MIN_RESPONSE_RATIO -- and the flag agreed, reading MANUAL. UNVERIFIED: whether a
    camera pinned in its OWN auto-exposure by the Windows Frame Server reports flags=1. None was
    available to test against, and it must not be reproduced by running the Camera app on the rig
    cameras. Treat the flag as a fast hint and keep cameras.secure_manual_exposure as the proof.

COM lives on the main thread. Everything here must be called from the thread that created the
objects -- these are apartment-threaded proxies, and a capture thread that touched one would be
a marshalling bug that surfaces as a random E_FAIL much later.
"""

from ctypes import POINTER, c_long, c_ulong

import comtypes
from comtypes import COMMETHOD, GUID, HRESULT, IUnknown
from comtypes.automation import VARIANT
from comtypes.persist import IPersist, IPropertyBag

# The system enumerator, the interface that hands out category enumerators, and the category
# that holds capture devices. Fixed by DirectShow; not ours to choose.
CLSID_SystemDeviceEnum = GUID("{62BE5D10-60EB-11d0-BD3B-00A0C911CE86}")
CLSID_VideoInputDeviceCategory = GUID("{860BB310-5D01-11d0-BD3B-00A0C911CE86}")
IID_IBaseFilter = GUID("{56A86895-0AD4-11CE-B03A-0020AF0BA770}")

# Set by BOTH control interfaces, and the whole reason this module exists.
FLAG_AUTO = 0x0001
FLAG_MANUAL = 0x0002

# THE EXTERNAL-TRIGGER SWITCH of the Arducam OV9281 UVC firmware. Property id 19 on
# IAMCameraControl is KSPROPERTY_CAMERACONTROL_AUTO_EXPOSURE_PRIORITY, the UVC
# CT_AE_PRIORITY_CONTROL. Windows draws it as the "Low Light Compensation" checkbox on the Camera
# Control tab, Arducam's application note calls it "low-brightness compensation", and Linux calls
# it exposure_dynamic_framerate. On these cameras 1 = one frame per rising edge on FSIN, 0 =
# free-running.
#
# MEASURED ON THIS RIG, 2026-09-11, two Arducam B0332 (OV9281), MJPG 1280x800, FSIN DISCONNECTED,
# on BOTH devices:
#   - Set(19, 1, FLAG_MANUAL) -> the readback is (1, 0) and the camera stops free-running. With
#     no pulses on FSIN, cap.read() still answers ok=True, about once a second, after a ~1000 ms
#     driver timeout, and the frame is ALL BLACK. Set(19, 0, FLAG_MANUAL) -> ~121 fps again,
#     immediately, with the media type unchanged.
#   - IT PERSISTS IN THE CAMERA across cap.release(), a rebind of this filter and a reopen of the
#     capture. So a camera left at 1 reopens at 1, and anything that opens one of these cameras
#     for ordinary free-running use must write 0 BEFORE opening it -- otherwise capture.open_all
#     waits out its settle timeout on a camera that is only delivering black frames at 1 Hz.
#   - Backlight Compensation is NOT this switch, whatever the application note says about the
#     name. Its range here is 0..2, default 1, and 0, 1 and 2 all leave the camera free-running
#     at ~121 fps; only the brightness changed.
#
# Deliberately NOT in CAMERA_CONTROL_PROPS: GetRange(19) succeeds but answers (0, 0, 0, 0, 0), so
# read() would carry a line saying "0 [0..0 step 0, default 0]", which says nothing at all.
AE_PRIORITY = 19

# Property ids, per interface. The two sets overlap numerically -- Exposure is CameraControl 4
# and Sharpness is VideoProcAmp 4 -- so a property is only ever named by the pair.
CAMERA_CONTROL_PROPS = (("Pan", 0), ("Tilt", 1), ("Roll", 2), ("Zoom", 3),
                        ("Exposure", 4), ("Iris", 5), ("Focus", 6))
VIDEO_PROC_AMP_PROPS = (("Brightness", 0), ("Contrast", 1), ("Hue", 2), ("Saturation", 3),
                        ("Sharpness", 4), ("Gamma", 5), ("ColorEnable", 6),
                        ("WhiteBalance", 7), ("BacklightCompensation", 8), ("Gain", 9))


class IMoniker(IPersist):
    """Only BindToObject and BindToStorage are declared for real.

    A COM interface is a vtable and nothing else, so every slot BEFORE the ones being called has
    to be accounted for or a call lands on the wrong function pointer. IMoniker derives from
    IPersistStream, which comtypes does not ship, so that interface's four slots are declared
    here as no-argument stubs: they hold the right positions and are never called.
    """

    _iid_ = GUID("{0000000F-0000-0000-C000-000000000046}")
    _methods_ = [
        COMMETHOD([], HRESULT, "IsDirty"),                      # IPersistStream, stub
        COMMETHOD([], HRESULT, "Load"),                         # stub
        COMMETHOD([], HRESULT, "Save"),                         # stub
        COMMETHOD([], HRESULT, "GetSizeMax"),                   # stub
        COMMETHOD([], HRESULT, "BindToObject",
                  (["in"], POINTER(IUnknown), "pbc"),
                  (["in"], POINTER(IUnknown), "pmkToLeft"),
                  (["in"], POINTER(GUID), "riidResult"),
                  (["out"], POINTER(POINTER(IUnknown)), "ppvResult")),
        COMMETHOD([], HRESULT, "BindToStorage",
                  (["in"], POINTER(IUnknown), "pbc"),
                  (["in"], POINTER(IUnknown), "pmkToLeft"),
                  (["in"], POINTER(GUID), "riid"),
                  (["out"], POINTER(POINTER(IUnknown)), "ppvObj")),
    ]


class IEnumMoniker(IUnknown):
    _iid_ = GUID("{00000102-0000-0000-C000-000000000046}")
    _methods_ = [
        COMMETHOD([], HRESULT, "Next",
                  (["in"], c_ulong, "celt"),
                  (["out"], POINTER(POINTER(IMoniker)), "rgelt"),
                  (["out"], POINTER(c_ulong), "pceltFetched")),
        COMMETHOD([], HRESULT, "Skip", (["in"], c_ulong, "celt")),
        COMMETHOD([], HRESULT, "Reset"),
        COMMETHOD([], HRESULT, "Clone",
                  (["out"], POINTER(POINTER(IUnknown)), "ppenum")),
    ]


class ICreateDevEnum(IUnknown):
    _iid_ = GUID("{29840822-5B84-11D0-BD3B-00A0C911CE86}")
    _methods_ = [
        COMMETHOD([], HRESULT, "CreateClassEnumerator",
                  (["in"], POINTER(GUID), "clsidDeviceClass"),
                  (["out"], POINTER(POINTER(IEnumMoniker)), "ppEnumMoniker"),
                  (["in"], c_ulong, "dwFlags")),
    ]


def _control_methods():
    """The three methods both control interfaces have, in this vtable order.

    IAMCameraControl and IAMVideoProcAmp are separate interfaces with different IIDs and
    different property ids, but identical method signatures, so the declaration is written once.
    A function rather than a shared list or a common base class: comtypes builds the vtable from
    _methods_ per class and requires an _iid_ on whichever class carries them, so each interface
    gets its own fresh copy. GetRange doubles as the capability test: S_OK for a property the
    device really has, an error for one it does not, which comtypes raises as COMError.
    """
    return [
        COMMETHOD([], HRESULT, "GetRange",
                  (["in"], c_long, "Property"),
                  (["out"], POINTER(c_long), "pMin"),
                  (["out"], POINTER(c_long), "pMax"),
                  (["out"], POINTER(c_long), "pSteppingDelta"),
                  (["out"], POINTER(c_long), "pDefault"),
                  (["out"], POINTER(c_long), "pCapsFlags")),
        COMMETHOD([], HRESULT, "Set",
                  (["in"], c_long, "Property"),
                  (["in"], c_long, "lValue"),
                  (["in"], c_long, "Flags")),
        COMMETHOD([], HRESULT, "Get",
                  (["in"], c_long, "Property"),
                  (["out"], POINTER(c_long), "lValue"),
                  (["out"], POINTER(c_long), "pFlags")),
    ]


class IAMCameraControl(IUnknown):
    _iid_ = GUID("{C6E13370-30AC-11d0-A18C-00A0C9118956}")
    _methods_ = _control_methods()


class IAMVideoProcAmp(IUnknown):
    _iid_ = GUID("{C6E13360-30AC-11d0-A18C-00A0C9118956}")
    _methods_ = _control_methods()


# CameraControl first, so Exposure -- the one property this rig lives or dies by -- is the first
# line of any readout built from this.
_INTERFACES = (("CameraControl", IAMCameraControl, CAMERA_CONTROL_PROPS),
               ("VideoProcAmp", IAMVideoProcAmp, VIDEO_PROC_AMP_PROPS))


def _ensure_com():
    """CoInitialize this thread, harmlessly, however many times it is called.

    comtypes initialises the main thread on import, so this is normally a no-op that bumps a
    reference count. It is called anyway because "already initialised" is free and a missing
    initialisation is a CO_E_NOTINITIALIZED on the first real call.
    """
    try:
        comtypes.CoInitialize()
    except OSError:
        pass                # already initialised on this thread, in another apartment model


def _moniker_string(moniker, name):
    """One IPropertyBag property off a moniker, or None if it has no such property.

    BindToStorage rather than BindToObject: the property bag is the moniker's own metadata, out
    of the registry, so nothing is instantiated and the device itself is not touched.
    """
    try:
        bag = moniker.BindToStorage(None, None, IPropertyBag._iid_).QueryInterface(IPropertyBag)
        return str(bag.Read(name, VARIANT(), None))
    except (comtypes.COMError, ValueError, TypeError):
        return None


def _monikers():
    """Every video input device moniker, in enumeration order. [] when there are none.

    CreateClassEnumerator answers S_FALSE and a NULL enumerator for an empty category -- not an
    error, so comtypes does not raise, so the NULL has to be looked for.
    """
    _ensure_com()
    dev_enum = comtypes.CoCreateInstance(CLSID_SystemDeviceEnum, interface=ICreateDevEnum)
    enum = dev_enum.CreateClassEnumerator(CLSID_VideoInputDeviceCategory, 0)
    if not enum:
        return []

    out = []
    while True:
        moniker, fetched = enum.Next(1)
        if not fetched or not moniker:
            return out
        out.append(moniker)


def list_devices():
    """[{index, friendlyName, devicePath}] for every video input device, in OpenCV's order.

    index is 0-based, as cv2.VideoCapture(index, cv2.CAP_DSHOW) counts; the app's device N is
    index N-1. See the module docstring. Opens nothing -- only the registry metadata is read --
    so this is safe to call while the cameras are streaming.
    """
    return [{"index": k,
             "friendlyName": _moniker_string(m, "FriendlyName"),
             "devicePath": _moniker_string(m, "DevicePath")}
            for k, m in enumerate(_monikers())]


class Controls:
    """The control interfaces of ONE device, held open so read() can be called again and again.

    Constructing this instantiates the device's capture filter (BindToObject). That is not a
    capture graph -- no pins are connected and no frames flow -- but it is a reference on the
    driver, so close() exists and is worth calling.

    Measured on this rig: both the reads AND the one write coexist with an OpenCV capture on the
    same device, in either order. It is not a claim on the device in the sense that makes
    cameras.CameraBusyError.

    Raises RuntimeError when there is no such device, and comtypes.COMError when the bind itself
    fails. Neither gets a class of its own: the only useful response to either is to say the
    readout is unavailable and show the reason.
    """

    def __init__(self, index):
        _ensure_com()
        monikers = _monikers()
        if not 0 <= index < len(monikers):
            raise RuntimeError(f"no video input device at index {index}; "
                               f"{len(monikers)} device(s) enumerated")
        self.index = index
        self.friendly_name = _moniker_string(monikers[index], "FriendlyName")
        self._filter = monikers[index].BindToObject(None, None, IID_IBaseFilter)
        self._controls = []
        for label, interface, props in _INTERFACES:
            try:
                self._controls.append((label, self._filter.QueryInterface(interface), props))
            except comtypes.COMError:
                # A device can expose one interface and not the other -- a camera with no
                # mechanical controls has no IAMCameraControl. Not an error; just fewer lines.
                pass
        if not self._controls:
            self.close()
            raise RuntimeError(f"device at index {index} exposes neither IAMCameraControl nor "
                               f"IAMVideoProcAmp")

    def read(self):
        """Every property the device actually supports, with its value, range and mode.

        A property counts as supported when GetRange returns S_OK. One that does not is dropped
        rather than reported as zero: an OV9281 has no Pan, and a line reading "Pan 0" would be
        a lie with a number on it.

        `auto` is True, False, or None when the driver sets neither flag -- which happens, and
        is worth showing as unknown rather than guessed as Manual.
        """
        out = []
        for label, control, props in self._controls:
            for name, prop_id in props:
                try:
                    lo, hi, step, default, _caps = control.GetRange(prop_id)
                except comtypes.COMError:
                    continue                        # the device does not have this property
                try:
                    value, flags = control.Get(prop_id)
                except comtypes.COMError:
                    value, flags = None, 0
                auto = None
                if flags & FLAG_AUTO:
                    auto = True
                elif flags & FLAG_MANUAL:
                    auto = False
                out.append({"interface": label, "name": name, "value": value,
                            "min": lo, "max": hi, "step": step, "default": default,
                            "auto": auto, "flags": flags})
        return out

    def get(self, interface, prop_id):
        """ONE property off ONE interface, as (value, flags). For ids read() does not carry.

        `interface` is the label -- "CameraControl" or "VideoProcAmp". read() is the readout of
        everything the device supports; this is the way to ask for a single id by number, which
        is what AE_PRIORITY needs, since it is deliberately not in the property lists.

        RuntimeError when the device does not expose that interface at all. A comtypes.COMError
        from the call itself is left to the caller: only the caller knows whether a property the
        driver refused is a failure or an answer.
        """
        control = next((c for label, c, _props in self._controls if label == interface), None)
        if control is None:
            raise RuntimeError(f"device at index {self.index} does not expose IAM{interface}")
        return control.Get(prop_id)

    def set(self, interface, prop_id, value, flags=FLAG_MANUAL):
        """Write ONE property. THE ONLY WRITE PATH IN THIS MODULE; nothing else here writes.

        `interface` is the label, as in get(). Returns nothing: the driver's answer is what a
        following get() reads, and cap.set()-style "it returned True" is not evidence.

        MEASURED ON THIS RIG, 2026-09-11: the write takes effect on a camera that ANOTHER
        thread's OpenCV capture graph is streaming from, through this separately bound filter,
        with the capture.CameraReader still running. No reader stop, no re-apply of the source
        and no format check are needed -- the media type was unchanged afterwards (still ~121 fps
        at 1280x800 after toggling AE_PRIORITY on and off again).

        RuntimeError when the device does not expose that interface; comtypes.COMError from the
        call itself propagates, for the reason get() gives.
        """
        control = next((c for label, c, _props in self._controls if label == interface), None)
        if control is None:
            raise RuntimeError(f"device at index {self.index} does not expose IAM{interface}")
        control.Set(prop_id, value, flags)

    def close(self):
        """Drop the interfaces and the filter. Idempotent; a second call does nothing."""
        self._controls = []
        self._filter = None
