# eye-calibration-recording

Aligns a pair of eye cameras, then records them — **one process**, driven by Spike2.

```
python eye-calibration-recording.py --session-id m123_run01 --seconds 30
```

| phase | Spike2 |
|---|---|
| 1. open both cameras, negotiate MJPG at full frame | — |
| 2. live alignment view; operator presses `a` | — |
| 3. put both cameras into FSIN trigger mode | — |
| 4. create the ready flag | sees it, drops a `SampleKey("S")` marker, **starts the pulse train** |
| 5. record one frame per FSIN pulse, or free-run + anchor if none come | sampling and pulsing |
| 6. write `session.json` | sees it and continues |

## Why one process

Each camera's strobe pin is wired to a Spike2 digital input, so Spike2 timestamps every
exposure. But the cameras strobe **continuously from the moment their graph starts** — through
format negotiation, through the two ~6 Hz restart segments, and through the whole alignment
session. On the bench that start-up is about 7 s per camera.

Keeping everything in one process means the cameras are **never released** between aligning and
recording: no second start-up, no second burst of settling pulses, and no need to hand the
accepted device order and rotation between two processes — they stay in memory.

## The Spike2 contract

`ProgStatus()` reports only `1` = running and `0` = terminated; it cannot carry the child's exit
code. Three file facts carry everything instead:

| fact | meaning |
|---|---|
| the app **creates** the ready flag | alignment accepted; recording begins now — and **line 1 is the session directory** |
| `ProgStatus()` hits 0 while the flag is absent | cancelled, or start-up failed → halt |
| `<out>/<session-id>/session.json` exists | the recording completed |

**The app never deletes the flag — Spike2 owns its lifetime, and must `FileDelete()` it before
`ProgRun`.** Not after: this process takes a couple of seconds just to import cv2, so anything it
cleared at start-up would be cleared long after a 0.1 s Spike2 poll loop had already seen the
stale file and charged ahead. The app warns on stderr if it finds a flag already there.

**Cancel needs no signal of its own.** Close the alignment window or press `Esc` and the app
exits without ever creating the flag, so `ProgStatus()` drops to `0` while Spike2 is still
polling. Same for any start-up failure.

`--session-id` is therefore **required from Spike2** — it is how Spike2 knows where to look for
`session.json`.

### Reading the session directory out of the flag

Line 1 of the ready flag is the full path of the directory the frames go to, alone on the line,
so Spike2 can save its own `.smrx` alongside them:

```
C:\Temp\test\m123_run01
eye-calibration-recording, 2026-08-18T14:31:07Z: alignment accepted, cameras streaming, waiting for Spike2.
LINE 1 above is the session directory -- Spike2 reads it to save alongside the frames.
Delete this file to start the recording.
```

```
var fh% := FileOpen(readyPath$, 8, 0);
var sessionDir$;
Read(sessionDir$);              ' line 1: the session directory
FileClose();
...
FileSaveAs(sessionDir$ + "\\eye_calibration.smrx", 0);
```

Two things to know. **The directory does not exist yet** when the flag is written — it is created
a moment later, when storage begins — so save at the *end* of the script, by which time the
recording has run and `session.json` is there to confirm it. And if your `--out` path ever
contains **spaces**, check that your Spike2 version's `Read()` returns the whole line rather than
the first token; keeping session paths space-free avoids the question entirely.

Anything machine-readable added to this file later must go on its own line **below** the path,
never above it.

### Spike2 script

```
const readyPath$ := "C:/Temp/eyecal_ready.flag";
const PY$        := "\"C:/Users/Public/RaymondLabCode2/.venv/Scripts/python.exe\"";
const APP$       := " \"C:/Users/Public/RaymondLabCode2/eye_sensor_calibration/recording_app/eye-calibration-recording.py\"";
const OUT$       := "C:/Temp/test";

var id$  := "m123_run01";
' the camera preset is detected from the devices unless --camera says otherwise
var cmd$ := PY$ + APP$ + " --seconds 10.0" +
            " --out \"" + OUT$ + "\" --session-id " + id$ +
            " --ready-flag \"" + readyPath$ + "\"";

FileDelete(readyPath$);         ' MUST be here, before ProgRun -- see above
SampleStart();

PrintLog("Starting: %s\n", cmd$);
var hProg% := ProgRun(cmd$, 1, 50, 5, 98, 90);
if hProg% < 0 then PrintLog("ProgRun failed (%d)\n", hProg%); halt endif;

' flag appears = accepted and recording; process ends first = cancelled
var st% := 1, fh% := -1;
repeat
    Yield(0.1);
    st% := ProgStatus(hProg%);
    fh% := FileStatus(readyPath$);
until (fh% <> -1) or (st% <> 1);

if st% <> 1 then PrintLog("Alignment cancelled - halting.\n"); halt endif;
SampleKey("S");

repeat
    Yield(0.05);
    st% := ProgStatus(hProg%);
until st% <> 1;

Yield(2);
SampleKey("s");
SampleStop();

if FileStatus(OUT$ + "/" + id$ + "/session.json") = -1 then
    PrintLog("Recording did not complete - halting.\n"); halt
endif;
PrintLog("Recording complete.\n");
```

`FileStatus()` is the right existence test: `-1` when the file is absent, non-negative when it is
there, and it does not open a view or consume the file.

### Where to put SampleStart()

The script above starts sampling **before** `ProgRun`, which is the simple arrangement: Spike2 is
definitely running by the time anything is recorded, and there is no race. The cost is that
camera start-up and the whole alignment session land in the Spike2 file as tens of thousands of
strobe pulses before the ones you care about. The anchor makes them harmless — it identifies the
recorded frames directly — but the file is bigger and busier.

If you want a **clean file** instead, move `SampleStart()` to after the flag appears and have
Spike2 delete the flag as a "go":

```
if st% <> 1 then PrintLog("Alignment cancelled - halting.\n"); halt endif;
SampleStart();
SampleKey("S");
FileDelete(readyPath$);         ' the "go"
```

and launch the app with `--handshake-wait-s 5`. It then waits for the flag to disappear before
storing its first frame, so nothing but the recording is in the file. It records anyway if
nobody answers within that window — a Spike2 that never replies must not be able to hang a
session.

## The strobe anchor

Spike2 timestamps every exposure, but nothing in the pulse train says *which pulse is stored
frame 0* — the camera has been free-running since long before recording began.

So the recording marks itself. At the start and again at the end, the exposure is raised by four
stops (16× the integration time) for 0.25 s. Two things happen, and either is enough:

- **in the Spike2 file** — those strobe pulses go conspicuously **wide**, because the strobe
  tracks exposure time. +4 stops from −10 is about 16 ms, which fits inside a 20 ms frame period,
  so the frame *rate* is unchanged and the timing being anchored is not perturbed.
- **in the stored data** — those same frames come out saturated. They are **stored, not
  discarded**, and their indices are written to `session.json`.

The anchor is therefore a frame index on one side and a pulse index on the other, with nothing
counted or assumed in between. And because there are two of them, **the pulse count between the
anchors must equal the stored frame count between them** — a complete drop check that no
queue-depth heuristic can give you.

Mapping a session:

1. Read `strobeAnchors.brightFrames.cam1.start` from `session.json` — a contiguous run of frame
   indices.
2. Find the matching run of wide pulses on that camera's TTL channel.
3. Align them; do the same at the end; confirm the two counts agree.
4. `ts.csv` carries per-frame host arrival times as an independent check — the inter-frame jitter
   is a fingerprint that can be cross-correlated against the pulse train. Measured against the
   strobe over ten sessions, `t_arrive_s` sits **4.8 ms sd** from the true pulse time: enough for
   ordering and rate, nowhere near enough for the sub-millisecond questions the TTL answers.
   The gap is USB transfer, driver buffering and OS scheduling, not clock error.

Set `"anchor": false` to turn it off. The app **proves at start-up** that an exposure change does
not re-negotiate the media type (`cameras.probe_anchor_safe`), and refuses to run if it does.

### Pairing

Do **not** pair by frame index. The two sensors free-run off independent oscillators with an
arbitrary phase offset of up to one frame period, plus drift; index *k* on one camera is not
simultaneous with index *k* on the other, and the two counts need not even match. Pair on
**strobe times** from Spike2. That is the whole point of having them.

## FSIN trigger mode

The anchor exists because nothing in a free-running pulse train says which pulse is frame 0.
**Under the FSIN trigger that question does not arise: stored frame `k` IS pulse `k`**, on every
camera, because the sensor exposes only when Spike2 tells it to. The whole mechanism is
`eyecal/trigger.py`, which also carries the measurements it rests on.

The switch is **one control write**, and OpenCV cannot reach it: `IAMCameraControl` property
**19**, `AUTO_EXPOSURE_PRIORITY` — the UVC `CT_AE_PRIORITY_CONTROL`, Windows' "Low Light
Compensation" checkbox, Linux's `exposure_dynamic_framerate`. It goes through a filter of the
app's own (`eyecal/dshow.py`), beside whatever OpenCV is doing, and it works while the capture
is streaming. Measured on this rig, 2026-09-11, on both Arducam B0332 (OV9281): `1` = one frame per
rising edge on FSIN, `0` = free-running at ~121 fps again, immediately.

### What happens after the operator presses `a`

1. **write** AE priority 1 on every device and read it straight back. Anything less than every
   device agreeing falls back — one camera pulsing and one free-running would produce two files
   whose frame indices mean different things.
2. **settle** a quarter second. Free-running frames keep landing for a moment after the write
   (measured: 8–16 per camera, 9 ms apart), and drained too early they would be stored as pulse
   0, 1, 2 and shift every index after them.
3. **drain** the queues. Everything in them predates the switch.
4. **create the ready flag.** This is what starts Spike2's pulse train, which is why the write
   comes first: a camera not yet in trigger mode would miss the first pulses.
5. **wait up to `trigger_wait_s`** (3 s) for a real frame from *every* camera. If none arrives,
   every camera goes back to free-run and the run records exactly as it always did, with the
   strobe anchor.
6. **record one frame per pulse.** It stops when the *train* does, not on the clock —
   `seconds + trigger_end_margin_s` is only a ceiling, so a Spike2 left pulsing cannot record
   until the disk fills.

**The camera remembers the setting** across a release, a reopen and the end of the process. So
free-run is written *before* anything is opened and again on **every** exit path, cancel and
crash included. Without that, the next ordinary recording gets one black frame a second and sits
out its whole settle timeout looking dead.

### Reading the `trigger` block in `session.json`

`mode` is the first thing to read, and it decides how the whole session is read:

| `mode` | what it means |
|---|---|
| `trigger` | stored frame `k` is pulse `k`. `strobeAnchors.enabled` is `false`: nothing to mark |
| `free-run` | an ordinary recording. Map it through the anchors, and read `reason` for why |

A fallback is **exit code 3** and a warning in `session.json`, because the operator must not
discover it by counting pulses in Spike2 afterwards. `reason` names it: the feature switched off,
`trigger write failed on device N`, or `no pulses within N s`. The rest of the block is
`aePriority` (the readback per device, and `aePriorityRestored` after a fallback), `pulseHz`,
`waitS`, `endMarginS`, `timeouts`, `framesDuringWait`, `lastRealFrameS` and `endReason`.

**The per-frame timestamp is Spike2's own TTL2 loopback.** The train Spike2 drives out to FSIN
is recorded back on TTL2, so every pulse carries a Spike2 timestamp — and pulse `k` is frame
`k`, with nothing counted or assumed in between. `ts.csv` still carries host arrival times;
they are for ordering and rate, not for timing.

**Timeout frames.** In trigger mode with no pulses the read does *not* fail — the driver answers
`ok=True` about once a second, after its ~1000 ms timeout, with an **all-zero** frame. Those are
counted as `timeouts` and never stored: one in the ledger would be a frame the sensor never took,
and it would shift every pulse after it. A really triggered frame cannot be all zero, because the
sensor's black-level pedestal sits near 28 counts at any exposure. One timeout at the end is how
a normal run *ends*; any others are one-second gaps in the train, and each is warned about.

Settings: `trigger` (on), `trigger_wait_s` (3 s), `trigger_end_margin_s` (5 s) and `pulse_hz`
(`null`). **Spike2 should pass `--pulse-hz`** — the rate it is actually driving. It is recorded,
and a delivered rate more than 10 percent away from it is a warning. `--no-trigger` records
free-running with the anchor, as before.

### Bench tools

```
python tests\fsin_trigger_check.py --seconds 5              records, then reports rates
python tests\fsin_trigger_check.py --seconds 5 --mode free  the control condition
python tests\fsin_trigger_view.py                           live: the panes ARE the pulses
```

`fsin_trigger_check.py` is the measurement: start the train, run it, and compare the delivered
rate against the pulse rate and against `--mode free`. `fsin_trigger_view.py` is the live view —
each pane goes black within 100 ms of the train stopping, so somebody can start and stop it and
watch. Both drive `eyecal/trigger.py`, so neither can disagree with the app, and both put every
camera back to free-run on the way out.

## Settings

`config.json` sits beside the script and is loaded automatically. The command line overrides it.
An unknown key is an error, not a shrug.

| key | default | notes |
|---|---|---|
| `camera` | `auto` | identifies the family (`elp` / `ov2311` / `ov9281`) from the geometries the device offers; set a preset name to force one |
| `devices` | `[1, 2]` | 1-based, as MATLAB's winvideo numbers them |
| `format`, `exposure` | `null` | override the preset |
| `seconds` | `30.0` | recording duration |
| `out` | `C:/Temp/test` | parent directory for sessions |
| `rotate180` | `false` | starting state; `r` changes it, and what is accepted rotates the **stored** frames |
| `align_preview_hz` | `25.0` | alignment redraw ceiling |
| `rec_preview_hz` | `10.0` | recording preview; `0` disables the window |
| `downsample` | `2` | alignment display pixel stride — never affects stored frames |
| `rec_downsample` | `4` | recording display pixel stride — **coupled to `rec_window_scale`**, see below |
| `rec_window_scale` | `0.5` | recording window as a fraction of the screen; alignment uses `0.92` |
| `handshake` | `true` | `false` skips the flag entirely (standalone testing) |
| `ready_flag` | `C:/Temp/eyecal_ready.flag` | |
| `handshake_wait_s` | `0.0` | `>0` waits that long for Spike2 to delete the flag before storing, then records regardless |
| `anchor` | `true` | the strobe bracket |
| `anchor_exposure` | `null` | `null` = normal + 4 stops |
| `anchor_hold_s` | `0.25` | |
| `trigger` | `true` | attempt FSIN trigger mode after accept; falls back if no pulses come |
| `trigger_wait_s` | `3.0` | how long to wait for the FIRST pulse on every camera |
| `trigger_end_margin_s` | `5.0` | how long past `seconds` a triggered run may go on |
| `pulse_hz` | `null` | the rate Spike2 is driving; pass `--pulse-hz`, and it is checked |

### Camera detection

`camera: auto` asks the cameras what they are. Each preset's native full frame is offered by one
family and no other, so the app asks every device for every preset's geometry and reads back what
the driver settled on. DirectShow does not synthesise a mode a camera does not have — it snaps to
the nearest one it does have, so an ELP asked for `1280x800` answers `1280x720` — and exactly one
preset reads back what it asked for. That is the identification. It costs about 5 s per camera,
all of it before Spike2 starts sampling.

Windows friendly names and USB VID/PID are deliberately **not** used: the two Arducam models
share a VID/PID, neither carries a serial number, and Windows caches the name against the port,
so moving a cable renames a camera. The limit is that this identifies the camera **family**, not
the unit — two of the same model cannot be told apart by anything software can see here. The
readbacks it decided on are kept in `session.json` under `requested.cameraDetection`.

### The two recording-preview settings are coupled

The recording preview is half-size on purpose — it is a progress indicator, not something to aim
a camera by, and it leaves room for the Spike2 display. But **shrinking the window without
shrinking the canvas makes it slower, not faster.**

`_letterbox` picks its interpolation by direction: `INTER_NEAREST` when scaling the canvas up,
`INTER_AREA` when scaling it down. Measured on this rig, `INTER_AREA` costs **5.75 ms per output
megapixel against `INTER_NEAREST`'s 0.78** — about 7× per pixel, which easily outweighs having
fewer pixels to write. So if the canvas ends up larger than the window, every refresh takes the
expensive path:

| `rec_downsample` | `rec_window_scale` | canvas → window | path | per refresh |
|---|---|---|---|---|
| 2 | 0.92 | 1600×600 → 1766×662 | upscale | 3.53 ms |
| 2 | 0.5 | 1600×600 → 883×331 | **downscale** | 4.78 ms |
| 4 | 0.92 | 800×300 → 1766×662 | upscale | 1.75 ms |
| **4** | **0.5** | 800×300 → 883×331 | upscale | **1.15 ms** |

Rule: keep the canvas smaller than the window. If you halve the window, at least halve the
canvas too. The window prints a note on stderr if the pairing is wrong, and
`tests/test_offline.py` asserts the shipped values are consistent for every preset.

None of this is a meaningful performance lever — the whole preview is ~3.5% of one core at
10 Hz. `rec_preview_hz` and `--no-preview` are the real ones. It is an ergonomics setting.

Alignment keys: `a` accept · `Esc` cancel · `c` colour · `r` rotate 180° · `s` swap sides ·
`h` hide overlay. Closing the window cancels. During recording, `q`/`Esc`/closing the preview
stops early — which still finalises the session and still counts as success.

## Output

```
<out>/<session-id>/
    frames/c1.bin, c1.json     flat frames + sidecar (shape, count, dtype)
    frames/c2.bin, c2.json
    first_frame_cam1.tif       lossless, double-clickable
    ts.csv                     cam, frame_idx, t_arrive_s, t_abs_posix, queue_depth,
                               queue_wait_ms, write_ms
    session.json               written LAST — its existence is the success signal
```

`c1` is whichever device the operator accepted as leftmost, recorded in
`session.json.acceptedDeviceOrder`. Frames are uncompressed and flat because PNG encoding
measures ~1164 ms per frame against ~1 ms for a raw write. Read them back with
`eyecal.session.read_frames(dir, cam)`, or four lines of numpy printed at the end of every run.

`ts.csv` columns — timestamps in seconds, durations in milliseconds, and the suffix says which:

| column | meaning |
|---|---|
| `t_arrive_s` | frame reached Python, relative to the start of recording |
| `t_abs_posix` | the same instant on the wall clock |
| `queue_depth` | frames waiting when this one came off the queue |
| `queue_wait_ms` | arrival → off the queue. The consumer falling behind |
| `write_ms` | off the queue → bytes down. The rotation, if any, and the file write |

`session.json.clock` pins `t_arrive_s` to the wall clock: `t0Perf`/`t0Posix` and `t1Perf`/`t1Posix`
bracket the run, and the per-frame fit (`absVsArrivePpm`, `absVsArriveResidualUs`,
`absVsArriveMaxStepUs`) says whether the OS moved the clock mid-run. A residual of a few
microseconds is a smooth slew and expected; a large one means a step, and any wall-clock alignment
for that session should be distrusted.

Sink cost is `queue_wait_ms + write_ms`. `warnings_for` gates on the **p95** and on the count of
frames over half a frame period — not on the maximum, because every run contains one 30–60 ms
outlier, so a peak-based gate fires on healthy and struggling sessions alike.

## Bench checks

Things only the real rig can answer. Worth doing once, at commissioning. `tools/strobe_timing.py`
does the measurement for items 1-2 from a recorded session and its `.smrx` -- point it at the
Spike2 file and read the report.

1. **Does the exposure change widen the strobe, or gap it?** Run a short session and look at the
   first and last 0.25 s on the TTL channels. Wide pulses = as designed. A gap means the driver
   restarts the graph — still a usable marker, but a few frames are lost inside the anchor
   window, so note it.
2. **Does a failed `cap.read()` still leave a pulse behind?** This is the only thing that can
   silently shift the mapping by one frame. Compare the pulse count between anchors against the
   stored frame count between them over a long run.
3. **Is there a trigger *input* on the camera board? Yes — and it works.** Measured on this rig,
   2026-09-11, on both Arducam B0332 (OV9281): FSIN driven from a Power1401 digital output gives
   one frame per rising edge, and the switch is `IAMCameraControl` property **19**
   (`AUTO_EXPOSURE_PRIORITY`), *not* Backlight Compensation, whatever the application note's name
   for it suggests — 0, 1 and 2 on that control all leave the camera free-running at ~121 fps and
   only change the brightness. Exposures are now simultaneous by construction, pulse count equals
   frame count, and the anchor is a fallback rather than the mechanism. See
   [FSIN trigger mode](#fsin-trigger-mode).

**If the app refuses with "exposure is not under manual control"**, read the means it prints. It
drives the exposure from −13 to −6 on every `AUTO_EXPOSURE` setting it knows, twice, with one
close and reopen in between, and lists what each attempt measured. Means pinned at the same middle
value on every line are the camera running its own auto-exposure; means near the black level
(below about 35 counts) are a dark scene instead — lens cap, lights off — which this check cannot
tell apart from auto. A real pin usually comes from the **Windows Camera app**: it, and anything
else that goes through Windows' Frame Server, can leave a camera in an auto mode these writes
cannot undo, and the state has appeared a minute or two *after* the app was closed. Do not use the
Camera app on the rig cameras. If it has been used, wait two minutes and retry, or disable and
re-enable the camera in Device Manager under Cameras.

### Use `tools/camera_check.py`, not the Windows Camera app

```
python tools\camera_check.py
```

That is what it is for. It shows both cameras live, side by side, and **writes nothing to any
camera control** — it opens them exactly as the app does, in the same load-bearing property
order, but hands the open a source dict holding only `FrameRate`, so each camera keeps whatever
state it was already in. Per camera it shows the device number and Windows' friendly name, the
negotiated media type and both frame rates, and then **every control the driver exposes** —
value, range, step, default, and whether it is on **AUTO**, in red. That last column is the one
thing the app itself cannot tell you: it comes from DirectShow's `IAMCameraControl` and
`IAMVideoProcAmp` through `eyecal/dshow.py`, and it is instant. Measured on this rig, both
OV9281s report `WhiteBalance` on AUTO and everything else on MANUAL out of the box.

The AUTO flag is a **hint, not the proof** — it says what the driver was last told. `m` runs the
app's real measurement (`cameras.secure_manual_exposure`) on demand and prints the evidence; that
verdict is the one that decides whether a session runs. `p` saves two PNGs: the canvas exactly as
shown, and the full-resolution frames tiled with nothing drawn on them. A camera that will not
open keeps its pane and shows the failure in red while the other one goes on running; `o` retries
it. Keys: `Esc` quit · `p` snapshot · `m` prove manual exposure · `o` reopen failed · `c` colour ·
`r` rotate · `s` swap sides · `h` hide overlay. `--snapshot PATH` and `--quit-after S` run it
unattended.

## Layout

```
eye-calibration-recording.py   entry point: CLI, phase sequencing, exit codes
eyecal/cameras.py              open + identify + format negotiation + verification, presets
eyecal/capture.py              CameraReader thread, open_all / close_all
eyecal/display.py              overlay, crosshair, tiling, letterbox, Window
eyecal/align.py                alignment phase -> accepted order and rotation
eyecal/record.py               recording phase: the strobe anchor, and one frame per FSIN pulse
eyecal/session.py              frame files, sidecars, ts.csv, session.json, report
eyecal/spike2.py               the ready flag
eyecal/dshow.py                DirectShow control readout: value, range and AUTO/MANUAL per property
eyecal/trigger.py              the FSIN trigger: the control write, the restore, the timeout frame
tests/test_offline.py          the pieces, without a camera
tests/test_endtoend.py         main() end to end, fake cameras, a fake Spike2 and a fake FSIN train
tests/fsin_trigger_check.py    bench check: record N seconds under the trigger and report rates
tests/fsin_trigger_view.py     bench check: watch the pulse train drive both cameras, live
tools/strobe_timing.py         measure a recorded session against its Spike2 file (bench check)
tools/session_video.py         render a session to calibration_video.mp4 with its timestamps on it
tools/camera_check.py          live camera state + snapshots; use INSTEAD of the Windows Camera app
```

`tests/` holds tests: they exercise the code and assert. `tools/` holds things you run against a
session after the fact; `strobe_timing.py` measures rather than asserts, and is kept out of
`tests/` so no test runner ever collects it.

Requires `opencv-python`, `numpy` and `comtypes`. The measured DirectShow behaviour this depends on was
verified on OpenCV 5.0.0 / Python 3.12; re-run both test files and one short real recording after
changing the interpreter.

## Exit codes

Useful from a terminal, but **not** what Spike2 reads:

| code | meaning |
|---|---|
| 0 | recorded and finalised |
| 1 | failed — camera would not open or could not be identified, wrong format, not streaming |
| 2 | bad command line |
| 3 | recorded, but a guard tripped — read `warnings` in `session.json` |
| 4 | cancelled by the operator at the alignment stage |
