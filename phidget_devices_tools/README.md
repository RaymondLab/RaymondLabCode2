# Phidget Live-Display GUIs

Two small GUIs for live-graphing Phidget sensors, built because the Phidget Control Panel cannot graph on macOS. Each GUI exists in two versions: a **single-file HTML app** (the primary version — the official phidget22 JavaScript library v3.25.1 is inlined, so the file is fully self-contained and works offline) and a **Python/matplotlib script** with equivalent features.

> These are prototypes / starting points — review, adapt, and verify before relying on them.

## The two GUIs

### 1. Gyro Live — `phidget_gyro_live.html` / `phidget_gyro_live.py`

**Hardware:** PhidgetSpatial Precision 3/3/3 (**MOT0110_0**), connected over USB.

Dual-trace rolling display for one gyroscope axis (z by default, selectable):

- **Top:** angular velocity (deg/s) — the primary, drift-free signal.
- **Bottom:** integrated angle (deg) — sanity check only; it drifts slowly, and a **Zero angle** button resets it.

On startup/connect the app averages ~2 s of resting data to auto-zero the gyro bias — **keep the sensor still for the first couple of seconds**. Extras: adjustable window length and data interval, average per-cycle **Min/Max** markers, a live **FFT** (wide + zoomed views with peak annotation), and CSV export of the data.

### 2. Light Live — `phidget_light_live.html` / `phidget_light_live.py`

**Hardware:** Light Phidget (**LUX1000_0**) on a VINT hub (e.g., a 1-port VINT hub).

Live illuminance display: lux vs. time in the top two-thirds, live stats (current/min/max/mean/std, sample rate, elapsed) plus a histogram of the displayed window in the bottom third. Deliberately **dark-mode** (near-black surfaces, dim amber trace) so the display adds as little light as possible to the room being measured and preserves dark adaptation. Extras: **log-scale** toggle for the lux axis, adjustable window length, device data interval (~100 ms–60 s) and illuminance change trigger (0 = report every interval), decimal precision, and CSV export.

## HTML versions — prerequisites & setup

The HTML files need no web server and no internet — just open the file in a browser. But the browser still needs a path to the hardware, and both paths require the **native Phidget drivers** to be installed.

> **Note on the embedded library:** each HTML file has the official **phidget22 JavaScript library v3.25.1** (release 1.25.20260408, ISC license) copied directly into it. That's what makes the files self-contained/offline, but it also means they are frozen at that version — they will **not** auto-update when Phidgets releases a newer library. If a future driver or firmware update breaks compatibility, replace the inlined library block in each file with the current version from [phidgets.com](https://www.phidgets.com/docs/Phidget22_API) (or the `phidget22` npm package).

### 1. Install the Phidget drivers

Download and install the Phidget22 drivers/Control Panel for your OS from [phidgets.com/docs/OS](https://www.phidgets.com/docs/Operating_System_Support). If the Phidget Control Panel already sees your device, the drivers are installed.

### 2. Open the HTML file and connect

Double-click `phidget_gyro_live.html` or `phidget_light_live.html` (or use File → Open in the browser), then choose one of the two connect buttons:

**Option A — "Connect USB" (direct WebUSB)**
- **Chrome or Edge only** (WebUSB is not available in Safari/Firefox).
- The device must **not be held by anything else**: quit the Phidget Control Panel, stop any Python scripts, and stop the Phidget Network Server first.
- Requires a PHIDUSB-stack device/hub.
- Click **Connect USB** and pick the device in the browser's chooser dialog.

**Option B — "Connect server" (works in any browser)**
1. On the machine the device is plugged into, open the Phidget Control Panel → **Network Server** tab → enable the **Phidget Web Server** (default port **8989**).
2. In the HTML page, enter the host (`localhost` if it's the same machine, otherwise that machine's hostname/IP) and port, then click **Connect server**.

### 3. Use it

- **Gyro:** hold the sensor still for the first ~2 s after connecting (bias auto-zero); changing the axis re-runs the auto-zero. Then use the toolbar for Zero angle / Min-Max / FFT / interval / window / Download CSV.
- **Light:** you can connect before plugging in the sensor — streaming starts the moment the LUX1000 attaches. Note that with a nonzero trigger (lx), a flat/quiet trace is expected behavior (readings are only reported when the value changes by at least the trigger amount).
- Kiosk/testing shortcut — auto-connect via URL parameters: `?autoconnect=net&host=localhost&port=8989`

## Python versions — prerequisites & setup

### 1. Install the Phidget drivers

Same as above — if the Phidget Control Panel sees your device, the native driver is already installed.

### 2. Install the Python packages

Requires a Python with Tk support (standard python.org and Anaconda installs include it; the scripts use the TkAgg matplotlib backend):

```bash
pip install phidget22 matplotlib
```

(numpy is pulled in by matplotlib.)

### 3. Run

Make sure nothing else is holding the device (quit the Phidget Control Panel and any browser session connected to it), then:

```bash
python phidget_gyro_live.py    # hold the sensor still for ~2 s at startup
python phidget_light_live.py
```

Each opens a maximized matplotlib window with a trimmed toolbar; the "configure subplots" button opens a custom dialog with layout sliders plus the live-view/device sliders (window length; and for the light script: trigger, interval, decimals). Defaults — axis, data interval, window length, CSV logging (`LOG_CSV`, off by default), colors, etc. — are constants in the **Config** section at the top of each script.
