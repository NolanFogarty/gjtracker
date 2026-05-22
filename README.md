# GJTRACKER

A browser-based Moon/Sun position tracker for EME (Earth-Moon-Earth) ham radio
operation. Python port of W7GJ's original Visual Basic 6 program.

## Setup

### 1. Install `uv`

`uv` is a fast Python package manager. It will handle downloading Python and
installing dependencies for you — you don't need Python pre-installed.

**macOS / Linux:**
```sh
curl -LsSf https://astral.sh/uv/install.sh | sh
```

**Windows (PowerShell):**
```powershell
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

After installation, open a new terminal so the `uv` command is on your `PATH`.
Verify with:

```sh
uv --version
```

Other install methods (Homebrew, winget, pip, etc.) are listed at
<https://docs.astral.sh/uv/getting-started/installation/>.

### 2. Get the script

Make sure `gjtracker.py` is on your machine. Any folder works.

### 3. Run it

From the folder containing `gjtracker.py`:

```sh
uv run gjtracker.py
```

The first run takes a few seconds while `uv` downloads Python (if needed) and
installs Flask into a cached environment. Subsequent runs start instantly.

Your browser should open automatically to <http://127.0.0.1:5000/>. If it
doesn't, open that URL manually.

## Using the app

1. Fill in **Home Station** latitude/longitude — or type a Maidenhead locator
   (e.g. `CN87xq`) and click **Locator → Lat/Lon**.
2. Optionally enable **DX Station** and do the same.
3. Pick a date range and a time increment.
4. Click **Calculate Tracking Schedule**.

The results table shows only rows where the object is above the configured
horizon at every enabled station. For Moon tracking with two stations, the
table includes the cross-polarization angle in the common window.

## Stopping the app

Press **Ctrl+C** in the terminal where it's running.

## Notes

- Everything runs locally — no internet connection is required after the first
  install of dependencies.
- The server only binds to `127.0.0.1`, so it's only reachable from your own
  machine.
- The original `CALL3.TXT` callsign database and `GJSKYTEMP.DAT` sky-temperature
  data files are **not** required for this port; the callsign lookup and
  degradation-vs-sky-noise features are omitted.
