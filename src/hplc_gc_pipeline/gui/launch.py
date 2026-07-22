"""Launch the Streamlit GUI as a subprocess.

`hplc gui [experiment_dir]` calls :func:`launch_gui`, which shells out to
``streamlit run app.py``. Streamlit runs its own web server (bound to
localhost) and opens a browser tab. Anything after ``--`` on the streamlit
command line is passed through to app.py as ``sys.argv``.

The server shuts down automatically when the last browser tab is closed.
"""

from __future__ import annotations

import os
import signal
import socket
import subprocess
import sys
import threading
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Optional

APP_PATH = Path(__file__).with_name("app.py")


def _skip_streamlit_onboarding() -> None:
    """Write an empty Streamlit credentials file if none exists.

    On first launch Streamlit otherwise prompts for an email on stdin, which
    hangs a non-interactive start. An empty email disables that prompt.
    """
    cred = Path.home() / ".streamlit" / "credentials.toml"
    if not cred.exists():
        cred.parent.mkdir(parents=True, exist_ok=True)
        cred.write_text('[general]\nemail = ""\n', encoding="utf-8")


def _is_port_in_use(port: int) -> bool:
    """Return True if something is already listening on *port*."""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        return s.connect_ex(("127.0.0.1", port)) == 0


def _free_port(port: int) -> None:
    """Kill any process currently listening on *port* (best-effort)."""
    try:
        result = subprocess.run(
            ["lsof", "-ti", f":{port}"],
            capture_output=True, text=True, timeout=5,
        )
        pids = [int(p) for p in result.stdout.strip().splitlines() if p.strip().isdigit()]
        for pid in pids:
            try:
                os.kill(pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
        if pids:
            time.sleep(1.5)
    except Exception:
        pass


def _active_sessions(port: int) -> int:
    """Return the current active-session count from Streamlit's metrics endpoint."""
    url = f"http://localhost:{port}/_stcore/metrics"
    try:
        resp = urllib.request.urlopen(url, timeout=2)
        for line in resp.read().decode().splitlines():
            if line.startswith("active_sessions "):
                return int(float(line.split()[-1]))
    except Exception:
        pass
    return -1  # server not yet up or unreachable


def _watch_and_terminate(proc: subprocess.Popen, port: int) -> None:
    """Daemon thread: terminates *proc* once all browser tabs are closed."""
    # Wait for the server to be up (up to 30 s)
    health = f"http://localhost:{port}/_stcore/health"
    for _ in range(60):
        try:
            urllib.request.urlopen(health, timeout=1)
            break
        except Exception:
            time.sleep(0.5)
    else:
        return  # server never came up

    ever_connected = False
    consecutive_zero = 0
    while proc.poll() is None:
        time.sleep(2)
        n = _active_sessions(port)
        if n > 0:
            ever_connected = True
            consecutive_zero = 0
        elif ever_connected and n == 0:
            consecutive_zero += 1
            if consecutive_zero >= 2:
                print("\nAll browser tabs closed — shutting down.", flush=True)
                proc.terminate()
                return
        else:
            # n == -1 (metrics unreachable) or not yet ever connected
            consecutive_zero = 0


def launch_gui(experiment_dir: Optional[Path] = None, port: int = 8501) -> int:
    """Start the Streamlit app. Returns the subprocess exit code."""
    try:
        import streamlit  # noqa: F401
    except ModuleNotFoundError:
        print(
            "The GUI needs the optional 'gui' dependencies. Install them with:\n"
            '    pip install -e ".[gui]"',
            file=sys.stderr,
        )
        return 1

    _skip_streamlit_onboarding()

    if _is_port_in_use(port):
        print(f"Port {port} is already in use — stopping previous instance…", flush=True)
        _free_port(port)

    cmd = [
        sys.executable, "-m", "streamlit", "run", str(APP_PATH),
        "--server.address", "localhost",
        "--server.port", str(port),
        "--browser.gatherUsageStats", "false",
    ]
    if experiment_dir is not None:
        cmd += ["--", str(Path(experiment_dir).expanduser().resolve())]

    print(f"Starting the HPLC GUI at http://localhost:{port}  (Ctrl-C to stop)")
    proc = subprocess.Popen(cmd)
    t = threading.Thread(target=_watch_and_terminate, args=(proc, port), daemon=True)
    t.start()
    try:
        proc.wait()
        return proc.returncode if proc.returncode is not None else 0
    except KeyboardInterrupt:
        proc.terminate()
        return 0
