"""Launch the Streamlit GUI as a subprocess.

`hplc gui [experiment_dir]` calls :func:`launch_gui`, which shells out to
``streamlit run app.py``. Streamlit runs its own web server (bound to
localhost) and opens a browser tab. Anything after ``--`` on the streamlit
command line is passed through to app.py as ``sys.argv``.
"""

from __future__ import annotations

import subprocess
import sys
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

    cmd = [
        sys.executable, "-m", "streamlit", "run", str(APP_PATH),
        "--server.address", "localhost",
        "--server.port", str(port),
        "--browser.gatherUsageStats", "false",
    ]
    if experiment_dir is not None:
        cmd += ["--", str(Path(experiment_dir).expanduser().resolve())]

    print(f"Starting the HPLC GUI at http://localhost:{port}  (Ctrl-C to stop)")
    try:
        return subprocess.call(cmd)
    except KeyboardInterrupt:
        return 0
