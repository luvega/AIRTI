import subprocess
import sys
from pathlib import Path


def test_tracked_background_panel_passes_integrity_check() -> None:
    panel = Path("data/reference/background_probes_v1.smi")

    result = subprocess.run(
        [sys.executable, "scripts/build_background_panel.py", "--check-only", str(panel)],
        capture_output=True,
        text=True,
        check=False,
        timeout=30,
    )

    assert result.returncode == 0, result.stderr
    assert "100 probes" in result.stdout
