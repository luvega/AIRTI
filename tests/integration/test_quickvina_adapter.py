import shutil
import subprocess

import pytest


@pytest.mark.skipif(shutil.which("qvina2") is None, reason="qvina2 binary not installed")
def test_quickvina_binary_exposes_help() -> None:
    result = subprocess.run(
        ["qvina2", "--help"], capture_output=True, text=True, check=False, timeout=30
    )

    assert result.returncode == 0
