import shutil
import subprocess

import pytest


@pytest.mark.skipif(shutil.which("boltz") is None, reason="Boltz-2 binary not installed")
def test_boltz_binary_exposes_predict_help() -> None:
    result = subprocess.run(
        ["boltz", "predict", "--help"],
        capture_output=True,
        text=True,
        check=False,
        timeout=30,
    )

    assert result.returncode == 0
