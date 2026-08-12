import shutil
import subprocess
from pathlib import Path

import pytest

from airti_tf.simulation.gromacs import parse_completed_ns


def test_fixture_log_reports_100ns() -> None:
    assert parse_completed_ns(Path("tests/fixtures/gromacs/md.log")) == pytest.approx(100.0)


@pytest.mark.skipif(shutil.which("gmx") is None, reason="GROMACS binary not installed")
def test_gromacs_binary_exposes_version() -> None:
    result = subprocess.run(
        ["gmx", "--version"], capture_output=True, text=True, check=False, timeout=30
    )

    assert result.returncode == 0
