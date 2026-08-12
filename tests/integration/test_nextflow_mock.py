import shutil
import subprocess
from pathlib import Path

import pytest


REQUIRED_MODULES = {
    "target_library.nf": "TARGET_LIBRARY",
    "ligand_prep.nf": "LIGAND_PREP",
    "screen.nf": "SCREEN",
    "refine.nf": "REFINE",
    "md.nf": "MD",
    "report.nf": "REPORT",
}


def test_workflow_declares_all_stage_modules() -> None:
    main = Path("workflow/main.nf")

    assert main.exists()
    source = main.read_text(encoding="utf-8")
    for filename, process_name in REQUIRED_MODULES.items():
        assert (Path("workflow/modules") / filename).exists()
        assert process_name in source


@pytest.mark.skipif(shutil.which("nextflow") is None, reason="Nextflow not installed")
def test_mock_workflow_reaches_report_and_resumes(tmp_path: Path) -> None:
    command = [
        "nextflow",
        "run",
        "workflow/main.nf",
        "-profile",
        "test",
        "--queries",
        str(Path("tests/fixtures/ligands.smi").resolve()),
        "--outdir",
        str(tmp_path / "delivery"),
        "-work-dir",
        str(tmp_path / "work"),
    ]
    first = subprocess.run(command, capture_output=True, text=True, timeout=120)
    second = subprocess.run(
        [*command, "-resume"], capture_output=True, text=True, timeout=120
    )

    assert first.returncode == 0, first.stdout + first.stderr
    assert second.returncode == 0, second.stdout + second.stderr
    assert (tmp_path / "delivery/final_report/report.md").exists()
    assert (tmp_path / "delivery/job_status.sqlite").exists()
    assert "cached:" in second.stdout.lower()

