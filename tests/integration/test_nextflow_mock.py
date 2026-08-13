import json
import os
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
NEXTFLOW = shutil.which("nextflow") or (
    str(Path(".tools/nextflow").resolve())
    if Path(".tools/nextflow").is_file()
    else None
)
NEXTFLOW_ENV = {
    **os.environ,
    "JAVA_HOME": str(Path(".tools/nextflow-env").resolve()),
    "NXF_HOME": str(Path(".tools/nextflow-home").resolve()),
    "NXF_VER": "24.10.4",
    "NXF_OFFLINE": "true",
}


def test_workflow_declares_all_stage_modules() -> None:
    main = Path("workflow/main.nf")

    assert main.exists()
    source = main.read_text(encoding="utf-8")
    for filename, process_name in REQUIRED_MODULES.items():
        assert (Path("workflow/modules") / filename).exists()
        assert process_name in source


@pytest.mark.skipif(NEXTFLOW is None, reason="Nextflow not installed")
def test_mock_workflow_reaches_report_and_resumes(tmp_path: Path) -> None:
    command = [
        str(NEXTFLOW),
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
    first = subprocess.run(
        command, capture_output=True, text=True, timeout=120, env=NEXTFLOW_ENV
    )
    second = subprocess.run(
        [*command, "-resume"],
        capture_output=True,
        text=True,
        timeout=120,
        env=NEXTFLOW_ENV,
    )

    assert first.returncode == 0, first.stdout + first.stderr
    assert second.returncode == 0, second.stdout + second.stderr
    assert (tmp_path / "delivery/final_report/report.md").exists()
    manifest = json.loads(
        (tmp_path / "delivery/final_report/report_manifest.json").read_text()
    )
    assert manifest["query_count"] == 2
    assert manifest["technical_success_rate"] == 1.0
    assert len(manifest["queries"]) == 2
    assert (tmp_path / "delivery/job_status.sqlite").exists()
    assert "cached:" in second.stdout.lower()
