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
    "evaluate_case.nf": "EVALUATE_CASE",
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


def test_production_modules_pass_manifests_with_their_asset_directories() -> None:
    ligand = Path("workflow/modules/ligand_prep.nf").read_text(encoding="utf-8")
    screen = Path("workflow/modules/screen.nf").read_text(encoding="utf-8")
    refine = Path("workflow/modules/refine.nf").read_text(encoding="utf-8")
    md = Path("workflow/modules/md.nf").read_text(encoding="utf-8")
    report = Path("workflow/modules/report.nf").read_text(encoding="utf-8")

    assert "prepared_ligands.jsonl" in ligand
    assert "--asset-dir prepared_ligands" in ligand
    assert "tuple path('prepared_ligands.jsonl'), path('prepared_ligands')" in ligand
    assert "--asset-dir screen_assets" in screen
    assert "tuple path('screened_candidates.jsonl'), path('screen_assets')" in screen
    assert "--asset-dir boltz_assets" in refine
    assert "tuple path('boltz_candidates.jsonl'), path('boltz_assets')" in refine
    assert "--asset-dir md_assets" in md
    assert "tuple path('md_candidates.jsonl'), path('md_assets')" in md
    assert "--candidates ${md_manifest}" in report
    assert "--resume" not in md
    assert "--protocol ${params.md_protocol}" in md


def test_gpu_stages_are_serialized_in_the_single_unified_image() -> None:
    config = Path("workflow/nextflow.config").read_text(encoding="utf-8")

    assert config.count("maxForks = 1") >= 2
    assert "md_protocol = 'production'" in config
    assert "production_image = 'airti-tf:0.2.0-gpu'" in config


def test_production_reference_defaults_use_the_documented_data_root() -> None:
    config = Path("workflow/nextflow.config").read_text(encoding="utf-8")

    assert (
        "target_manifest = "
        "'/data/airti-target-fishing/reference/human_canonical_targets.jsonl'"
        in config
    )
    assert "target_assets = '/data/airti-target-fishing/reference/targets'" in config
    assert "/data/airti/reference/" not in config


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
