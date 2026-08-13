import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest

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
SMOKE_BATCHES = (
    Path("data/benchmark/smoke_v1_batch_a.smi"),
    Path("data/benchmark/smoke_v1_batch_b.smi"),
)


def test_ten_case_smoke_is_partitioned_into_service_sized_batches() -> None:
    identifiers: list[str] = []
    for batch in SMOKE_BATCHES:
        rows = [line.split() for line in batch.read_text().splitlines() if line.strip()]
        assert 1 <= len(rows) <= 5
        identifiers.extend(row[1] for row in rows)

    assert len(identifiers) == 10
    assert len(set(identifiers)) == 10


@pytest.mark.skipif(NEXTFLOW is None, reason="Nextflow not installed")
def test_ten_case_orchestration_smoke_delivery_and_resume(tmp_path: Path) -> None:
    run_root = tmp_path / "orchestration-smoke"
    process = subprocess.run(
        ["./scripts/run_orchestration_smoke.sh", str(run_root)],
        capture_output=True,
        text=True,
        timeout=180,
        env=NEXTFLOW_ENV,
    )
    assert process.returncode == 0, process.stdout + process.stderr

    summary = json.loads((run_root / "orchestration_summary.json").read_text())
    assert summary["validation_scope"] == "orchestration_mock_only"
    assert summary["query_count"] == 10
    assert summary["batch_count"] == 2
    assert summary["max_queries_per_batch"] == 5
    assert summary["technical_success_rate"] >= 0.95
    assert summary["all_metrics_traceable"] is True
    assert summary["unsupported_targets_have_no_numeric_score"] is True
    assert summary["resumed_without_recomputation"] is True
    assert (run_root / "artifacts.sha256").stat().st_size > 0
