import hashlib
import json
from pathlib import Path

import pytest

from airti_tf.reporting.render import ArtifactIntegrityError, write_report_delivery


def context_for(root: Path) -> dict[str, object]:
    artifact = root / "source" / "coverage.csv"
    artifact.parent.mkdir(parents=True)
    artifact.write_text("status,count\nready,1\n", encoding="utf-8")
    digest = hashlib.sha256(artifact.read_bytes()).hexdigest()
    return {
        "project_id": "trace-demo",
        "generated_at": "2026-08-12T18:00:00+08:00",
        "query_ligands": ["query-1"],
        "conclusion": "结果支持一个候选靶点，建议实验复核。",
        "coverage": {
            "total": 1,
            "ready": 1,
            "unsupported": 0,
            "artifact_id": "coverage",
        },
        "top_targets": [],
        "artifacts": [
            {
                "artifact_id": "coverage",
                "task_id": "task-1",
                "path": str(artifact.relative_to(root)),
                "sha256": digest,
            }
        ],
        "limitations": ["纯计算证据需实验验证。"],
        "wet_lab_recommendations": ["进行直接结合实验。"],
    }


def test_delivery_manifest_validates_existing_artifacts(tmp_path: Path) -> None:
    context = context_for(tmp_path)
    output = write_report_delivery(context, tmp_path / "final_report", artifact_root=tmp_path)

    assert output.report_path.exists()
    assert output.manifest_path.exists()
    manifest = json.loads(output.manifest_path.read_text())
    assert manifest["release_status"] == "passed"
    assert manifest["artifacts"][0]["task_id"] == "task-1"
    for directory in ["tables", "evidence_cards", "structures", "pymol", "trajectories", "audit"]:
        assert (tmp_path / "final_report" / directory).is_dir()


def test_hash_mismatch_blocks_report_release(tmp_path: Path) -> None:
    context = context_for(tmp_path)
    artifacts = context["artifacts"]
    assert isinstance(artifacts, list)
    artifacts[0]["sha256"] = "0" * 64

    with pytest.raises(ArtifactIntegrityError):
        write_report_delivery(context, tmp_path / "final_report", artifact_root=tmp_path)

