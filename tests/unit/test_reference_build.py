import hashlib
import json
from pathlib import Path

from typer.testing import CliRunner

from airti_tf.cli import app
from airti_tf.manifest_io import write_jsonl_atomic
from airti_tf.sources.uniprot import UniProtRecord
from airti_tf.stages import TargetPocketRow
from airti_tf.targets.build import (
    ReferenceBuildContext,
    build_reference_targets,
    evaluate_reference_gate,
)


def _record(target_id: str, *, transmembrane_segments: int = 0) -> UniProtRecord:
    sequence = "MPEPTIDE"
    return UniProtRecord(
        uniprot_id=target_id,
        gene_primary=f"GENE{target_id[-1]}",
        taxonomy_id=9606,
        sequence=sequence,
        sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
        reviewed=True,
        release="2026_03",
        transmembrane_segments=transmembrane_segments,
    )


def _write_proteome(path: Path, records: list[UniProtRecord]) -> None:
    write_jsonl_atomic(path, [record.model_dump(mode="json") for record in records])


def _unsupported(record: UniProtRecord, _context: ReferenceBuildContext) -> list[TargetPocketRow]:
    return [
        TargetPocketRow(
            schema_version="1.1",
            target_id=record.uniprot_id,
            gene_symbol=record.gene_primary,
            family="unknown",
            status="unsupported",
            unsupported_reason="no_structure",
            environment="membrane" if record.transmembrane_segments else "soluble",
        )
    ]


def test_reference_build_is_full_coverage_and_resumable(tmp_path: Path) -> None:
    records = [_record("P00001"), _record("P00002", transmembrane_segments=6)]
    proteome = tmp_path / "proteome.jsonl"
    panel = tmp_path / "panel.smi"
    output = tmp_path / "targets.jsonl"
    _write_proteome(proteome, records)
    panel.write_text("CC probe\n", encoding="utf-8")
    calls: list[str] = []

    def builder(record: UniProtRecord, context: ReferenceBuildContext) -> list[TargetPocketRow]:
        calls.append(record.uniprot_id)
        return _unsupported(record, context)

    first = build_reference_targets(
        proteome_manifest=proteome,
        root=tmp_path,
        background_panel=panel,
        max_pockets=3,
        workers=2,
        resume=True,
        output_manifest=output,
        target_builder=builder,
    )
    second = build_reference_targets(
        proteome_manifest=proteome,
        root=tmp_path,
        background_panel=panel,
        max_pockets=3,
        workers=2,
        resume=True,
        output_manifest=output,
        target_builder=lambda *_args: (_ for _ in ()).throw(AssertionError("rebuilt")),
    )

    rows = [json.loads(line) for line in output.read_text().splitlines()]
    assert calls == ["P00001", "P00002"]
    assert first.target_count == 2
    assert first.row_count == 2
    assert first.unsupported_target_count == 2
    assert first.resumed_target_count == 0
    assert second.resumed_target_count == 2
    assert {row["target_id"] for row in rows} == {"P00001", "P00002"}
    assert next(row for row in rows if row["target_id"] == "P00002")["environment"] == "membrane"


def test_reference_build_records_per_target_failure_without_dropping_target(tmp_path: Path) -> None:
    records = [_record("P00001"), _record("P00002")]
    proteome = tmp_path / "proteome.jsonl"
    panel = tmp_path / "panel.smi"
    output = tmp_path / "targets.jsonl"
    _write_proteome(proteome, records)
    panel.write_text("CC probe\n", encoding="utf-8")

    def builder(record: UniProtRecord, context: ReferenceBuildContext) -> list[TargetPocketRow]:
        if record.uniprot_id == "P00002":
            raise RuntimeError("synthetic failure")
        return _unsupported(record, context)

    summary = build_reference_targets(
        proteome_manifest=proteome,
        root=tmp_path,
        background_panel=panel,
        max_pockets=3,
        workers=1,
        resume=False,
        output_manifest=output,
        target_builder=builder,
    )

    rows = [json.loads(line) for line in output.read_text().splitlines()]
    failed = next(row for row in rows if row["target_id"] == "P00002")
    checkpoint = json.loads(
        (tmp_path / "build-state/P00002.json").read_text()
    )
    assert summary.failed_target_count == 1
    assert failed["status"] == "failed"
    assert failed["unsupported_reason"] == "build_error:RuntimeError"
    assert checkpoint["error_detail"] == "RuntimeError: synthetic failure"


def test_evidence_gate_requires_exact_coverage_no_failures_and_gold_ready(tmp_path: Path) -> None:
    records = [_record("P00001"), _record("P00002")]
    proteome = tmp_path / "proteome.jsonl"
    targets = tmp_path / "targets.jsonl"
    _write_proteome(proteome, records)
    write_jsonl_atomic(
        targets,
        [
            _unsupported(records[0], ReferenceBuildContext(root=tmp_path, background_panel=tmp_path / "x", max_pockets=3))[0].model_dump(mode="json"),
            TargetPocketRow(
                schema_version="1.1",
                target_id="P00002",
                gene_symbol="GENE2",
                family="unknown",
                status="failed",
                unsupported_reason="build_error:RuntimeError",
            ).model_dump(mode="json"),
        ],
    )

    gate = evaluate_reference_gate(
        proteome_manifest=proteome,
        target_manifest=targets,
        expected_target_count=2,
        gold_target_ids=["P00001"],
    )

    assert gate.passed is False
    assert "failed_targets:1" in gate.violations
    assert "gold_not_ready:P00001" in gate.violations


def test_build_reference_cli_forwards_production_controls(monkeypatch, tmp_path: Path) -> None:
    proteome = tmp_path / "proteome.jsonl"
    panel = tmp_path / "panel.smi"
    root = tmp_path / "reference"
    output = root / "targets.jsonl"
    proteome.write_text("{}\n", encoding="utf-8")
    panel.write_text("CC probe\n", encoding="utf-8")
    observed: dict[str, object] = {}

    class Summary:
        failed_target_count = 0

        @staticmethod
        def model_dump_json(indent: int) -> str:
            assert indent == 2
            return '{"target_count": 2}'

    def fake_build(**kwargs):
        observed.update(kwargs)
        return Summary()

    monkeypatch.setattr("airti_tf.targets.build.build_reference_targets", fake_build)
    result = CliRunner().invoke(
        app,
        [
            "targets",
            "build-reference",
            "--proteome",
            str(proteome),
            "--root",
            str(root),
            "--background-panel",
            str(panel),
            "--max-pockets",
            "3",
            "--workers",
            "64",
            "--calibration-workers",
            "8",
            "--resume",
            "--output",
            str(output),
        ],
    )

    assert result.exit_code == 0, result.stdout
    assert observed["max_pockets"] == 3
    assert observed["workers"] == 64
    assert observed["calibration_workers"] == 8
    assert observed["resume"] is True
    assert observed["output_manifest"] == output


def test_gate_reference_cli_writes_gate_and_exits_nonzero_on_failure(
    tmp_path: Path,
) -> None:
    records = [_record("P00001"), _record("P00002")]
    proteome = tmp_path / "proteome.jsonl"
    targets = tmp_path / "targets.jsonl"
    output = tmp_path / "gate.json"
    _write_proteome(proteome, records)
    write_jsonl_atomic(
        targets,
        [
            _unsupported(
                record,
                ReferenceBuildContext(
                    root=tmp_path,
                    background_panel=tmp_path / "x",
                    max_pockets=3,
                ),
            )[0].model_dump(mode="json")
            for record in records
        ],
    )

    result = CliRunner().invoke(
        app,
        [
            "targets",
            "gate-reference",
            "--proteome",
            str(proteome),
            "--targets",
            str(targets),
            "--expected-target-count",
            "2",
            "--gold-target-id",
            "P00001",
            "--output",
            str(output),
        ],
    )

    assert result.exit_code == 4
    payload = json.loads(output.read_text())
    assert payload["passed"] is False
    assert payload["violations"] == ["gold_not_ready:P00001"]
