"""Resumable, full-coverage construction of the human target reference library."""

from __future__ import annotations

import hashlib
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from pydantic import BaseModel, ConfigDict, Field

from airti_tf.manifest_io import (
    content_sha256,
    read_jsonl,
    write_artifact,
    write_jsonl_atomic,
)
from airti_tf.sources.uniprot import UniProtRecord
from airti_tf.stages import TargetPocketRow


class ReferenceBuildContext(BaseModel):
    """Frozen controls and paths passed to one target builder."""

    model_config = ConfigDict(extra="forbid", arbitrary_types_allowed=True)

    root: Path
    background_panel: Path
    max_pockets: int = Field(gt=0, le=10)
    calibration_workers: int = Field(default=1, gt=0, le=64)


TargetBuilder = Callable[
    [UniProtRecord, ReferenceBuildContext], list[TargetPocketRow]
]


class ReferenceBuildSummary(BaseModel):
    model_config = ConfigDict(extra="forbid")

    manifest_path: Path
    manifest_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    target_count: int = Field(gt=0)
    row_count: int = Field(gt=0)
    ready_target_count: int = Field(ge=0)
    unsupported_target_count: int = Field(ge=0)
    failed_target_count: int = Field(ge=0)
    resumed_target_count: int = Field(ge=0)


class ReferenceGateResult(BaseModel):
    model_config = ConfigDict(extra="forbid")

    passed: bool
    expected_target_count: int = Field(gt=0)
    observed_target_count: int = Field(ge=0)
    ready_target_count: int = Field(ge=0)
    unsupported_target_count: int = Field(ge=0)
    failed_target_count: int = Field(ge=0)
    gold_ready: list[str]
    violations: list[str]


class _TargetCheckpoint(BaseModel):
    model_config = ConfigDict(extra="forbid")

    schema_version: str = "1.0"
    target_id: str
    build_fingerprint: str = Field(pattern=r"^[0-9a-f]{64}$")
    rows: list[dict[str, object]]
    error_detail: str | None = None


def _default_target_builder(
    record: UniProtRecord, context: ReferenceBuildContext
) -> list[TargetPocketRow]:
    from airti_tf.targets.production import build_production_target

    return build_production_target(record, context)


def _fingerprint(
    record: UniProtRecord,
    *,
    panel_sha256: str,
    max_pockets: int,
    calibration_workers: int,
) -> str:
    return content_sha256(
        {
            "builder_schema": "1.3",
            "target_id": record.uniprot_id,
            "sequence_sha256": record.sequence_sha256,
            "release": record.release,
            "transmembrane_segments": record.transmembrane_segments,
            "membrane_associated": record.membrane_associated,
            "background_panel_sha256": panel_sha256,
            "max_pockets": max_pockets,
            "calibration_workers": calibration_workers,
        }
    )


def _failed_row(record: UniProtRecord, error: Exception) -> TargetPocketRow:
    return TargetPocketRow(
        schema_version="1.1",
        target_id=record.uniprot_id,
        gene_symbol=record.gene_primary,
        family="unknown",
        status="failed",
        unsupported_reason=f"build_error:{type(error).__name__}",
        environment=(
            "membrane"
            if record.transmembrane_segments or record.membrane_associated
            else "soluble"
        ),
    )


def _validate_target_rows(
    record: UniProtRecord,
    rows: list[TargetPocketRow],
    *,
    max_pockets: int,
) -> list[TargetPocketRow]:
    if not rows:
        raise ValueError("target builder returned no rows")
    if len(rows) > max_pockets:
        raise ValueError("target builder exceeded max_pockets")
    if any(row.target_id != record.uniprot_id for row in rows):
        raise ValueError("target builder returned a row for another accession")
    statuses = {row.status for row in rows}
    if len(statuses) != 1:
        raise ValueError("target builder returned mixed target statuses")
    if next(iter(statuses)) != "ready" and len(rows) != 1:
        raise ValueError("non-ready targets must have exactly one row")
    return rows


def build_reference_targets(
    *,
    proteome_manifest: Path,
    root: Path,
    background_panel: Path,
    max_pockets: int,
    workers: int,
    resume: bool,
    output_manifest: Path,
    calibration_workers: int = 1,
    target_builder: TargetBuilder = _default_target_builder,
) -> ReferenceBuildSummary:
    """Build every canonical target, checkpointing each accession atomically."""
    if workers < 1:
        raise ValueError("workers must be positive")
    if not 1 <= calibration_workers <= 64:
        raise ValueError("calibration_workers must be between 1 and 64")
    if not background_panel.is_file():
        raise ValueError("background probe panel is missing")
    records = [
        UniProtRecord.model_validate(row) for row in read_jsonl(proteome_manifest)
    ]
    if not records:
        raise ValueError("canonical proteome manifest is empty")
    if any(
        record.taxonomy_id != 9606 or not record.reviewed for record in records
    ):
        raise ValueError("reference build requires reviewed human records only")
    accessions = [record.uniprot_id for record in records]
    if len(accessions) != len(set(accessions)):
        raise ValueError("canonical proteome contains duplicate accessions")

    root.mkdir(parents=True, exist_ok=True)
    checkpoint_dir = root / "build-state"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    panel_sha256 = hashlib.sha256(background_panel.read_bytes()).hexdigest()
    context = ReferenceBuildContext(
        root=root.resolve(),
        background_panel=background_panel.resolve(),
        max_pockets=max_pockets,
        calibration_workers=calibration_workers,
    )
    completed: dict[str, list[TargetPocketRow]] = {}
    pending: list[tuple[UniProtRecord, str, Path]] = []
    resumed_target_count = 0
    for record in records:
        fingerprint = _fingerprint(
            record,
            panel_sha256=panel_sha256,
            max_pockets=max_pockets,
            calibration_workers=calibration_workers,
        )
        checkpoint_path = checkpoint_dir / f"{record.uniprot_id}.json"
        if resume and checkpoint_path.is_file():
            try:
                checkpoint = _TargetCheckpoint.model_validate_json(
                    checkpoint_path.read_text(encoding="utf-8")
                )
                rows = [
                    TargetPocketRow.model_validate(row) for row in checkpoint.rows
                ]
                if (
                    checkpoint.target_id == record.uniprot_id
                    and checkpoint.build_fingerprint == fingerprint
                ):
                    completed[record.uniprot_id] = _validate_target_rows(
                        record, rows, max_pockets=max_pockets
                    )
                    resumed_target_count += 1
                    continue
            except (OSError, ValueError):
                pass
        pending.append((record, fingerprint, checkpoint_path))

    def build_one(
        record: UniProtRecord, fingerprint: str, checkpoint_path: Path
    ) -> tuple[str, list[TargetPocketRow]]:
        error_detail: str | None = None
        try:
            rows = _validate_target_rows(
                record,
                target_builder(record, context),
                max_pockets=max_pockets,
            )
        except Exception as error:
            rows = [_failed_row(record, error)]
            error_detail = f"{type(error).__name__}: {error}"[:4000]
        checkpoint = _TargetCheckpoint(
            target_id=record.uniprot_id,
            build_fingerprint=fingerprint,
            rows=[row.model_dump(mode="json") for row in rows],
            error_detail=error_detail,
        )
        write_artifact(checkpoint_path, checkpoint.model_dump(mode="json"))
        return record.uniprot_id, rows

    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(build_one, record, fingerprint, checkpoint_path): record
            for record, fingerprint, checkpoint_path in pending
        }
        for future in as_completed(futures):
            target_id, rows = future.result()
            completed[target_id] = rows

    ordered_rows = [
        row
        for record in sorted(records, key=lambda item: item.uniprot_id)
        for row in sorted(
            completed[record.uniprot_id],
            key=lambda item: item.pocket_id or "",
        )
    ]
    manifest_sha256 = write_jsonl_atomic(
        output_manifest,
        [row.model_dump(mode="json") for row in ordered_rows],
    )
    status_by_target = {
        target_id: rows[0].status for target_id, rows in completed.items()
    }
    return ReferenceBuildSummary(
        manifest_path=output_manifest,
        manifest_sha256=manifest_sha256,
        target_count=len(records),
        row_count=len(ordered_rows),
        ready_target_count=sum(
            status == "ready" for status in status_by_target.values()
        ),
        unsupported_target_count=sum(
            status == "unsupported" for status in status_by_target.values()
        ),
        failed_target_count=sum(
            status == "failed" for status in status_by_target.values()
        ),
        resumed_target_count=resumed_target_count,
    )


def evaluate_reference_gate(
    *,
    proteome_manifest: Path,
    target_manifest: Path,
    expected_target_count: int,
    gold_target_ids: list[str],
) -> ReferenceGateResult:
    """Apply the frozen all-target, zero-failure and gold-readiness gate."""
    proteome = [
        UniProtRecord.model_validate(row) for row in read_jsonl(proteome_manifest)
    ]
    targets = [
        TargetPocketRow.model_validate(row) for row in read_jsonl(target_manifest)
    ]
    canonical_ids = {record.uniprot_id for record in proteome}
    observed_ids = {row.target_id for row in targets}
    status_by_target: dict[str, str] = {}
    violations: list[str] = []
    for row in targets:
        previous = status_by_target.setdefault(row.target_id, row.status)
        if previous != row.status:
            violations.append(f"mixed_status:{row.target_id}")
    if len(proteome) != expected_target_count:
        violations.append(
            f"proteome_count:{len(proteome)}!={expected_target_count}"
        )
    if observed_ids != canonical_ids:
        violations.append("canonical_coverage_mismatch")
    failed_count = sum(
        status == "failed" for status in status_by_target.values()
    )
    if failed_count:
        violations.append(f"failed_targets:{failed_count}")
    gold_ready = sorted(
        target_id
        for target_id in gold_target_ids
        if status_by_target.get(target_id) == "ready"
    )
    for target_id in sorted(set(gold_target_ids) - set(gold_ready)):
        violations.append(f"gold_not_ready:{target_id}")
    return ReferenceGateResult(
        passed=not violations,
        expected_target_count=expected_target_count,
        observed_target_count=len(status_by_target),
        ready_target_count=sum(
            status == "ready" for status in status_by_target.values()
        ),
        unsupported_target_count=sum(
            status == "unsupported" for status in status_by_target.values()
        ),
        failed_target_count=failed_count,
        gold_ready=gold_ready,
        violations=violations,
    )
