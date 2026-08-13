"""File-based production stage contracts used by the Nextflow CLI bridge."""

from __future__ import annotations

import hashlib
import math
import os
import shutil
import statistics
import subprocess
from collections.abc import Callable
from datetime import UTC, datetime
from pathlib import Path
from typing import Literal

import yaml
from pydantic import BaseModel, ConfigDict, Field, model_validator

from airti_tf.ligands.prepare import prepare_ligand, validate_query_batch
from airti_tf.manifest_io import read_jsonl, write_bytes_atomic, write_jsonl_atomic
from airti_tf.pockets.receptor import DockingBox
from airti_tf.ranking.consensus import TargetEvidence, rank_targets
from airti_tf.reporting.render import write_report_delivery
from airti_tf.refinement.boltz2 import (
    BoltzJob,
    BoltzSeedResult,
    InsufficientBoltzSeedsError,
    build_boltz_yaml,
    run_boltz_seed,
    summarize_boltz_seeds,
)
from airti_tf.screening.calibration import (
    BackgroundDistribution,
    ScreenHit,
    empirical_percentile,
    route_screen_candidates,
)
from airti_tf.screening.quickvina import (
    DockingJob,
    DockingSeedResult,
    InsufficientSeedsError,
    run_docking_seed,
    summarize_seeds,
)
from airti_tf.simulation.gromacs import plan_md_replicas
from airti_tf.state import StateStore

Profile = Literal["local", "production"]
PDBQTPreparer = Callable[[Path, Path], None]
DockingRunner = Callable[[DockingJob, int], DockingSeedResult]
BoltzRunner = Callable[[BoltzJob, int], BoltzSeedResult]


class LigandBundleSummary(BaseModel):
    model_config = ConfigDict(extra="forbid")

    query_count: int = Field(gt=0, le=5)
    state_count: int = Field(ge=0)
    failed_query_count: int = Field(ge=0)
    manifest_path: Path
    manifest_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")


class PreparedLigandRow(BaseModel):
    model_config = ConfigDict(extra="forbid")

    schema_version: Literal["1.0"] = "1.0"
    query_id: str
    ligand_id: str
    ligand_state_id: str | None
    canonical_smiles: str | None
    atom_count: int | None = Field(default=None, gt=0)
    formal_charge: int | None
    status: Literal["succeeded", "failed"]
    error_code: str | None
    uncertainty_flags: list[str]
    sdf_path: Path | None
    pdbqt_path: Path | None

    @model_validator(mode="after")
    def validate_state(self) -> "PreparedLigandRow":
        assets = (
            self.ligand_state_id,
            self.canonical_smiles,
            self.atom_count,
            self.formal_charge,
            self.sdf_path,
            self.pdbqt_path,
        )
        if self.status == "succeeded" and any(value is None for value in assets):
            raise ValueError("successful prepared ligand row lacks required assets")
        if self.status == "failed" and not self.error_code:
            raise ValueError("failed prepared ligand row requires error_code")
        return self


class TargetPocketRow(BaseModel):
    model_config = ConfigDict(extra="forbid")

    schema_version: Literal["1.0"] = "1.0"
    target_id: str
    gene_symbol: str | None = None
    family: str
    status: Literal["ready", "unsupported", "failed"]
    unsupported_reason: str | None = None
    sequence: str | None = Field(default=None, pattern=r"^[A-Z]+$")
    sequence_sha256: str | None = Field(default=None, pattern=r"^[0-9a-f]{64}$")
    model_sequence: str | None = Field(default=None, pattern=r"^[A-Z]+$")
    model_sequence_sha256: str | None = Field(
        default=None, pattern=r"^[0-9a-f]{64}$"
    )
    model_sequence_start: int | None = Field(default=None, gt=0)
    model_sequence_end: int | None = Field(default=None, gt=0)
    structure_quality: float | None = Field(default=None, ge=0, le=1)
    structure_id: str | None = None
    structure_source: Literal["pdb", "alphafold"] | None = None
    structure_path: Path | None = None
    calibration_path: Path | None = None
    pocket_id: str | None = None
    receptor_pdbqt_path: Path | None = None
    box: DockingBox | None = None
    background_affinities: list[float] = Field(default_factory=list)
    msa_path: Path | None = None
    msa_database_version: str | None = None
    pocket_residues: list[int] = Field(default_factory=list)
    model_pocket_residues: list[int] = Field(default_factory=list)

    @model_validator(mode="after")
    def validate_state(self) -> "TargetPocketRow":
        required = (
            self.sequence,
            self.sequence_sha256,
            self.model_sequence,
            self.model_sequence_sha256,
            self.model_sequence_start,
            self.model_sequence_end,
            self.structure_quality,
            self.structure_id,
            self.structure_source,
            self.structure_path,
            self.calibration_path,
            self.pocket_id,
            self.receptor_pdbqt_path,
            self.box,
            self.msa_path,
            self.msa_database_version,
        )
        if self.status == "ready":
            if (
                any(value is None for value in required)
                or not self.pocket_residues
                or not self.model_pocket_residues
            ):
                raise ValueError("ready target pocket row lacks required production assets")
            assert self.sequence is not None
            assert self.model_sequence is not None
            assert self.model_sequence_sha256 is not None
            assert self.model_sequence_start is not None
            assert self.model_sequence_end is not None
            if self.model_sequence_end > len(self.sequence):
                raise ValueError("model sequence coordinates exceed canonical sequence")
            expected_model = self.sequence[
                self.model_sequence_start - 1 : self.model_sequence_end
            ]
            if expected_model != self.model_sequence:
                raise ValueError("model sequence does not match canonical coordinates")
            if hashlib.sha256(self.model_sequence.encode()).hexdigest() != (
                self.model_sequence_sha256
            ):
                raise ValueError("model sequence SHA-256 mismatch")
            if any(
                residue < 1 or residue > len(self.model_sequence)
                for residue in self.model_pocket_residues
            ):
                raise ValueError("model-local pocket residue is out of range")
            expected_canonical_pocket = [
                self.model_sequence_start + residue - 1
                for residue in self.model_pocket_residues
            ]
            if expected_canonical_pocket != self.pocket_residues:
                raise ValueError("canonical and model-local pocket residues disagree")
            BackgroundDistribution(
                pocket_id=str(self.pocket_id), affinities=self.background_affinities
            )
        elif self.status == "unsupported" and not self.unsupported_reason:
            raise ValueError("unsupported target row requires unsupported_reason")
        return self


class ScreenBundleSummary(BaseModel):
    model_config = ConfigDict(extra="forbid")

    query_count: int = Field(gt=0, le=5)
    ready_target_count: int = Field(ge=0)
    unsupported_target_count: int = Field(ge=0)
    failed_target_count: int = Field(ge=0)
    successful_docking_job_count: int = Field(ge=0)
    failed_docking_job_count: int = Field(ge=0)
    candidate_count: int = Field(ge=0)
    manifest_path: Path
    manifest_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")


class TargetCoverageRecord(BaseModel):
    model_config = ConfigDict(extra="forbid")

    total: int = Field(ge=0)
    ready: int = Field(ge=0)
    unsupported: int = Field(ge=0)
    failed: int = Field(ge=0)


class ScreenCandidateRow(BaseModel):
    model_config = ConfigDict(extra="forbid")

    schema_version: Literal["1.0"] = "1.0"
    query_id: str
    ligand_id: str
    ligand_state_id: str
    ligand_smiles: str
    ligand_atom_count: int = Field(gt=0, le=128)
    ligand_formal_charge: int
    target_id: str
    gene_symbol: str | None = None
    family: str
    pocket_id: str
    pocket_residues: list[int] = Field(min_length=1)
    model_pocket_residues: list[int] = Field(min_length=1)
    sequence: str = Field(pattern=r"^[A-Z]+$")
    sequence_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    model_sequence: str = Field(pattern=r"^[A-Z]+$")
    model_sequence_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    model_sequence_start: int = Field(gt=0)
    model_sequence_end: int = Field(gt=0)
    msa_path: Path
    msa_database_version: str
    structure_quality: float = Field(ge=0, le=1)
    affinity_median: float
    calibrated_score: float = Field(ge=0, le=1)
    seed_range: float = Field(ge=0)
    pose_consistency: float = Field(ge=0, le=1)
    seed_success_count: int = Field(ge=2, le=3)
    selection_reason: str
    screen_rank: int = Field(gt=0)
    target_coverage: TargetCoverageRecord


class BoltzRefinementSummary(BaseModel):
    model_config = ConfigDict(extra="forbid")

    query_count: int = Field(gt=0, le=5)
    selected_candidate_count: int = Field(ge=0)
    succeeded_candidate_count: int = Field(ge=0)
    failed_candidate_count: int = Field(ge=0)
    manifest_path: Path
    manifest_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")


class BoltzCandidateRow(ScreenCandidateRow):
    boltz_status: Literal["succeeded", "failed"]
    boltz_error_code: str | None
    boltz_seed_errors: dict[str, str | None]
    boltz_seed_success_count: int = Field(ge=0, le=3)
    boltz_confidence_median: float | None = Field(default=None, ge=0, le=1)
    boltz_ligand_iptm_median: float | None = Field(default=None, ge=0, le=1)
    boltz_affinity_probability_median: float | None = Field(default=None, ge=0, le=1)
    boltz_affinity_pred_value_median: float | None = None
    boltz_pocket_constraint_median: float | None = Field(default=None, ge=0, le=1)
    boltz_confidence_range: float | None = Field(default=None, ge=0)
    boltz_score: float | None = Field(default=None, ge=0, le=1)
    boltz_structure_path: Path | None = None

    @model_validator(mode="after")
    def validate_boltz_state(self) -> "BoltzCandidateRow":
        if self.boltz_status == "succeeded" and (
            self.boltz_score is None or self.boltz_structure_path is None
        ):
            raise ValueError("successful Boltz candidate lacks score or structure")
        if self.boltz_status == "failed" and not self.boltz_error_code:
            raise ValueError("failed Boltz candidate requires error_code")
        return self


class MDStageResult(BaseModel):
    model_config = ConfigDict(extra="forbid")

    replica: int = Field(gt=0)
    status: Literal["stable", "unstable", "failed"]
    completed_ns: float = Field(ge=0)
    md_score: float | None = Field(default=None, ge=0, le=1)
    error_code: str | None
    trajectory_path: Path | None = None
    checkpoint_path: Path | None = None

    @model_validator(mode="after")
    def validate_md_state(self) -> "MDStageResult":
        if self.status in {"stable", "unstable"} and (
            self.md_score is None
            or self.trajectory_path is None
            or self.checkpoint_path is None
        ):
            raise ValueError("successful MD replica lacks score or required artifacts")
        if self.status == "failed" and not self.error_code:
            raise ValueError("failed MD replica requires error_code")
        return self


MDRunner = Callable[[BoltzCandidateRow, Path, int], MDStageResult]


class MDBundleSummary(BaseModel):
    model_config = ConfigDict(extra="forbid")

    query_count: int = Field(gt=0, le=5)
    selected_candidate_count: int = Field(ge=0)
    planned_replica_count: int = Field(ge=0)
    succeeded_candidate_count: int = Field(ge=0)
    failed_candidate_count: int = Field(ge=0)
    manifest_path: Path
    manifest_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")


class MDCandidateRow(BoltzCandidateRow):
    md_rank: int = Field(gt=0)
    md_status: Literal["stable", "unstable", "failed"]
    md_error_code: str | None
    md_replica_success_count: int = Field(ge=0, le=3)
    completed_ns: float = Field(ge=0)
    md_score: float | None = Field(default=None, ge=0, le=1)
    md_replicas: list[MDStageResult]

    @model_validator(mode="after")
    def validate_final_state(self) -> "MDCandidateRow":
        if self.md_status in {"stable", "unstable"} and self.md_score is None:
            raise ValueError("successful MD candidate lacks md_score")
        if self.md_status == "failed" and not self.md_error_code:
            raise ValueError("failed MD candidate requires error_code")
        return self


class ReportBundleSummary(BaseModel):
    model_config = ConfigDict(extra="forbid")

    project_id: str
    query_count: int = Field(gt=0, le=5)
    reported_target_count: int = Field(ge=0)
    coverage: TargetCoverageRecord
    report_path: Path
    manifest_path: Path
    state_db: Path
    report_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")


def _parse_smiles_queries(path: Path, *, max_molecules: int) -> list[tuple[str, str]]:
    if not 1 <= max_molecules <= 5:
        raise ValueError("max_molecules must be between 1 and 5")
    queries: list[tuple[str, str]] = []
    for line_number, raw_line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) != 2:
            raise ValueError(
                f"SMILES input line {line_number} must contain exactly SMILES and query_id"
            )
        smiles, query_id = fields
        queries.append((query_id, smiles))
    if len(queries) > max_molecules:
        raise ValueError(f"query batch must contain 1 to {max_molecules} molecules")
    validate_query_batch([smiles for _, smiles in queries])
    identifiers = [query_id for query_id, _ in queries]
    if len(set(identifiers)) != len(identifiers):
        raise ValueError("query identifiers must be unique")
    return queries


def _run_meeko(input_sdf: Path, output_pdbqt: Path) -> None:
    process = subprocess.run(
        [
            "mk_prepare_ligand.py",
            "--mol",
            str(input_sdf),
            "--out",
            str(output_pdbqt),
            "--add_index_map",
        ],
        check=False,
        capture_output=True,
        text=True,
        timeout=300,
    )
    if process.returncode != 0:
        raise RuntimeError(
            "Meeko ligand preparation failed: "
            + (process.stderr.strip() or process.stdout.strip() or "unknown error")
        )
    if not output_pdbqt.is_file() or output_pdbqt.stat().st_size == 0:
        raise RuntimeError("Meeko ligand preparation did not create a non-empty PDBQT")


def _portable_path(path: Path, *, relative_to: Path) -> str:
    return os.path.relpath(path.resolve(), start=relative_to.resolve())


def prepare_ligand_bundle(
    input_path: Path,
    *,
    output_manifest: Path,
    asset_dir: Path,
    profile: Profile,
    max_molecules: int,
    pdbqt_preparer: PDBQTPreparer = _run_meeko,
) -> LigandBundleSummary:
    """Prepare a service-sized SMILES batch and write a portable state manifest."""
    queries = _parse_smiles_queries(input_path, max_molecules=max_molecules)
    asset_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    failed_queries = 0
    for query_id, smiles in queries:
        result = prepare_ligand(smiles, profile=profile)
        if result.status == "failed":
            failed_queries += 1
            rows.append(
                {
                    "schema_version": "1.0",
                    "query_id": query_id,
                    "ligand_id": result.ligand_id,
                    "ligand_state_id": None,
                    "canonical_smiles": None,
                    "atom_count": None,
                    "formal_charge": None,
                    "status": "failed",
                    "error_code": result.error_code,
                    "uncertainty_flags": result.uncertainty_flags,
                    "sdf_path": None,
                    "pdbqt_path": None,
                }
            )
            continue
        for state in result.states:
            stem = state.ligand_state_id
            sdf_path = asset_dir / f"{stem}.sdf"
            pdbqt_path = asset_dir / f"{stem}.pdbqt"
            sdf_path.write_text(state.mol_block.rstrip() + "\n$$$$\n", encoding="utf-8")
            pdbqt_preparer(sdf_path, pdbqt_path)
            rows.append(
                {
                    "schema_version": "1.0",
                    "query_id": query_id,
                    "ligand_id": result.ligand_id,
                    "ligand_state_id": state.ligand_state_id,
                    "canonical_smiles": state.canonical_smiles,
                    "atom_count": state.atom_count,
                    "formal_charge": state.formal_charge,
                    "status": "succeeded",
                    "error_code": None,
                    "uncertainty_flags": result.uncertainty_flags,
                    "sdf_path": _portable_path(sdf_path, relative_to=output_manifest.parent),
                    "pdbqt_path": _portable_path(
                        pdbqt_path, relative_to=output_manifest.parent
                    ),
                }
            )
    manifest_sha256 = write_jsonl_atomic(output_manifest, rows)
    return LigandBundleSummary(
        query_count=len(queries),
        state_count=sum(row["status"] == "succeeded" for row in rows),
        failed_query_count=failed_queries,
        manifest_path=output_manifest,
        manifest_sha256=manifest_sha256,
    )


def _default_docking_runner(job: DockingJob, seed: int) -> DockingSeedResult:
    return run_docking_seed(job, seed=seed)


def _resolve_manifest_path(path: Path, *, manifest: Path) -> Path:
    return path if path.is_absolute() else manifest.parent / path


def _pdbqt_coordinates(path: Path) -> list[tuple[float, float, float]]:
    coordinates: list[tuple[float, float, float]] = []
    if not path.is_file():
        return coordinates
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            coordinates.append(
                (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            )
        except ValueError:
            continue
    return coordinates


def _pose_consistency(job: DockingJob, successful_seeds: list[int]) -> float:
    poses = [
        _pdbqt_coordinates(job.output_dir / f"{job.job_id}.seed{seed}.poses.pdbqt")
        for seed in successful_seeds
    ]
    if len(poses) < 2 or any(not pose for pose in poses):
        return 0.0
    if len({len(pose) for pose in poses}) != 1:
        return 0.0
    pairwise: list[float] = []
    for first_index, first in enumerate(poses):
        for second in poses[first_index + 1 :]:
            squared = sum(
                (left[axis] - right[axis]) ** 2
                for left, right in zip(first, second, strict=True)
                for axis in range(3)
            )
            pairwise.append(math.sqrt(squared / len(first)))
    median_rmsd = sorted(pairwise)[len(pairwise) // 2]
    return round(max(0.0, min(1.0, 1.0 - median_rmsd / 4.0)), 6)


def screen_ligand_bundle(
    ligand_manifest: Path,
    target_manifest: Path,
    *,
    output_manifest: Path,
    asset_dir: Path,
    top_n: int,
    docking_runner: DockingRunner = _default_docking_runner,
) -> ScreenBundleSummary:
    """Dock all prepared states against ready pockets and route per-query hits."""
    ligands = [PreparedLigandRow.model_validate(row) for row in read_jsonl(ligand_manifest)]
    if any(row.status == "failed" for row in ligands):
        raise ValueError("ligand manifest contains failed queries; screening is blocked")
    targets = [TargetPocketRow.model_validate(row) for row in read_jsonl(target_manifest)]
    ready_targets = [row for row in targets if row.status == "ready"]
    query_ids = sorted({row.query_id for row in ligands})
    if not query_ids:
        raise ValueError("ligand manifest contains no prepared states")
    if not ready_targets:
        raise ValueError("target manifest contains no ready pockets")

    asset_dir.mkdir(parents=True, exist_ok=True)
    hits_by_query: dict[str, list[ScreenHit]] = {query_id: [] for query_id in query_ids}
    source_by_key: dict[tuple[str, str, str, str], tuple[PreparedLigandRow, TargetPocketRow, int]] = {}
    successful_jobs = 0
    failed_jobs = 0
    for ligand in ligands:
        assert ligand.ligand_state_id is not None
        assert ligand.pdbqt_path is not None
        for target in ready_targets:
            assert target.pocket_id is not None
            assert target.receptor_pdbqt_path is not None
            assert target.box is not None
            job_material = (
                f"{ligand.query_id}|{ligand.ligand_state_id}|"
                f"{target.target_id}|{target.pocket_id}"
            )
            job_id = hashlib.sha256(job_material.encode()).hexdigest()[:20]
            job = DockingJob(
                job_id=job_id,
                receptor_pdbqt=_resolve_manifest_path(
                    target.receptor_pdbqt_path, manifest=target_manifest
                ),
                ligand_pdbqt=_resolve_manifest_path(
                    ligand.pdbqt_path, manifest=ligand_manifest
                ),
                box=target.box,
                output_dir=asset_dir / job_id,
            )
            seed_records = [docking_runner(job, seed) for seed in (11, 29, 47)]
            try:
                summary = summarize_seeds(seed_records)
            except InsufficientSeedsError:
                failed_jobs += 1
                continue
            successful_jobs += 1
            background = BackgroundDistribution(
                pocket_id=target.pocket_id,
                affinities=target.background_affinities,
            )
            hit = ScreenHit(
                target_id=target.target_id,
                family=target.family,
                pocket_id=target.pocket_id,
                ligand_state_id=ligand.ligand_state_id,
                affinity_median=summary.affinity_median,
                calibrated_score=empirical_percentile(
                    query=summary.affinity_median,
                    background=background.affinities,
                ),
                seed_range=summary.affinity_range,
                pose_consistency=_pose_consistency(job, summary.successful_seeds),
            )
            hits_by_query[ligand.query_id].append(hit)
            source_by_key[
                (
                    ligand.query_id,
                    target.target_id,
                    target.pocket_id,
                    ligand.ligand_state_id,
                )
            ] = (ligand, target, summary.seed_success_count)

    rows: list[dict[str, object]] = []
    target_statuses = {row.target_id: row.status for row in targets}
    target_coverage = TargetCoverageRecord(
        total=len(target_statuses),
        ready=sum(status == "ready" for status in target_statuses.values()),
        unsupported=sum(
            status == "unsupported" for status in target_statuses.values()
        ),
        failed=sum(status == "failed" for status in target_statuses.values()),
    )
    for query_id in query_ids:
        routed = route_screen_candidates(hits_by_query[query_id], top_n=top_n)
        for rank, candidate in enumerate(routed, 1):
            ligand, target, seed_success_count = source_by_key[
                (
                    query_id,
                    candidate.target_id,
                    candidate.best_pocket_id,
                    candidate.best_state_id,
                )
            ]
            assert ligand.canonical_smiles is not None
            assert ligand.atom_count is not None
            assert ligand.formal_charge is not None
            assert target.sequence is not None
            assert target.sequence_sha256 is not None
            assert target.model_sequence is not None
            assert target.model_sequence_sha256 is not None
            assert target.model_sequence_start is not None
            assert target.model_sequence_end is not None
            assert target.structure_quality is not None
            assert target.msa_path is not None
            assert target.msa_database_version is not None
            source_msa = _resolve_manifest_path(
                target.msa_path, manifest=target_manifest
            )
            msa_dir = asset_dir / "msa"
            msa_dir.mkdir(exist_ok=True)
            staged_msa = msa_dir / (
                f"{target.target_id}.{target.model_sequence_sha256[:12]}.a3m"
            )
            if not staged_msa.is_file():
                shutil.copyfile(source_msa, staged_msa)
            rows.append(
                {
                    "schema_version": "1.0",
                    "query_id": query_id,
                    "ligand_id": ligand.ligand_id,
                    "ligand_state_id": candidate.best_state_id,
                    "ligand_smiles": ligand.canonical_smiles,
                    "ligand_atom_count": ligand.atom_count,
                    "ligand_formal_charge": ligand.formal_charge,
                    "target_id": candidate.target_id,
                    "gene_symbol": target.gene_symbol,
                    "family": candidate.family,
                    "pocket_id": candidate.best_pocket_id,
                    "pocket_residues": target.pocket_residues,
                    "model_pocket_residues": target.model_pocket_residues,
                    "sequence": target.sequence,
                    "sequence_sha256": target.sequence_sha256,
                    "model_sequence": target.model_sequence,
                    "model_sequence_sha256": target.model_sequence_sha256,
                    "model_sequence_start": target.model_sequence_start,
                    "model_sequence_end": target.model_sequence_end,
                    "msa_path": _portable_path(
                        staged_msa, relative_to=output_manifest.parent
                    ),
                    "msa_database_version": target.msa_database_version,
                    "structure_quality": target.structure_quality,
                    "affinity_median": candidate.affinity_median,
                    "calibrated_score": candidate.calibrated_score,
                    "seed_range": candidate.seed_range,
                    "pose_consistency": candidate.pose_consistency,
                    "seed_success_count": seed_success_count,
                    "selection_reason": candidate.selection_reason,
                    "screen_rank": rank,
                    "target_coverage": target_coverage.model_dump(mode="json"),
                }
            )
    digest = write_jsonl_atomic(output_manifest, rows)
    return ScreenBundleSummary(
        query_count=len(query_ids),
        ready_target_count=target_coverage.ready,
        unsupported_target_count=target_coverage.unsupported,
        failed_target_count=target_coverage.failed,
        successful_docking_job_count=successful_jobs,
        failed_docking_job_count=failed_jobs,
        candidate_count=len(rows),
        manifest_path=output_manifest,
        manifest_sha256=digest,
    )


def _default_boltz_runner(job: BoltzJob, seed: int) -> BoltzSeedResult:
    return run_boltz_seed(job, seed=seed)


def refine_boltz_bundle(
    screen_manifest: Path,
    *,
    output_manifest: Path,
    asset_dir: Path,
    profile: Profile,
    top_n: int,
    cache_path: Path | None = None,
    boltz_runner: BoltzRunner = _default_boltz_runner,
) -> BoltzRefinementSummary:
    """Refine each query's routed candidates with a three-seed Boltz consensus."""
    candidates = [
        ScreenCandidateRow.model_validate(row) for row in read_jsonl(screen_manifest)
    ]
    query_ids = sorted({row.query_id for row in candidates})
    if not query_ids:
        raise ValueError("screen manifest contains no candidates")
    selected: list[ScreenCandidateRow] = []
    for query_id in query_ids:
        per_query = sorted(
            (row for row in candidates if row.query_id == query_id),
            key=lambda row: (row.screen_rank, row.target_id, row.pocket_id),
        )
        selected.extend(per_query[:top_n])

    asset_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    succeeded_count = 0
    failed_count = 0
    for candidate in selected:
        job_material = (
            f"{candidate.query_id}|{candidate.ligand_state_id}|"
            f"{candidate.target_id}|{candidate.pocket_id}"
        )
        job_id = hashlib.sha256(job_material.encode()).hexdigest()[:20]
        seed_results: list[BoltzSeedResult] = []
        for seed in (11, 29, 47):
            seed_dir = asset_dir / job_id / f"seed-{seed}"
            seed_dir.mkdir(parents=True, exist_ok=True)
            job = BoltzJob(
                job_id=job_id,
                target_id=candidate.target_id,
                sequence=candidate.model_sequence,
                sequence_sha256=candidate.model_sequence_sha256,
                msa_path=_resolve_manifest_path(
                    candidate.msa_path, manifest=screen_manifest
                ),
                msa_database_version=candidate.msa_database_version,
                ligand_state_id=candidate.ligand_state_id,
                ligand_smiles=candidate.ligand_smiles,
                ligand_atom_count=candidate.ligand_atom_count,
                pocket_residues=candidate.model_pocket_residues,
                input_yaml=seed_dir / f"{job_id}.yaml",
                output_dir=seed_dir / "output",
                cache_path=cache_path,
            )
            payload = build_boltz_yaml(job, profile=profile)
            write_bytes_atomic(
                job.input_yaml,
                yaml.safe_dump(payload, sort_keys=False).encode("utf-8"),
            )
            seed_results.append(boltz_runner(job, seed))

        base = candidate.model_dump(mode="json")
        errors = {
            str(result.seed): result.error_code
            for result in seed_results
            if result.status == "failed"
        }
        try:
            summary = summarize_boltz_seeds(seed_results)
        except InsufficientBoltzSeedsError:
            failed_count += 1
            rows.append(
                {
                    **base,
                    "boltz_status": "failed",
                    "boltz_error_code": "insufficient_successful_seeds",
                    "boltz_seed_errors": errors,
                    "boltz_seed_success_count": sum(
                        result.status == "succeeded" and result.severe_clash is False
                        for result in seed_results
                    ),
                    "boltz_confidence_median": None,
                    "boltz_ligand_iptm_median": None,
                    "boltz_affinity_probability_median": None,
                    "boltz_affinity_pred_value_median": None,
                    "boltz_pocket_constraint_median": None,
                    "boltz_confidence_range": None,
                    "boltz_score": None,
                    "boltz_structure_path": None,
                }
            )
            continue

        succeeded_count += 1
        usable = [
            result
            for result in seed_results
            if result.status == "succeeded"
            and result.severe_clash is False
            and result.structure_path is not None
        ]
        representative = sorted(
            usable,
            key=lambda result: (
                -(result.confidence_score or 0.0),
                result.seed,
            ),
        )[0]
        assert representative.structure_path is not None
        boltz_score = round(
            (
                summary.confidence_median
                + summary.ligand_iptm_median
                + summary.affinity_probability_median
                + summary.pocket_constraint_median
            )
            / 4,
            6,
        )
        rows.append(
            {
                **base,
                "boltz_status": "succeeded",
                "boltz_error_code": None,
                "boltz_seed_errors": errors,
                "boltz_seed_success_count": summary.seed_success_count,
                "boltz_confidence_median": summary.confidence_median,
                "boltz_ligand_iptm_median": summary.ligand_iptm_median,
                "boltz_affinity_probability_median": (
                    summary.affinity_probability_median
                ),
                "boltz_affinity_pred_value_median": (
                    summary.affinity_pred_value_median
                ),
                "boltz_pocket_constraint_median": summary.pocket_constraint_median,
                "boltz_confidence_range": summary.confidence_range,
                "boltz_score": boltz_score,
                "boltz_structure_path": _portable_path(
                    representative.structure_path,
                    relative_to=output_manifest.parent,
                ),
            }
        )
    digest = write_jsonl_atomic(output_manifest, rows)
    return BoltzRefinementSummary(
        query_count=len(query_ids),
        selected_candidate_count=len(selected),
        succeeded_candidate_count=succeeded_count,
        failed_candidate_count=failed_count,
        manifest_path=output_manifest,
        manifest_sha256=digest,
    )


def _unimplemented_md_runner(
    candidate: BoltzCandidateRow,
    run_dir: Path,
    replica: int,
) -> MDStageResult:
    del candidate, run_dir
    return MDStageResult(
        replica=replica,
        status="failed",
        completed_ns=0.0,
        md_score=None,
        error_code="md_system_builder_unimplemented",
    )


def run_md_bundle(
    boltz_manifest: Path,
    *,
    output_manifest: Path,
    asset_dir: Path,
    top_n: int,
    md_runner: MDRunner = _unimplemented_md_runner,
) -> MDBundleSummary:
    """Route Boltz candidates into Top10/Top3 MD replicas and aggregate evidence."""
    candidates = [
        BoltzCandidateRow.model_validate(row) for row in read_jsonl(boltz_manifest)
    ]
    query_ids = sorted({row.query_id for row in candidates})
    if not query_ids:
        raise ValueError("Boltz manifest contains no candidates")
    asset_dir.mkdir(parents=True, exist_ok=True)
    selected: list[BoltzCandidateRow] = []
    rank_by_key: dict[tuple[str, str], int] = {}
    for query_id in query_ids:
        eligible = [
            row
            for row in candidates
            if row.query_id == query_id and row.boltz_status == "succeeded"
        ]
        evidence = [
            TargetEvidence(
                target_id=row.target_id,
                ligand_id=row.ligand_id,
                status="ready",
                vina_score=row.calibrated_score,
                docking_consistency=row.pose_consistency,
                structure_quality=row.structure_quality,
                boltz_score=row.boltz_score,
                successful_seeds=row.boltz_seed_success_count,
                boltz_seed_spread=row.boltz_confidence_range,
                heavy_atom_count=row.ligand_atom_count,
            )
            for row in eligible
        ]
        ranked = rank_targets(evidence, stage="boltz").ranked[:top_n]
        lookup = {row.target_id: row for row in eligible}
        for rank, item in enumerate(ranked, 1):
            selected.append(lookup[item.target_id])
            rank_by_key[(query_id, item.target_id)] = rank

    replica_results: dict[tuple[str, str], list[MDStageResult]] = {}
    planned_count = 0
    for query_id in query_ids:
        per_query = sorted(
            (row for row in selected if row.query_id == query_id),
            key=lambda row: rank_by_key[(query_id, row.target_id)],
        )
        candidate_by_target = {row.target_id: row for row in per_query}
        for plan in plan_md_replicas([row.target_id for row in per_query]):
            candidate = candidate_by_target[plan.target_id]
            run_material = f"{query_id}|{candidate.target_id}|{plan.replica}"
            run_id = hashlib.sha256(run_material.encode()).hexdigest()[:20]
            result = md_runner(
                candidate,
                asset_dir / run_id,
                plan.replica,
            )
            replica_results.setdefault((query_id, candidate.target_id), []).append(
                result
            )
            planned_count += 1

    rows: list[dict[str, object]] = []
    succeeded_count = 0
    failed_count = 0
    for candidate in selected:
        key = (candidate.query_id, candidate.target_id)
        results = sorted(replica_results[key], key=lambda item: item.replica)
        valid = [result for result in results if result.status in {"stable", "unstable"}]
        required = 2 if len(results) == 3 else 1
        replica_payload = []
        for result in results:
            payload = result.model_dump(mode="json")
            for path_field in ("trajectory_path", "checkpoint_path"):
                raw_path = getattr(result, path_field)
                payload[path_field] = (
                    _portable_path(raw_path, relative_to=output_manifest.parent)
                    if raw_path is not None
                    else None
                )
            replica_payload.append(payload)
        base = candidate.model_dump(mode="json")
        if len(valid) < required:
            failed_count += 1
            rows.append(
                {
                    **base,
                    "md_rank": rank_by_key[key],
                    "md_status": "failed",
                    "md_error_code": "insufficient_successful_replicas",
                    "md_replica_success_count": len(valid),
                    "completed_ns": max(
                        (result.completed_ns for result in results), default=0.0
                    ),
                    "md_score": None,
                    "md_replicas": replica_payload,
                }
            )
            continue
        succeeded_count += 1
        md_scores = [result.md_score for result in valid if result.md_score is not None]
        rows.append(
            {
                **base,
                "md_rank": rank_by_key[key],
                "md_status": (
                    "unstable" if any(result.status == "unstable" for result in valid) else "stable"
                ),
                "md_error_code": None,
                "md_replica_success_count": len(valid),
                "completed_ns": min(result.completed_ns for result in valid),
                "md_score": statistics.median(md_scores),
                "md_replicas": replica_payload,
            }
        )
    digest = write_jsonl_atomic(output_manifest, rows)
    return MDBundleSummary(
        query_count=len(query_ids),
        selected_candidate_count=len(selected),
        planned_replica_count=planned_count,
        succeeded_candidate_count=succeeded_count,
        failed_candidate_count=failed_count,
        manifest_path=output_manifest,
        manifest_sha256=digest,
    )


def render_report_bundle(
    md_manifest: Path,
    *,
    output_dir: Path,
    state_db: Path,
    project_id: str | None,
) -> ReportBundleSummary:
    """Rank final candidates and render a traceable computation-only report."""
    candidates = [
        MDCandidateRow.model_validate(row) for row in read_jsonl(md_manifest)
    ]
    query_ids = sorted({row.query_id for row in candidates})
    if not query_ids:
        raise ValueError("MD manifest contains no candidates")
    coverage_values = {
        tuple(row.target_coverage.model_dump().values()) for row in candidates
    }
    if len(coverage_values) != 1:
        raise ValueError("MD candidates carry inconsistent target coverage")
    coverage = candidates[0].target_coverage
    input_sha256 = hashlib.sha256(md_manifest.read_bytes()).hexdigest()
    resolved_project_id = project_id or f"AIRTI-{input_sha256[:12].upper()}"
    store = StateStore(state_db)
    task_id = store.register(
        "render-report", input_sha256, run_id=resolved_project_id
    )

    artifact_id = "md-candidates"
    top_targets: list[dict[str, object]] = []
    for query_id in query_ids:
        per_query = [row for row in candidates if row.query_id == query_id]
        ranked = rank_targets(
            [
                TargetEvidence(
                    target_id=row.target_id,
                    ligand_id=row.ligand_id,
                    status="ready",
                    vina_score=row.calibrated_score,
                    docking_consistency=row.pose_consistency,
                    structure_quality=row.structure_quality,
                    boltz_score=row.boltz_score,
                    md_score=row.md_score,
                    md_status=row.md_status,
                    successful_seeds=row.boltz_seed_success_count,
                    boltz_seed_spread=row.boltz_confidence_range,
                    heavy_atom_count=row.ligand_atom_count,
                )
                for row in per_query
            ],
            stage="final",
        ).ranked[:5]
        for item in ranked:
            top_targets.append(
                {
                    "rank": item.rank,
                    "query_id": query_id,
                    "target_id": item.target_id,
                    "gene_symbol": next(
                        (
                            row.gene_symbol
                            for row in per_query
                            if row.target_id == item.target_id and row.gene_symbol
                        ),
                        item.target_id,
                    ),
                    "priority": item.priority,
                    "evidence_tier": item.evidence_tier,
                    "vina_score": item.vina_score,
                    "boltz_score": item.boltz_score,
                    "md_score": item.md_score,
                    "uncertainty_flags": item.uncertainty_flags,
                    "artifact_id": artifact_id,
                }
            )

    context = {
        "project_id": resolved_project_id,
        "generated_at": datetime.now(UTC).isoformat(),
        "query_ligands": query_ids,
        "conclusion": (
            "结果形成结构与动力学计算证据支持的候选优先级，"
            "建议开展正交实验验证。"
        ),
        "coverage": {
            **coverage.model_dump(mode="json"),
            "artifact_id": artifact_id,
        },
        "top_targets": top_targets,
        "artifacts": [
            {
                "artifact_id": artifact_id,
                "task_id": task_id,
                "path": _portable_path(
                    md_manifest, relative_to=md_manifest.parent
                ),
                "sha256": input_sha256,
            }
        ],
        "limitations": [
            "本报告仅汇总计算证据，不能替代直接结合、细胞内靶点占有或表型因果实验。",
            "unsupported 与 failed 蛋白不参与数值排序，其数量作为覆盖率单独报告。",
            "最终排序受结构质量、口袋定义、构象采样和打分函数适用域限制。",
        ],
        "wet_lab_recommendations": [
            "对每个查询分子的 Top 5 候选开展 CETSA/TPP 与 SPR/MST 等正交验证。",
            "对通过结合验证的候选开展 CRISPR/回补和表型救援实验。",
        ],
    }

    if not store.claim(task_id):
        report_path = output_dir / "report.md"
        manifest_path = output_dir / "report_manifest.json"
        if not report_path.is_file() or not manifest_path.is_file():
            raise RuntimeError("report task already claimed but delivery is incomplete")
        report_sha256 = hashlib.sha256(report_path.read_bytes()).hexdigest()
        return ReportBundleSummary(
            project_id=resolved_project_id,
            query_count=len(query_ids),
            reported_target_count=len(top_targets),
            coverage=coverage,
            report_path=report_path,
            manifest_path=manifest_path,
            state_db=state_db,
            report_sha256=report_sha256,
        )

    try:
        delivery = write_report_delivery(
            context,
            output_dir,
            artifact_root=md_manifest.parent,
        )
        store.transition(task_id, "succeeded", output_hash=delivery.report_sha256)
        store.register_artifact(
            task_id,
            sha256=input_sha256,
            path=str(md_manifest.resolve()),
        )
        store.register_artifact(
            task_id,
            sha256=delivery.report_sha256,
            path=str(delivery.report_path.resolve()),
        )
    except Exception as error:
        store.transition(
            task_id,
            "failed",
            error_code=type(error).__name__,
        )
        raise
    return ReportBundleSummary(
        project_id=resolved_project_id,
        query_count=len(query_ids),
        reported_target_count=len(top_targets),
        coverage=coverage,
        report_path=delivery.report_path,
        manifest_path=delivery.manifest_path,
        state_db=state_db,
        report_sha256=delivery.report_sha256,
    )
