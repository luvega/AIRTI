"""Versioned records exchanged between pipeline stages."""

from __future__ import annotations

from enum import StrEnum
from pathlib import Path
from typing import Annotated, Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator

Sha256 = Annotated[str, Field(pattern=r"^[0-9a-f]{64}$")]


class TargetStatus(StrEnum):
    READY = "ready"
    UNSUPPORTED = "unsupported"
    FAILED = "failed"


class StageStatus(StrEnum):
    PENDING = "pending"
    RUNNING = "running"
    SUCCEEDED = "succeeded"
    FAILED = "failed"
    SKIPPED = "skipped"


class ProvenanceRecord(BaseModel):
    model_config = ConfigDict(extra="forbid")

    schema_version: Literal["1.0"] = "1.0"
    record_id: str = Field(min_length=1)
    input_sha256: Sha256
    tool_version: str = Field(min_length=1)


class TargetRecord(ProvenanceRecord):
    uniprot_id: str = Field(min_length=1)
    sequence: str = Field(pattern=r"^[A-Z]+$")
    status: TargetStatus
    unsupported_reason: str | None = None
    calibrated_score: float | None = None

    @model_validator(mode="after")
    def preserve_missingness(self) -> "TargetRecord":
        if self.status != TargetStatus.READY and self.calibrated_score is not None:
            raise ValueError("unsupported or failed targets cannot receive a numeric score")
        if self.status == TargetStatus.UNSUPPORTED and not self.unsupported_reason:
            raise ValueError("unsupported target requires unsupported_reason")
        return self


class LigandRecord(ProvenanceRecord):
    ligand_id: str
    canonical_smiles: str
    status: StageStatus
    error_code: str | None = None


class PocketRecord(ProvenanceRecord):
    pocket_id: str
    target_id: str
    status: StageStatus
    unsupported_reason: str | None = None


class DockingRecord(ProvenanceRecord):
    docking_id: str
    ligand_state_id: str
    pocket_id: str
    seed: int
    status: StageStatus
    affinity_kcal_mol: float | None = None
    error_code: str | None = None


class BoltzRecord(ProvenanceRecord):
    boltz_id: str
    ligand_state_id: str
    target_id: str
    seed: int
    status: StageStatus
    affinity_score: float | None = None
    confidence_score: float | None = None
    error_code: str | None = None


class MDRecord(ProvenanceRecord):
    md_id: str
    ligand_state_id: str
    target_id: str
    replica: int = Field(gt=0)
    status: StageStatus
    completed_ns: float | None = Field(default=None, ge=0)
    md_score: float | None = None
    error_code: str | None = None


class RankedTarget(ProvenanceRecord):
    target_id: str
    ligand_id: str
    rank: int = Field(gt=0)
    status: StageStatus
    candidate_priority: float = Field(ge=0, le=1)
    evidence_tier: str = "partial"


class ArtifactRecord(ProvenanceRecord):
    artifact_id: str
    task_id: str
    path: Path
    sha256: Sha256
    status: StageStatus
