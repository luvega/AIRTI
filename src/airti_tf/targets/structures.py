"""Deterministic target-structure quality gates and selection."""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator

StructureSource = Literal["pdb", "alphafold"]
UnsupportedReason = Literal[
    "no_structure",
    "low_coverage",
    "low_confidence",
    "sequence_mismatch",
    "unsupported_chemistry",
]


class StructureCandidate(BaseModel):
    model_config = ConfigDict(extra="forbid")

    structure_id: str
    source: StructureSource
    coverage: float = Field(ge=0, le=1)
    sequence_identity: float = Field(default=1.0, ge=0, le=1)
    mainchain_missing_fraction: float = Field(default=0.0, ge=0, le=1)
    resolution: float | None = Field(default=None, gt=0)
    has_ligand: bool = False
    confidence: float | None = Field(default=None, ge=0, le=1)
    pae_supported: bool = True
    unsupported_chemistry: bool = False


class StructureSelection(BaseModel):
    model_config = ConfigDict(extra="forbid")

    status: Literal["ready", "unsupported"]
    structure_id: str | None = None
    source: StructureSource | None = None
    score: float | None = Field(default=None, ge=0, le=1)
    unsupported_reason: UnsupportedReason | None = None

    @model_validator(mode="after")
    def validate_state(self) -> "StructureSelection":
        if self.status == "ready" and (
            self.structure_id is None or self.source is None or self.score is None
        ):
            raise ValueError("ready structure selection requires structure and score")
        if self.status == "unsupported" and (
            self.unsupported_reason is None or self.score is not None
        ):
            raise ValueError("unsupported structure requires reason and no score")
        return self


def _eligibility(
    candidate: StructureCandidate,
) -> tuple[int, float] | UnsupportedReason:
    if candidate.unsupported_chemistry:
        return "unsupported_chemistry"
    if candidate.sequence_identity < 0.98:
        return "sequence_mismatch"
    if candidate.coverage < 0.70 or candidate.mainchain_missing_fraction > 0.10:
        return "low_coverage"

    if candidate.source == "pdb":
        if candidate.resolution is None:
            return "low_confidence"
        if candidate.has_ligand and candidate.resolution <= 3.0:
            quality = min(1.0, 0.55 * candidate.coverage + 0.45 * (1 - candidate.resolution / 10))
            return (0, quality)
        if not candidate.has_ligand and candidate.resolution <= 2.8:
            quality = min(1.0, 0.50 * candidate.coverage + 0.50 * (1 - candidate.resolution / 10))
            return (1, quality)
        return "low_confidence"

    if (
        candidate.confidence is None
        or candidate.confidence < 0.70
        or not candidate.pae_supported
    ):
        return "low_confidence"
    quality = min(1.0, 0.55 * candidate.coverage + 0.45 * candidate.confidence)
    return (2, quality)


def choose_structure(candidates: list[StructureCandidate]) -> StructureSelection:
    """Choose one usable structure without dropping uncomputable targets."""
    if not candidates:
        return StructureSelection(status="unsupported", unsupported_reason="no_structure")

    eligible: list[tuple[int, float, StructureCandidate]] = []
    rejected_reasons: list[UnsupportedReason] = []
    for candidate in candidates:
        result = _eligibility(candidate)
        if isinstance(result, str):
            rejected_reasons.append(result)
        else:
            priority, quality = result
            eligible.append((priority, quality, candidate))

    if eligible:
        priority, quality, selected = sorted(
            eligible,
            key=lambda item: (item[0], -item[1], item[2].structure_id),
        )[0]
        del priority
        return StructureSelection(
            status="ready",
            structure_id=selected.structure_id,
            source=selected.source,
            score=quality,
        )

    reason_order: tuple[UnsupportedReason, ...] = (
        "unsupported_chemistry",
        "sequence_mismatch",
        "low_coverage",
        "low_confidence",
        "no_structure",
    )
    reason = next(item for item in reason_order if item in rejected_reasons)
    return StructureSelection(status="unsupported", unsupported_reason=reason)
