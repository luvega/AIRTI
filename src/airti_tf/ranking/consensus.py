"""Frozen, stage-aware evidence ranking for reverse target fishing."""

from __future__ import annotations

from collections.abc import Iterable
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator

Stage = Literal["screen", "boltz", "final"]


class FrozenWeightError(TypeError):
    """Raised when scientific weights are mutated after freezing."""


class StageWeights:
    """Minimal immutable weight record with an explicit mutation error."""

    __slots__ = ("_frozen", "boltz", "docking_consistency", "md", "structure_quality", "vina")
    _frozen: bool
    vina: float
    docking_consistency: float
    structure_quality: float
    boltz: float
    md: float

    def __init__(
        self,
        *,
        vina: float,
        docking_consistency: float = 0.0,
        structure_quality: float,
        boltz: float = 0.0,
        md: float = 0.0,
    ) -> None:
        object.__setattr__(self, "_frozen", False)
        object.__setattr__(self, "vina", vina)
        object.__setattr__(self, "docking_consistency", docking_consistency)
        object.__setattr__(self, "structure_quality", structure_quality)
        object.__setattr__(self, "boltz", boltz)
        object.__setattr__(self, "md", md)
        object.__setattr__(self, "_frozen", True)

    def __setattr__(self, name: str, value: object) -> None:
        if getattr(self, "_frozen", False):
            raise FrozenWeightError(f"ranking weight {name} is frozen")
        object.__setattr__(self, name, value)

    def as_dict(self) -> dict[str, float]:
        return {
            "vina": self.vina,
            "boltz": self.boltz,
            "md": self.md,
            "docking_consistency": self.docking_consistency,
            "structure_quality": self.structure_quality,
        }


class FrozenWeights:
    __slots__ = ("boltz", "final", "risk_penalty_max", "screen", "version")
    screen: StageWeights
    boltz: StageWeights
    final: StageWeights
    risk_penalty_max: float
    version: str

    def __init__(self) -> None:
        object.__setattr__(
            self,
            "screen",
            StageWeights(
                vina=0.65, docking_consistency=0.20, structure_quality=0.15
            ),
        )
        object.__setattr__(
            self,
            "boltz",
            StageWeights(
                vina=0.35,
                boltz=0.40,
                docking_consistency=0.10,
                structure_quality=0.15,
            ),
        )
        object.__setattr__(
            self,
            "final",
            StageWeights(vina=0.25, boltz=0.30, md=0.30, structure_quality=0.15),
        )
        object.__setattr__(self, "risk_penalty_max", 0.15)
        object.__setattr__(self, "version", "retrospective-v1")

    def __setattr__(self, name: str, value: object) -> None:
        raise FrozenWeightError(f"ranking configuration {name} is frozen")


FROZEN_WEIGHTS = FrozenWeights()


class TargetEvidence(BaseModel):
    """Target-level evidence after pocket and ligand-state aggregation."""

    model_config = ConfigDict(extra="forbid")

    target_id: str
    ligand_id: str
    status: Literal["ready", "unsupported", "failed"]
    vina_score: float | None = None
    docking_consistency: float | None = Field(default=None, ge=0, le=1)
    structure_quality: float | None = Field(default=None, ge=0, le=1)
    boltz_score: float | None = None
    md_score: float | None = Field(default=None, ge=0, le=1)
    md_status: Literal["stable", "unstable", "failed"] | None = None
    successful_seeds: int = Field(default=0, ge=0)
    structure_low_confidence: bool = False
    boltz_seed_spread: float | None = Field(default=None, ge=0)
    severe_clash: bool = False
    heavy_atom_count: int = Field(default=0, ge=0)

    @model_validator(mode="after")
    def preserve_unsupported_missingness(self) -> "TargetEvidence":
        if self.status == "unsupported" and any(
            value is not None
            for value in (
                self.vina_score,
                self.docking_consistency,
                self.structure_quality,
                self.boltz_score,
                self.md_score,
            )
        ):
            raise ValueError("unsupported targets cannot carry numeric evidence")
        return self


class Coverage(BaseModel):
    model_config = ConfigDict(extra="forbid")

    total: int
    ready: int
    unsupported: int
    failed: int


class RankedEvidence(BaseModel):
    model_config = ConfigDict(extra="forbid")

    target_id: str
    ligand_id: str
    rank: int = Field(gt=0)
    priority: float = Field(ge=0, le=1)
    evidence_tier: Literal["screen_only", "two_stage", "full_three_stage", "partial_computational"]
    raw_scores: dict[str, float | None]
    normalized_scores: dict[str, float | None]
    risk_penalty: float = Field(ge=0, le=0.15)
    uncertainty_flags: list[str]
    successful_seeds: int
    rank_reason: str
    weight_version: str
    vina_score: float | None
    boltz_score: float | None
    md_score: float | None
    structure_quality: float | None
    source: TargetEvidence


class RankingResult(BaseModel):
    model_config = ConfigDict(extra="forbid")

    stage: Stage
    ranked: list[RankedEvidence]
    full_evidence: list[RankedEvidence]
    partial_evidence: list[RankedEvidence]
    coverage: Coverage
    weight_version: str


def _percentiles(values: Iterable[float | None]) -> list[float | None]:
    material = list(values)
    observed = [value for value in material if value is not None]
    if not observed:
        return [None] * len(material)
    # Ties receive the same empirical CDF percentile. Higher evidence is better.
    return [
        None
        if value is None
        else sum(candidate <= value for candidate in observed) / len(observed)
        for value in material
    ]


def _risk(evidence: TargetEvidence, stage: Stage) -> tuple[float, list[str]]:
    flags: list[str] = []
    if evidence.structure_low_confidence:
        flags.append("structure_low_confidence")
    if evidence.boltz_seed_spread is not None and evidence.boltz_seed_spread > 0.25:
        flags.append("boltz_seed_disagreement")
    if evidence.severe_clash:
        flags.append("severe_clash")
    if evidence.heavy_atom_count > 56:
        flags.append("large_ligand")
    if stage in {"boltz", "final"} and evidence.boltz_score is None:
        flags.append("missing_boltz")
    if stage == "final" and evidence.md_score is None:
        flags.append("missing_md")
    # Every declared uncertainty has equal, transparent influence; cap globally.
    return min(FROZEN_WEIGHTS.risk_penalty_max, 0.03 * len(flags)), flags


def _tier(evidence: TargetEvidence, stage: Stage) -> str:
    if stage == "screen":
        return "screen_only"
    if stage == "boltz" and evidence.boltz_score is not None:
        return "two_stage"
    if (
        stage == "final"
        and evidence.boltz_score is not None
        and evidence.md_score is not None
        and evidence.md_status in {"stable", "unstable"}
    ):
        return "full_three_stage"
    return "partial_computational"


def rank_targets(records: list[TargetEvidence], *, stage: Stage) -> RankingResult:
    """Rank ready targets while preserving unsupported/failed coverage counts."""
    coverage = Coverage(
        total=len(records),
        ready=sum(record.status == "ready" for record in records),
        unsupported=sum(record.status == "unsupported" for record in records),
        failed=sum(record.status == "failed" for record in records),
    )
    ready = [record for record in records if record.status == "ready"]
    fields = ["vina", "boltz", "md", "docking_consistency", "structure_quality"]
    raw_by_field: dict[str, list[float | None]] = {
        "vina": [record.vina_score for record in ready],
        "boltz": [record.boltz_score for record in ready],
        "md": [record.md_score for record in ready],
        "docking_consistency": [record.docking_consistency for record in ready],
        "structure_quality": [record.structure_quality for record in ready],
    }
    normalized_by_field = {
        field: _percentiles(raw_by_field[field]) for field in fields
    }
    weights = getattr(FROZEN_WEIGHTS, stage).as_dict()

    provisional: list[RankedEvidence] = []
    for index, evidence in enumerate(ready):
        raw = {field: raw_by_field[field][index] for field in fields}
        normalized = {field: normalized_by_field[field][index] for field in fields}
        available = {
            field: value
            for field, value in normalized.items()
            if value is not None and weights[field] > 0
        }
        denominator = sum(weights[field] for field in available)
        base = (
            sum(weights[field] * value for field, value in available.items())
            / denominator
            if denominator
            else 0.0
        )
        penalty, flags = _risk(evidence, stage)
        tier = _tier(evidence, stage)
        provisional.append(
            RankedEvidence(
                target_id=evidence.target_id,
                ligand_id=evidence.ligand_id,
                rank=1,
                priority=max(0.0, min(1.0, round(base - penalty, 6))),
                evidence_tier=tier,  # type: ignore[arg-type]
                raw_scores=raw,
                normalized_scores=normalized,
                risk_penalty=penalty,
                uncertainty_flags=flags,
                successful_seeds=evidence.successful_seeds,
                rank_reason=(
                    f"stage={stage}; tier={tier}; available="
                    + ",".join(sorted(available))
                    + (f"; risks={','.join(flags)}" if flags else "; risks=none")
                ),
                weight_version=FROZEN_WEIGHTS.version,
                vina_score=evidence.vina_score,
                boltz_score=evidence.boltz_score,
                md_score=evidence.md_score,
                structure_quality=evidence.structure_quality,
                source=evidence,
            )
        )

    tier_order = {
        "full_three_stage": 3,
        "two_stage": 2,
        "screen_only": 1,
        "partial_computational": 0,
    }
    ordered = sorted(
        provisional,
        key=lambda item: (
            -tier_order[item.evidence_tier],
            -item.priority,
            -item.successful_seeds,
            -(item.structure_quality if item.structure_quality is not None else -1.0),
            item.target_id,
        ),
    )
    ranked = [item.model_copy(update={"rank": rank}) for rank, item in enumerate(ordered, 1)]
    full = [item for item in ranked if item.evidence_tier == "full_three_stage"]
    partial = [item for item in ranked if item.evidence_tier != "full_three_stage"]
    return RankingResult(
        stage=stage,
        ranked=ranked,
        full_evidence=full,
        partial_evidence=partial,
        coverage=coverage,
        weight_version=FROZEN_WEIGHTS.version,
    )
