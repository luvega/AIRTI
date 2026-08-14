"""Blinded, tier-aware evaluation for retrospective target-fishing cases."""

from __future__ import annotations

from pathlib import Path
from typing import Literal, cast

import yaml
from pydantic import BaseModel, ConfigDict, Field, model_validator

from airti_tf.manifest_io import read_jsonl
from airti_tf.ranking.consensus import TargetEvidence, rank_targets

AnchorTier = Literal["gold", "silver"]
EvaluationStage = Literal["screen", "boltz", "final"]
MDStatus = Literal["stable", "unstable", "failed"]


class PositiveAnchor(BaseModel):
    model_config = ConfigDict(extra="forbid")

    target_id: str = Field(min_length=1)
    gene_symbol: str = Field(min_length=1)
    tier: AnchorTier


class CaseDefinition(BaseModel):
    model_config = ConfigDict(extra="forbid")

    schema_version: Literal["1.0"] = "1.0"
    case_id: str = Field(min_length=1)
    query_id: str = Field(min_length=1)
    query_inchikey: str = Field(pattern=r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
    positive_anchors: list[PositiveAnchor] = Field(min_length=1)
    novelty_exclusions: list[str] = Field(default_factory=list)
    upstream_visible_target_ids: list[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def prevent_anchor_leakage(self) -> "CaseDefinition":
        anchors = {anchor.target_id for anchor in self.positive_anchors}
        leaked = anchors.intersection(self.upstream_visible_target_ids)
        if leaked:
            raise ValueError(f"anchor leakage into upstream stages: {sorted(leaked)}")
        if len(anchors) != len(self.positive_anchors):
            raise ValueError("positive anchor target IDs must be unique")
        if not anchors.issubset(self.novelty_exclusions):
            raise ValueError("all positive anchors must be novelty exclusions")
        return self


class TierMetrics(BaseModel):
    model_config = ConfigDict(extra="forbid")

    anchor_count: int = Field(ge=0)
    recall_at: dict[int, float]
    mean_reciprocal_rank: float = Field(ge=0, le=1)


class StageEvaluation(BaseModel):
    model_config = ConfigDict(extra="forbid")

    anchor_ranks: dict[str, int | None]
    tiers: dict[AnchorTier, TierMetrics]


class CaseEvaluation(BaseModel):
    model_config = ConfigDict(extra="forbid")

    schema_version: Literal["1.0"] = "1.0"
    case_id: str
    query_id: str
    stages: dict[EvaluationStage, StageEvaluation]
    exploratory_targets: list[str]
    interpretation: Literal["computational_candidates_only"] = (
        "computational_candidates_only"
    )


def load_case_definition(path: Path) -> CaseDefinition:
    """Load a strict case file and reject target labels visible upstream."""
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("case definition must be a YAML object")
    return CaseDefinition.model_validate(payload)


def _as_float(value: object) -> float:
    if not isinstance(value, (str, bytes, int, float)):
        raise ValueError(f"expected numeric manifest value, received {value!r}")
    return float(value)


def _as_int(value: object) -> int:
    if not isinstance(value, (str, bytes, int, float)):
        raise ValueError(f"expected integer manifest value, received {value!r}")
    return int(value)


def _evidence(row: dict[str, object], *, include_md: bool) -> TargetEvidence:
    return TargetEvidence(
        target_id=str(row["target_id"]),
        ligand_id=str(row.get("ligand_id", "case-query")),
        status="ready",
        vina_score=_as_float(row["calibrated_score"]),
        docking_consistency=_as_float(row.get("pose_consistency", 0.0)),
        structure_quality=_as_float(row["structure_quality"]),
        boltz_score=(
            _as_float(row["boltz_score"])
            if row.get("boltz_score") is not None
            else None
        ),
        md_score=(
            _as_float(row["md_score"])
            if include_md and row.get("md_score") is not None
            else None
        ),
        md_status=(
            cast(MDStatus, str(row["md_status"])) if include_md else None
        ),
        successful_seeds=_as_int(row.get("boltz_seed_success_count", 0)),
        boltz_seed_spread=(
            _as_float(row["boltz_confidence_range"])
            if row.get("boltz_confidence_range") is not None
            else None
        ),
        severe_clash=bool(row.get("severe_clash", False)),
        heavy_atom_count=_as_int(row.get("ligand_atom_count", 0)),
    )


def _rank_map(
    rows: list[dict[str, object]], *, stage: EvaluationStage, query_id: str
) -> tuple[dict[str, int], list[str]]:
    selected = [row for row in rows if row.get("query_id") == query_id]
    if stage == "screen":
        ranked = sorted(
            selected,
            key=lambda row: (_as_int(row["screen_rank"]), str(row["target_id"])),
        )
        return (
            {
                str(row["target_id"]): _as_int(row["screen_rank"])
                for row in ranked
            },
            [str(row["target_id"]) for row in ranked],
        )
    if stage == "boltz":
        selected = [row for row in selected if row.get("boltz_status") == "succeeded"]
        result = rank_targets(
            [_evidence(row, include_md=False) for row in selected], stage="boltz"
        )
    else:
        selected = [
            row
            for row in selected
            if row.get("md_status") in {"stable", "unstable"}
        ]
        result = rank_targets(
            [_evidence(row, include_md=True) for row in selected], stage="final"
        )
    ordered = [item.target_id for item in result.ranked]
    return ({target_id: rank for rank, target_id in enumerate(ordered, 1)}, ordered)


def _stage_evaluation(
    case: CaseDefinition, ranks: dict[str, int]
) -> StageEvaluation:
    anchor_ranks = {
        anchor.target_id: ranks.get(anchor.target_id)
        for anchor in case.positive_anchors
    }
    tiers: dict[AnchorTier, TierMetrics] = {}
    for tier in ("gold", "silver"):
        anchors = [anchor for anchor in case.positive_anchors if anchor.tier == tier]
        observed = [anchor_ranks[anchor.target_id] for anchor in anchors]
        recall_at = {
            cutoff: (
                sum(rank is not None and rank <= cutoff for rank in observed)
                / len(anchors)
                if anchors
                else 0.0
            )
            for cutoff in (10, 50, 100)
        }
        reciprocal = sum(1 / rank for rank in observed if rank is not None)
        tiers[tier] = TierMetrics(
            anchor_count=len(anchors),
            recall_at=recall_at,
            mean_reciprocal_rank=reciprocal / len(anchors) if anchors else 0.0,
        )
    return StageEvaluation(anchor_ranks=anchor_ranks, tiers=tiers)


def evaluate_case_manifests(
    *,
    case_path: Path,
    screen_manifest: Path,
    boltz_manifest: Path,
    md_manifest: Path,
) -> CaseEvaluation:
    """Evaluate frozen stage outputs without exposing anchors to upstream stages."""
    case = load_case_definition(case_path)
    stage_rows = {
        "screen": read_jsonl(screen_manifest),
        "boltz": read_jsonl(boltz_manifest),
        "final": read_jsonl(md_manifest),
    }
    evaluations: dict[EvaluationStage, StageEvaluation] = {}
    final_rows: list[dict[str, object]] = []
    final_order: list[str] = []
    for stage in ("screen", "boltz", "final"):
        ranks, order = _rank_map(
            stage_rows[stage], stage=stage, query_id=case.query_id
        )
        evaluations[stage] = _stage_evaluation(case, ranks)
        if stage == "final":
            final_rows = stage_rows[stage]
            final_order = order
    by_target = {
        str(row["target_id"]): row
        for row in final_rows
        if row.get("query_id") == case.query_id
    }
    excluded = set(case.novelty_exclusions)
    exploratory = [
        target_id
        for target_id in final_order
        if target_id not in excluded
        and by_target[target_id].get("md_status") == "stable"
    ][:5]
    return CaseEvaluation(
        case_id=case.case_id,
        query_id=case.query_id,
        stages=evaluations,
        exploratory_targets=exploratory,
    )
