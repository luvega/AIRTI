"""Per-pocket empirical calibration and family-diverse target routing."""

from __future__ import annotations

import math
from collections import Counter, defaultdict

from pydantic import BaseModel, ConfigDict, Field, model_validator


class InsufficientBackgroundError(RuntimeError):
    """Raised when a pocket lacks the required valid probe distribution."""


class BackgroundDistribution(BaseModel):
    model_config = ConfigDict(extra="forbid")

    pocket_id: str
    affinities: list[float]

    @model_validator(mode="after")
    def require_valid_distribution(self) -> "BackgroundDistribution":
        valid = [value for value in self.affinities if math.isfinite(value)]
        if len(valid) < 95:
            raise InsufficientBackgroundError(
                f"background distribution has {len(valid)} valid values; 95 required"
            )
        self.affinities = valid
        return self


class ScreenHit(BaseModel):
    model_config = ConfigDict(extra="forbid")

    target_id: str
    family: str
    pocket_id: str
    ligand_state_id: str
    affinity_median: float
    calibrated_score: float = Field(ge=0, le=1)
    seed_range: float = Field(ge=0)
    pose_consistency: float = Field(ge=0, le=1)


class RoutedCandidate(BaseModel):
    model_config = ConfigDict(extra="forbid")

    target_id: str
    family: str
    best_pocket_id: str
    best_state_id: str
    calibrated_score: float = Field(ge=0, le=1)
    second_best_score: float | None = Field(default=None, ge=0, le=1)
    affinity_median: float
    seed_range: float
    pose_consistency: float
    selection_reason: str


def empirical_percentile(*, query: float, background: list[float]) -> float:
    """Return a smoothed percentile where more-negative affinity ranks higher."""
    valid = [value for value in background if math.isfinite(value)]
    if not valid:
        raise ValueError("background distribution cannot be empty")
    return (sum(value >= query for value in valid) + 1) / (len(valid) + 1)


def _aggregate_hits(records: list[ScreenHit]) -> list[RoutedCandidate]:
    by_target: dict[str, list[ScreenHit]] = defaultdict(list)
    for record in records:
        by_target[record.target_id].append(record)

    aggregated: list[RoutedCandidate] = []
    for target_id, hits in by_target.items():
        ordered = sorted(
            hits,
            key=lambda hit: (
                -hit.calibrated_score,
                hit.seed_range,
                -hit.pose_consistency,
                hit.pocket_id,
                hit.ligand_state_id,
            ),
        )
        best = ordered[0]
        aggregated.append(
            RoutedCandidate(
                target_id=target_id,
                family=best.family,
                best_pocket_id=best.pocket_id,
                best_state_id=best.ligand_state_id,
                calibrated_score=best.calibrated_score,
                second_best_score=(
                    ordered[1].calibrated_score if len(ordered) > 1 else None
                ),
                affinity_median=best.affinity_median,
                seed_range=best.seed_range,
                pose_consistency=best.pose_consistency,
                selection_reason="unassigned",
            )
        )
    return sorted(
        aggregated,
        key=lambda item: (
            -item.calibrated_score,
            item.seed_range,
            -item.pose_consistency,
            item.target_id,
        ),
    )


def route_screen_candidates(
    records: list[ScreenHit],
    *,
    top_n: int,
    primary_n: int = 240,
    family_cap: int = 15,
) -> list[RoutedCandidate]:
    """Aggregate pockets/states and apply deterministic family diversification."""
    if top_n <= 0 or primary_n <= 0 or family_cap <= 0:
        raise ValueError("routing limits must be positive")
    ordered = _aggregate_hits(records)
    target_primary = min(primary_n, top_n)
    selected: list[RoutedCandidate] = []
    selected_ids: set[str] = set()
    family_counts: Counter[str] = Counter()

    def add(candidate: RoutedCandidate, reason: str) -> None:
        selected.append(candidate.model_copy(update={"selection_reason": reason}))
        selected_ids.add(candidate.target_id)
        family_counts[candidate.family] += 1

    for candidate in ordered:
        if len(selected) >= target_primary:
            break
        if family_counts[candidate.family] < family_cap:
            add(candidate, "global_primary")

    covered_families = set(family_counts)
    for candidate in ordered:
        if len(selected) >= top_n:
            break
        if candidate.target_id in selected_ids or candidate.family in covered_families:
            continue
        add(candidate, "family_diversity")
        covered_families.add(candidate.family)

    for candidate in ordered:
        if len(selected) >= top_n:
            break
        if (
            candidate.target_id not in selected_ids
            and family_counts[candidate.family] < family_cap
        ):
            add(candidate, "score_fill")

    return selected

