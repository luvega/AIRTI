"""Leakage-resistant retrieval benchmarks and immutable release decisions."""

from __future__ import annotations

import json
import random
from pathlib import Path
from typing import Literal

import yaml
from pydantic import BaseModel, ConfigDict, Field

from airti_tf.manifest_io import content_sha256, write_bytes_atomic

DatasetRole = Literal[
    "smoke", "retrospective_train", "retrospective_validation", "blind"
]


class DataLeakageError(ValueError):
    """Raised when a protected benchmark split influences model selection."""


class ConfidenceInterval(BaseModel):
    model_config = ConfigDict(extra="forbid")

    estimate: float = Field(ge=0, le=1)
    lower: float = Field(ge=0, le=1)
    upper: float = Field(ge=0, le=1)
    iterations: int = Field(gt=0)
    seed: int


class FrozenWeightArtifact(BaseModel):
    model_config = ConfigDict(extra="forbid")

    path: Path
    sha256: str = Field(pattern=r"^[0-9a-f]{64}$")


class ReleaseMetrics(BaseModel):
    model_config = ConfigDict(extra="forbid")

    blind_success_at_100: float = Field(ge=0, le=1)
    technical_success_rate: float = Field(ge=0, le=1)
    boltz_success_at_k: float = Field(ge=0, le=1)
    vina_success_at_k: float = Field(ge=0, le=1)
    successful_target_families: int = Field(ge=0)
    median_top20_jaccard: float = Field(ge=0, le=1)
    failures_with_error_codes: bool
    report_integrity_passed: bool


class ReleaseCriterion(BaseModel):
    model_config = ConfigDict(extra="forbid")

    criterion: str
    passed: bool
    observed: float | int | bool
    requirement: str


class ReleaseDecision(BaseModel):
    model_config = ConfigDict(extra="forbid")

    schema_version: Literal["1.0"] = "1.0"
    status: Literal["pass", "fail"]
    dataset_role: DatasetRole
    criteria: list[ReleaseCriterion]
    frozen_weights_sha256: str | None = None


def _validate_inputs(
    truth: dict[str, set[str]], ranked: dict[str, list[str]], k: int
) -> None:
    if k <= 0:
        raise ValueError("k must be positive")
    if not truth:
        raise ValueError("truth cannot be empty")
    if any(not targets for targets in truth.values()):
        raise ValueError("truth sets cannot be empty")
    extra = set(ranked) - set(truth)
    if extra:
        raise ValueError(f"rankings contain ligands absent from truth: {sorted(extra)}")


def success_at_k(
    truth: dict[str, set[str]], ranked: dict[str, list[str]], *, k: int
) -> float:
    """Fraction of ligands with any eligible direct human target in Top-k."""
    _validate_inputs(truth, ranked, k)
    hits = [
        bool(truth[ligand_id].intersection(ranked.get(ligand_id, [])[:k]))
        for ligand_id in truth
    ]
    return sum(hits) / len(hits)


def recall_at_k(
    truth: dict[str, set[str]], ranked: dict[str, list[str]], *, k: int
) -> float:
    """Macro-average target-level recall at k for multi-target compounds."""
    _validate_inputs(truth, ranked, k)
    recalls = [
        len(truth[ligand_id].intersection(ranked.get(ligand_id, [])[:k]))
        / len(truth[ligand_id])
        for ligand_id in truth
    ]
    return sum(recalls) / len(recalls)


def reciprocal_rank(
    truth: dict[str, set[str]], ranked: dict[str, list[str]]
) -> float:
    """Mean reciprocal rank of the first eligible direct human target."""
    _validate_inputs(truth, ranked, k=1)
    values: list[float] = []
    for ligand_id, eligible in truth.items():
        first = next(
            (
                rank
                for rank, target_id in enumerate(ranked.get(ligand_id, []), start=1)
                if target_id in eligible
            ),
            None,
        )
        values.append(0.0 if first is None else 1.0 / first)
    return sum(values) / len(values)


def top_k_jaccard(first: list[str], second: list[str], *, k: int) -> float:
    """Set stability of two ranked lists, restricted to Top-k."""
    if k <= 0:
        raise ValueError("k must be positive")
    first_set = set(first[:k])
    second_set = set(second[:k])
    union = first_set | second_set
    return 1.0 if not union else len(first_set & second_set) / len(union)


def bootstrap_success_at_k(
    truth: dict[str, set[str]],
    ranked: dict[str, list[str]],
    *,
    k: int,
    iterations: int = 10_000,
    seed: int = 20260812,
) -> ConfidenceInterval:
    """Ligand-level percentile bootstrap with a recorded random seed."""
    _validate_inputs(truth, ranked, k)
    if iterations <= 0:
        raise ValueError("iterations must be positive")
    ligand_ids = sorted(truth)
    hit_by_ligand = {
        ligand_id: float(
            bool(truth[ligand_id].intersection(ranked.get(ligand_id, [])[:k]))
        )
        for ligand_id in ligand_ids
    }
    generator = random.Random(seed)
    bootstrapped = sorted(
        sum(hit_by_ligand[generator.choice(ligand_ids)] for _ in ligand_ids)
        / len(ligand_ids)
        for _ in range(iterations)
    )
    lower_index = int(0.025 * (iterations - 1))
    upper_index = int(0.975 * (iterations - 1))
    return ConfidenceInterval(
        estimate=sum(hit_by_ligand.values()) / len(ligand_ids),
        lower=bootstrapped[lower_index],
        upper=bootstrapped[upper_index],
        iterations=iterations,
        seed=seed,
    )


def fit_weights(*, dataset_role: DatasetRole) -> dict[str, float]:
    """Guard the only split authorized to fit scientific ranking weights."""
    if dataset_role != "retrospective_train":
        raise DataLeakageError(
            f"weights may only be fit on retrospective_train, not {dataset_role}"
        )
    return {
        "vina": 0.25,
        "boltz": 0.30,
        "md": 0.30,
        "structure_quality": 0.15,
    }


def freeze_weights(weights: dict[str, object], path: Path) -> FrozenWeightArtifact:
    """Write canonical frozen weights and return the content identity."""
    raw = yaml.safe_dump(
        weights, allow_unicode=True, sort_keys=True, default_flow_style=False
    ).encode("utf-8")
    digest = write_bytes_atomic(path, raw)
    return FrozenWeightArtifact(path=path, sha256=digest)


def _criterion(
    name: str, observed: float | int | bool, passed: bool, requirement: str
) -> ReleaseCriterion:
    return ReleaseCriterion(
        criterion=name,
        observed=observed,
        passed=passed,
        requirement=requirement,
    )


def evaluate_release(
    metrics: ReleaseMetrics,
    *,
    dataset_role: DatasetRole = "retrospective_validation",
    frozen_weights_path: Path | None = None,
) -> ReleaseDecision:
    """Evaluate every pre-registered production release criterion."""
    frozen_digest: str | None = None
    if dataset_role == "blind":
        if frozen_weights_path is None or not frozen_weights_path.is_file():
            raise DataLeakageError("blind evaluation requires an existing frozen weights file")
        frozen_digest = content_sha256(frozen_weights_path.read_bytes())

    criteria = [
        _criterion(
            "blind_success_at_100",
            metrics.blind_success_at_100,
            metrics.blind_success_at_100 >= 0.30,
            ">= 0.30",
        ),
        _criterion(
            "technical_success_rate",
            metrics.technical_success_rate,
            metrics.technical_success_rate >= 0.95,
            ">= 0.95",
        ),
        _criterion(
            "boltz_not_worse_than_vina",
            metrics.boltz_success_at_k - metrics.vina_success_at_k,
            metrics.boltz_success_at_k >= metrics.vina_success_at_k,
            "Boltz-2 Success@k >= QuickVina2 Success@k",
        ),
        _criterion(
            "target_family_coverage",
            metrics.successful_target_families,
            metrics.successful_target_families >= 3,
            ">= 3 families",
        ),
        _criterion(
            "top20_stability",
            metrics.median_top20_jaccard,
            metrics.median_top20_jaccard >= 0.70,
            "median Top20 Jaccard >= 0.70",
        ),
        _criterion(
            "failure_auditability",
            metrics.failures_with_error_codes,
            metrics.failures_with_error_codes,
            "all failures have error_code and rerunnable task_id",
        ),
        _criterion(
            "report_integrity",
            metrics.report_integrity_passed,
            metrics.report_integrity_passed,
            "all report integrity checks passed",
        ),
    ]
    return ReleaseDecision(
        status="pass" if all(item.passed for item in criteria) else "fail",
        dataset_role=dataset_role,
        criteria=criteria,
        frozen_weights_sha256=frozen_digest,
    )


def write_release_decision(decision: ReleaseDecision, path: Path) -> Path:
    """Write a checksum-bearing decision artifact atomically."""
    payload = decision.model_dump(mode="json")
    payload["decision_sha256"] = content_sha256(payload)
    raw = (json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n").encode(
        "utf-8"
    )
    write_bytes_atomic(path, raw)
    return path

