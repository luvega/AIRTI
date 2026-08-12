"""Boltz-2 YAML, execution, output parsing, and multi-seed consensus."""

from __future__ import annotations

import json
import math
import re
import statistics
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator


class MissingMSAError(FileNotFoundError):
    """Raised when production inference lacks its pinned MSA."""


class InsufficientBoltzSeedsError(RuntimeError):
    """Raised when fewer than two non-clashing seeds complete."""


class BoltzJob(BaseModel):
    model_config = ConfigDict(extra="forbid")

    job_id: str
    target_id: str
    sequence: str = Field(pattern=r"^[A-Z]+$")
    sequence_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    msa_path: Path
    msa_database_version: str
    ligand_state_id: str
    ligand_smiles: str
    ligand_atom_count: int = Field(gt=0, le=128)
    pocket_residues: list[int] = Field(min_length=1)
    input_yaml: Path
    output_dir: Path


class BoltzSeedResult(BaseModel):
    model_config = ConfigDict(extra="forbid")

    seed: int
    status: Literal["succeeded", "failed"]
    confidence_score: float | None = Field(default=None, ge=0, le=1)
    ligand_iptm: float | None = Field(default=None, ge=0, le=1)
    affinity_probability: float | None = Field(default=None, ge=0, le=1)
    affinity_pred_value: float | None = None
    pocket_constraint_fraction: float | None = Field(default=None, ge=0, le=1)
    severe_clash: bool | None = None
    structure_path: Path | None = None
    error_code: str | None = None

    @model_validator(mode="after")
    def validate_state(self) -> "BoltzSeedResult":
        required = (
            self.confidence_score,
            self.ligand_iptm,
            self.affinity_probability,
            self.affinity_pred_value,
            self.pocket_constraint_fraction,
            self.severe_clash,
        )
        if self.status == "succeeded" and any(value is None for value in required):
            raise ValueError("successful Boltz result lacks a required metric")
        if self.status == "failed" and not self.error_code:
            raise ValueError("failed Boltz result requires error_code")
        return self


class BoltzSummary(BaseModel):
    model_config = ConfigDict(extra="forbid")

    seed_success_count: int = Field(ge=2, le=3)
    confidence_median: float = Field(ge=0, le=1)
    ligand_iptm_median: float = Field(ge=0, le=1)
    affinity_probability_median: float = Field(ge=0, le=1)
    affinity_pred_value_median: float
    pocket_constraint_median: float = Field(ge=0, le=1)
    confidence_range: float = Field(ge=0)
    successful_seeds: list[int]


def msa_cache_path(
    cache_root: Path,
    *,
    uniprot_id: str,
    sequence_sha256: str,
    database_version: str,
) -> Path:
    """Return an MSA cache path tied to accession, sequence, and databases."""
    safe_database = re.sub(r"[^A-Za-z0-9_.-]+", "-", database_version)
    return cache_root / f"{uniprot_id}.{sequence_sha256[:12]}.{safe_database}.a3m"


def build_boltz_yaml(
    job: BoltzJob, *, profile: Literal["local", "production"] = "production"
) -> dict[str, Any]:
    """Build the official Boltz YAML schema for one protein-small molecule pair."""
    if profile == "production" and not job.msa_path.is_file():
        raise MissingMSAError(job.msa_path)
    msa_value = str(job.msa_path) if job.msa_path.is_file() else "empty"
    return {
        "version": 1,
        "sequences": [
            {
                "protein": {
                    "id": "A",
                    "sequence": job.sequence,
                    "msa": msa_value,
                }
            },
            {"ligand": {"id": "B", "smiles": job.ligand_smiles}},
        ],
        "constraints": [
            {
                "pocket": {
                    "binder": "B",
                    "contacts": [["A", residue] for residue in job.pocket_residues],
                    "max_distance": 6.0,
                    "force": True,
                }
            }
        ],
        "properties": [{"affinity": {"binder": "B"}}],
    }


def build_boltz_command(job: BoltzJob, *, seed: int) -> list[str]:
    """Build a deterministic three-sample Boltz-2 prediction command."""
    return [
        "boltz",
        "predict",
        str(job.input_yaml),
        "--out_dir",
        str(job.output_dir),
        "--model",
        "boltz2",
        "--diffusion_samples",
        "3",
        "--diffusion_samples_affinity",
        "3",
        "--max_parallel_samples",
        "1",
        "--seed",
        str(seed),
        "--use_potentials",
    ]


def _load_json(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"JSON output is not an object: {path}")
    return payload


def _finite_float(payload: dict[str, Any], key: str) -> float:
    value = float(payload[key])
    if not math.isfinite(value):
        raise ValueError(f"non-finite Boltz metric: {key}")
    return value


def parse_boltz_output(
    output_dir: Path, *, input_stem: str, seed: int
) -> BoltzSeedResult:
    """Parse official confidence/affinity JSON plus AIRTI structural QC."""
    prediction_dir = output_dir / "predictions" / input_stem
    confidence_path = prediction_dir / f"confidence_{input_stem}_model_0.json"
    affinity_path = prediction_dir / f"affinity_{input_stem}.json"
    quality_path = prediction_dir / "airti_quality.json"
    structure_path = prediction_dir / f"{input_stem}_model_0.cif"
    required_paths = (confidence_path, affinity_path, quality_path, structure_path)
    if not all(path.is_file() for path in required_paths):
        return BoltzSeedResult(seed=seed, status="failed", error_code="output_missing")
    try:
        confidence = _load_json(confidence_path)
        affinity = _load_json(affinity_path)
        quality = _load_json(quality_path)
        result = BoltzSeedResult(
            seed=seed,
            status="succeeded",
            confidence_score=_finite_float(confidence, "confidence_score"),
            ligand_iptm=_finite_float(confidence, "ligand_iptm"),
            affinity_probability=_finite_float(
                affinity, "affinity_probability_binary"
            ),
            affinity_pred_value=_finite_float(affinity, "affinity_pred_value"),
            pocket_constraint_fraction=_finite_float(
                quality, "pocket_constraint_fraction"
            ),
            severe_clash=bool(quality["severe_clash"]),
            structure_path=structure_path,
        )
    except (KeyError, TypeError, ValueError, json.JSONDecodeError):
        return BoltzSeedResult(seed=seed, status="failed", error_code="nan_output")
    if result.pocket_constraint_fraction is not None and result.pocket_constraint_fraction < 0.5:
        return BoltzSeedResult(
            seed=seed, status="failed", error_code="constraint_violation"
        )
    return result


def summarize_boltz_seeds(records: list[BoltzSeedResult]) -> BoltzSummary:
    """Require two non-clashing seeds and return a robust median consensus."""
    successful = [
        record
        for record in records
        if record.status == "succeeded" and record.severe_clash is False
    ]
    if len(successful) < 2:
        raise InsufficientBoltzSeedsError(
            f"only {len(successful)} of {len(records)} non-clashing seeds succeeded"
        )

    def values(name: str) -> list[float]:
        extracted = [getattr(record, name) for record in successful]
        return [float(value) for value in extracted if value is not None]

    confidences = values("confidence_score")
    return BoltzSummary(
        seed_success_count=len(successful),
        confidence_median=statistics.median(confidences),
        ligand_iptm_median=statistics.median(values("ligand_iptm")),
        affinity_probability_median=statistics.median(
            values("affinity_probability")
        ),
        affinity_pred_value_median=statistics.median(values("affinity_pred_value")),
        pocket_constraint_median=statistics.median(
            values("pocket_constraint_fraction")
        ),
        confidence_range=max(confidences) - min(confidences),
        successful_seeds=sorted(record.seed for record in successful),
    )


def classify_boltz_failure(stderr: str) -> str:
    """Map Boltz failures to stable codes used by retry policy."""
    lowered = stderr.lower()
    if "out of memory" in lowered or "cuda oom" in lowered:
        return "cuda_oom"
    if "yaml" in lowered:
        return "invalid_yaml"
    if "129 atoms" in lowered or "too large" in lowered:
        return "ligand_too_large"
    if "nan" in lowered or "non-finite" in lowered:
        return "nan_output"
    if "constraint" in lowered:
        return "constraint_violation"
    if "msa" in lowered:
        return "missing_msa"
    return "nonzero_exit"
