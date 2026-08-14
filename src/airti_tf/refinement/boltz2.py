"""Boltz-2 YAML, execution, output parsing, and multi-seed consensus."""

from __future__ import annotations

import json
import math
import re
import statistics
import subprocess
from collections.abc import Callable
from pathlib import Path
from typing import Any, Literal

from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from pydantic import BaseModel, ConfigDict, Field, model_validator

from airti_tf.manifest_io import write_artifact, write_bytes_atomic


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
    cofactors: list[str] = Field(default_factory=list)
    pocket_residues: list[int] = Field(min_length=1)
    input_yaml: Path
    output_dir: Path
    cache_path: Path | None = None


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


class BoltzStructureQuality(BaseModel):
    model_config = ConfigDict(extra="forbid")

    pocket_constraint_fraction: float = Field(ge=0, le=1)
    severe_clash: bool
    minimum_interatomic_distance_a: float = Field(ge=0)
    evaluated_pocket_residue_count: int = Field(gt=0)


class BoltzExecutionResult(BaseModel):
    model_config = ConfigDict(extra="forbid")

    return_code: int
    stdout: str
    stderr: str
    timed_out: bool


BoltzExecutor = Callable[[list[str], int], BoltzExecutionResult]


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
    cofactor_sequences = [
        {"ligand": {"id": chr(ord("C") + index), "ccd": ccd_id}}
        for index, ccd_id in enumerate(job.cofactors)
    ]
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
            *cofactor_sequences,
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
    command = [
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
        "--num_workers",
        "0",
        "--seed",
        str(seed),
        "--use_potentials",
    ]
    if job.cache_path is not None:
        command.extend(["--cache", str(job.cache_path)])
    return command


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
    if (
        "no locator available" in lowered
        or "no writable cache directories" in lowered
        or "permission denied" in lowered and "cache" in lowered
    ):
        return "runtime_environment"
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
    if re.search(
        r"(?:msa|a3m).{0,80}(?:not found|no such file|does not exist|missing)",
        lowered,
        flags=re.DOTALL,
    ):
        return "missing_msa"
    return "nonzero_exit"


def _to_text(value: str | bytes | None) -> str:
    if value is None:
        return ""
    return value.decode(errors="replace") if isinstance(value, bytes) else value


def _subprocess_executor(
    command: list[str], timeout_seconds: int
) -> BoltzExecutionResult:
    try:
        process = subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
        )
    except subprocess.TimeoutExpired as error:
        return BoltzExecutionResult(
            return_code=124,
            stdout=_to_text(error.stdout),
            stderr=_to_text(error.stderr),
            timed_out=True,
        )
    return BoltzExecutionResult(
        return_code=process.returncode,
        stdout=process.stdout,
        stderr=process.stderr,
        timed_out=False,
    )


def run_boltz_seed(
    job: BoltzJob,
    *,
    seed: int,
    executor: BoltzExecutor = _subprocess_executor,
    timeout_seconds: int = 7_200,
) -> BoltzSeedResult:
    """Execute one Boltz seed, generate structural QC, and parse its outputs."""
    job.output_dir.mkdir(parents=True, exist_ok=True)
    command = build_boltz_command(job, seed=seed)
    result = executor(command, timeout_seconds)
    stdout_sha256 = write_bytes_atomic(
        job.output_dir / "boltz.stdout.log", result.stdout.encode("utf-8")
    )
    stderr_sha256 = write_bytes_atomic(
        job.output_dir / "boltz.stderr.log", result.stderr.encode("utf-8")
    )
    write_artifact(
        job.output_dir / "boltz.execution.json",
        {
            "command": command,
            "return_code": result.return_code,
            "seed": seed,
            "stderr_sha256": stderr_sha256,
            "stdout_sha256": stdout_sha256,
            "timed_out": result.timed_out,
        },
    )
    if result.timed_out:
        return BoltzSeedResult(seed=seed, status="failed", error_code="timeout")
    if result.return_code != 0:
        return BoltzSeedResult(
            seed=seed,
            status="failed",
            error_code=classify_boltz_failure(result.stderr or result.stdout),
        )
    input_stem = job.input_yaml.stem
    result_root = job.output_dir / f"boltz_results_{input_stem}"
    if not result_root.is_dir() and (job.output_dir / "predictions").is_dir():
        result_root = job.output_dir
    prediction_dir = result_root / "predictions" / input_stem
    structure_path = prediction_dir / f"{input_stem}_model_0.cif"
    if not structure_path.is_file():
        return BoltzSeedResult(seed=seed, status="failed", error_code="output_missing")
    quality_path = prediction_dir / "airti_quality.json"
    if not quality_path.is_file():
        try:
            write_boltz_quality(
                structure_path,
                pocket_residues=job.pocket_residues,
            )
        except (OSError, TypeError, ValueError):
            return BoltzSeedResult(
                seed=seed,
                status="failed",
                error_code="structural_qc_failed",
            )
    return parse_boltz_output(result_root, input_stem=input_stem, seed=seed)


def _as_list(value: object) -> list[str]:
    if isinstance(value, list):
        return [str(item) for item in value]
    return [str(value)]


def assess_boltz_structure(
    structure_path: Path,
    *,
    pocket_residues: list[int],
    contact_distance_a: float = 6.0,
    clash_distance_a: float = 1.2,
) -> BoltzStructureQuality:
    """Measure predicted pocket retention and protein-ligand atomic clashes."""
    if not pocket_residues:
        raise ValueError("at least one pocket residue is required for Boltz QC")
    payload = MMCIF2Dict(str(structure_path))
    required = (
        "_atom_site.auth_asym_id",
        "_atom_site.auth_seq_id",
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z",
    )
    missing = [key for key in required if key not in payload]
    if missing:
        raise ValueError(f"Boltz mmCIF lacks atom-site fields: {missing}")
    chains = _as_list(payload["_atom_site.auth_asym_id"])
    residues = _as_list(payload["_atom_site.auth_seq_id"])
    xs = _as_list(payload["_atom_site.Cartn_x"])
    ys = _as_list(payload["_atom_site.Cartn_y"])
    zs = _as_list(payload["_atom_site.Cartn_z"])
    lengths = {len(chains), len(residues), len(xs), len(ys), len(zs)}
    if len(lengths) != 1:
        raise ValueError("Boltz mmCIF atom-site columns have inconsistent lengths")

    protein_by_residue: dict[int, list[tuple[float, float, float]]] = {}
    ligand_atoms: list[tuple[float, float, float]] = []
    protein_atoms: list[tuple[float, float, float]] = []
    for chain, raw_residue, raw_x, raw_y, raw_z in zip(
        chains, residues, xs, ys, zs, strict=True
    ):
        coordinate = (float(raw_x), float(raw_y), float(raw_z))
        if chain == "B":
            ligand_atoms.append(coordinate)
            continue
        if chain != "A":
            continue
        protein_atoms.append(coordinate)
        try:
            residue = int(raw_residue)
        except ValueError:
            continue
        protein_by_residue.setdefault(residue, []).append(coordinate)
    if not protein_atoms or not ligand_atoms:
        raise ValueError("Boltz mmCIF must contain protein chain A and ligand chain B")

    def distance(left: tuple[float, float, float], right: tuple[float, float, float]) -> float:
        return math.sqrt(sum((a - b) ** 2 for a, b in zip(left, right, strict=True)))

    minimum = min(distance(protein, ligand) for protein in protein_atoms for ligand in ligand_atoms)
    unique_pocket_residues = sorted(set(pocket_residues))
    contacted = 0
    for residue in unique_pocket_residues:
        atoms = protein_by_residue.get(residue, [])
        if any(
            distance(protein, ligand) <= contact_distance_a
            for protein in atoms
            for ligand in ligand_atoms
        ):
            contacted += 1
    return BoltzStructureQuality(
        pocket_constraint_fraction=contacted / len(unique_pocket_residues),
        severe_clash=minimum < clash_distance_a,
        minimum_interatomic_distance_a=minimum,
        evaluated_pocket_residue_count=len(unique_pocket_residues),
    )


def write_boltz_quality(
    structure_path: Path,
    *,
    pocket_residues: list[int],
) -> Path:
    """Write the structural QC sidecar consumed by ``parse_boltz_output``."""
    quality = assess_boltz_structure(
        structure_path,
        pocket_residues=pocket_residues,
    )
    output = structure_path.parent / "airti_quality.json"
    write_artifact(output, quality.model_dump(mode="json"))
    return output
