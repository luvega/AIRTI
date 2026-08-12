"""Trajectory completeness checks and conservative MD evidence scoring."""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, ConfigDict, Field


class TrajectoryMetrics(BaseModel):
    model_config = ConfigDict(extra="forbid")

    frames: int = Field(ge=0)
    expected_frames: int = Field(gt=0)
    completed_ns: float = Field(ge=0)
    energy_continuous: bool
    pbc_processed: bool
    protein_ca_rmsd_nm: float = Field(ge=0)
    ligand_pocket_rmsd_nm: float = Field(ge=0)
    pocket_rmsf_nm: float = Field(ge=0)
    hbond_occupancy: float = Field(ge=0, le=1)
    hydrophobic_occupancy: float = Field(ge=0, le=1)
    salt_bridge_occupancy: float = Field(ge=0, le=1)
    ligand_center_distance_nm: float = Field(ge=0)
    contact_atom_count: float = Field(ge=0)
    late_unbound_fraction: float = Field(ge=0, le=1)
    mmgbsa_kcal_mol: float | None = None


class TrajectoryAnalysis(BaseModel):
    model_config = ConfigDict(extra="forbid")

    status: Literal["stable", "unstable", "failed"]
    md_score: float | None = Field(default=None, ge=0, le=1)
    error_code: str | None = None
    flags: list[str] = Field(default_factory=list)
    mmgbsa_kcal_mol: float | None = None
    mmgbsa_role: Literal["interpretive_only"] = "interpretive_only"


def _clipped_inverse(value: float, scale: float = 1.0) -> float:
    return max(0.0, min(1.0, 1.0 - value / scale))


def analyze_trajectory(metrics: TrajectoryMetrics) -> TrajectoryAnalysis:
    """QC and score a completed trajectory without treating MM/GBSA as truth."""
    if metrics.completed_ns < 95.0 or metrics.frames / metrics.expected_frames < 0.95:
        return TrajectoryAnalysis(
            status="failed",
            error_code="trajectory_incomplete",
            mmgbsa_kcal_mol=metrics.mmgbsa_kcal_mol,
        )
    if not metrics.energy_continuous:
        return TrajectoryAnalysis(
            status="failed",
            error_code="energy_discontinuity",
            mmgbsa_kcal_mol=metrics.mmgbsa_kcal_mol,
        )
    if not metrics.pbc_processed:
        return TrajectoryAnalysis(
            status="failed",
            error_code="pbc_processing_missing",
            mmgbsa_kcal_mol=metrics.mmgbsa_kcal_mol,
        )

    structure = sum(
        [
            _clipped_inverse(metrics.protein_ca_rmsd_nm),
            _clipped_inverse(metrics.ligand_pocket_rmsd_nm),
            _clipped_inverse(metrics.pocket_rmsf_nm),
        ]
    ) / 3.0
    contacts = sum(
        [
            metrics.hbond_occupancy,
            metrics.hydrophobic_occupancy,
            metrics.salt_bridge_occupancy,
        ]
    ) / 3.0
    retention = 1.0 - metrics.late_unbound_fraction
    score = round(0.300 * structure + 0.395 * contacts + 0.305 * retention, 3)

    flags: list[str] = []
    if metrics.late_unbound_fraction >= 0.80:
        flags.append("ligand_departed_pocket")
    if metrics.contact_atom_count < 3:
        flags.append("low_contact_count")

    return TrajectoryAnalysis(
        status="unstable" if flags else "stable",
        md_score=score,
        flags=flags,
        mmgbsa_kcal_mol=metrics.mmgbsa_kcal_mol,
    )

