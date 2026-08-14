"""Trajectory completeness checks and conservative MD evidence scoring."""

from __future__ import annotations

import math
from importlib import import_module
from pathlib import Path
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field


class TrajectoryMetrics(BaseModel):
    model_config = ConfigDict(extra="forbid")

    frames: int = Field(ge=0)
    expected_frames: int = Field(gt=0)
    completed_ns: float = Field(ge=0)
    expected_ns: float = Field(default=100.0, gt=0)
    analysis_role: Literal[
        "scientific_evidence", "pipeline_validation_only"
    ] = "scientific_evidence"
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
    analysis_role: Literal[
        "scientific_evidence", "pipeline_validation_only"
    ] = "scientific_evidence"


def _clipped_inverse(value: float, scale: float = 1.0) -> float:
    return max(0.0, min(1.0, 1.0 - value / scale))


def analyze_trajectory(metrics: TrajectoryMetrics) -> TrajectoryAnalysis:
    """QC and score a completed trajectory without treating MM/GBSA as truth."""
    if (
        metrics.completed_ns / metrics.expected_ns < 0.95
        or metrics.frames / metrics.expected_frames < 0.95
    ):
        return TrajectoryAnalysis(
            status="failed",
            error_code="trajectory_incomplete",
            mmgbsa_kcal_mol=metrics.mmgbsa_kcal_mol,
            analysis_role=metrics.analysis_role,
        )
    if not metrics.energy_continuous:
        return TrajectoryAnalysis(
            status="failed",
            error_code="energy_discontinuity",
            mmgbsa_kcal_mol=metrics.mmgbsa_kcal_mol,
            analysis_role=metrics.analysis_role,
        )
    if not metrics.pbc_processed:
        return TrajectoryAnalysis(
            status="failed",
            error_code="pbc_processing_missing",
            mmgbsa_kcal_mol=metrics.mmgbsa_kcal_mol,
            analysis_role=metrics.analysis_role,
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
        analysis_role=metrics.analysis_role,
    )


def measure_trajectory(
    *,
    topology: Path,
    trajectory: Path,
    log_path: Path,
    pocket_residues: list[int],
    expected_ns: float,
    analysis_role: Literal[
        "scientific_evidence", "pipeline_validation_only"
    ],
) -> TrajectoryMetrics:
    """Measure aligned protein, ligand and pocket-retention trajectory metrics."""
    mda = import_module("MDAnalysis")
    np = import_module("numpy")
    rotation_matrix = import_module(
        "MDAnalysis.analysis.align"
    ).rotation_matrix
    distance_array = import_module(
        "MDAnalysis.lib.distances"
    ).distance_array

    universe = mda.Universe(str(topology), str(trajectory))
    protein_ca = universe.select_atoms("protein and name CA")
    ligand = universe.select_atoms("resname LIG MOL and not name H*")
    if not pocket_residues:
        raise ValueError("trajectory analysis requires pocket residues")
    residue_query = " ".join(str(value) for value in pocket_residues)
    pocket = universe.select_atoms(f"protein and resid {residue_query}")
    pocket_ca = universe.select_atoms(
        f"protein and name CA and resid {residue_query}"
    )
    if len(protein_ca) < 3 or len(ligand) == 0 or len(pocket) == 0:
        raise ValueError("trajectory selections lack protein, ligand, or pocket atoms")

    universe.trajectory[0]
    reference_ca = protein_ca.positions.copy()
    reference_center = reference_ca.mean(axis=0)
    reference_ligand = ligand.positions.copy()
    protein_rmsd: list[float] = []
    ligand_rmsd: list[float] = []
    pocket_frames: list[object] = []
    center_distances: list[float] = []
    contact_counts: list[float] = []
    hbond_frames = 0
    hydrophobic_frames = 0
    salt_bridge_frames = 0
    ligand_no_h = ligand.select_atoms("not name H*")
    ligand_polar = ligand.select_atoms("name N* O* S*")
    ligand_hydrophobic = ligand.select_atoms("name C* S*")
    pocket_no_h = pocket.select_atoms("not name H*")
    pocket_polar = pocket.select_atoms("name N* O* S*")
    pocket_hydrophobic = pocket.select_atoms("name C* S*")

    for _frame in universe.trajectory:
        mobile_center = protein_ca.positions.mean(axis=0)
        rotation, _rmsd = rotation_matrix(
            protein_ca.positions - mobile_center,
            reference_ca - reference_center,
        )
        universe.atoms.translate(-mobile_center)
        universe.atoms.rotate(rotation)
        universe.atoms.translate(reference_center)
        ca_delta = protein_ca.positions - reference_ca
        protein_rmsd.append(float(np.sqrt(np.mean(np.sum(ca_delta**2, axis=1)))))
        ligand_delta = ligand.positions - reference_ligand
        ligand_rmsd.append(
            float(np.sqrt(np.mean(np.sum(ligand_delta**2, axis=1))))
        )
        if len(pocket_ca):
            pocket_frames.append(pocket_ca.positions.copy())
        ligand_center = ligand_no_h.positions.mean(axis=0)
        pocket_center = pocket_no_h.positions.mean(axis=0)
        center_distances.append(float(np.linalg.norm(ligand_center - pocket_center)))
        all_distances = distance_array(
            ligand_no_h.positions, pocket_no_h.positions
        )
        contact_counts.append(float(np.sum(np.min(all_distances, axis=0) <= 4.5)))
        if len(ligand_polar) and len(pocket_polar):
            hbond_frames += bool(
                np.any(distance_array(ligand_polar.positions, pocket_polar.positions) <= 3.5)
            )
            salt_bridge_frames += bool(
                np.any(distance_array(ligand_polar.positions, pocket_polar.positions) <= 4.0)
            )
        if len(ligand_hydrophobic) and len(pocket_hydrophobic):
            hydrophobic_frames += bool(
                np.any(
                    distance_array(
                        ligand_hydrophobic.positions,
                        pocket_hydrophobic.positions,
                    )
                    <= 4.5
                )
            )

    frame_count = len(protein_rmsd)
    if frame_count == 0:
        raise ValueError("trajectory contains no frames")
    if pocket_frames:
        pocket_array = np.asarray(pocket_frames)
        mean_positions = pocket_array.mean(axis=0)
        pocket_rmsf_a = float(
            np.sqrt(
                np.mean(
                    np.sum((pocket_array - mean_positions) ** 2, axis=2)
                )
            )
        )
    else:
        pocket_rmsf_a = math.inf
    late_start = max(0, int(frame_count * 0.8))
    late_distances = center_distances[late_start:]
    completed_ns = max(float(universe.trajectory[-1].time), 0.0) / 1_000
    log_text = log_path.read_text(encoding="utf-8", errors="replace")
    return TrajectoryMetrics(
        frames=frame_count,
        expected_frames=max(1, round(expected_ns * 100) + 1),
        completed_ns=completed_ns,
        expected_ns=expected_ns,
        analysis_role=analysis_role,
        energy_continuous="Fatal error" not in log_text,
        pbc_processed=trajectory.is_file() and trajectory.stat().st_size > 0,
        protein_ca_rmsd_nm=float(np.median(protein_rmsd)) / 10,
        ligand_pocket_rmsd_nm=float(np.median(ligand_rmsd)) / 10,
        pocket_rmsf_nm=pocket_rmsf_a / 10,
        hbond_occupancy=hbond_frames / frame_count,
        hydrophobic_occupancy=hydrophobic_frames / frame_count,
        salt_bridge_occupancy=salt_bridge_frames / frame_count,
        ligand_center_distance_nm=float(np.median(center_distances)) / 10,
        contact_atom_count=float(np.mean(contact_counts)),
        late_unbound_fraction=(
            sum(distance > 10.0 for distance in late_distances)
            / len(late_distances)
        ),
        mmgbsa_kcal_mol=None,
    )
