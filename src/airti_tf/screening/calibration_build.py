"""Build per-pocket empirical docking background distributions."""

from __future__ import annotations

import csv
import hashlib
import statistics
import subprocess
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from pydantic import BaseModel, ConfigDict, Field

from airti_tf.ligands.prepare import prepare_ligand
from airti_tf.manifest_io import write_artifact, write_bytes_atomic
from airti_tf.pockets.receptor import DockingBox
from airti_tf.screening.quickvina import (
    DockingJob,
    DockingSeedResult,
    InsufficientSeedsError,
    run_docking_seed,
    summarize_seeds,
)

PDBQTPreparer = Callable[[Path, Path], None]
DockingRunner = Callable[[DockingJob, int], DockingSeedResult]


class PocketCalibrationSummary(BaseModel):
    model_config = ConfigDict(extra="forbid")

    pocket_id: str
    expected_probe_count: int = Field(gt=0)
    successful_probe_count: int = Field(ge=0)
    failed_probe_count: int = Field(ge=0)
    output_path: Path
    output_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")


class _ProbeResult(BaseModel):
    model_config = ConfigDict(extra="forbid")

    source_id: str
    status: str
    affinity_kcal_mol: float | None = None
    state_count: int = Field(ge=0)
    successful_state_count: int = Field(ge=0)
    error_code: str | None = None


def _run_meeko(input_sdf: Path, output_pdbqt: Path) -> None:
    result = subprocess.run(
        [
            "mk_prepare_ligand.py",
            "--mol",
            str(input_sdf),
            "--out",
            str(output_pdbqt),
            "--add_index_map",
        ],
        check=False,
        capture_output=True,
        text=True,
        timeout=300,
    )
    if result.returncode != 0:
        detail = result.stderr.strip() or result.stdout.strip() or "unknown error"
        raise RuntimeError(f"Meeko background-probe preparation failed: {detail}")
    if not output_pdbqt.is_file() or output_pdbqt.stat().st_size == 0:
        raise RuntimeError("Meeko did not create a background-probe PDBQT")


def _run_seed(job: DockingJob, seed: int) -> DockingSeedResult:
    return run_docking_seed(job, seed=seed)


def _read_panel(path: Path, *, expected_probe_count: int) -> list[tuple[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"source_id", "canonical_smiles"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError("background panel lacks source_id or canonical_smiles")
    if len(rows) != expected_probe_count:
        raise ValueError(
            f"background panel has {len(rows)} probes; {expected_probe_count} required"
        )
    probes = [
        (str(row["source_id"]).strip(), str(row["canonical_smiles"]).strip())
        for row in rows
    ]
    if any(not source_id or not smiles for source_id, smiles in probes):
        raise ValueError("background panel contains empty IDs or structures")
    if len({source_id for source_id, _ in probes}) != len(probes):
        raise ValueError("background panel source IDs must be unique")
    return probes


def calibrate_pocket_background(
    panel_tsv: Path,
    *,
    receptor_pdbqt: Path,
    box: DockingBox,
    pocket_id: str,
    output: Path,
    asset_dir: Path,
    expected_probe_count: int = 100,
    minimum_successful_probes: int = 95,
    workers: int = 8,
    pdbqt_preparer: PDBQTPreparer = _run_meeko,
    docking_runner: DockingRunner = _run_seed,
) -> PocketCalibrationSummary:
    """Dock a fixed panel and retain one best-state three-seed value per probe."""
    if not receptor_pdbqt.is_file():
        raise ValueError(f"receptor PDBQT does not exist: {receptor_pdbqt}")
    if not 1 <= minimum_successful_probes <= expected_probe_count:
        raise ValueError("invalid minimum successful background-probe count")
    if workers < 1:
        raise ValueError("workers must be positive")
    probes = _read_panel(panel_tsv, expected_probe_count=expected_probe_count)
    asset_dir.mkdir(parents=True, exist_ok=True)
    # Dimorphite-DL/RDKit preparation is intentionally serialized. Some of its
    # chemistry tables are process-global and concurrent enumeration can yield
    # non-reproducible state counts.
    prepared_by_source = {
        source_id: prepare_ligand(smiles, profile="production")
        for source_id, smiles in probes
    }
    pocket_key = hashlib.sha256(pocket_id.encode()).hexdigest()[:16]

    def calibrate_probe(probe: tuple[str, str]) -> _ProbeResult:
        source_id, _smiles = probe
        prepared = prepared_by_source[source_id]
        if prepared.status != "succeeded":
            return _ProbeResult(
                source_id=source_id,
                status="failed",
                state_count=0,
                successful_state_count=0,
                error_code=prepared.error_code or "ligand_preparation_failed",
            )
        probe_key = hashlib.sha256(source_id.encode()).hexdigest()[:16]
        prepared_dir = asset_dir / "prepared" / probe_key
        prepared_dir.mkdir(parents=True, exist_ok=True)
        docking_dir = asset_dir / "pockets" / pocket_key / probe_key
        state_affinities: list[float] = []
        for state_index, state in enumerate(prepared.states, 1):
            sdf = prepared_dir / f"state-{state_index}.sdf"
            pdbqt = prepared_dir / f"state-{state_index}.pdbqt"
            write_bytes_atomic(
                sdf,
                (state.mol_block.rstrip() + "\n$$$$\n").encode("utf-8"),
            )
            try:
                if not pdbqt.is_file() or pdbqt.stat().st_size == 0:
                    pdbqt_preparer(sdf, pdbqt)
            except (OSError, RuntimeError, subprocess.SubprocessError):
                continue
            job = DockingJob(
                job_id=f"{pocket_key}-{probe_key}-state-{state_index}",
                receptor_pdbqt=receptor_pdbqt,
                ligand_pdbqt=pdbqt,
                box=box,
                output_dir=docking_dir / f"state-{state_index}",
                # The executor owns calibration parallelism; without this,
                # every QuickVina process tries to claim all visible CPUs.
                cpu=1,
            )
            seed_results = [docking_runner(job, seed) for seed in (11, 29, 47)]
            try:
                state_affinities.append(summarize_seeds(seed_results).affinity_median)
            except InsufficientSeedsError:
                continue
        if not state_affinities:
            return _ProbeResult(
                source_id=source_id,
                status="failed",
                state_count=len(prepared.states),
                successful_state_count=0,
                error_code="no_successful_docking_state",
            )
        return _ProbeResult(
            source_id=source_id,
            status="succeeded",
            affinity_kcal_mol=min(state_affinities),
            state_count=len(prepared.states),
            successful_state_count=len(state_affinities),
        )

    with ThreadPoolExecutor(max_workers=workers) as executor:
        results = list(executor.map(calibrate_probe, probes))
    successful = [
        result for result in results if result.status == "succeeded"
    ]
    if len(successful) < minimum_successful_probes:
        raise ValueError(
            f"only {len(successful)} successful background probes; "
            f"{minimum_successful_probes} required"
        )
    affinities = [
        result.affinity_kcal_mol
        for result in successful
        if result.affinity_kcal_mol is not None
    ]
    payload = {
        "schema_version": "1.0",
        "pocket_id": pocket_id,
        "panel_path": str(panel_tsv.resolve()),
        "panel_sha256": hashlib.sha256(panel_tsv.read_bytes()).hexdigest(),
        "receptor_sha256": hashlib.sha256(receptor_pdbqt.read_bytes()).hexdigest(),
        "seeds": [11, 29, 47],
        "aggregation": "best_state_of_three_seed_medians",
        "expected_probe_count": expected_probe_count,
        "successful_probe_count": len(successful),
        "failed_probe_count": len(results) - len(successful),
        "background_affinities": affinities,
        "affinity_median": statistics.median(affinities),
        "probes": [result.model_dump(mode="json") for result in results],
    }
    written = write_artifact(output, payload)
    return PocketCalibrationSummary(
        pocket_id=pocket_id,
        expected_probe_count=expected_probe_count,
        successful_probe_count=len(successful),
        failed_probe_count=len(results) - len(successful),
        output_path=written.path,
        output_sha256=written.sha256,
    )
