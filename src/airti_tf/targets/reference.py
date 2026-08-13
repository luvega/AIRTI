"""Compile and validate a frozen human target-library release bundle."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field

from airti_tf.manifest_io import (
    read_jsonl,
    write_artifact,
    write_bytes_atomic,
    write_jsonl_atomic,
)
from airti_tf.pockets.receptor import DockingBox
from airti_tf.sources.uniprot import UniProtRecord
from airti_tf.stages import TargetPocketRow


class ReferenceBundleSummary(BaseModel):
    model_config = ConfigDict(extra="forbid")

    manifest_path: Path
    manifest_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    target_count: int = Field(gt=0)
    ready_target_count: int = Field(ge=0)
    unsupported_target_count: int = Field(ge=0)
    failed_target_count: int = Field(ge=0)
    artifact_count: int = Field(gt=0)


class PilotTargetManifestSummary(BaseModel):
    model_config = ConfigDict(extra="forbid")

    manifest_path: Path
    manifest_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    target_count: int = Field(gt=0)
    ready_target_count: int = Field(ge=0)
    unsupported_target_count: int = Field(ge=0)


def _relative_asset(path: Path, *, root: Path) -> tuple[Path, str]:
    if path.is_absolute():
        resolved = path.resolve()
    else:
        resolved = (root / path).resolve()
    try:
        relative = resolved.relative_to(root.resolve())
    except ValueError as error:
        raise ValueError(f"reference artifact escapes root: {path}") from error
    if not resolved.is_file() or resolved.stat().st_size == 0:
        raise ValueError(f"reference artifact is missing or empty: {relative}")
    return resolved, relative.as_posix()


def _artifact_record(path: Path, *, root: Path) -> dict[str, object]:
    resolved, relative = _relative_asset(path, root=root)
    raw = resolved.read_bytes()
    return {
        "path": relative,
        "sha256": hashlib.sha256(raw).hexdigest(),
        "size_bytes": len(raw),
    }


def compile_reference_bundle(
    *,
    root: Path,
    proteome_manifest: Path,
    target_manifest: Path,
    output_manifest: Path,
) -> ReferenceBundleSummary:
    """Validate canonical coverage and hash every runtime-required artifact."""
    root = root.resolve()
    proteome_path, _ = _relative_asset(proteome_manifest, root=root)
    target_path, _ = _relative_asset(target_manifest, root=root)
    output_path = output_manifest.resolve()
    try:
        output_path.relative_to(root)
    except ValueError as error:
        raise ValueError("reference manifest output must be inside reference root") from error

    proteome = [
        UniProtRecord.model_validate(row) for row in read_jsonl(proteome_path)
    ]
    targets = [
        TargetPocketRow.model_validate(row) for row in read_jsonl(target_path)
    ]
    if not proteome:
        raise ValueError("canonical proteome manifest is empty")
    if not targets:
        raise ValueError("target pocket manifest is empty")
    releases = {record.release for record in proteome}
    if len(releases) != 1:
        raise ValueError("canonical proteome contains mixed UniProt releases")
    if any(record.taxonomy_id != 9606 for record in proteome):
        raise ValueError("reference bundle must contain human taxonomy 9606 only")
    if any(not record.reviewed for record in proteome):
        raise ValueError("production canonical proteome must contain reviewed entries only")
    proteome_by_id = {record.uniprot_id: record for record in proteome}
    if len(proteome_by_id) != len(proteome):
        raise ValueError("canonical proteome contains duplicate accessions")
    target_ids = {row.target_id for row in targets}
    if target_ids != set(proteome_by_id):
        missing = sorted(set(proteome_by_id) - target_ids)
        extra = sorted(target_ids - set(proteome_by_id))
        raise ValueError(
            "canonical target coverage mismatch: "
            f"missing={missing[:10]}, extra={extra[:10]}"
        )

    status_by_target: dict[str, str] = {}
    runtime_assets: set[Path] = {proteome_path, target_path}
    for target in targets:
        previous = status_by_target.setdefault(target.target_id, target.status)
        if previous != target.status:
            raise ValueError(f"mixed statuses for target {target.target_id}")
        canonical = proteome_by_id[target.target_id]
        if target.status != "ready":
            continue
        if (
            target.sequence != canonical.sequence
            or target.sequence_sha256 != canonical.sequence_sha256
        ):
            raise ValueError(f"sequence mismatch for ready target {target.target_id}")
        assert target.model_sequence is not None
        assert target.model_sequence_sha256 is not None
        assert target.model_sequence_start is not None
        assert target.model_sequence_end is not None
        expected_model_sequence = canonical.sequence[
            target.model_sequence_start - 1 : target.model_sequence_end
        ]
        if (
            target.model_sequence != expected_model_sequence
            or target.model_sequence_sha256
            != hashlib.sha256(expected_model_sequence.encode()).hexdigest()
        ):
            raise ValueError(
                f"model sequence mismatch for ready target {target.target_id}"
            )
        assert target.receptor_pdbqt_path is not None
        assert target.msa_path is not None
        assert target.structure_path is not None
        assert target.calibration_path is not None
        if any(
            path.is_absolute()
            for path in (
                target.receptor_pdbqt_path,
                target.msa_path,
                target.structure_path,
                target.calibration_path,
            )
        ):
            raise ValueError("target runtime asset paths must be relative to reference root")
        receptor, _ = _relative_asset(target.receptor_pdbqt_path, root=root)
        msa, _ = _relative_asset(target.msa_path, root=root)
        structure, _ = _relative_asset(target.structure_path, root=root)
        calibration, _ = _relative_asset(target.calibration_path, root=root)
        runtime_assets.update((receptor, msa, structure, calibration))

    artifacts = [
        _artifact_record(path, root=root)
        for path in sorted(runtime_assets, key=lambda item: str(item))
    ]
    ready_target_count = sum(
        status == "ready" for status in status_by_target.values()
    )
    unsupported_target_count = sum(
        status == "unsupported" for status in status_by_target.values()
    )
    failed_target_count = sum(
        status == "failed" for status in status_by_target.values()
    )
    payload = {
        "schema_version": "1.0",
        "organism": "Homo sapiens",
        "taxonomy_id": 9606,
        "proteome_id": "UP000005640",
        "uniprot_release": next(iter(releases)),
        "target_count": len(status_by_target),
        "ready_target_count": ready_target_count,
        "unsupported_target_count": unsupported_target_count,
        "failed_target_count": failed_target_count,
        "artifacts": artifacts,
    }
    written = write_artifact(output_path, payload)
    return ReferenceBundleSummary(
        manifest_path=written.path,
        manifest_sha256=written.sha256,
        target_count=len(status_by_target),
        ready_target_count=ready_target_count,
        unsupported_target_count=unsupported_target_count,
        failed_target_count=failed_target_count,
        artifact_count=len(artifacts),
    )


def _distance(
    left: tuple[float, float, float], right: tuple[float, float, float]
) -> float:
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(left, right, strict=True)))


def _pocket_residues_from_ligand(
    structure_pdb: Path,
    *,
    target_id: str,
    ligand_resname: str,
    contact_distance_a: float = 6.0,
) -> list[int]:
    text = structure_pdb.read_text(encoding="utf-8")
    mappings: list[tuple[str, int, int, int]] = []
    for line in text.splitlines():
        if not line.startswith("DBREF "):
            continue
        fields = line.split()
        if len(fields) >= 10 and fields[6] == target_id:
            mappings.append(
                (fields[2], int(fields[3]), int(fields[4]), int(fields[8]))
            )

    def canonical_residue(chain: str, residue: int) -> int:
        for mapped_chain, pdb_start, pdb_end, uniprot_start in mappings:
            if chain == mapped_chain and pdb_start <= residue <= pdb_end:
                return uniprot_start + residue - pdb_start
        if mappings:
            raise ValueError(
                f"PDB residue {chain}:{residue} is outside the UniProt DBREF mapping"
            )
        return residue

    protein_atoms: list[tuple[int, tuple[float, float, float]]] = []
    ligand_atoms: list[tuple[float, float, float]] = []
    for line in text.splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        try:
            coordinate = (
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54]),
            )
        except ValueError:
            continue
        residue_name = line[17:20].strip()
        if line.startswith("HETATM") and residue_name == ligand_resname:
            ligand_atoms.append(coordinate)
        elif line.startswith("ATOM  "):
            try:
                residue_number = int(line[22:26])
            except ValueError:
                continue
            protein_atoms.append(
                (canonical_residue(line[21:22], residue_number), coordinate)
            )
    if not ligand_atoms:
        raise ValueError(f"structure contains no ligand residue {ligand_resname}")
    residues = sorted(
        {
            residue
            for residue, protein in protein_atoms
            if any(
                _distance(protein, ligand) <= contact_distance_a
                for ligand in ligand_atoms
            )
        }
    )
    if not residues:
        raise ValueError("co-crystal ligand has no contacting protein residues")
    return residues


def _model_sequence_bounds(structure_pdb: Path, *, target_id: str) -> tuple[int, int]:
    mappings: list[tuple[int, int]] = []
    for line in structure_pdb.read_text(encoding="utf-8").splitlines():
        if not line.startswith("DBREF "):
            continue
        fields = line.split()
        if len(fields) >= 10 and fields[6] == target_id:
            mappings.append((int(fields[8]), int(fields[9])))
    if not mappings:
        raise ValueError("pilot PDB contains no UniProt DBREF sequence mapping")
    start = min(item[0] for item in mappings)
    end = max(item[1] for item in mappings)
    if start < 1 or end < start:
        raise ValueError("pilot PDB contains an invalid UniProt DBREF range")
    return start, end


def build_pilot_target_manifest(
    *,
    root: Path,
    proteome_manifest: Path,
    target_id: str,
    family: str,
    structure_id: str,
    structure_source: Literal["pdb", "alphafold"],
    structure_pdb: Path,
    ligand_resname: str,
    receptor_pdbqt: Path,
    calibration_json: Path,
    box: DockingBox | dict[str, object],
    structure_quality: float,
    msa_database_version: str,
    output_manifest: Path,
) -> PilotTargetManifestSummary:
    """Build an honest full-coverage manifest with selected pilot targets ready."""
    root = root.resolve()
    proteome_path, _ = _relative_asset(proteome_manifest, root=root)
    structure_path, structure_relative = _relative_asset(structure_pdb, root=root)
    _receptor_path, receptor_relative = _relative_asset(receptor_pdbqt, root=root)
    calibration_path, calibration_relative = _relative_asset(
        calibration_json, root=root
    )
    proteome = [
        UniProtRecord.model_validate(row) for row in read_jsonl(proteome_path)
    ]
    by_id = {record.uniprot_id: record for record in proteome}
    if target_id not in by_id:
        raise ValueError(f"pilot target is absent from canonical proteome: {target_id}")
    calibration = json.loads(calibration_path.read_text(encoding="utf-8"))
    affinities = calibration.get("background_affinities")
    successful_count = calibration.get("successful_probe_count")
    if (
        not isinstance(affinities, list)
        or len(affinities) < 95
        or not isinstance(successful_count, int)
        or successful_count < 95
    ):
        raise ValueError("pilot pocket calibration has fewer than 95 valid probes")
    validated_box = DockingBox.model_validate(box)
    pocket_residues = _pocket_residues_from_ligand(
        structure_path,
        target_id=target_id,
        ligand_resname=ligand_resname,
    )
    target = by_id[target_id]
    model_sequence_start, model_sequence_end = _model_sequence_bounds(
        structure_path, target_id=target_id
    )
    if model_sequence_end > len(target.sequence):
        raise ValueError("pilot PDB DBREF range exceeds the canonical sequence")
    model_sequence = target.sequence[
        model_sequence_start - 1 : model_sequence_end
    ]
    model_sequence_sha256 = hashlib.sha256(model_sequence.encode()).hexdigest()
    model_pocket_residues = [
        residue - model_sequence_start + 1 for residue in pocket_residues
    ]
    if any(
        residue < 1 or residue > len(model_sequence)
        for residue in model_pocket_residues
    ):
        raise ValueError("pilot pocket residues fall outside the model sequence")
    msa_path = root / "targets" / target_id / f"{target_id}.single_sequence.a3m"
    write_bytes_atomic(
        msa_path,
        f">{target_id}:{model_sequence_start}-{model_sequence_end}\n"
        f"{model_sequence}\n".encode("utf-8"),
    )
    msa_relative = msa_path.relative_to(root).as_posix()
    rows: list[dict[str, object]] = []
    for record in sorted(proteome, key=lambda item: item.uniprot_id):
        if record.uniprot_id != target_id:
            rows.append(
                {
                    "schema_version": "1.0",
                    "target_id": record.uniprot_id,
                    "gene_symbol": record.gene_primary,
                    "family": "unknown",
                    "status": "unsupported",
                    "unsupported_reason": "pilot_not_built",
                }
            )
            continue
        rows.append(
            {
                "schema_version": "1.0",
                "target_id": record.uniprot_id,
                "gene_symbol": record.gene_primary,
                "family": family,
                "status": "ready",
                "unsupported_reason": None,
                "sequence": record.sequence,
                "sequence_sha256": record.sequence_sha256,
                "model_sequence": model_sequence,
                "model_sequence_sha256": model_sequence_sha256,
                "model_sequence_start": model_sequence_start,
                "model_sequence_end": model_sequence_end,
                "structure_quality": structure_quality,
                "structure_id": structure_id,
                "structure_source": structure_source,
                "structure_path": structure_relative,
                "calibration_path": calibration_relative,
                "pocket_id": f"{target_id}:{structure_id}:{ligand_resname}",
                "receptor_pdbqt_path": receptor_relative,
                "box": validated_box.model_dump(mode="json"),
                "background_affinities": [float(value) for value in affinities],
                "msa_path": msa_relative,
                "msa_database_version": msa_database_version,
                "pocket_residues": pocket_residues,
                "model_pocket_residues": model_pocket_residues,
            }
        )
    digest = write_jsonl_atomic(output_manifest, rows)
    return PilotTargetManifestSummary(
        manifest_path=output_manifest,
        manifest_sha256=digest,
        target_count=len(rows),
        ready_target_count=1,
        unsupported_target_count=len(rows) - 1,
    )
