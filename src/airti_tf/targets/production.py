"""Production target acquisition, pocketing, receptor preparation and calibration."""

from __future__ import annotations

import hashlib
import importlib
import json
import math
import subprocess
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal, cast
from urllib.error import HTTPError, URLError
from urllib.request import urlopen

from airti_tf.manifest_io import write_artifact, write_bytes_atomic
from airti_tf.pockets.fpocket import (
    PocketQC,
    build_fpocket_command,
    parse_fpocket,
    select_qualified_pockets,
)
from airti_tf.pockets.receptor import (
    AtomCoordinate,
    BoxTooLargeError,
    DockingBox,
    build_box,
    build_meeko_command,
)
from airti_tf.screening.calibration_build import calibrate_pocket_background
from airti_tf.sources.uniprot import UniProtRecord
from airti_tf.stages import CofactorRecord, TargetPocketRow
from airti_tf.targets.structures import StructureCandidate, choose_structure

PDBe_BEST_STRUCTURES = "https://www.ebi.ac.uk/pdbe/api/mappings/best_structures"
RCSB_PDB_DOWNLOAD = "https://files.rcsb.org/download"
RCSB_LIGAND_DOWNLOAD = "https://files.rcsb.org/ligands/download"
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction"
OPM_PDB_DOWNLOAD = "https://opm-assets.storage.googleapis.com/pdb"
_DOWNLOAD_SEMAPHORE = threading.BoundedSemaphore(1)


@dataclass(frozen=True)
class _ProductionStructure:
    candidate: StructureCandidate
    chain_id: str
    model_start: int
    model_end: int
    pdb_start: int
    download_url: str


def _read_url(url: str, *, optional: bool = False) -> bytes | None:
    for attempt in range(5):
        try:
            # The EBI endpoints intermittently reset concurrent TLS sessions.
            # Serialize only network I/O; pocketing and docking remain parallel.
            with _DOWNLOAD_SEMAPHORE:
                with urlopen(url, timeout=120) as response:  # noqa: S310
                    return cast(bytes, response.read())
        except HTTPError as error:
            if optional and error.code == 404:
                return None
            if error.code not in {429, 500, 502, 503, 504} or attempt == 4:
                raise
        except (URLError, TimeoutError):
            if attempt == 4:
                raise
        time.sleep(0.5 * (2**attempt))
    raise RuntimeError("unreachable HTTP retry state")


def _read_json_url(url: str, *, optional: bool = False) -> Any:
    raw = _read_url(url, optional=optional)
    if raw is None:
        return None
    return json.loads(raw)


def _method(value: str) -> Literal["x-ray", "electron_microscopy"] | None:
    normalized = value.lower()
    if "x-ray" in normalized:
        return "x-ray"
    if "electron" in normalized:
        return "electron_microscopy"
    return None


def _pdbe_structures(record: UniProtRecord) -> list[_ProductionStructure]:
    payload = _read_json_url(
        f"{PDBe_BEST_STRUCTURES}/{record.uniprot_id}", optional=True
    )
    if not isinstance(payload, dict):
        return []
    rows = payload.get(record.uniprot_id) or payload.get(record.uniprot_id.lower())
    if not isinstance(rows, list):
        return []
    structures: list[_ProductionStructure] = []
    for row in rows:
        if not isinstance(row, dict):
            continue
        method = _method(str(row.get("experimental_method", "")))
        resolution = row.get("resolution")
        if method is None or not isinstance(resolution, (int, float)):
            continue
        coverage = float(row.get("coverage", 0.0))
        model_start = int(row.get("unp_start", 1))
        model_end = int(row.get("unp_end", len(record.sequence)))
        pdb_start = int(row.get("start", 1))
        pdb_id = str(row.get("pdb_id", "")).lower()
        chain_id = str(row.get("chain_id", ""))
        if (
            not pdb_id
            or not chain_id
            or model_start < 1
            or model_end > len(record.sequence)
        ):
            continue
        structures.append(
            _ProductionStructure(
                candidate=StructureCandidate(
                    structure_id=pdb_id.upper(),
                    source="pdb",
                    experimental_method=method,
                    coverage=coverage,
                    mainchain_missing_fraction=max(0.0, 1.0 - coverage),
                    resolution=float(resolution),
                    has_ligand=False,
                ),
                chain_id=chain_id,
                model_start=model_start,
                model_end=model_end,
                pdb_start=pdb_start,
                download_url=f"{RCSB_PDB_DOWNLOAD}/{pdb_id}.pdb",
            )
        )
    return structures


def _alphafold_structure(record: UniProtRecord) -> _ProductionStructure | None:
    payload = _read_json_url(
        f"{ALPHAFOLD_API}/{record.uniprot_id}", optional=True
    )
    if not isinstance(payload, list) or not payload or not isinstance(payload[0], dict):
        return None
    row = payload[0]
    pdb_url = row.get("pdbUrl")
    confidence = row.get("globalMetricValue")
    if not isinstance(pdb_url, str) or not isinstance(confidence, (int, float)):
        return None
    model_start = int(row.get("uniprotStart", 1))
    model_end = int(row.get("uniprotEnd", len(record.sequence)))
    if model_start < 1 or model_end > len(record.sequence):
        return None
    return _ProductionStructure(
        candidate=StructureCandidate(
            structure_id=str(row.get("entryId", f"AF-{record.uniprot_id}-F1")),
            source="alphafold",
            coverage=(model_end - model_start + 1) / len(record.sequence),
            confidence=float(confidence) / 100,
            pae_supported=bool(row.get("paeDocUrl")),
        ),
        chain_id=str(row.get("chainId", "A")),
        model_start=model_start,
        model_end=model_end,
        pdb_start=1,
        download_url=pdb_url,
    )


def _select_structure(record: UniProtRecord) -> _ProductionStructure | None:
    structures = _pdbe_structures(record)
    selection = choose_structure([item.candidate for item in structures])
    if selection.status != "ready":
        alphafold = _alphafold_structure(record)
        if alphafold is not None:
            structures.append(alphafold)
            selection = choose_structure([item.candidate for item in structures])
    if selection.status != "ready":
        return None
    return next(
        item
        for item in structures
        if item.candidate.structure_id == selection.structure_id
        and item.candidate.source == selection.source
    )


def _clean_structure(
    raw: bytes, *, chain_id: str, keep_hem: bool
) -> tuple[bytes, bool]:
    text = raw.decode("utf-8")
    records: list[str] = []
    retained_hem = False
    for line in text.splitlines():
        if line.startswith(("HEADER", "TITLE ", "EXPDTA", "DBREF ", "REMARK 465")):
            records.append(line)
        elif (
            line.startswith("ATOM  ")
            and line[21:22] == chain_id
            and line[16:17] in {" ", "A"}
        ):
            records.append(f"{line[:16]} {line[17:]}")
        elif (
            keep_hem
            and line.startswith("HETATM")
            and line[17:20].strip() == "HEM"
            and line[21:22] in {chain_id, " "}
            and line[16:17] in {" ", "A"}
        ):
            records.append(f"{line[:16]} {line[17:]}")
            retained_hem = True
    if not any(line.startswith("ATOM  ") for line in records):
        raise ValueError(f"selected structure contains no atoms for chain {chain_id}")
    records.extend(["TER", "END"])
    return ("\n".join(records) + "\n").encode(), retained_hem


def _patch_heme_cif_charges(raw: str) -> str:
    """Patch the three CCD charge fields required by Meeko chemtempgen.

    The pristine CCD file remains beside this generated derivative.  Refusing
    unexpected row counts or source charges makes an upstream CCD change
    visible instead of silently creating a chemically different template.
    """
    expected = {"FE": "?", "O2A": "-1", "O2D": "-1"}
    observed = {atom_id: 0 for atom_id in expected}
    output: list[str] = []
    for line in raw.splitlines():
        fields = line.split()
        atom_id = (
            fields[1]
            if len(fields) >= 5
            and fields[0] == "HEM"
            and fields[1] == fields[2]
            else None
        )
        if atom_id in expected:
            if fields[4] != expected[atom_id]:
                raise ValueError(
                    f"unexpected HEM {atom_id} formal charge: {fields[4]}"
                )
            observed[atom_id] += 1
            fields[4] = "0"
            line = " ".join(fields)
        output.append(line)
    invalid = [atom_id for atom_id, count in observed.items() if count != 1]
    if invalid:
        raise ValueError(
            "HEM CCD must contain exactly one expected atom row for: "
            + ", ".join(invalid)
        )
    return "\n".join(output) + "\n"


def _write_meeko_heme_template(*, cif_path: Path, output_path: Path) -> Path:
    """Generate a chemically connected HEM template with Meeko itself."""
    patched_path = output_path.with_name("HEM_meeko_source.cif")
    patched = _patch_heme_cif_charges(cif_path.read_text(encoding="utf-8"))
    patched_sha256 = write_bytes_atomic(patched_path, patched.encode())

    chemtempgen = importlib.import_module("meeko.chemtempgen")
    component = chemtempgen.ChemicalComponent.from_cif(str(patched_path), "HEM")
    component = component.make_canonical(
        chemtempgen.acidic_proton_loc_canonical
    )
    component.ResidueTemplate_check()
    chemtempgen.export_chem_templates_to_json(
        [component], json_fname=str(output_path)
    )
    if not output_path.is_file() or output_path.stat().st_size == 0:
        raise RuntimeError("Meeko did not produce the HEM residue template")
    write_artifact(
        output_path.with_name("HEM_meeko_template.provenance.json"),
        {
            "schema_version": "1.0",
            "source": str(cif_path.name),
            "source_sha256": hashlib.sha256(cif_path.read_bytes()).hexdigest(),
            "patched_source": str(patched_path.name),
            "patched_source_sha256": patched_sha256,
            "transformations": {
                "FE.formal_charge": ["?", "0"],
                "O2A.formal_charge": ["-1", "0"],
                "O2D.formal_charge": ["-1", "0"],
            },
            "purpose": "Meeko docking receptor template only; not an MD parameter set",
        },
    )
    return output_path


def _run(command: list[str], *, cwd: Path, timeout: int) -> None:
    result = subprocess.run(
        command,
        cwd=cwd,
        check=False,
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    if result.returncode != 0:
        detail = result.stderr.strip() or result.stdout.strip() or "unknown error"
        raise RuntimeError(f"command failed ({command[0]}): {detail}")


def _pocket_atoms(path: Path) -> list[AtomCoordinate]:
    atoms: list[AtomCoordinate] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        try:
            atoms.append(
                AtomCoordinate(
                    coord=(
                        float(line[30:38]),
                        float(line[38:46]),
                        float(line[46:54]),
                    )
                )
            )
        except ValueError:
            continue
    return atoms


def _dbref_mappings(
    structure: Path, *, target_id: str, chain_id: str
) -> list[tuple[int, int, int]]:
    mappings: list[tuple[int, int, int]] = []
    for line in structure.read_text(encoding="utf-8").splitlines():
        if not line.startswith("DBREF "):
            continue
        fields = line.split()
        if len(fields) >= 10 and fields[2] == chain_id and fields[6] == target_id:
            mappings.append((int(fields[3]), int(fields[4]), int(fields[8])))
    return mappings


def _pocket_residues(
    pocket: Path,
    *,
    mappings: list[tuple[int, int, int]],
    fallback_pdb_start: int,
    fallback_uniprot_start: int,
) -> list[int]:
    residues: set[int] = set()
    for line in pocket.read_text(encoding="utf-8").splitlines():
        if not line.startswith("ATOM  "):
            continue
        try:
            pdb_residue = int(line[22:26])
        except ValueError:
            continue
        canonical = None
        for pdb_start, pdb_end, uniprot_start in mappings:
            if pdb_start <= pdb_residue <= pdb_end:
                canonical = uniprot_start + pdb_residue - pdb_start
                break
        if canonical is None:
            canonical = fallback_uniprot_start + pdb_residue - fallback_pdb_start
        residues.add(canonical)
    return sorted(residue for residue in residues if residue > 0)


def _mean_pocket_plddt(path: Path) -> float:
    values: list[float] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("ATOM  "):
            try:
                values.append(float(line[60:66]))
            except ValueError:
                continue
    return sum(values) / len(values) if values else 0.0


def _select_dockable_pockets(
    qualified: list[PocketQC],
    *,
    pocket_output: Path,
    max_pockets: int,
) -> list[tuple[PocketQC, Path, DockingBox]]:
    """Keep the best pockets representable by an untruncated Vina box."""
    selected: list[tuple[PocketQC, Path, DockingBox]] = []
    for pocket_qc in qualified:
        pocket_pdb = (
            pocket_output
            / "pockets"
            / f"pocket{pocket_qc.pocket.rank}_atm.pdb"
        )
        try:
            box = build_box(
                _pocket_atoms(pocket_pdb),
                padding_a=5.0,
                min_size_a=18.0,
                max_size_a=32.0,
            )
        except BoxTooLargeError:
            continue
        selected.append((pocket_qc, pocket_pdb, box))
        if len(selected) == max_pockets:
            break
    return selected


def _unsupported(record: UniProtRecord, reason: str) -> list[TargetPocketRow]:
    return [
        TargetPocketRow(
            schema_version="1.1",
            target_id=record.uniprot_id,
            gene_symbol=record.gene_primary,
            family=record.protein_family or "unclassified",
            status="unsupported",
            unsupported_reason=reason,
            environment=(
                "membrane"
                if record.transmembrane_segments or record.membrane_associated
                else "soluble"
            ),
        )
    ]


def build_production_target(
    record: UniProtRecord, context: Any
) -> list[TargetPocketRow]:
    """Build up to ``max_pockets`` calibrated rows for one canonical target."""
    selected = _select_structure(record)
    if selected is None:
        return _unsupported(record, "no_qualified_structure")
    target_dir = Path(context.root) / "targets" / record.uniprot_id
    target_dir.mkdir(parents=True, exist_ok=True)
    raw_url = selected.download_url
    orientation_source: str | None = None
    environment: Literal["soluble", "soluble_construct", "membrane"] = "soluble"
    if record.transmembrane_segments or record.membrane_associated:
        if selected.candidate.source != "pdb":
            return _unsupported(record, "membrane_orientation_unavailable")
        opm_url = (
            f"{OPM_PDB_DOWNLOAD}/{selected.candidate.structure_id.lower()}.pdb"
        )
        opm_raw = _read_url(opm_url, optional=True)
        if opm_raw is None:
            return _unsupported(record, "membrane_orientation_unavailable")
        raw = opm_raw
        orientation_source = opm_url
        environment = "membrane"
    else:
        downloaded = _read_url(raw_url)
        if downloaded is None:
            raise RuntimeError("structure download unexpectedly returned no data")
        raw = downloaded
        if selected.model_start > 1 or selected.model_end < len(record.sequence):
            environment = "soluble_construct"
    raw_path = target_dir / "source.pdb"
    write_bytes_atomic(raw_path, raw)
    cleaned, retained_hem = _clean_structure(
        raw, chain_id=selected.chain_id, keep_hem=True
    )
    structure_path = target_dir / "structure.pdb"
    structure_hash = write_bytes_atomic(structure_path, cleaned)

    _run(build_fpocket_command(structure_path), cwd=target_dir, timeout=1800)
    pocket_output = target_dir / f"{structure_path.stem}_out"
    candidates = parse_fpocket(
        pocket_output,
        target_id=record.uniprot_id,
        structure_hash=structure_hash,
    )
    if selected.candidate.source == "alphafold":
        candidates = [
            candidate.model_copy(
                update={
                    "mean_plddt": _mean_pocket_plddt(
                        pocket_output
                        / "pockets"
                        / f"pocket{candidate.rank}_atm.pdb"
                    )
                }
            )
            for candidate in candidates
        ]
    qualified = select_qualified_pockets(
        candidates,
        limit=len(candidates),
        structure_source=selected.candidate.source,
    )
    if not qualified:
        return _unsupported(record, "no_qualified_pocket")
    dockable = _select_dockable_pockets(
        qualified,
        pocket_output=pocket_output,
        max_pockets=int(context.max_pockets),
    )
    if not dockable:
        return _unsupported(record, "no_dockable_pocket")

    templates: dict[str, Path] = {}
    cofactors: list[CofactorRecord] = []
    if retained_hem:
        heme_cif = target_dir / "HEM.cif"
        if not heme_cif.is_file():
            heme_raw = _read_url(f"{RCSB_LIGAND_DOWNLOAD}/HEM.cif")
            if heme_raw is None:
                raise RuntimeError("HEM template download returned no data")
            write_bytes_atomic(heme_cif, heme_raw)
        templates["HEM"] = _write_meeko_heme_template(
            cif_path=heme_cif,
            output_path=target_dir / "HEM_meeko_template.json",
        )
        cofactors.append(
            CofactorRecord(
                ccd_id="HEM",
                parameter_id="p450-ferric-thiolate-v1",
            )
        )

    model_sequence = record.sequence[
        selected.model_start - 1 : selected.model_end
    ]
    model_sequence_sha256 = hashlib.sha256(model_sequence.encode()).hexdigest()
    msa = target_dir / f"{record.uniprot_id}.single_sequence.a3m"
    write_bytes_atomic(
        msa,
        (
            f">{record.uniprot_id}:{selected.model_start}-{selected.model_end}\n"
            f"{model_sequence}\n"
        ).encode(),
    )
    mappings = _dbref_mappings(
        structure_path,
        target_id=record.uniprot_id,
        chain_id=selected.chain_id,
    )
    rows: list[TargetPocketRow] = []
    for pocket_qc, pocket_pdb, box in dockable:
        pocket = pocket_qc.pocket
        pocket_dir = target_dir / f"pocket-{pocket.rank}"
        pocket_dir.mkdir(exist_ok=True)
        receptor_prefix = pocket_dir / "receptor"
        _run(
            build_meeko_command(
                input_pdb=structure_path,
                output_prefix=receptor_prefix,
                box=box,
                residue_templates=templates,
            ),
            cwd=pocket_dir,
            timeout=1800,
        )
        receptor_pdbqt = receptor_prefix.with_suffix(".pdbqt")
        if not receptor_pdbqt.is_file() or receptor_pdbqt.stat().st_size == 0:
            raise RuntimeError("Meeko did not produce a receptor PDBQT")
        calibration_path = pocket_dir / "calibration.json"
        calibrate_pocket_background(
            Path(context.background_panel),
            receptor_pdbqt=receptor_pdbqt,
            box=box,
            pocket_id=pocket.pocket_id,
            output=calibration_path,
            asset_dir=pocket_dir / "calibration-assets",
            expected_probe_count=100,
            minimum_successful_probes=95,
            workers=int(context.calibration_workers),
        )
        calibration = json.loads(calibration_path.read_text(encoding="utf-8"))
        canonical_residues = _pocket_residues(
            pocket_pdb,
            mappings=mappings,
            fallback_pdb_start=selected.pdb_start,
            fallback_uniprot_start=selected.model_start,
        )
        canonical_residues = [
            residue
            for residue in canonical_residues
            if selected.model_start <= residue <= selected.model_end
        ]
        model_residues = [
            residue - selected.model_start + 1
            for residue in canonical_residues
        ]
        if not model_residues:
            raise ValueError("qualified pocket has no mapped model residues")
        rows.append(
            TargetPocketRow(
                schema_version="1.1",
                target_id=record.uniprot_id,
                gene_symbol=record.gene_primary,
                family=record.protein_family or "unclassified",
                status="ready",
                environment=environment,
                orientation_source=orientation_source,
                cofactors=cofactors,
                sequence=record.sequence,
                sequence_sha256=record.sequence_sha256,
                model_sequence=model_sequence,
                model_sequence_sha256=model_sequence_sha256,
                model_sequence_start=selected.model_start,
                model_sequence_end=selected.model_end,
                structure_quality=choose_structure(
                    [selected.candidate]
                ).score,
                structure_id=selected.candidate.structure_id,
                structure_source=selected.candidate.source,
                structure_path=structure_path.relative_to(context.root),
                calibration_path=calibration_path.relative_to(context.root),
                pocket_id=pocket.pocket_id,
                receptor_pdbqt_path=receptor_pdbqt.relative_to(context.root),
                box=box,
                background_affinities=[
                    float(value)
                    for value in calibration["background_affinities"]
                    if math.isfinite(float(value))
                ],
                msa_path=msa.relative_to(context.root),
                msa_database_version="single-sequence-v1",
                pocket_residues=canonical_residues,
                model_pocket_residues=model_residues,
            )
        )
    return rows
