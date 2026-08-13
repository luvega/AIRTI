import hashlib
import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

from airti_tf.cli import app
from airti_tf.runtime import _reference_manifest_status
from airti_tf.targets.reference import compile_reference_bundle
from airti_tf.targets.reference import build_pilot_target_manifest


def _sha256(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def _write_reference_fixture(root: Path) -> tuple[Path, Path]:
    root.mkdir()
    assets = root / "targets"
    assets.mkdir()
    sequence = "MPEPTIDE"
    proteome = root / "human_canonical_proteome.jsonl"
    proteome.write_text(
        json.dumps(
            {
                "schema_version": "1.0",
                "uniprot_id": "P00533",
                "gene_primary": "EGFR",
                "taxonomy_id": 9606,
                "sequence": sequence,
                "sequence_sha256": _sha256(sequence),
                "reviewed": True,
                "release": "2026_03",
                "isoform_aliases": [],
            }
        )
        + "\n"
        + json.dumps(
            {
                "schema_version": "1.0",
                "uniprot_id": "P0UNSP",
                "gene_primary": None,
                "taxonomy_id": 9606,
                "sequence": "MUNSUPPORTED",
                "sequence_sha256": _sha256("MUNSUPPORTED"),
                "reviewed": True,
                "release": "2026_03",
                "isoform_aliases": [],
            }
        )
        + "\n",
        encoding="utf-8",
    )
    (assets / "P00533.pdbqt").write_text("ATOM\n", encoding="utf-8")
    (assets / "P00533.a3m").write_text(">P00533\nMPEPTIDE\n", encoding="utf-8")
    (assets / "P00533.pdb").write_text("ATOM\n", encoding="utf-8")
    (assets / "P00533.calibration.json").write_text("{}\n", encoding="utf-8")
    targets = root / "human_canonical_targets.jsonl"
    targets.write_text(
        json.dumps(
            {
                "schema_version": "1.0",
                "target_id": "P00533",
                "gene_symbol": "EGFR",
                "family": "kinase",
                "status": "ready",
                "unsupported_reason": None,
                "sequence": sequence,
                "sequence_sha256": _sha256(sequence),
                "model_sequence": sequence,
                "model_sequence_sha256": _sha256(sequence),
                "model_sequence_start": 1,
                "model_sequence_end": len(sequence),
                "structure_quality": 0.9,
                "structure_id": "1M17",
                "structure_source": "pdb",
                "structure_path": "targets/P00533.pdb",
                "calibration_path": "targets/P00533.calibration.json",
                "pocket_id": "P00533:p1",
                "receptor_pdbqt_path": "targets/P00533.pdbqt",
                "box": {
                    "center": [0.0, 0.0, 0.0],
                    "size": [18.0, 18.0, 18.0],
                },
                "background_affinities": [
                    -5.0 - index / 100 for index in range(100)
                ],
                "msa_path": "targets/P00533.a3m",
                "msa_database_version": "uniref30-2025-06",
                "pocket_residues": [1, 2, 3],
                "model_pocket_residues": [1, 2, 3],
            }
        )
        + "\n"
        + json.dumps(
            {
                "schema_version": "1.0",
                "target_id": "P0UNSP",
                "gene_symbol": None,
                "family": "unknown",
                "status": "unsupported",
                "unsupported_reason": "no_structure",
            }
        )
        + "\n",
        encoding="utf-8",
    )
    return proteome, targets


def test_compile_reference_bundle_hashes_all_runtime_assets(tmp_path: Path) -> None:
    root = tmp_path / "reference"
    proteome, targets = _write_reference_fixture(root)

    summary = compile_reference_bundle(
        root=root,
        proteome_manifest=proteome,
        target_manifest=targets,
        output_manifest=root / "reference_manifest.json",
    )

    payload = json.loads(summary.manifest_path.read_text(encoding="utf-8"))
    assert summary.target_count == 2
    assert summary.ready_target_count == 1
    assert summary.unsupported_target_count == 1
    assert summary.failed_target_count == 0
    assert payload["taxonomy_id"] == 9606
    assert payload["uniprot_release"] == "2026_03"
    assert {item["path"] for item in payload["artifacts"]} == {
        "human_canonical_proteome.jsonl",
        "human_canonical_targets.jsonl",
        "targets/P00533.a3m",
        "targets/P00533.calibration.json",
        "targets/P00533.pdb",
        "targets/P00533.pdbqt",
    }
    assert _reference_manifest_status(summary.manifest_path) == "valid"


def test_reference_bundle_rejects_missing_canonical_target(tmp_path: Path) -> None:
    root = tmp_path / "reference"
    proteome, targets = _write_reference_fixture(root)
    targets.write_text(targets.read_text().splitlines()[0] + "\n", encoding="utf-8")

    with pytest.raises(ValueError, match="canonical target coverage mismatch"):
        compile_reference_bundle(
            root=root,
            proteome_manifest=proteome,
            target_manifest=targets,
            output_manifest=root / "reference_manifest.json",
        )


def test_reference_status_detects_artifact_tampering(tmp_path: Path) -> None:
    root = tmp_path / "reference"
    proteome, targets = _write_reference_fixture(root)
    summary = compile_reference_bundle(
        root=root,
        proteome_manifest=proteome,
        target_manifest=targets,
        output_manifest=root / "reference_manifest.json",
    )
    (root / "targets/P00533.a3m").write_text("tampered\n", encoding="utf-8")

    assert _reference_manifest_status(summary.manifest_path) == "invalid"


def test_compile_reference_cli_writes_machine_readable_summary(tmp_path: Path) -> None:
    root = tmp_path / "reference"
    proteome, targets = _write_reference_fixture(root)
    output = root / "reference_manifest.json"

    result = CliRunner().invoke(
        app,
        [
            "targets",
            "compile-reference",
            "--root",
            str(root),
            "--proteome",
            str(proteome),
            "--targets",
            str(targets),
            "--output",
            str(output),
        ],
    )

    assert result.exit_code == 0, result.stdout
    assert output.is_file()
    assert '"target_count": 2' in result.stdout


def test_build_pilot_manifest_preserves_full_canonical_coverage(tmp_path: Path) -> None:
    root = tmp_path / "reference"
    proteome, _targets = _write_reference_fixture(root)
    structure = root / "targets/1M17_ligand.pdb"
    structure.write_text(
        "DBREF  TEST A    1     2  UNP    P00533   EGFR_HUMAN       2      3\n"
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C\n"
        "ATOM      2  CA  GLY A   2      20.000  20.000  20.000  1.00 20.00           C\n"
        "HETATM    3  C1  LIG A 999       1.000   0.000   0.000  1.00 20.00           C\n",
        encoding="utf-8",
    )
    receptor = root / "targets/P00533.pdbqt"
    calibration = root / "targets/calibration.json"
    calibration.write_text(
        json.dumps(
            {
                "successful_probe_count": 100,
                "background_affinities": [-5.0 - index / 100 for index in range(100)],
            }
        ),
        encoding="utf-8",
    )

    summary = build_pilot_target_manifest(
        root=root,
        proteome_manifest=proteome,
        target_id="P00533",
        family="kinase",
        structure_id="1M17",
        structure_source="pdb",
        structure_pdb=structure,
        ligand_resname="LIG",
        receptor_pdbqt=receptor,
        calibration_json=calibration,
        box={"center": [0, 0, 0], "size": [18, 18, 18]},
        structure_quality=0.9,
        msa_database_version="single-sequence-pilot",
        output_manifest=root / "pilot_targets.jsonl",
    )

    rows = [json.loads(line) for line in summary.manifest_path.read_text().splitlines()]
    ready = next(row for row in rows if row["target_id"] == "P00533")
    unsupported = next(row for row in rows if row["target_id"] == "P0UNSP")
    assert summary.target_count == 2
    assert ready["status"] == "ready"
    assert ready["gene_symbol"] == "EGFR"
    assert ready["pocket_residues"] == [2]
    assert ready["model_pocket_residues"] == [1]
    assert ready["model_sequence_start"] == 2
    assert ready["model_sequence_end"] == 3
    assert ready["model_sequence"] == "PE"
    assert (root / ready["msa_path"]).is_file()
    assert unsupported["status"] == "unsupported"
    assert unsupported["unsupported_reason"] == "pilot_not_built"
