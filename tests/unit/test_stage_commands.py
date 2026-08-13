import hashlib
import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

from airti_tf.cli import app
from airti_tf.refinement.boltz2 import BoltzSeedResult
from airti_tf.screening.quickvina import DockingSeedResult
from airti_tf.stages import (
    BoltzRefinementSummary,
    MDBundleSummary,
    MDStageResult,
    ReportBundleSummary,
    ScreenBundleSummary,
    prepare_ligand_bundle,
    refine_boltz_bundle,
    render_report_bundle,
    run_md_bundle,
    screen_ligand_bundle,
)
from airti_tf.state import StateStore


def test_cli_exposes_all_production_stage_commands() -> None:
    result = CliRunner().invoke(app, ["--help"])

    assert result.exit_code == 0
    for command in (
        "prepare-ligands",
        "screen",
        "refine-boltz",
        "run-md",
        "render-report",
    ):
        assert command in result.stdout


def test_prepare_ligand_bundle_writes_portable_manifest_and_pdbqt(
    tmp_path: Path,
) -> None:
    queries = tmp_path / "queries.smi"
    queries.write_text(
        "CC(=O)Oc1ccccc1C(=O)O aspirin\n"
        "Cn1c(=O)c2c(ncn2C)n(C)c1=O caffeine\n",
        encoding="utf-8",
    )
    manifest = tmp_path / "prepared_ligands.jsonl"
    assets = tmp_path / "prepared_ligands"

    def prepare_pdbqt(input_sdf: Path, output_pdbqt: Path) -> None:
        assert input_sdf.is_file()
        output_pdbqt.write_text("REMARK prepared\n", encoding="utf-8")

    summary = prepare_ligand_bundle(
        queries,
        output_manifest=manifest,
        asset_dir=assets,
        profile="production",
        max_molecules=5,
        pdbqt_preparer=prepare_pdbqt,
    )

    rows = [json.loads(line) for line in manifest.read_text().splitlines()]
    assert summary.query_count == 2
    assert summary.failed_query_count == 0
    assert {row["query_id"] for row in rows} == {"aspirin", "caffeine"}
    assert all(row["status"] == "succeeded" for row in rows)
    assert all(not Path(row["pdbqt_path"]).is_absolute() for row in rows)
    assert all((manifest.parent / row["pdbqt_path"]).is_file() for row in rows)
    assert all((manifest.parent / row["sdf_path"]).is_file() for row in rows)


def test_prepare_ligand_bundle_rejects_duplicate_query_ids(tmp_path: Path) -> None:
    queries = tmp_path / "queries.smi"
    queries.write_text("CCO duplicate\nCCN duplicate\n", encoding="utf-8")

    with pytest.raises(ValueError, match="query identifiers must be unique"):
        prepare_ligand_bundle(
            queries,
            output_manifest=tmp_path / "prepared.jsonl",
            asset_dir=tmp_path / "assets",
            profile="local",
            max_molecules=5,
            pdbqt_preparer=lambda _input, _output: None,
        )


def test_screen_bundle_calibrates_ready_targets_and_preserves_coverage(
    tmp_path: Path,
) -> None:
    ligand_assets = tmp_path / "ligands"
    target_assets = tmp_path / "targets"
    ligand_assets.mkdir()
    target_assets.mkdir()
    ligand_pdbqt = ligand_assets / "state-1.pdbqt"
    receptor_pdbqt = target_assets / "P00533.pdbqt"
    target_msa = target_assets / "P00533.a3m"
    target_structure = target_assets / "P00533.pdb"
    target_calibration = target_assets / "P00533.calibration.json"
    ligand_pdbqt.write_text("ATOM      1  C   LIG A   1       0.0 0.0 0.0\n")
    receptor_pdbqt.write_text("ATOM      1  C   REC A   1       0.0 0.0 0.0\n")
    target_msa.write_text(">P00533\nMPEPTIDE\n", encoding="utf-8")
    target_structure.write_text("ATOM\n", encoding="utf-8")
    target_calibration.write_text("{}\n", encoding="utf-8")
    ligand_manifest = tmp_path / "prepared.jsonl"
    ligand_manifest.write_text(
        json.dumps(
            {
                "schema_version": "1.0",
                "query_id": "aspirin",
                "ligand_id": "ligand-1",
                "ligand_state_id": "state-1",
                "canonical_smiles": "CC(=O)Oc1ccccc1C(=O)O",
                "atom_count": 13,
                "formal_charge": 0,
                "status": "succeeded",
                "error_code": None,
                "uncertainty_flags": [],
                "sdf_path": "ligands/state-1.sdf",
                "pdbqt_path": "ligands/state-1.pdbqt",
            }
        )
        + "\n"
    )
    target_manifest = tmp_path / "targets.jsonl"
    ready = {
        "schema_version": "1.0",
        "target_id": "P00533",
        "family": "kinase",
        "status": "ready",
        "unsupported_reason": None,
        "sequence": "MPEPTIDE",
        "sequence_sha256": "a" * 64,
        "model_sequence": "MPEPTIDE",
        "model_sequence_sha256": hashlib.sha256(b"MPEPTIDE").hexdigest(),
        "model_sequence_start": 1,
        "model_sequence_end": 8,
        "structure_quality": 0.9,
        "structure_id": "1M17",
        "structure_source": "pdb",
        "structure_path": "targets/P00533.pdb",
        "calibration_path": "targets/P00533.calibration.json",
        "pocket_id": "P00533:p1",
        "receptor_pdbqt_path": "targets/P00533.pdbqt",
        "box": {"center": [0.0, 0.0, 0.0], "size": [18.0, 18.0, 18.0]},
        "background_affinities": [-5.0 - index / 100 for index in range(100)],
        "msa_path": "targets/P00533.a3m",
        "msa_database_version": "pilot-v1",
        "pocket_residues": [1, 2, 3],
        "model_pocket_residues": [1, 2, 3],
    }
    unsupported = {
        "schema_version": "1.0",
        "target_id": "P0UNSP",
        "family": "unknown",
        "status": "unsupported",
        "unsupported_reason": "no_structure",
    }
    target_manifest.write_text(
        json.dumps(ready) + "\n" + json.dumps(unsupported) + "\n",
        encoding="utf-8",
    )

    def run_seed(job: object, seed: int) -> DockingSeedResult:
        del job
        affinity = {11: -8.0, 29: -8.2, 47: -7.9}[seed]
        return DockingSeedResult(
            seed=seed,
            status="succeeded",
            affinity_kcal_mol=affinity,
            pose_count=1,
        )

    output = tmp_path / "screened.jsonl"
    summary = screen_ligand_bundle(
        ligand_manifest,
        target_manifest,
        output_manifest=output,
        asset_dir=tmp_path / "screen-assets",
        top_n=300,
        docking_runner=run_seed,
    )

    rows = [json.loads(line) for line in output.read_text().splitlines()]
    assert summary.query_count == 1
    assert summary.ready_target_count == 1
    assert summary.unsupported_target_count == 1
    assert summary.failed_docking_job_count == 0
    assert len(rows) == 1
    assert rows[0]["query_id"] == "aspirin"
    assert rows[0]["target_id"] == "P00533"
    assert rows[0]["ligand_formal_charge"] == 0
    assert rows[0]["seed_success_count"] == 3
    assert 0.0 <= rows[0]["calibrated_score"] <= 1.0
    assert rows[0]["target_coverage"] == {
        "total": 2,
        "ready": 1,
        "unsupported": 1,
        "failed": 0,
    }
    assert Path(rows[0]["msa_path"]).parts[0] == "screen-assets"
    assert (output.parent / rows[0]["msa_path"]).read_text() == target_msa.read_text()
    assert "P0UNSP" not in {row["target_id"] for row in rows}


def test_screen_cli_forwards_file_contract(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    ligands = tmp_path / "ligands.jsonl"
    targets = tmp_path / "targets.jsonl"
    ligands.write_text("{}\n")
    targets.write_text("{}\n")
    output = tmp_path / "screened.jsonl"
    observed: dict[str, object] = {}

    def fake_screen(
        ligand_manifest: Path,
        target_manifest: Path,
        *,
        output_manifest: Path,
        asset_dir: Path,
        top_n: int,
    ) -> ScreenBundleSummary:
        observed.update(
            {
                "ligands": ligand_manifest,
                "targets": target_manifest,
                "output": output_manifest,
                "asset_dir": asset_dir,
                "top_n": top_n,
            }
        )
        output_manifest.write_text("", encoding="utf-8")
        return ScreenBundleSummary(
            query_count=1,
            ready_target_count=1,
            unsupported_target_count=0,
            failed_target_count=0,
            successful_docking_job_count=1,
            failed_docking_job_count=0,
            candidate_count=1,
            manifest_path=output_manifest,
            manifest_sha256="a" * 64,
        )

    monkeypatch.setattr("airti_tf.stages.screen_ligand_bundle", fake_screen)
    result = CliRunner().invoke(
        app,
        [
            "screen",
            "--ligands",
            str(ligands),
            "--targets",
            str(targets),
            "--output",
            str(output),
            "--asset-dir",
            str(tmp_path / "screen-assets"),
            "--top-n",
            "25",
        ],
    )

    assert result.exit_code == 0, result.stdout
    assert observed["ligands"] == ligands
    assert observed["targets"] == targets
    assert observed["output"] == output
    assert observed["top_n"] == 25


def test_refine_boltz_bundle_writes_three_seed_consensus_and_structure(
    tmp_path: Path,
) -> None:
    msa = tmp_path / "P00533.a3m"
    msa.write_text(">query\nMPEPTIDE\n", encoding="utf-8")
    screened = tmp_path / "screened.jsonl"
    screened.write_text(
        json.dumps(
            {
                "schema_version": "1.0",
                "query_id": "aspirin",
                "ligand_id": "ligand-1",
                "ligand_state_id": "state-1",
                "ligand_smiles": "CC(=O)Oc1ccccc1C(=O)O",
                "ligand_atom_count": 13,
                "ligand_formal_charge": 0,
                "target_id": "P00533",
                "family": "kinase",
                "pocket_id": "P00533:p1",
                "pocket_residues": [1, 2, 3],
                "model_pocket_residues": [1, 2, 3],
                "sequence": "MPEPTIDE",
                "sequence_sha256": "a" * 64,
                "model_sequence": "MPEPTIDE",
                "model_sequence_sha256": hashlib.sha256(b"MPEPTIDE").hexdigest(),
                "model_sequence_start": 1,
                "model_sequence_end": 8,
                "msa_path": "P00533.a3m",
                "msa_database_version": "pilot-v1",
                "structure_quality": 0.9,
                "affinity_median": -8.1,
                "calibrated_score": 0.95,
                "seed_range": 0.3,
                "pose_consistency": 0.8,
                "seed_success_count": 3,
                "selection_reason": "global_primary",
                "screen_rank": 1,
                "target_coverage": {
                    "total": 2,
                    "ready": 1,
                    "unsupported": 1,
                    "failed": 0,
                },
            }
        )
        + "\n",
        encoding="utf-8",
    )

    def run_seed(job: object, seed: int) -> BoltzSeedResult:
        output_dir = Path(getattr(job, "output_dir"))
        structure = output_dir / f"model-{seed}.cif"
        structure.parent.mkdir(parents=True, exist_ok=True)
        structure.write_text("data_model\n", encoding="utf-8")
        return BoltzSeedResult(
            seed=seed,
            status="succeeded",
            confidence_score={11: 0.7, 29: 0.8, 47: 0.9}[seed],
            ligand_iptm=0.75,
            affinity_probability=0.85,
            affinity_pred_value=-1.2,
            pocket_constraint_fraction=0.9,
            severe_clash=False,
            structure_path=structure,
        )

    output = tmp_path / "refined.jsonl"
    summary = refine_boltz_bundle(
        screened,
        output_manifest=output,
        asset_dir=tmp_path / "boltz-assets",
        profile="production",
        top_n=30,
        boltz_runner=run_seed,
    )

    rows = [json.loads(line) for line in output.read_text().splitlines()]
    assert summary.query_count == 1
    assert summary.selected_candidate_count == 1
    assert summary.succeeded_candidate_count == 1
    assert summary.failed_candidate_count == 0
    assert len(rows) == 1
    assert rows[0]["boltz_status"] == "succeeded"
    assert rows[0]["boltz_seed_success_count"] == 3
    assert rows[0]["boltz_confidence_median"] == pytest.approx(0.8)
    assert 0.0 <= rows[0]["boltz_score"] <= 1.0
    assert not Path(rows[0]["boltz_structure_path"]).is_absolute()
    assert (output.parent / rows[0]["boltz_structure_path"]).is_file()


def test_refine_boltz_cli_forwards_file_contract(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    candidates = tmp_path / "screened.jsonl"
    candidates.write_text("{}\n")
    output = tmp_path / "refined.jsonl"
    observed: dict[str, object] = {}

    def fake_refine(
        screen_manifest: Path,
        *,
        output_manifest: Path,
        asset_dir: Path,
        profile: str,
        top_n: int,
        cache_path: Path,
    ) -> BoltzRefinementSummary:
        observed.update(
            {
                "candidates": screen_manifest,
                "output": output_manifest,
                "asset_dir": asset_dir,
                "profile": profile,
                "top_n": top_n,
                "cache_path": cache_path,
            }
        )
        output_manifest.write_text("", encoding="utf-8")
        return BoltzRefinementSummary(
            query_count=1,
            selected_candidate_count=1,
            succeeded_candidate_count=1,
            failed_candidate_count=0,
            manifest_path=output_manifest,
            manifest_sha256="a" * 64,
        )

    monkeypatch.setattr("airti_tf.stages.refine_boltz_bundle", fake_refine)
    result = CliRunner().invoke(
        app,
        [
            "refine-boltz",
            "--candidates",
            str(candidates),
            "--output",
            str(output),
            "--asset-dir",
            str(tmp_path / "boltz-assets"),
            "--profile",
            "production",
            "--top-n",
            "20",
            "--cache",
            str(tmp_path / "boltz-cache"),
        ],
    )

    assert result.exit_code == 0, result.stdout
    assert observed["candidates"] == candidates
    assert observed["profile"] == "production"
    assert observed["top_n"] == 20
    assert observed["cache_path"] == tmp_path / "boltz-cache"


def test_run_md_bundle_routes_top_candidate_to_three_replicas(tmp_path: Path) -> None:
    structure = tmp_path / "model.cif"
    structure.write_text("data_model\n", encoding="utf-8")
    refined = tmp_path / "refined.jsonl"
    refined.write_text(
        json.dumps(
            {
                "schema_version": "1.0",
                "query_id": "aspirin",
                "ligand_id": "ligand-1",
                "ligand_state_id": "state-1",
                "ligand_smiles": "CC(=O)Oc1ccccc1C(=O)O",
                "ligand_atom_count": 13,
                "ligand_formal_charge": 0,
                "target_id": "P00533",
                "family": "kinase",
                "pocket_id": "P00533:p1",
                "pocket_residues": [1, 2, 3],
                "model_pocket_residues": [1, 2, 3],
                "sequence": "MPEPTIDE",
                "sequence_sha256": "a" * 64,
                "model_sequence": "MPEPTIDE",
                "model_sequence_sha256": hashlib.sha256(b"MPEPTIDE").hexdigest(),
                "model_sequence_start": 1,
                "model_sequence_end": 8,
                "msa_path": "P00533.a3m",
                "msa_database_version": "pilot-v1",
                "structure_quality": 0.9,
                "affinity_median": -8.1,
                "calibrated_score": 0.95,
                "seed_range": 0.3,
                "pose_consistency": 0.8,
                "seed_success_count": 3,
                "selection_reason": "global_primary",
                "screen_rank": 1,
                "target_coverage": {
                    "total": 2,
                    "ready": 1,
                    "unsupported": 1,
                    "failed": 0,
                },
                "boltz_status": "succeeded",
                "boltz_error_code": None,
                "boltz_seed_errors": {},
                "boltz_seed_success_count": 3,
                "boltz_confidence_median": 0.8,
                "boltz_ligand_iptm_median": 0.75,
                "boltz_affinity_probability_median": 0.85,
                "boltz_affinity_pred_value_median": -1.2,
                "boltz_pocket_constraint_median": 0.9,
                "boltz_confidence_range": 0.2,
                "boltz_score": 0.825,
                "boltz_structure_path": "model.cif",
            }
        )
        + "\n",
        encoding="utf-8",
    )

    def run_replica(candidate: object, run_dir: Path, replica: int) -> MDStageResult:
        del candidate
        run_dir.mkdir(parents=True, exist_ok=True)
        trajectory = run_dir / "md.xtc"
        trajectory.write_bytes(f"replica-{replica}".encode())
        checkpoint = run_dir / "md.cpt"
        checkpoint.write_bytes(b"checkpoint")
        return MDStageResult(
            replica=replica,
            status="stable",
            completed_ns=100.0,
            md_score=0.8 + replica / 100,
            error_code=None,
            trajectory_path=trajectory,
            checkpoint_path=checkpoint,
        )

    output = tmp_path / "md.jsonl"
    summary = run_md_bundle(
        refined,
        output_manifest=output,
        asset_dir=tmp_path / "md-assets",
        top_n=10,
        md_runner=run_replica,
    )

    rows = [json.loads(line) for line in output.read_text().splitlines()]
    assert summary.query_count == 1
    assert summary.selected_candidate_count == 1
    assert summary.planned_replica_count == 3
    assert summary.succeeded_candidate_count == 1
    assert len(rows) == 1
    assert rows[0]["md_status"] == "stable"
    assert rows[0]["md_replica_success_count"] == 3
    assert rows[0]["completed_ns"] == 100.0
    assert rows[0]["md_score"] == pytest.approx(0.82)
    assert len(rows[0]["md_replicas"]) == 3


def test_run_md_cli_forwards_file_contract(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    candidates = tmp_path / "refined.jsonl"
    candidates.write_text("{}\n")
    output = tmp_path / "md.jsonl"
    observed: dict[str, object] = {}

    def fake_md(
        boltz_manifest: Path,
        *,
        output_manifest: Path,
        asset_dir: Path,
        top_n: int,
    ) -> MDBundleSummary:
        observed.update(
            {
                "candidates": boltz_manifest,
                "output": output_manifest,
                "asset_dir": asset_dir,
                "top_n": top_n,
            }
        )
        output_manifest.write_text("", encoding="utf-8")
        return MDBundleSummary(
            query_count=1,
            selected_candidate_count=1,
            planned_replica_count=3,
            succeeded_candidate_count=1,
            failed_candidate_count=0,
            manifest_path=output_manifest,
            manifest_sha256="a" * 64,
        )

    monkeypatch.setattr("airti_tf.stages.run_md_bundle", fake_md)
    result = CliRunner().invoke(
        app,
        [
            "run-md",
            "--candidates",
            str(candidates),
            "--output",
            str(output),
            "--asset-dir",
            str(tmp_path / "md-assets"),
            "--top-n",
            "8",
        ],
    )

    assert result.exit_code == 0, result.stdout
    assert observed["candidates"] == candidates
    assert observed["top_n"] == 8


def test_render_report_bundle_uses_carried_proteome_coverage(tmp_path: Path) -> None:
    candidates = tmp_path / "md.jsonl"
    candidates.write_text(
        json.dumps(
            {
                "schema_version": "1.0",
                "query_id": "aspirin",
                "ligand_id": "ligand-1",
                "ligand_state_id": "state-1",
                "ligand_smiles": "CC(=O)Oc1ccccc1C(=O)O",
                "ligand_atom_count": 13,
                "ligand_formal_charge": 0,
                "target_id": "P00533",
                "family": "kinase",
                "pocket_id": "P00533:p1",
                "pocket_residues": [1, 2, 3],
                "model_pocket_residues": [1, 2, 3],
                "sequence": "MPEPTIDE",
                "sequence_sha256": "a" * 64,
                "model_sequence": "MPEPTIDE",
                "model_sequence_sha256": hashlib.sha256(b"MPEPTIDE").hexdigest(),
                "model_sequence_start": 1,
                "model_sequence_end": 8,
                "msa_path": "P00533.a3m",
                "msa_database_version": "pilot-v1",
                "structure_quality": 0.9,
                "affinity_median": -8.1,
                "calibrated_score": 0.95,
                "seed_range": 0.3,
                "pose_consistency": 0.8,
                "seed_success_count": 3,
                "selection_reason": "global_primary",
                "screen_rank": 1,
                "target_coverage": {
                    "total": 20435,
                    "ready": 16800,
                    "unsupported": 3600,
                    "failed": 35,
                },
                "boltz_status": "succeeded",
                "boltz_error_code": None,
                "boltz_seed_errors": {},
                "boltz_seed_success_count": 3,
                "boltz_confidence_median": 0.8,
                "boltz_ligand_iptm_median": 0.75,
                "boltz_affinity_probability_median": 0.85,
                "boltz_affinity_pred_value_median": -1.2,
                "boltz_pocket_constraint_median": 0.9,
                "boltz_confidence_range": 0.2,
                "boltz_score": 0.825,
                "boltz_structure_path": "model.cif",
                "md_rank": 1,
                "md_status": "stable",
                "md_error_code": None,
                "md_replica_success_count": 3,
                "completed_ns": 100.0,
                "md_score": 0.82,
                "md_replicas": [],
            }
        )
        + "\n",
        encoding="utf-8",
    )
    output_dir = tmp_path / "final_report"
    state_db = tmp_path / "job_status.sqlite"

    summary = render_report_bundle(
        candidates,
        output_dir=output_dir,
        state_db=state_db,
        project_id="AIRTI-PILOT-001",
    )

    report = (output_dir / "report.md").read_text(encoding="utf-8")
    assert summary.query_count == 1
    assert summary.reported_target_count == 1
    assert summary.coverage.total == 20435
    assert "aspirin" in report
    assert "P00533" in report
    assert "20435" in report
    assert "failed 靶点" in report
    assert StateStore(state_db).task_count() == 1
    assert StateStore(state_db).artifact_count() == 2


def test_render_report_cli_forwards_file_contract(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    candidates = tmp_path / "md.jsonl"
    candidates.write_text("{}\n")
    output_dir = tmp_path / "report"
    state_db = tmp_path / "state.sqlite"
    observed: dict[str, object] = {}

    def fake_report(
        md_manifest: Path,
        *,
        output_dir: Path,
        state_db: Path,
        project_id: str | None,
    ) -> ReportBundleSummary:
        observed.update(
            {
                "candidates": md_manifest,
                "output_dir": output_dir,
                "state_db": state_db,
                "project_id": project_id,
            }
        )
        output_dir.mkdir(parents=True)
        report_path = output_dir / "report.md"
        manifest_path = output_dir / "report_manifest.json"
        report_path.write_text("# report\n")
        manifest_path.write_text("{}\n")
        return ReportBundleSummary(
            project_id="AIRTI-PILOT-001",
            query_count=1,
            reported_target_count=1,
            coverage={"total": 2, "ready": 1, "unsupported": 1, "failed": 0},
            report_path=report_path,
            manifest_path=manifest_path,
            state_db=state_db,
            report_sha256="a" * 64,
        )

    monkeypatch.setattr("airti_tf.stages.render_report_bundle", fake_report)
    result = CliRunner().invoke(
        app,
        [
            "render-report",
            "--candidates",
            str(candidates),
            "--output",
            str(output_dir),
            "--state-db",
            str(state_db),
            "--project-id",
            "AIRTI-PILOT-001",
        ],
    )

    assert result.exit_code == 0, result.stdout
    assert observed["candidates"] == candidates
    assert observed["output_dir"] == output_dir
    assert observed["state_db"] == state_db
    assert observed["project_id"] == "AIRTI-PILOT-001"
