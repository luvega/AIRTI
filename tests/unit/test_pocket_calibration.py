import csv
import json
from pathlib import Path

import pytest

from airti_tf.pockets.receptor import DockingBox
from airti_tf.screening.calibration_build import calibrate_pocket_background
from airti_tf.screening.quickvina import DockingSeedResult


def _write_panel(path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "source_id",
                "canonical_smiles",
                "molecular_weight",
                "clogp",
                "hbd",
                "hba",
                "rotatable_bonds",
                "formal_charge",
                "mw_stratum",
            ),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {
                "source_id": "probe-1",
                "canonical_smiles": "CC(=O)Oc1ccccc1C(=O)O",
                "molecular_weight": "180.2",
                "clogp": "1.2",
                "hbd": "1",
                "hba": "4",
                "rotatable_bonds": "2",
                "formal_charge": "0",
                "mw_stratum": "0",
            }
        )
        writer.writerow(
            {
                "source_id": "probe-2",
                "canonical_smiles": "Cn1c(=O)c2c(ncn2C)n(C)c1=O",
                "molecular_weight": "194.2",
                "clogp": "-0.1",
                "hbd": "0",
                "hba": "6",
                "rotatable_bonds": "0",
                "formal_charge": "0",
                "mw_stratum": "0",
            }
        )


def test_calibration_uses_best_state_and_requires_successful_probes(
    tmp_path: Path,
) -> None:
    panel = tmp_path / "panel.tsv"
    _write_panel(panel)
    receptor = tmp_path / "receptor.pdbqt"
    receptor.write_text("ATOM\n", encoding="utf-8")
    prepared_paths: list[Path] = []

    def prepare_pdbqt(_input_sdf: Path, output_pdbqt: Path) -> None:
        prepared_paths.append(output_pdbqt)
        output_pdbqt.write_text("ATOM\n", encoding="utf-8")

    def run_seed(job: object, seed: int) -> DockingSeedResult:
        state_index = int(Path(getattr(job, "ligand_pdbqt")).stem.rsplit("-", 1)[-1])
        return DockingSeedResult(
            seed=seed,
            status="succeeded",
            affinity_kcal_mol=-6.0 - state_index - seed / 1000,
            pose_count=1,
        )

    output = tmp_path / "calibration.json"
    summary = calibrate_pocket_background(
        panel,
        receptor_pdbqt=receptor,
        box=DockingBox(center=(0, 0, 0), size=(18, 18, 18)),
        pocket_id="P00533:p1",
        output=output,
        asset_dir=tmp_path / "assets",
        expected_probe_count=2,
        minimum_successful_probes=2,
        workers=2,
        pdbqt_preparer=prepare_pdbqt,
        docking_runner=run_seed,
    )

    payload = json.loads(output.read_text(encoding="utf-8"))
    assert summary.successful_probe_count == 2
    assert summary.failed_probe_count == 0
    assert len(payload["background_affinities"]) == 2
    assert all(record["state_count"] >= 1 for record in payload["probes"])
    assert payload["panel_sha256"]

    first_prepare_count = len(prepared_paths)
    calibrate_pocket_background(
        panel,
        receptor_pdbqt=receptor,
        box=DockingBox(center=(1, 1, 1), size=(18, 18, 18)),
        pocket_id="P00533:p2",
        output=tmp_path / "calibration-2.json",
        asset_dir=tmp_path / "assets",
        expected_probe_count=2,
        minimum_successful_probes=2,
        workers=2,
        pdbqt_preparer=prepare_pdbqt,
        docking_runner=run_seed,
    )
    assert len(prepared_paths) == first_prepare_count


def test_calibration_fails_closed_below_minimum_probe_count(tmp_path: Path) -> None:
    panel = tmp_path / "panel.tsv"
    _write_panel(panel)
    receptor = tmp_path / "receptor.pdbqt"
    receptor.write_text("ATOM\n", encoding="utf-8")

    def fail_seed(_job: object, seed: int) -> DockingSeedResult:
        return DockingSeedResult(seed=seed, status="failed", error_code="failed")

    with pytest.raises(ValueError, match="successful background probes"):
        calibrate_pocket_background(
            panel,
            receptor_pdbqt=receptor,
            box=DockingBox(center=(0, 0, 0), size=(18, 18, 18)),
            pocket_id="P00533:p1",
            output=tmp_path / "calibration.json",
            asset_dir=tmp_path / "assets",
            expected_probe_count=2,
            minimum_successful_probes=2,
            workers=1,
            pdbqt_preparer=lambda _input, output: output.write_text("ATOM\n"),
            docking_runner=fail_seed,
        )
