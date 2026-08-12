import pytest

from airti_tf.simulation.analysis import TrajectoryMetrics, analyze_trajectory


def metrics(**overrides: object) -> TrajectoryMetrics:
    payload: dict[str, object] = {
        "frames": 50_001,
        "expected_frames": 50_001,
        "completed_ns": 100.0,
        "energy_continuous": True,
        "pbc_processed": True,
        "protein_ca_rmsd_nm": 0.22,
        "ligand_pocket_rmsd_nm": 0.18,
        "pocket_rmsf_nm": 0.15,
        "hbond_occupancy": 0.65,
        "hydrophobic_occupancy": 0.72,
        "salt_bridge_occupancy": 0.3,
        "ligand_center_distance_nm": 0.4,
        "contact_atom_count": 12.0,
        "late_unbound_fraction": 0.05,
        "mmgbsa_kcal_mol": -32.0,
    }
    payload.update(overrides)
    return TrajectoryMetrics.model_validate(payload)


def test_unfinished_trajectory_has_no_final_md_score() -> None:
    result = analyze_trajectory(
        metrics(frames=20_000, expected_frames=50_001, completed_ns=40.0)
    )

    assert result.status == "failed"
    assert result.md_score is None
    assert result.error_code == "trajectory_incomplete"


def test_complete_bound_trajectory_receives_score() -> None:
    result = analyze_trajectory(metrics())

    assert result.status == "stable"
    assert result.md_score == pytest.approx(0.755, abs=0.001)
    assert result.mmgbsa_kcal_mol == -32.0
    assert result.mmgbsa_role == "interpretive_only"


def test_ligand_persistently_leaving_pocket_is_unstable_but_preserved() -> None:
    result = analyze_trajectory(metrics(late_unbound_fraction=0.85))

    assert result.status == "unstable"
    assert result.md_score is not None
    assert "ligand_departed_pocket" in result.flags


def test_energy_discontinuity_is_technical_failure() -> None:
    result = analyze_trajectory(metrics(energy_continuous=False))

    assert result.status == "failed"
    assert result.error_code == "energy_discontinuity"
    assert result.md_score is None

