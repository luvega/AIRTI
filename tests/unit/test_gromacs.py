from pathlib import Path

import pytest

from airti_tf.simulation.gromacs import (
    build_grompp_command,
    build_mdrun_command,
    build_system_commands,
    plan_md_replicas,
    render_production_mdp,
    validate_parameterization,
)


def test_md_protocol_is_exactly_100_ns() -> None:
    mdp = render_production_mdp()

    assert mdp["dt"] == pytest.approx(0.002)
    assert mdp["nsteps"] == 50_000_000
    assert mdp["nstxout-compressed"] == 5_000
    assert mdp["nstenergy"] == 5_000
    assert mdp["ref-t"] == 300
    assert mdp["ref-p"] == 1.0


def test_resume_uses_checkpoint_when_present(tmp_path: Path) -> None:
    checkpoint = tmp_path / "md.cpt"
    checkpoint.write_bytes(b"checkpoint")

    command = build_mdrun_command(tmp_path)

    assert command[command.index("-cpi") + 1] == str(checkpoint)
    assert "-append" in command


def test_fresh_run_does_not_reference_missing_checkpoint(tmp_path: Path) -> None:
    command = build_mdrun_command(tmp_path)

    assert "-cpi" not in command
    assert "-append" not in command


def test_grompp_uses_npt_coordinates_and_fixed_topology(tmp_path: Path) -> None:
    command = build_grompp_command(tmp_path)

    assert command == [
        "gmx",
        "grompp",
        "-f",
        str(tmp_path / "md.mdp"),
        "-c",
        str(tmp_path / "npt.gro"),
        "-p",
        str(tmp_path / "topol.top"),
        "-o",
        str(tmp_path / "md.tpr"),
    ]


def test_system_build_uses_dodecahedron_tip3p_and_015m_salt() -> None:
    commands = build_system_commands(Path("complex.pdb"), Path("run"))
    flattened = [token for command in commands for token in command]

    assert "dodecahedron" in flattened
    assert "1.0" in flattened
    assert "tip3p" in flattened
    assert "0.15" in flattened


def test_top3_receive_two_additional_independent_replicas() -> None:
    targets = [f"P{index:05}" for index in range(1, 11)]
    plan = plan_md_replicas(targets)

    assert len(plan) == 16
    assert [item.replica for item in plan if item.target_id == "P00001"] == [1, 2, 3]
    assert [item.replica for item in plan if item.target_id == "P00004"] == [1]
    assert len({item.velocity_seed for item in plan if item.target_id == "P00001"}) == 3


def test_parameterization_failure_has_no_fallback_forcefield() -> None:
    result = validate_parameterization(antechamber_ok=False, parmed_ok=True)

    assert result.status == "failed"
    assert result.error_code == "ligand_parameterization_failed"
    assert result.forcefield is None

