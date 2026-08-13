import subprocess
from pathlib import Path

import pytest

from airti_tf.simulation.gromacs import (
    build_grompp_command,
    build_md_system,
    build_mdrun_command,
    build_system_commands,
    extract_boltz_complex,
    plan_md_replicas,
    prepare_ligand_for_amber,
    render_production_mdp,
    render_system_inputs,
    validate_parameterization,
)


MINIMAL_BOLTZ_CIF = """data_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 N N ALA Axp . 0.0 0.0 0.0 1.0 90.0 1 ALA A N 1
ATOM 2 C CA ALA Axp . 1.4 0.0 0.0 1.0 90.0 1 ALA A CA 1
HETATM 3 C C1 LIG Bxp . 2.0 2.0 2.0 1.0 90.0 1 LIG B C1 1
HETATM 4 O O1 LIG Bxp . 3.2 2.0 2.0 1.0 90.0 1 LIG B O1 1
#
"""


def test_md_protocol_is_exactly_100_ns() -> None:
    mdp = render_production_mdp()

    assert mdp["dt"] == pytest.approx(0.002)
    assert mdp["nsteps"] == 50_000_000
    assert mdp["nstxout-compressed"] == 5_000
    assert mdp["nstenergy"] == 5_000
    assert mdp["ref-t"] == 300
    assert mdp["ref-p"] == 1.0
    assert mdp["tc-grps"] == "System"
    assert mdp["tau-t"] == pytest.approx(0.1)
    assert mdp["tau-p"] == pytest.approx(2.0)
    assert mdp["pbc"] == "xyz"


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


def test_system_build_uses_real_ambertools_parmed_and_gromacs_commands(
    tmp_path: Path,
) -> None:
    render_system_inputs(tmp_path, velocity_seed=12345, ligand_formal_charge=-1)
    commands = build_system_commands(Path("complex.cif"), tmp_path)
    flattened = [token for command in commands for token in command]

    assert commands[0][0] == "antechamber"
    assert "-c" in commands[0] and "bcc" in commands[0]
    assert "-at" in commands[0] and "gaff2" in commands[0]
    assert any(command[0] == "parmchk2" for command in commands)
    assert any(command[0] == "pdb4amber" for command in commands)
    assert any(command[:2] == ["tleap", "-f"] for command in commands)
    assert any(
        command[0] == "python" and command[1].endswith("convert_parmed.py")
        for command in commands
    )
    assert "--source-complex" not in flattened
    assert "0.15" in flattened
    assert "leaprc.protein.ff19SB" in (tmp_path / "build.leap.in").read_text()
    assert "leaprc.gaff2" in (tmp_path / "build.leap.in").read_text()
    assert "leaprc.water.tip3p" in (tmp_path / "build.leap.in").read_text()
    assert "addIonsRand COM Na+ 1" in (tmp_path / "build.leap.in").read_text()
    assert "addIonsRand COM Cl- 1" in (tmp_path / "build.leap.in").read_text()
    assert "format=\"gromacs\"" in (tmp_path / "convert_parmed.py").read_text()
    assert "[ WAT ]" in (tmp_path / "write_water_index.py").read_text()
    genion = next(command for command in commands if command[:2] == ["gmx", "genion"])
    assert genion[genion.index("-n") + 1].endswith("water.ndx")
    assert genion[genion.index("-pname") + 1] == "NA"
    assert genion[genion.index("-nname") + 1] == "CL"
    assert '"Na+": "NA"' in (tmp_path / "convert_parmed.py").read_text()
    assert '"Cl-": "CL"' in (tmp_path / "convert_parmed.py").read_text()
    assert "ligand_formal_charge = -1" in (
        tmp_path / "convert_parmed.py"
    ).read_text()
    assert "charge_normalization.json" in (
        tmp_path / "convert_parmed.py"
    ).read_text()
    assert "gen-seed = 12345" in (tmp_path / "nvt.mdp").read_text()


def test_boltz_complex_is_split_into_protein_and_single_ligand(tmp_path: Path) -> None:
    complex_cif = tmp_path / "complex.cif"
    complex_cif.write_text(MINIMAL_BOLTZ_CIF, encoding="utf-8")

    summary = extract_boltz_complex(complex_cif, output_dir=tmp_path / "split")

    protein = summary.protein_pdb.read_text()
    ligand = summary.ligand_pdb.read_text()
    assert summary.protein_atom_count == 2
    assert summary.ligand_atom_count == 2
    assert protein.startswith("ATOM")
    assert "HETATM" in ligand
    assert " LIG " in ligand

    ligand_sdf = prepare_ligand_for_amber(
        summary.ligand_pdb,
        ligand_smiles="CO",
        ligand_formal_charge=0,
        output_sdf=tmp_path / "split/ligand.sdf",
    )
    assert ligand_sdf.is_file()
    assert " V2000" in ligand_sdf.read_text()


def test_system_builder_fails_closed_on_ligand_parameterization(
    tmp_path: Path,
) -> None:
    complex_cif = tmp_path / "complex.cif"
    complex_cif.write_text(MINIMAL_BOLTZ_CIF, encoding="utf-8")

    def fail_antechamber(
        command: list[str], cwd: Path, stdin: str | None
    ) -> subprocess.CompletedProcess[str]:
        assert command[0] == "antechamber"
        assert cwd.is_absolute()
        assert stdin is None
        return subprocess.CompletedProcess(command, returncode=1, stdout="", stderr="bad")

    result = build_md_system(
        complex_cif,
        run_dir=tmp_path / "run",
        ligand_smiles="CO",
        ligand_formal_charge=0,
        velocity_seed=12345,
        command_runner=fail_antechamber,
    )

    assert result.status == "failed"
    assert result.error_code == "ligand_parameterization_failed"
    assert result.completed_step == "prepare_ligand_coordinates"


def test_system_builder_requires_every_equilibration_artifact(tmp_path: Path) -> None:
    complex_cif = tmp_path / "complex.cif"
    complex_cif.write_text(MINIMAL_BOLTZ_CIF, encoding="utf-8")
    observed_stdin: list[str | None] = []

    def create_expected_outputs(
        command: list[str], cwd: Path, stdin: str | None
    ) -> subprocess.CompletedProcess[str]:
        observed_stdin.append(stdin)
        if command[0] == "tleap":
            for filename in ("solvated.prmtop", "solvated.inpcrd"):
                (cwd / filename).write_text("generated\n")
        elif command[0] == "python" and command[1].endswith("convert_parmed.py"):
            for filename in (
                "topol.top",
                "solvated.gro",
                "charge_normalization.json",
            ):
                (cwd / filename).write_text("generated\n")
        elif command[0] == "python":
            (cwd / "water.ndx").write_text("[ WAT ]\n1 2 3\n")
        elif command[:2] == ["gmx", "mdrun"]:
            stem = Path(command[command.index("-deffnm") + 1])
            stem.with_suffix(".gro").write_text("generated\n")
            if stem.name in {"nvt", "npt"}:
                stem.with_suffix(".cpt").write_text("generated\n")
        else:
            output = Path(command[command.index("-o") + 1])
            output.write_text("generated\n")
        return subprocess.CompletedProcess(command, returncode=0, stdout="ok", stderr="")

    result = build_md_system(
        complex_cif,
        run_dir=tmp_path / "run",
        ligand_smiles="CO",
        ligand_formal_charge=0,
        velocity_seed=12345,
        command_runner=create_expected_outputs,
    )

    assert result.status == "succeeded"
    assert result.completed_step == "equilibrate_npt"
    assert result.forcefield == "ff19SB+GAFF2/TIP3P"
    assert "WAT\n" in observed_stdin
    assert (tmp_path / "run/system_build_status.json").is_file()


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
