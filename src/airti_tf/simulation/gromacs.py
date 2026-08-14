"""Auditable AmberTools/ParmEd/GROMACS molecular-dynamics protocol."""

from __future__ import annotations

import hashlib
import json
import re
import subprocess
from collections.abc import Callable
from pathlib import Path
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field

from airti_tf.manifest_io import write_artifact, write_bytes_atomic
from airti_tf.simulation.analysis import (
    TrajectoryAnalysis,
    analyze_trajectory,
    measure_trajectory,
)

CommandRunner = Callable[
    [list[str], Path, str | None], subprocess.CompletedProcess[str]
]
MDProtocol = Literal["smoke", "production"]
MDSystemKind = Literal["soluble", "membrane"]


class MDReplicaPlan(BaseModel):
    """One independently seeded production trajectory."""

    model_config = ConfigDict(extra="forbid")

    target_id: str
    replica: int
    velocity_seed: int


class ParameterizationResult(BaseModel):
    """Strict ligand/system parameterization gate."""

    model_config = ConfigDict(extra="forbid")

    status: Literal["succeeded", "failed"]
    error_code: str | None = None
    forcefield: str | None = None


class ComplexComponents(BaseModel):
    """Protein and ligand coordinates extracted from a Boltz complex."""

    model_config = ConfigDict(extra="forbid")

    protein_pdb: Path
    ligand_pdb: Path
    cofactor_pdb: Path | None = None
    protein_atom_count: int
    ligand_atom_count: int
    cofactor_atom_count: int = 0


class MDSystemBuildResult(BaseModel):
    """Fail-closed result of system construction and equilibration."""

    model_config = ConfigDict(extra="forbid")

    status: Literal["succeeded", "failed"]
    error_code: str | None = None
    completed_step: str | None = None
    forcefield: str | None = None
    run_dir: Path
    protein_atom_count: int = 0
    ligand_atom_count: int = 0
    command_logs: list[Path] = Field(default_factory=list)


class MDTrajectoryResult(BaseModel):
    """Production/smoke trajectory execution and conservative analysis result."""

    model_config = ConfigDict(extra="forbid")

    status: Literal["stable", "unstable", "failed"]
    completed_ns: float = Field(ge=0)
    md_score: float | None = Field(default=None, ge=0, le=1)
    error_code: str | None = None
    trajectory_path: Path | None = None
    checkpoint_path: Path | None = None
    metrics_path: Path | None = None


def _column(
    data: dict[str, object], key: str, *, count: int, default: str
) -> list[str]:
    raw = data.get(key)
    if raw is None:
        return [default] * count
    values = raw if isinstance(raw, list) else [raw]
    if len(values) != count:
        raise ValueError(f"inconsistent mmCIF atom-site column: {key}")
    return [str(value) for value in values]


def _pdb_atom_line(
    *,
    group: str,
    serial: int,
    atom_name: str,
    residue_name: str,
    chain: str,
    residue_number: int,
    x: float,
    y: float,
    z: float,
    occupancy: float,
    b_factor: float,
    element: str,
) -> str:
    normalized_atom = atom_name[:4]
    atom_field = (
        normalized_atom.ljust(4)
        if len(normalized_atom) == 4
        else f" {normalized_atom:<3}"
    )
    return (
        f"{group:<6}{serial:>5} {atom_field} {residue_name[:3]:>3} "
        f"{chain[:1]:1}{residue_number:>4}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}"
        f"{occupancy:>6.2f}{b_factor:>6.2f}          {element[:2]:>2}\n"
    )


def extract_boltz_complex(
    complex_cif: Path,
    *,
    output_dir: Path,
) -> ComplexComponents:
    """Split the fixed Boltz chain-A/chain-B complex into Amber inputs."""
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict

    if not complex_cif.is_file() or complex_cif.stat().st_size == 0:
        raise ValueError(f"Boltz complex is missing or empty: {complex_cif}")
    raw_data = MMCIF2Dict(str(complex_cif))
    data: dict[str, object] = dict(raw_data)
    groups_raw = data.get("_atom_site.group_PDB")
    groups = groups_raw if isinstance(groups_raw, list) else [groups_raw]
    if not groups or groups == [None]:
        raise ValueError("Boltz complex contains no atom_site records")
    count = len(groups)
    group_values = [str(value) for value in groups]
    elements = _column(data, "_atom_site.type_symbol", count=count, default="")
    atom_names = _column(data, "_atom_site.label_atom_id", count=count, default="X")
    residues = _column(data, "_atom_site.label_comp_id", count=count, default="UNK")
    label_chains = _column(data, "_atom_site.label_asym_id", count=count, default="")
    auth_chains = _column(data, "_atom_site.auth_asym_id", count=count, default="")
    label_sequence_ids = _column(
        data, "_atom_site.label_seq_id", count=count, default="."
    )
    auth_sequence_ids = _column(
        data, "_atom_site.auth_seq_id", count=count, default="1"
    )
    xs = _column(data, "_atom_site.Cartn_x", count=count, default="nan")
    ys = _column(data, "_atom_site.Cartn_y", count=count, default="nan")
    zs = _column(data, "_atom_site.Cartn_z", count=count, default="nan")
    occupancies = _column(data, "_atom_site.occupancy", count=count, default="1")
    b_factors = _column(
        data, "_atom_site.B_iso_or_equiv", count=count, default="0"
    )
    model_numbers = _column(
        data, "_atom_site.pdbx_PDB_model_num", count=count, default="1"
    )

    protein_lines: list[str] = []
    ligand_lines: list[str] = []
    cofactor_lines: list[str] = []
    ligand_components: set[str] = set()
    waters = {"HOH", "WAT", "TIP3", "SOL"}
    for index in range(count):
        if model_numbers[index] != "1":
            continue
        group = group_values[index]
        chain = (
            auth_chains[index]
            if auth_chains[index] not in {"", ".", "?"}
            else label_chains[index]
        )
        is_protein = group == "ATOM" and chain == "A"
        is_ligand = (
            group == "HETATM"
            and chain == "B"
            and residues[index].upper() not in waters
        )
        is_cofactor = (
            group == "HETATM"
            and chain not in {"A", "B"}
            and residues[index].upper() not in waters
        )
        if not is_protein and not is_ligand and not is_cofactor:
            continue
        try:
            residue_number = (
                int(
                    float(
                        label_sequence_ids[index]
                        if label_sequence_ids[index] not in {"", ".", "?"}
                        else auth_sequence_ids[index]
                    )
                )
                if is_protein
                else 1
            )
            coordinate = (float(xs[index]), float(ys[index]), float(zs[index]))
            occupancy = float(occupancies[index])
            b_factor = float(b_factors[index])
        except ValueError as error:
            raise ValueError("Boltz complex contains invalid atom coordinates") from error
        output_group = "ATOM" if is_protein else "HETATM"
        output_chain = "A" if is_protein else ("B" if is_ligand else "C")
        output_residue = (
            residues[index]
            if is_protein or is_cofactor
            else "LIG"
        )
        line = _pdb_atom_line(
            group=output_group,
            serial=len(protein_lines) + len(ligand_lines) + 1,
            atom_name=atom_names[index],
            residue_name=output_residue,
            chain=output_chain,
            residue_number=residue_number,
            x=coordinate[0],
            y=coordinate[1],
            z=coordinate[2],
            occupancy=occupancy,
            b_factor=b_factor,
            element=elements[index],
        )
        if is_protein:
            protein_lines.append(line)
        elif is_ligand:
            ligand_components.add(residues[index])
            ligand_lines.append(line)
        else:
            cofactor_lines.append(line)
    if not protein_lines:
        raise ValueError("Boltz complex contains no protein atoms in chain A")
    if not ligand_lines:
        raise ValueError("Boltz complex contains no ligand atoms in chain B")
    if len(ligand_components) != 1:
        raise ValueError("Boltz complex must contain exactly one ligand component")

    output_dir.mkdir(parents=True, exist_ok=True)
    protein_path = output_dir / "protein.pdb"
    ligand_path = output_dir / "ligand.pdb"
    cofactor_path = output_dir / "cofactor.pdb"
    write_bytes_atomic(
        protein_path, ("".join(protein_lines) + "TER\nEND\n").encode("utf-8")
    )
    write_bytes_atomic(ligand_path, ("".join(ligand_lines) + "END\n").encode("utf-8"))
    if cofactor_lines:
        write_bytes_atomic(
            cofactor_path,
            ("".join(cofactor_lines) + "END\n").encode("utf-8"),
        )
    return ComplexComponents(
        protein_pdb=protein_path,
        ligand_pdb=ligand_path,
        cofactor_pdb=cofactor_path if cofactor_lines else None,
        protein_atom_count=len(protein_lines),
        ligand_atom_count=len(ligand_lines),
        cofactor_atom_count=len(cofactor_lines),
    )


def align_boltz_complex_to_reference(
    complex_cif: Path,
    *,
    reference_pdb: Path,
    output_cif: Path,
) -> Path:
    """Rigidly align the Boltz protein and all ligands to an OPM reference."""
    from Bio.Align import PairwiseAligner
    from Bio.PDB import MMCIFIO, MMCIFParser, PDBParser, Superimposer
    from Bio.SeqUtils import seq1

    mobile_structure = MMCIFParser(QUIET=True).get_structure(
        "mobile", str(complex_cif)
    )
    reference_structure = PDBParser(QUIET=True).get_structure(
        "reference", str(reference_pdb)
    )

    def ca_atoms(structure: object) -> tuple[str, list[object]]:
        residues = []
        for residue in structure.get_residues():
            if "CA" in residue and residue.id[0] == " ":
                residues.append(residue)
        sequence = "".join(seq1(residue.resname, undef_code="X") for residue in residues)
        return sequence, [residue["CA"] for residue in residues]

    reference_sequence, reference_atoms = ca_atoms(reference_structure)
    mobile_sequence, mobile_atoms = ca_atoms(mobile_structure)
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(reference_sequence, mobile_sequence)[0]
    paired_reference = []
    paired_mobile = []
    for reference_index, mobile_index in zip(
        alignment.indices[0], alignment.indices[1], strict=True
    ):
        if reference_index >= 0 and mobile_index >= 0:
            paired_reference.append(reference_atoms[reference_index])
            paired_mobile.append(mobile_atoms[mobile_index])
    if len(paired_reference) < 3:
        raise ValueError("fewer than three aligned C-alpha atoms for membrane orientation")
    superimposer = Superimposer()
    superimposer.set_atoms(paired_reference, paired_mobile)
    superimposer.apply(list(mobile_structure.get_atoms()))
    output_cif.parent.mkdir(parents=True, exist_ok=True)
    writer = MMCIFIO()
    writer.set_structure(mobile_structure)
    writer.save(str(output_cif))
    if not output_cif.is_file() or output_cif.stat().st_size == 0:
        raise RuntimeError("aligned Boltz complex was not written")
    return output_cif


def prepare_ligand_for_amber(
    ligand_pdb: Path,
    *,
    ligand_smiles: str,
    ligand_formal_charge: int,
    output_sdf: Path,
) -> Path:
    """Restore SMILES bond orders and explicit hydrogens on Boltz coordinates."""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    template = Chem.MolFromSmiles(ligand_smiles)
    observed = Chem.MolFromPDBFile(
        str(ligand_pdb),
        sanitize=False,
        removeHs=True,
        proximityBonding=True,
    )
    if template is None:
        raise ValueError("ligand SMILES cannot be parsed for MD parameterization")
    if observed is None:
        raise ValueError("Boltz ligand coordinates cannot be parsed")
    if template.GetNumHeavyAtoms() != observed.GetNumHeavyAtoms():
        raise ValueError("Boltz ligand atom count does not match the ligand SMILES")
    try:
        assigned = AllChem.AssignBondOrdersFromTemplate(
            template, observed
        )
        Chem.SanitizeMol(assigned)
    except (RuntimeError, ValueError) as error:
        raise ValueError("Boltz ligand topology does not match the ligand SMILES") from error
    if Chem.GetFormalCharge(assigned) != ligand_formal_charge:
        raise ValueError("ligand formal charge does not match the prepared state")
    with_hydrogens = Chem.AddHs(assigned, addCoords=True)
    mol_block = Chem.MolToMolBlock(with_hydrogens)
    write_bytes_atomic(
        output_sdf, (mol_block.rstrip() + "\n$$$$\n").encode("utf-8")
    )
    return output_sdf


def render_md_mdp(
    *, protocol: MDProtocol, system_kind: MDSystemKind
) -> dict[str, str | int | float]:
    """Return the frozen 1 ns smoke or 100 ns production protocol.

    GROMACS uses ps for ``dt``. Therefore 50,000,000 steps at 0.002 ps
    equal exactly 100,000 ps (100 ns). Coordinates and energies are written
    every 5,000 steps, i.e. every 10 ps.
    """
    settings: dict[str, str | int | float] = {
        "integrator": "md",
        "dt": 0.002,
        "nsteps": 500_000 if protocol == "smoke" else 50_000_000,
        "nstxout-compressed": 5_000,
        "nstenergy": 5_000,
        "continuation": "yes",
        "constraints": "h-bonds",
        "constraint-algorithm": "lincs",
        "cutoff-scheme": "Verlet",
        "coulombtype": "PME",
        "rcoulomb": 1.0,
        "vdwtype": "Cut-off",
        "rvdw": 1.0,
        "tcoupl": "V-rescale",
        "tc-grps": "System",
        "tau-t": 0.1,
        "ref-t": 300,
        "pcoupl": "C-rescale",
        "tau-p": 2.0,
        "ref-p": 1.0,
        "compressibility": 4.5e-5,
        "gen-vel": "no",
        "pbc": "xyz",
    }
    if system_kind == "membrane":
        settings.update(
            {
                "pcoupltype": "semiisotropic",
                "ref-p": "1.0 1.0",
                "compressibility": "4.5e-5 4.5e-5",
            }
        )
    return settings


def render_production_mdp() -> dict[str, str | int | float]:
    """Return the fixed soluble 100 ns, 300 K, 1 bar protocol."""
    return render_md_mdp(protocol="production", system_kind="soluble")


def build_membrane_command(
    *,
    complex_pdb: Path,
    output_pdb: Path,
    ligand_frcmod: Path,
    ligand_lib: Path,
    cofactor_parameters: list[tuple[Path, Path]] | None = None,
) -> list[str]:
    """Build the pinned POPC/cholesterol PACKMOL-Memgen command."""
    required = (complex_pdb, ligand_frcmod, ligand_lib)
    if any(not path.is_file() or path.stat().st_size == 0 for path in required):
        raise ValueError("membrane builder input is missing or empty")
    command = [
        "packmol-memgen",
        "-p",
        str(complex_pdb),
        "-o",
        str(output_pdb),
        "--preoriented",
        "--keepligs",
        "-l",
        "POPC:CHL1",
        "-r",
        "4:1",
        "--salt",
        "--salt_c",
        "Na+",
        "--salt_a",
        "Cl-",
        "--saltcon",
        "0.15",
        "--parametrize",
        "--ffprot",
        "ff19SB",
        "--fflip",
        "lipid21",
        "--ffwat",
        "tip3p",
        "--gaff2",
        "--ligand_param",
        f"{ligand_frcmod}:{ligand_lib}",
        "--notprotonate",
        "--noprogress",
        "--overwrite",
    ]
    for frcmod, library in cofactor_parameters or []:
        if any(
            not path.is_file() or path.stat().st_size == 0
            for path in (frcmod, library)
        ):
            raise ValueError("cofactor parameter adapter is missing or empty")
        command.extend(["--ligand_param", f"{frcmod}:{library}"])
    return command


def render_membrane_inputs(
    run_dir: Path,
    *,
    velocity_seed: int,
    protocol: MDProtocol,
    cofactor_parameter_ids: list[str],
    cofactor_parameter_root: Path,
) -> bool:
    """Write membrane MDPs and a fail-closed cofactor-adapter preflight."""
    if not 1 <= velocity_seed <= 2_147_483_647:
        raise ValueError("velocity seed must be a positive signed 32-bit integer")
    run_dir.mkdir(parents=True, exist_ok=True)
    adapters: list[dict[str, object]] = []
    for parameter_id in cofactor_parameter_ids:
        adapter_dir = cofactor_parameter_root / parameter_id
        manifest = adapter_dir / "adapter.json"
        frcmod = adapter_dir / "cofactor.frcmod"
        library = adapter_dir / "cofactor.lib"
        missing = [
            str(path)
            for path in (manifest, frcmod, library)
            if not path.is_file() or path.stat().st_size == 0
        ]
        violations: list[str] = []
        if not missing:
            try:
                payload = json.loads(manifest.read_text(encoding="utf-8"))
                if payload.get("parameter_id") != parameter_id:
                    violations.append("parameter_id_mismatch")
                expected_hashes = {
                    frcmod: payload.get("frcmod_sha256"),
                    library: payload.get("library_sha256"),
                }
                for path, expected_hash in expected_hashes.items():
                    observed_hash = hashlib.sha256(path.read_bytes()).hexdigest()
                    if expected_hash != observed_hash:
                        violations.append(f"sha256_mismatch:{path.name}")
                if not payload.get("source") or not payload.get("chemical_state"):
                    violations.append("provenance_incomplete")
            except (OSError, json.JSONDecodeError, AttributeError):
                violations.append("adapter_manifest_invalid")
        adapters.append(
            {
                "parameter_id": parameter_id,
                "adapter_dir": str(adapter_dir),
                "frcmod": str(frcmod),
                "library": str(library),
                "missing": missing,
                "violations": violations,
                "passed": not missing and not violations,
            }
        )
    passed = all(bool(adapter["passed"]) for adapter in adapters)
    write_artifact(
        run_dir / "cofactor_preflight.json",
        {"schema_version": "1.0", "passed": passed, "adapters": adapters},
    )
    common: dict[str, str | int | float] = {
        "integrator": "md",
        "dt": 0.002,
        "nstxout-compressed": 5_000,
        "nstenergy": 1_000,
        "constraints": "h-bonds",
        "constraint-algorithm": "lincs",
        "cutoff-scheme": "Verlet",
        "coulombtype": "PME",
        "rcoulomb": 1.0,
        "vdwtype": "Cut-off",
        "rvdw": 1.0,
        "tcoupl": "V-rescale",
        "tc-grps": "System",
        "tau-t": 0.1,
        "ref-t": 300,
        "pbc": "xyz",
    }
    nvt: dict[str, str | int | float] = {
        **common,
        "nsteps": 250_000,
        "continuation": "no",
        "pcoupl": "no",
        "gen-vel": "yes",
        "gen-temp": 300,
        "gen-seed": velocity_seed,
    }
    npt: dict[str, str | int | float] = {
        **common,
        "nsteps": 2_500_000,
        "continuation": "yes",
        "pcoupl": "C-rescale",
        "pcoupltype": "semiisotropic",
        "tau-p": 2.0,
        "ref-p": "1.0 1.0",
        "compressibility": "4.5e-5 4.5e-5",
        "gen-vel": "no",
    }
    minimization: dict[str, str | int | float] = {
        "integrator": "steep",
        "emtol": 1000.0,
        "emstep": 0.01,
        "nsteps": 50_000,
        "cutoff-scheme": "Verlet",
        "coulombtype": "PME",
        "rcoulomb": 1.0,
        "vdwtype": "Cut-off",
        "rvdw": 1.0,
        "constraints": "h-bonds",
        "pbc": "xyz",
    }
    mdp_files: dict[str, dict[str, str | int | float]] = {
        "min.mdp": minimization,
        "nvt.mdp": nvt,
        "npt.mdp": npt,
        "md.mdp": render_md_mdp(protocol=protocol, system_kind="membrane"),
    }
    for filename, settings in mdp_files.items():
        write_bytes_atomic(run_dir / filename, _mdp_text(settings).encode())
    return passed


def _combine_complex_pdb(components: ComplexComponents, output: Path) -> Path:
    parts: list[str] = []
    for path in (
        components.protein_pdb,
        components.ligand_pdb,
        components.cofactor_pdb,
    ):
        if path is None:
            continue
        parts.extend(
            line
            for line in path.read_text(encoding="utf-8").splitlines()
            if line not in {"END", "TER"}
        )
    write_bytes_atomic(output, ("\n".join(parts) + "\nTER\nEND\n").encode())
    return output


def _cofactor_parameter_files(
    parameter_ids: list[str], *, root: Path
) -> tuple[list[tuple[Path, Path]], list[str]]:
    parameters: list[tuple[Path, Path]] = []
    leap_lines: list[str] = []
    for parameter_id in parameter_ids:
        adapter_dir = root / parameter_id
        manifest = json.loads(
            (adapter_dir / "adapter.json").read_text(encoding="utf-8")
        )
        if manifest.get("parameter_id") != parameter_id:
            raise ValueError("cofactor adapter identity mismatch")
        raw_lines = manifest.get("leap_lines", [])
        if not isinstance(raw_lines, list) or not all(
            isinstance(line, str) for line in raw_lines
        ):
            raise ValueError("cofactor adapter leap_lines must be strings")
        parameters.append(
            (
                adapter_dir / "cofactor.frcmod",
                adapter_dir / "cofactor.lib",
            )
        )
        leap_lines.extend(raw_lines)
    return parameters, leap_lines


def build_membrane_md_system(
    complex_cif: Path,
    *,
    run_dir: Path,
    ligand_smiles: str,
    ligand_formal_charge: int,
    velocity_seed: int,
    protocol: MDProtocol,
    cofactor_parameter_ids: list[str],
    cofactor_parameter_root: Path,
    command_runner: CommandRunner | None = None,
) -> MDSystemBuildResult:
    """Build and equilibrate a preoriented POPC/CHL1 membrane system."""
    run_dir = run_dir.resolve()
    run_dir.mkdir(parents=True, exist_ok=True)
    executor = command_runner or _run_system_command
    try:
        components = extract_boltz_complex(complex_cif, output_dir=run_dir)
        prepare_ligand_for_amber(
            components.ligand_pdb,
            ligand_smiles=ligand_smiles,
            ligand_formal_charge=ligand_formal_charge,
            output_sdf=run_dir / "ligand.sdf",
        )
        adapters_ok = render_membrane_inputs(
            run_dir,
            velocity_seed=velocity_seed,
            protocol=protocol,
            cofactor_parameter_ids=cofactor_parameter_ids,
            cofactor_parameter_root=cofactor_parameter_root,
        )
        if not adapters_ok:
            result = MDSystemBuildResult(
                status="failed",
                error_code="cofactor_parameter_adapter_missing",
                completed_step="cofactor_preflight",
                run_dir=run_dir,
                protein_atom_count=components.protein_atom_count,
                ligand_atom_count=components.ligand_atom_count,
            )
            _write_build_status(run_dir, result)
            return result
        cofactor_parameters, leap_lines = _cofactor_parameter_files(
            cofactor_parameter_ids, root=cofactor_parameter_root
        )
        complex_pdb = _combine_complex_pdb(
            components, run_dir / "oriented_complex.pdb"
        )
        write_bytes_atomic(
            run_dir / "ligand.lib.in",
            (
                "source leaprc.gaff2\n"
                "LIG = loadMol2 ligand.mol2\n"
                "saveOff LIG ligand.lib\n"
                "quit\n"
            ).encode(),
        )
    except (OSError, RuntimeError, ValueError, json.JSONDecodeError) as error:
        result = MDSystemBuildResult(
            status="failed",
            error_code="invalid_boltz_structure",
            run_dir=run_dir,
        )
        write_bytes_atomic(run_dir / "membrane_setup.stderr.log", str(error).encode())
        _write_build_status(run_dir, result)
        return result

    commands: list[tuple[str, str, list[str], tuple[str, ...]]] = [
        (
            "parameterize_ligand",
            "ligand_parameterization_failed",
            [
                "antechamber",
                "-i",
                str(run_dir / "ligand.sdf"),
                "-fi",
                "sdf",
                "-o",
                str(run_dir / "ligand.mol2"),
                "-fo",
                "mol2",
                "-c",
                "bcc",
                "-at",
                "gaff2",
                "-nc",
                str(ligand_formal_charge),
                "-s",
                "2",
            ],
            ("ligand.mol2",),
        ),
        (
            "check_ligand_parameters",
            "ligand_parameterization_failed",
            [
                "parmchk2",
                "-i",
                str(run_dir / "ligand.mol2"),
                "-f",
                "mol2",
                "-o",
                str(run_dir / "ligand.frcmod"),
                "-s",
                "gaff2",
            ],
            ("ligand.frcmod",),
        ),
        (
            "write_ligand_library",
            "ligand_parameterization_failed",
            ["tleap", "-f", str(run_dir / "ligand.lib.in")],
            ("ligand.lib",),
        ),
    ]
    logs: list[Path] = []
    completed_step = "prepare_ligand_coordinates"
    for index, (step, error_code, command, expected) in enumerate(commands, 1):
        execution = executor(command, run_dir, None)
        log = run_dir / f"{index:02d}-{step}.log"
        write_bytes_atomic(
            log,
            (execution.stdout + "\n" + execution.stderr).encode(errors="replace"),
        )
        logs.append(log)
        if execution.returncode != 0 or not all(
            (run_dir / filename).is_file()
            and (run_dir / filename).stat().st_size > 0
            for filename in expected
        ):
            result = MDSystemBuildResult(
                status="failed",
                error_code=error_code,
                completed_step=completed_step,
                run_dir=run_dir,
                protein_atom_count=components.protein_atom_count,
                ligand_atom_count=components.ligand_atom_count,
                command_logs=logs,
            )
            _write_build_status(run_dir, result)
            return result
        completed_step = step

    membrane_pdb = run_dir / "membrane.pdb"
    membrane_execution: subprocess.CompletedProcess[str] | None = None
    try:
        membrane_command = build_membrane_command(
            complex_pdb=complex_pdb,
            output_pdb=membrane_pdb,
            ligand_frcmod=run_dir / "ligand.frcmod",
            ligand_lib=run_dir / "ligand.lib",
            cofactor_parameters=cofactor_parameters,
        )
        for line in leap_lines:
            membrane_command.extend(["--leapline", line])
        membrane_execution = executor(
            membrane_command, run_dir, None
        )
    except (OSError, subprocess.SubprocessError, ValueError) as error:
        membrane_execution = None
        detail = str(error)
    else:
        assert membrane_execution is not None
        detail = membrane_execution.stdout + "\n" + membrane_execution.stderr
    log = run_dir / "04-pack_membrane.log"
    write_bytes_atomic(log, detail.encode(errors="replace"))
    logs.append(log)
    topology_candidates = sorted(run_dir.glob("*_lipid.top"))
    coordinate_candidates = sorted(run_dir.glob("*_lipid.crd"))
    if (
        membrane_execution is None
        or membrane_execution.returncode != 0
        or len(topology_candidates) != 1
        or len(coordinate_candidates) != 1
    ):
        result = MDSystemBuildResult(
            status="failed",
            error_code="membrane_assembly_failed",
            completed_step=completed_step,
            run_dir=run_dir,
            protein_atom_count=components.protein_atom_count,
            ligand_atom_count=components.ligand_atom_count,
            command_logs=logs,
        )
        _write_build_status(run_dir, result)
        return result
    completed_step = "pack_membrane"

    conversion = f'''from pathlib import Path
import parmed as pmd
root = Path(__file__).resolve().parent
system = pmd.load_file({str(topology_candidates[0])!r}, xyz={str(coordinate_candidates[0])!r})
system.save(str(root / "topol.top"), format="gromacs", overwrite=True)
system.save(str(root / "initial.gro"), format="gro", overwrite=True)
'''
    write_bytes_atomic(run_dir / "convert_membrane.py", conversion.encode())
    equilibration: list[tuple[str, str, list[str], tuple[str, ...]]] = [
        (
            "convert_topology",
            "topology_conversion_failed",
            ["python", str(run_dir / "convert_membrane.py")],
            ("topol.top", "initial.gro"),
        ),
        (
            "prepare_minimization",
            "minimization_failed",
            [
                "gmx",
                "grompp",
                "-f",
                str(run_dir / "min.mdp"),
                "-c",
                str(run_dir / "initial.gro"),
                "-p",
                str(run_dir / "topol.top"),
                "-o",
                str(run_dir / "min.tpr"),
            ],
            ("min.tpr",),
        ),
        (
            "minimize",
            "minimization_failed",
            ["gmx", "mdrun", "-deffnm", str(run_dir / "min")],
            ("min.gro",),
        ),
        (
            "prepare_nvt",
            "equilibration_failed",
            [
                "gmx",
                "grompp",
                "-f",
                str(run_dir / "nvt.mdp"),
                "-c",
                str(run_dir / "min.gro"),
                "-r",
                str(run_dir / "min.gro"),
                "-p",
                str(run_dir / "topol.top"),
                "-o",
                str(run_dir / "nvt.tpr"),
            ],
            ("nvt.tpr",),
        ),
        (
            "equilibrate_nvt",
            "equilibration_failed",
            ["gmx", "mdrun", "-deffnm", str(run_dir / "nvt")],
            ("nvt.gro", "nvt.cpt"),
        ),
        (
            "prepare_npt",
            "equilibration_failed",
            [
                "gmx",
                "grompp",
                "-f",
                str(run_dir / "npt.mdp"),
                "-c",
                str(run_dir / "nvt.gro"),
                "-r",
                str(run_dir / "nvt.gro"),
                "-t",
                str(run_dir / "nvt.cpt"),
                "-p",
                str(run_dir / "topol.top"),
                "-o",
                str(run_dir / "npt.tpr"),
            ],
            ("npt.tpr",),
        ),
        (
            "equilibrate_npt",
            "equilibration_failed",
            ["gmx", "mdrun", "-deffnm", str(run_dir / "npt")],
            ("npt.gro", "npt.cpt"),
        ),
    ]
    for offset, (step, error_code, command, expected) in enumerate(
        equilibration, 5
    ):
        execution = executor(command, run_dir, None)
        log = run_dir / f"{offset:02d}-{step}.log"
        write_bytes_atomic(
            log,
            (execution.stdout + "\n" + execution.stderr).encode(errors="replace"),
        )
        logs.append(log)
        if execution.returncode != 0 or not all(
            (run_dir / filename).is_file()
            and (run_dir / filename).stat().st_size > 0
            for filename in expected
        ):
            result = MDSystemBuildResult(
                status="failed",
                error_code=error_code,
                completed_step=completed_step,
                run_dir=run_dir,
                protein_atom_count=components.protein_atom_count,
                ligand_atom_count=components.ligand_atom_count,
                command_logs=logs,
            )
            _write_build_status(run_dir, result)
            return result
        completed_step = step
    result = MDSystemBuildResult(
        status="succeeded",
        completed_step=completed_step,
        forcefield="ff19SB+Lipid21+GAFF2/TIP3P",
        run_dir=run_dir,
        protein_atom_count=components.protein_atom_count,
        ligand_atom_count=components.ligand_atom_count,
        command_logs=logs,
    )
    _write_build_status(run_dir, result)
    return result


def _mdp_text(settings: dict[str, str | int | float]) -> str:
    return "".join(f"{key} = {value}\n" for key, value in settings.items())


def render_system_inputs(
    run_dir: Path,
    *,
    velocity_seed: int,
    ligand_formal_charge: int = 0,
    protocol: MDProtocol = "production",
) -> None:
    """Write the frozen ff19SB/GAFF2/TIP3P build and equilibration inputs."""
    if not 1 <= velocity_seed <= 2_147_483_647:
        raise ValueError("velocity seed must be a positive signed 32-bit integer")
    run_dir.mkdir(parents=True, exist_ok=True)
    leap_input = """source leaprc.protein.ff19SB
source leaprc.gaff2
source leaprc.water.tip3p
LIG = loadMol2 ligand.mol2
loadAmberParams ligand.frcmod
PROT = loadPdb protein.amber.pdb
COM = combine { PROT LIG }
check COM
solvateOct COM TIP3PBOX 10.0
addIonsRand COM Na+ 1
addIonsRand COM Cl- 1
addIonsRand COM Na+ 0
addIonsRand COM Cl- 0
saveAmberParm COM solvated.prmtop solvated.inpcrd
savePdb COM solvated.pdb
quit
"""
    parmed_script = """import json
from pathlib import Path

import parmed as pmd

root = Path(__file__).resolve().parent
ligand_formal_charge = __LIGAND_FORMAL_CHARGE__
system = pmd.load_file(
    str(root / "solvated.prmtop"),
    xyz=str(root / "solvated.inpcrd"),
)
ligand_residues = [residue for residue in system.residues if residue.name == "MOL"]
if len(ligand_residues) != 1:
    raise ValueError("expected exactly one MOL ligand residue")
ligand_atoms = list(ligand_residues[0].atoms)
observed_charge = sum(float(atom.charge) for atom in ligand_atoms)
charge_correction = ligand_formal_charge - observed_charge
if abs(charge_correction) > 0.05:
    raise ValueError(
        "ligand charge residual exceeds the 0.05 e normalization limit: "
        f"{charge_correction:.8f}"
    )
per_atom_correction = charge_correction / len(ligand_atoms)
for atom in ligand_atoms:
    atom.charge += per_atom_correction
for residue in system.residues:
    original_name = residue.name
    if original_name not in {"Na+", "Cl-"}:
        continue
    normalized_name = {"Na+": "NA", "Cl-": "CL"}[original_name]
    residue.name = normalized_name
    for atom in residue.atoms:
        atom.name = normalized_name
(root / "charge_normalization.json").write_text(
    json.dumps(
        {
            "atom_count": len(ligand_atoms),
            "correction_e": charge_correction,
            "observed_charge_e": observed_charge,
            "per_atom_correction_e": per_atom_correction,
            "target_formal_charge_e": ligand_formal_charge,
        },
        sort_keys=True,
        separators=(",", ":"),
    )
    + "\\n",
    encoding="utf-8",
)
system.save(str(root / "topol.top"), format="gromacs", overwrite=True)
system.save(str(root / "solvated.gro"), format="gro", overwrite=True)
"""
    parmed_script = parmed_script.replace(
        "__LIGAND_FORMAL_CHARGE__", str(ligand_formal_charge)
    )
    water_index_script = """from pathlib import Path

root = Path(__file__).resolve().parent
lines = (root / "solvated.gro").read_text(encoding="utf-8").splitlines()
atom_count = int(lines[1].strip())
atom_lines = lines[2 : 2 + atom_count]
water_atoms = [
    index
    for index, line in enumerate(atom_lines, 1)
    if line[5:10].strip() == "WAT"
]
if not water_atoms:
    raise ValueError("ParmEd GROMACS coordinates contain no WAT solvent")
chunks = [
    " ".join(str(value) for value in water_atoms[start : start + 15])
    for start in range(0, len(water_atoms), 15)
]
(root / "water.ndx").write_text(
    "[ WAT ]\\n" + "\\n".join(chunks) + "\\n",
    encoding="utf-8",
)
"""
    inputs: dict[str, str] = {
        "build.leap.in": leap_input,
        "convert_parmed.py": parmed_script,
        "write_water_index.py": water_index_script,
        "ions.mdp": _mdp_text(
            {
                "integrator": "steep",
                "emtol": 1000.0,
                "nsteps": 500,
                "cutoff-scheme": "Verlet",
                "coulombtype": "Cut-off",
                "rcoulomb": 1.0,
                "rvdw": 1.0,
                "pbc": "xyz",
            }
        ),
        "min.mdp": _mdp_text(
            {
                "integrator": "steep",
                "emtol": 1000.0,
                "emstep": 0.01,
                "nsteps": 50_000,
                "cutoff-scheme": "Verlet",
                "coulombtype": "PME",
                "rcoulomb": 1.0,
                "vdwtype": "Cut-off",
                "rvdw": 1.0,
                "constraints": "h-bonds",
                "pbc": "xyz",
            }
        ),
        "nvt.mdp": _mdp_text(
            {
                "integrator": "md",
                "dt": 0.002,
                "nsteps": 50_000,
                "nstxout-compressed": 5_000,
                "nstenergy": 1_000,
                "continuation": "no",
                "constraints": "h-bonds",
                "constraint-algorithm": "lincs",
                "cutoff-scheme": "Verlet",
                "coulombtype": "PME",
                "rcoulomb": 1.0,
                "vdwtype": "Cut-off",
                "rvdw": 1.0,
                "tcoupl": "V-rescale",
                "tc-grps": "System",
                "tau-t": 0.1,
                "ref-t": 300,
                "pcoupl": "no",
                "gen-vel": "yes",
                "gen-temp": 300,
                "gen-seed": velocity_seed,
                "pbc": "xyz",
            }
        ),
        "npt.mdp": _mdp_text(
            {
                "integrator": "md",
                "dt": 0.002,
                "nsteps": 250_000,
                "nstxout-compressed": 5_000,
                "nstenergy": 1_000,
                "continuation": "yes",
                "constraints": "h-bonds",
                "constraint-algorithm": "lincs",
                "cutoff-scheme": "Verlet",
                "coulombtype": "PME",
                "rcoulomb": 1.0,
                "vdwtype": "Cut-off",
                "rvdw": 1.0,
                "tcoupl": "V-rescale",
                "tc-grps": "System",
                "tau-t": 0.1,
                "ref-t": 300,
                "pcoupl": "C-rescale",
                "tau-p": 2.0,
                "ref-p": 1.0,
                "compressibility": 4.5e-5,
                "gen-vel": "no",
                "pbc": "xyz",
            }
        ),
        "md.mdp": _mdp_text(
            render_md_mdp(protocol=protocol, system_kind="soluble")
        ),
    }
    for filename, content in inputs.items():
        write_bytes_atomic(run_dir / filename, content.encode("utf-8"))


def build_grompp_command(run_dir: Path) -> list[str]:
    """Build the non-interactive production preprocessing command."""
    return [
        "gmx",
        "grompp",
        "-f",
        str(run_dir / "md.mdp"),
        "-c",
        str(run_dir / "npt.gro"),
        "-p",
        str(run_dir / "topol.top"),
        "-o",
        str(run_dir / "md.tpr"),
    ]


def build_mdrun_command(run_dir: Path) -> list[str]:
    """Build a checkpoint-aware production command."""
    command = ["gmx", "mdrun", "-deffnm", str(run_dir / "md")]
    checkpoint = run_dir / "md.cpt"
    if checkpoint.exists():
        command.extend(["-cpi", str(checkpoint), "-append"])
    return command


def build_system_commands(
    complex_pdb: Path,
    run_dir: Path,
    *,
    ligand_formal_charge: int = 0,
) -> list[list[str]]:
    """Build the strict AmberTools-to-GROMACS system/equilibration sequence."""
    del complex_pdb
    return [
        [
            "antechamber",
            "-i",
            str(run_dir / "ligand.sdf"),
            "-fi",
            "sdf",
            "-o",
            str(run_dir / "ligand.mol2"),
            "-fo",
            "mol2",
            "-c",
            "bcc",
            "-at",
            "gaff2",
            "-nc",
            str(ligand_formal_charge),
            "-s",
            "2",
        ],
        [
            "parmchk2",
            "-i",
            str(run_dir / "ligand.mol2"),
            "-f",
            "mol2",
            "-o",
            str(run_dir / "ligand.frcmod"),
            "-s",
            "gaff2",
        ],
        [
            "pdb4amber",
            "-i",
            str(run_dir / "protein.pdb"),
            "-o",
            str(run_dir / "protein.amber.pdb"),
            "--nohyd",
        ],
        ["tleap", "-f", str(run_dir / "build.leap.in")],
        ["python", str(run_dir / "convert_parmed.py")],
        [
            "gmx",
            "grompp",
            "-f",
            str(run_dir / "ions.mdp"),
            "-c",
            str(run_dir / "solvated.gro"),
            "-p",
            str(run_dir / "topol.top"),
            "-o",
            str(run_dir / "ions.tpr"),
        ],
        ["python", str(run_dir / "write_water_index.py")],
        [
            "gmx",
            "genion",
            "-s",
            str(run_dir / "ions.tpr"),
            "-n",
            str(run_dir / "water.ndx"),
            "-p",
            str(run_dir / "topol.top"),
            "-o",
            str(run_dir / "ions.gro"),
            "-pname",
            "NA",
            "-nname",
            "CL",
            "-conc",
            "0.15",
            "-neutral",
        ],
        [
            "gmx",
            "grompp",
            "-f",
            str(run_dir / "min.mdp"),
            "-c",
            str(run_dir / "ions.gro"),
            "-p",
            str(run_dir / "topol.top"),
            "-o",
            str(run_dir / "min.tpr"),
        ],
        ["gmx", "mdrun", "-deffnm", str(run_dir / "min")],
        [
            "gmx",
            "grompp",
            "-f",
            str(run_dir / "nvt.mdp"),
            "-c",
            str(run_dir / "min.gro"),
            "-r",
            str(run_dir / "min.gro"),
            "-p",
            str(run_dir / "topol.top"),
            "-o",
            str(run_dir / "nvt.tpr"),
        ],
        ["gmx", "mdrun", "-deffnm", str(run_dir / "nvt")],
        [
            "gmx",
            "grompp",
            "-f",
            str(run_dir / "npt.mdp"),
            "-c",
            str(run_dir / "nvt.gro"),
            "-r",
            str(run_dir / "nvt.gro"),
            "-t",
            str(run_dir / "nvt.cpt"),
            "-p",
            str(run_dir / "topol.top"),
            "-o",
            str(run_dir / "npt.tpr"),
        ],
        ["gmx", "mdrun", "-deffnm", str(run_dir / "npt")],
    ]


def _run_system_command(
    command: list[str], cwd: Path, stdin: str | None
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        cwd=cwd,
        input=stdin,
        capture_output=True,
        text=True,
        check=False,
        timeout=21_600,
    )


def _write_build_status(run_dir: Path, result: MDSystemBuildResult) -> None:
    write_artifact(run_dir / "system_build_status.json", result.model_dump(mode="json"))


def build_md_system(
    complex_cif: Path,
    *,
    run_dir: Path,
    ligand_smiles: str,
    ligand_formal_charge: int,
    velocity_seed: int,
    protocol: MDProtocol = "production",
    command_runner: CommandRunner = _run_system_command,
) -> MDSystemBuildResult:
    """Build, salt, minimize and equilibrate one Boltz protein–ligand complex."""
    complex_cif = complex_cif.resolve()
    run_dir = run_dir.resolve()
    run_dir.mkdir(parents=True, exist_ok=True)
    try:
        components = extract_boltz_complex(complex_cif, output_dir=run_dir)
        prepare_ligand_for_amber(
            components.ligand_pdb,
            ligand_smiles=ligand_smiles,
            ligand_formal_charge=ligand_formal_charge,
            output_sdf=run_dir / "ligand.sdf",
        )
        render_system_inputs(
            run_dir,
            velocity_seed=velocity_seed,
            ligand_formal_charge=ligand_formal_charge,
            protocol=protocol,
        )
    except (OSError, RuntimeError, ValueError) as error:
        result = MDSystemBuildResult(
            status="failed",
            error_code="invalid_boltz_structure",
            completed_step=None,
            run_dir=run_dir,
        )
        write_bytes_atomic(run_dir / "extract_complex.stderr.log", str(error).encode())
        _write_build_status(run_dir, result)
        return result

    steps: list[tuple[str, str, tuple[str, ...]]] = [
        ("parameterize_ligand", "ligand_parameterization_failed", ("ligand.mol2",)),
        ("check_ligand_parameters", "ligand_parameterization_failed", ("ligand.frcmod",)),
        ("prepare_protein", "protein_preparation_failed", ("protein.amber.pdb",)),
        (
            "assemble_and_solvate",
            "system_assembly_failed",
            ("solvated.prmtop", "solvated.inpcrd"),
        ),
        (
            "convert_topology",
            "topology_conversion_failed",
            ("topol.top", "solvated.gro", "charge_normalization.json"),
        ),
        ("prepare_ions", "ionization_failed", ("ions.tpr",)),
        ("index_water", "ionization_failed", ("water.ndx",)),
        ("add_ions", "ionization_failed", ("ions.gro",)),
        ("prepare_minimization", "minimization_failed", ("min.tpr",)),
        ("minimize", "minimization_failed", ("min.gro",)),
        ("prepare_nvt", "equilibration_failed", ("nvt.tpr",)),
        ("equilibrate_nvt", "equilibration_failed", ("nvt.gro", "nvt.cpt")),
        ("prepare_npt", "equilibration_failed", ("npt.tpr",)),
        ("equilibrate_npt", "equilibration_failed", ("npt.gro", "npt.cpt")),
    ]
    commands = build_system_commands(
        complex_cif,
        run_dir,
        ligand_formal_charge=ligand_formal_charge,
    )
    if len(commands) != len(steps):
        raise RuntimeError("internal MD system-build command contract mismatch")

    completed_step = "prepare_ligand_coordinates"
    logs: list[Path] = []
    for index, (command, step) in enumerate(zip(commands, steps, strict=True), 1):
        step_name, error_code, expected_files = step
        stdin = "WAT\n" if command[:2] == ["gmx", "genion"] else None
        log_path = run_dir / f"{index:02d}-{step_name}.log"
        try:
            execution = command_runner(command, run_dir, stdin)
            log_text = (
                "$ "
                + " ".join(command)
                + "\n\n[stdout]\n"
                + execution.stdout
                + "\n[stderr]\n"
                + execution.stderr
            )
            write_bytes_atomic(log_path, log_text.encode("utf-8", errors="replace"))
            logs.append(log_path)
        except (OSError, subprocess.SubprocessError) as error:
            write_bytes_atomic(log_path, str(error).encode("utf-8", errors="replace"))
            logs.append(log_path)
            execution = None
        files_valid = all(
            (run_dir / filename).is_file()
            and (run_dir / filename).stat().st_size > 0
            for filename in expected_files
        )
        if execution is None or execution.returncode != 0 or not files_valid:
            result = MDSystemBuildResult(
                status="failed",
                error_code=error_code,
                completed_step=completed_step,
                run_dir=run_dir,
                protein_atom_count=components.protein_atom_count,
                ligand_atom_count=components.ligand_atom_count,
                command_logs=logs,
            )
            _write_build_status(run_dir, result)
            return result
        completed_step = step_name

    result = MDSystemBuildResult(
        status="succeeded",
        completed_step=completed_step,
        forcefield="ff19SB+GAFF2/TIP3P",
        run_dir=run_dir,
        protein_atom_count=components.protein_atom_count,
        ligand_atom_count=components.ligand_atom_count,
        command_logs=logs,
    )
    _write_build_status(run_dir, result)
    return result


def _stable_seed(target_id: str, replica: int) -> int:
    material = f"airti-md-v1:{target_id}:{replica}".encode()
    # GROMACS expects a positive signed 32-bit seed.
    return int.from_bytes(hashlib.sha256(material).digest()[:4], "big") % 2_147_483_646 + 1


def plan_md_replicas(target_ids: list[str]) -> list[MDReplicaPlan]:
    """Plan Top10 MD, with three independent replicas for Top3."""
    plans: list[MDReplicaPlan] = []
    for rank, target_id in enumerate(target_ids[:10], start=1):
        replica_count = 3 if rank <= 3 else 1
        for replica in range(1, replica_count + 1):
            plans.append(
                MDReplicaPlan(
                    target_id=target_id,
                    replica=replica,
                    velocity_seed=_stable_seed(target_id, replica),
                )
            )
    return plans


def validate_parameterization(
    *, antechamber_ok: bool, parmed_ok: bool
) -> ParameterizationResult:
    """Reject a failed parameterization instead of changing force fields."""
    if not antechamber_ok:
        return ParameterizationResult(
            status="failed", error_code="ligand_parameterization_failed"
        )
    if not parmed_ok:
        return ParameterizationResult(
            status="failed", error_code="topology_conversion_failed"
        )
    return ParameterizationResult(status="succeeded", forcefield="ff19SB+GAFF2/TIP3P")


def parse_completed_ns(log_path: Path) -> float:
    """Read the last GROMACS time value (ps) and return completed ns."""
    text = log_path.read_text(encoding="utf-8", errors="replace")
    values = re.findall(r"^\s*Time\s*[:=]\s*([0-9]+(?:\.[0-9]+)?)\s*(?:ps)?\s*$", text, re.MULTILINE)
    if not values:
        values = re.findall(r"^\s*Step\s+Time\s*\n\s*\d+\s+([0-9]+(?:\.[0-9]+)?)", text, re.MULTILINE)
    if not values:
        raise ValueError(f"no completed simulation time found in {log_path}")
    return max(float(value) for value in values) / 1_000.0


def run_md_trajectory(
    *,
    run_dir: Path,
    protocol: MDProtocol,
    pocket_residues: list[int],
    command_runner: CommandRunner = _run_system_command,
) -> MDTrajectoryResult:
    """Run checkpoint-aware MD, process PBC, and compute trajectory evidence."""
    expected_ns = 1.0 if protocol == "smoke" else 100.0
    analysis_role: Literal[
        "scientific_evidence", "pipeline_validation_only"
    ] = (
        "pipeline_validation_only"
        if protocol == "smoke"
        else "scientific_evidence"
    )
    commands = [
        build_grompp_command(run_dir),
        build_mdrun_command(run_dir),
    ]
    for index, command in enumerate(commands, 1):
        result = command_runner(command, run_dir, None)
        write_bytes_atomic(
            run_dir / f"trajectory-{index:02d}.log",
            (result.stdout + "\n" + result.stderr).encode(errors="replace"),
        )
        if result.returncode != 0:
            return MDTrajectoryResult(
                status="failed",
                completed_ns=0,
                error_code=(
                    "production_preprocessing_failed"
                    if index == 1
                    else "production_run_failed"
                ),
            )
    trajectory = run_dir / "md.xtc"
    checkpoint = run_dir / "md.cpt"
    log = run_dir / "md.log"
    if not all(
        path.is_file() and path.stat().st_size > 0
        for path in (trajectory, checkpoint, log)
    ):
        return MDTrajectoryResult(
            status="failed",
            completed_ns=0,
            error_code="trajectory_artifact_missing",
        )
    pbc_trajectory = run_dir / "md_pbc.xtc"
    pbc_command = [
        "gmx",
        "trjconv",
        "-s",
        str(run_dir / "md.tpr"),
        "-f",
        str(trajectory),
        "-o",
        str(pbc_trajectory),
        "-pbc",
        "mol",
        "-center",
    ]
    pbc_result = command_runner(pbc_command, run_dir, "Protein\nSystem\n")
    write_bytes_atomic(
        run_dir / "trajectory-03-pbc.log",
        (pbc_result.stdout + "\n" + pbc_result.stderr).encode(errors="replace"),
    )
    if (
        pbc_result.returncode != 0
        or not pbc_trajectory.is_file()
        or pbc_trajectory.stat().st_size == 0
    ):
        return MDTrajectoryResult(
            status="failed",
            completed_ns=parse_completed_ns(log),
            error_code="pbc_processing_failed",
            trajectory_path=trajectory,
            checkpoint_path=checkpoint,
        )
    try:
        metrics = measure_trajectory(
            topology=run_dir / "md.tpr",
            trajectory=pbc_trajectory,
            log_path=log,
            pocket_residues=pocket_residues,
            expected_ns=expected_ns,
            analysis_role=analysis_role,
        )
        analysis: TrajectoryAnalysis = analyze_trajectory(metrics)
    except (OSError, RuntimeError, TypeError, ValueError) as error:
        write_bytes_atomic(
            run_dir / "trajectory_analysis.stderr.log",
            str(error).encode(errors="replace"),
        )
        return MDTrajectoryResult(
            status="failed",
            completed_ns=parse_completed_ns(log),
            error_code="trajectory_analysis_failed",
            trajectory_path=pbc_trajectory,
            checkpoint_path=checkpoint,
        )
    metrics_path = run_dir / "trajectory_metrics.json"
    write_artifact(metrics_path, metrics.model_dump(mode="json"))
    write_artifact(
        run_dir / "trajectory_analysis.json",
        analysis.model_dump(mode="json"),
    )
    return MDTrajectoryResult(
        status=analysis.status,
        completed_ns=metrics.completed_ns,
        md_score=analysis.md_score,
        error_code=analysis.error_code,
        trajectory_path=pbc_trajectory,
        checkpoint_path=checkpoint,
        metrics_path=metrics_path,
    )
