"""Auditable AmberTools/ParmEd/GROMACS molecular-dynamics protocol."""

from __future__ import annotations

import hashlib
import re
import subprocess
from collections.abc import Callable
from pathlib import Path
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field

from airti_tf.manifest_io import write_artifact, write_bytes_atomic

CommandRunner = Callable[
    [list[str], Path, str | None], subprocess.CompletedProcess[str]
]


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
    protein_atom_count: int
    ligand_atom_count: int


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
    raw_data = MMCIF2Dict(str(complex_cif))  # type: ignore[no-untyped-call]
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
        if not is_protein and not is_ligand:
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
        output_chain = "A" if is_protein else "B"
        output_residue = residues[index] if is_protein else "LIG"
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
        else:
            ligand_components.add(residues[index])
            ligand_lines.append(line)
    if not protein_lines:
        raise ValueError("Boltz complex contains no protein atoms in chain A")
    if not ligand_lines:
        raise ValueError("Boltz complex contains no ligand atoms in chain B")
    if len(ligand_components) != 1:
        raise ValueError("Boltz complex must contain exactly one ligand component")

    output_dir.mkdir(parents=True, exist_ok=True)
    protein_path = output_dir / "protein.pdb"
    ligand_path = output_dir / "ligand.pdb"
    write_bytes_atomic(
        protein_path, ("".join(protein_lines) + "TER\nEND\n").encode("utf-8")
    )
    write_bytes_atomic(ligand_path, ("".join(ligand_lines) + "END\n").encode("utf-8"))
    return ComplexComponents(
        protein_pdb=protein_path,
        ligand_pdb=ligand_path,
        protein_atom_count=len(protein_lines),
        ligand_atom_count=len(ligand_lines),
    )


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
        assigned = AllChem.AssignBondOrdersFromTemplate(  # type: ignore[no-untyped-call]
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


def render_production_mdp() -> dict[str, str | int | float]:
    """Return the fixed 100 ns, 300 K, 1 bar production protocol.

    GROMACS uses ps for ``dt``. Therefore 50,000,000 steps at 0.002 ps
    equal exactly 100,000 ps (100 ns). Coordinates and energies are written
    every 5,000 steps, i.e. every 10 ps.
    """
    return {
        "integrator": "md",
        "dt": 0.002,
        "nsteps": 50_000_000,
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


def _mdp_text(settings: dict[str, str | int | float]) -> str:
    return "".join(f"{key} = {value}\n" for key, value in settings.items())


def render_system_inputs(
    run_dir: Path,
    *,
    velocity_seed: int,
    ligand_formal_charge: int = 0,
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
        "md.mdp": _mdp_text(render_production_mdp()),
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
