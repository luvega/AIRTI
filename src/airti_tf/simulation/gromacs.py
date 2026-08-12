"""Auditable AmberTools/ParmEd/GROMACS molecular-dynamics protocol."""

from __future__ import annotations

import hashlib
import re
from pathlib import Path
from typing import Literal

from pydantic import BaseModel, ConfigDict


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
        "ref-t": 300,
        "pcoupl": "C-rescale",
        "ref-p": 1.0,
        "compressibility": 4.5e-5,
        "gen-vel": "no",
    }


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


def build_system_commands(complex_pdb: Path, run_dir: Path) -> list[list[str]]:
    """Describe the strict Amber-to-GROMACS system-build sequence.

    The companion LEaP input must load ``leaprc.protein.ff19SB``,
    ``leaprc.gaff2`` and TIP3P. ParmEd performs the topology conversion; no
    silent force-field substitution is allowed.
    """
    return [
        [
            "tleap",
            "-f",
            str(run_dir / "build.leap.in"),
            "--source-complex",
            str(complex_pdb),
        ],
        [
            "parmed",
            "-p",
            str(run_dir / "complex.prmtop"),
            "-c",
            str(run_dir / "complex.inpcrd"),
            "-i",
            str(run_dir / "convert.parmed.in"),
        ],
        [
            "gmx",
            "editconf",
            "-f",
            str(run_dir / "amber.gro"),
            "-o",
            str(run_dir / "boxed.gro"),
            "-bt",
            "dodecahedron",
            "-d",
            "1.0",
        ],
        [
            "gmx",
            "solvate",
            "-cp",
            str(run_dir / "boxed.gro"),
            "-cs",
            "tip3p",
            "-p",
            str(run_dir / "topol.top"),
            "-o",
            str(run_dir / "solvated.gro"),
        ],
        [
            "gmx",
            "genion",
            "-s",
            str(run_dir / "ions.tpr"),
            "-p",
            str(run_dir / "topol.top"),
            "-o",
            str(run_dir / "ions.gro"),
            "-conc",
            "0.15",
            "-neutral",
        ],
    ]


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

