"""Stable fpocket command, parser, and pocket quality gates."""

from __future__ import annotations

import hashlib
import re
from pathlib import Path
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator

PocketReason = Literal[
    "volume_too_small",
    "too_few_residues",
    "low_confidence",
    "inaccessible",
    "backbone_clash",
]


class PocketCandidate(BaseModel):
    model_config = ConfigDict(extra="forbid")

    pocket_id: str
    target_id: str
    rank: int = Field(gt=0)
    volume_a3: float = Field(ge=0)
    druggability: float = Field(ge=0, le=1)
    fpocket_score: float
    residue_count: int = Field(ge=0)
    mean_plddt: float = Field(default=100.0, ge=0, le=100)
    exposed_fraction: float = Field(default=1.0, ge=0, le=1)
    severe_backbone_clash: bool = False
    known_ligand_overlap: bool = False


class PocketQC(BaseModel):
    model_config = ConfigDict(extra="forbid")

    pocket: PocketCandidate
    status: Literal["ready", "unsupported"]
    unsupported_reason: PocketReason | None = None
    score: float | None = Field(default=None, ge=0, le=1)

    @model_validator(mode="after")
    def validate_state(self) -> "PocketQC":
        if self.status == "ready" and self.score is None:
            raise ValueError("ready pocket requires score")
        if self.status == "unsupported" and (
            self.unsupported_reason is None or self.score is not None
        ):
            raise ValueError("unsupported pocket requires reason and no score")
        return self


def build_fpocket_command(input_pdb: Path) -> list[str]:
    """Build the documented fpocket invocation."""
    return ["fpocket", "-f", str(input_pdb)]


def _field(block: str, name: str) -> float:
    match = re.search(rf"^\s*{re.escape(name)}\s*:\s*([-+0-9.eE]+)\s*$", block, re.MULTILINE)
    if match is None:
        raise ValueError(f"fpocket output lacks field: {name}")
    return float(match.group(1))


def _residue_count(pocket_pdb: Path) -> int:
    residues: set[tuple[str, str]] = set()
    if not pocket_pdb.is_file():
        return 0
    with pocket_pdb.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            fields = line.split()
            if len(fields) >= 6:
                residues.add((fields[4], fields[5]))
    return len(residues)


def parse_fpocket(
    output_dir: Path,
    *,
    target_id: str,
    structure_hash: str | None = None,
) -> list[PocketCandidate]:
    """Parse fpocket output into stable target/pocket records."""
    info_files = sorted(output_dir.glob("*_info.txt"))
    if len(info_files) != 1:
        raise ValueError(f"expected one fpocket info file in {output_dir}")
    info_path = info_files[0]
    raw = info_path.read_bytes()
    stable_structure_hash = structure_hash or hashlib.sha256(raw).hexdigest()
    text = raw.decode("utf-8")
    blocks = re.findall(
        r"^Pocket\s+(\d+)\s*:\s*(.*?)(?=^Pocket\s+\d+\s*:|\Z)",
        text,
        flags=re.MULTILINE | re.DOTALL,
    )
    candidates: list[PocketCandidate] = []
    for raw_rank, block in blocks:
        rank = int(raw_rank)
        pocket_pdb = output_dir / "pockets" / f"pocket{rank}_atm.pdb"
        candidates.append(
            PocketCandidate(
                pocket_id=f"{target_id}:{stable_structure_hash}:{rank}",
                target_id=target_id,
                rank=rank,
                volume_a3=_field(block, "Volume"),
                druggability=_field(block, "Druggability Score"),
                fpocket_score=_field(block, "Score"),
                residue_count=_residue_count(pocket_pdb),
            )
        )
    return sorted(candidates, key=lambda pocket: pocket.rank)


def qc_pocket(
    pocket: PocketCandidate,
    *,
    structure_source: Literal["pdb", "alphafold"] = "pdb",
) -> PocketQC:
    """Apply hard pocket quality gates without generating zero scores."""
    if pocket.volume_a3 < 80:
        return PocketQC(
            pocket=pocket, status="unsupported", unsupported_reason="volume_too_small"
        )
    if pocket.residue_count < 6:
        return PocketQC(
            pocket=pocket, status="unsupported", unsupported_reason="too_few_residues"
        )
    if structure_source == "alphafold" and pocket.mean_plddt < 70:
        return PocketQC(
            pocket=pocket, status="unsupported", unsupported_reason="low_confidence"
        )
    if pocket.exposed_fraction < 0.05:
        return PocketQC(
            pocket=pocket, status="unsupported", unsupported_reason="inaccessible"
        )
    if pocket.severe_backbone_clash:
        return PocketQC(
            pocket=pocket, status="unsupported", unsupported_reason="backbone_clash"
        )
    confidence = pocket.mean_plddt / 100
    volume_score = min(pocket.volume_a3 / 400, 1.0)
    score = (
        0.45 * pocket.druggability
        + 0.25 * volume_score
        + 0.20 * confidence
        + 0.10 * pocket.exposed_fraction
    )
    return PocketQC(pocket=pocket, status="ready", score=min(score, 1.0))


def select_qualified_pockets(
    candidates: list[PocketCandidate],
    *,
    limit: int = 5,
    structure_source: Literal["pdb", "alphafold"] = "pdb",
) -> list[PocketQC]:
    """Return up to five qualified pockets with ligand-overlap priority."""
    qualified = [
        result
        for candidate in candidates
        if (result := qc_pocket(candidate, structure_source=structure_source)).status == "ready"
    ]
    return sorted(
        qualified,
        key=lambda result: (
            not result.pocket.known_ligand_overlap,
            -(result.score or 0),
            result.pocket.rank,
        ),
    )[:limit]
