"""Deterministic molecular-weight-stratified background probe selection."""

from __future__ import annotations

import csv
import io
from collections import Counter
from pathlib import Path

from pydantic import BaseModel, ConfigDict
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, Lipinski, rdFingerprintGenerator
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.Chem.MolStandardize import rdMolStandardize

from airti_tf.manifest_io import write_bytes_atomic


class BackgroundCandidate(BaseModel):
    model_config = ConfigDict(extra="forbid")

    source_id: str
    canonical_smiles: str
    molecular_weight: float
    clogp: float
    hbd: int
    hba: int
    rotatable_bonds: int
    formal_charge: int
    fingerprint: bytes


class BackgroundProbe(BaseModel):
    model_config = ConfigDict(extra="forbid")

    source_id: str
    canonical_smiles: str
    molecular_weight: float
    clogp: float
    hbd: int
    hba: int
    rotatable_bonds: int
    formal_charge: int
    mw_stratum: int


_FINGERPRINT_GENERATOR = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
_PAINS_PARAMETERS = FilterCatalogParams()
_PAINS_PARAMETERS.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
_PAINS_CATALOG = FilterCatalog(_PAINS_PARAMETERS)
_REACTIVE_PATTERNS = tuple(
    Chem.MolFromSmarts(pattern)
    for pattern in (
        "[CX3](=[OX1])[F,Cl,Br,I]",
        "[SX4](=[OX1])(=[OX1])[F,Cl,Br,I]",
        "[NX2]=[CX2]=[OX1]",
        "[NX1]#[NX2+][NX1-]",
    )
)
_ALLOWED_ATOMIC_NUMBERS = {5, 6, 7, 8, 9, 15, 16, 17, 35, 53}


def _tanimoto(first: bytes, second: bytes) -> float:
    first_bits = int.from_bytes(first)
    second_bits = int.from_bytes(second)
    intersection = (first_bits & second_bits).bit_count()
    union = (first_bits | second_bits).bit_count()
    return intersection / union if union else 1.0


def build_chembl_candidate(source_id: str, smiles: str) -> BackgroundCandidate | None:
    """Standardize and filter one ChEMBL small molecule for background use."""
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return None
    parent = rdMolStandardize.FragmentParent(rdMolStandardize.Cleanup(molecule))
    if not any(atom.GetAtomicNum() == 6 for atom in parent.GetAtoms()):
        return None
    if any(atom.GetAtomicNum() not in _ALLOWED_ATOMIC_NUMBERS for atom in parent.GetAtoms()):
        return None
    molecular_weight = float(Descriptors.MolWt(parent))
    if not 100 <= molecular_weight <= 900:
        return None
    formal_charge = int(Chem.GetFormalCharge(parent))
    if not -2 <= formal_charge <= 2:
        return None
    if any(
        label == "?"
        for _, label in Chem.FindMolChiralCenters(
            parent, includeUnassigned=True, useLegacyImplementation=False
        )
    ):
        return None
    if _PAINS_CATALOG.HasMatch(parent):
        return None
    if any(
        pattern is not None and parent.HasSubstructMatch(pattern)
        for pattern in _REACTIVE_PATTERNS
    ):
        return None
    canonical = str(Chem.MolToSmiles(parent, canonical=True, isomericSmiles=True))
    fingerprint = _FINGERPRINT_GENERATOR.GetFingerprint(parent)
    return BackgroundCandidate(
        source_id=source_id,
        canonical_smiles=canonical,
        molecular_weight=molecular_weight,
        clogp=float(Descriptors.MolLogP(parent)),
        hbd=int(Lipinski.NumHDonors(parent)),
        hba=int(Lipinski.NumHAcceptors(parent)),
        rotatable_bonds=int(Lipinski.NumRotatableBonds(parent)),
        formal_charge=formal_charge,
        fingerprint=bytes(DataStructs.BitVectToBinaryText(fingerprint)),
    )


def _maxmin_pick(
    candidates: list[BackgroundCandidate], count: int, *, seed: int
) -> list[BackgroundCandidate]:
    ordered = sorted(candidates, key=lambda item: item.source_id)
    if len(ordered) < count:
        raise ValueError(f"stratum has {len(ordered)} candidates; {count} required")
    first_index = seed % len(ordered)
    selected = [ordered[first_index]]
    remaining = ordered[:first_index] + ordered[first_index + 1 :]
    while len(selected) < count:
        ranked = sorted(
            remaining,
            key=lambda candidate: (
                -min(
                    1 - _tanimoto(candidate.fingerprint, chosen.fingerprint)
                    for chosen in selected
                ),
                candidate.source_id,
            ),
        )
        chosen = ranked[0]
        selected.append(chosen)
        remaining.remove(chosen)
    return selected


def select_background_panel(
    candidates: list[BackgroundCandidate],
    *,
    panel_size: int,
    strata: int,
    seed: int,
) -> list[BackgroundProbe]:
    """Select an equal MaxMin sample from molecular-weight quantile strata."""
    if panel_size <= 0 or strata <= 0 or panel_size % strata:
        raise ValueError("panel_size must be positive and divisible by strata")
    ordered = sorted(candidates, key=lambda item: (item.molecular_weight, item.source_id))
    if len(ordered) < panel_size:
        raise ValueError("candidate pool is smaller than requested panel")
    per_stratum = panel_size // strata
    panel: list[BackgroundProbe] = []
    for stratum in range(strata):
        start = len(ordered) * stratum // strata
        end = len(ordered) * (stratum + 1) // strata
        selected = _maxmin_pick(
            ordered[start:end], per_stratum, seed=seed + stratum
        )
        panel.extend(
            BackgroundProbe(
                source_id=item.source_id,
                canonical_smiles=item.canonical_smiles,
                molecular_weight=item.molecular_weight,
                clogp=item.clogp,
                hbd=item.hbd,
                hba=item.hba,
                rotatable_bonds=item.rotatable_bonds,
                formal_charge=item.formal_charge,
                mw_stratum=stratum,
            )
            for item in selected
        )
    return panel


_PANEL_FIELDS = (
    "source_id",
    "canonical_smiles",
    "molecular_weight",
    "clogp",
    "hbd",
    "hba",
    "rotatable_bonds",
    "formal_charge",
    "mw_stratum",
)


def write_panel_tsv(path: Path, panel: list[BackgroundProbe]) -> str:
    """Write the immutable background panel in a transparent tabular format."""
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=_PANEL_FIELDS, delimiter="\t")
    writer.writeheader()
    for probe in panel:
        row = probe.model_dump()
        row["molecular_weight"] = f"{probe.molecular_weight:.4f}"
        row["clogp"] = f"{probe.clogp:.4f}"
        writer.writerow(row)
    return write_bytes_atomic(path, buffer.getvalue().encode("utf-8"))


def validate_panel_tsv(path: Path) -> dict[str, int]:
    """Fail unless a panel contains 100 unique probes in ten equal strata."""
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 100:
        raise ValueError(f"background panel has {len(rows)} probes; 100 required")
    if len({row["source_id"] for row in rows}) != 100:
        raise ValueError("background source IDs must be unique")
    if len({row["canonical_smiles"] for row in rows}) != 100:
        raise ValueError("background structures must be unique")
    counts = Counter(int(row["mw_stratum"]) for row in rows)
    if counts != Counter({stratum: 10 for stratum in range(10)}):
        raise ValueError(f"background strata are unbalanced: {dict(counts)}")
    return {"probe_count": 100, "strata": 10, "per_stratum": 10}
