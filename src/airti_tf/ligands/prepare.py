"""Deterministic ligand standardization, protomer enumeration, and 3D preparation."""

from __future__ import annotations

import hashlib
from typing import Literal

from dimorphite_dl import protonate_smiles
from pydantic import BaseModel, ConfigDict, Field
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

Profile = Literal["local", "production"]


class LigandState(BaseModel):
    model_config = ConfigDict(extra="forbid")

    ligand_state_id: str
    canonical_smiles: str
    ph_min: float
    ph_max: float
    atom_count: int = Field(gt=0)
    formal_charge: int
    mol_block: str


class LigandPreparationResult(BaseModel):
    model_config = ConfigDict(extra="forbid")

    ligand_id: str
    original_smiles: str
    status: Literal["succeeded", "failed"]
    fragment_count: int = 0
    states: list[LigandState] = Field(default_factory=list)
    error_code: str | None = None
    uncertainty_flags: list[str] = Field(default_factory=list)


def validate_query_batch(smiles: list[str]) -> list[str]:
    """Validate the service-level query cardinality and uniqueness contract."""
    if not 1 <= len(smiles) <= 5:
        raise ValueError("query batch must contain 1 to 5 molecules")
    normalized = [item.strip() for item in smiles]
    if any(not item for item in normalized):
        raise ValueError("query SMILES cannot be empty")
    if len(set(normalized)) != len(normalized):
        raise ValueError("query molecules must be unique")
    return normalized


def _failure(
    smiles: str,
    ligand_id: str,
    error_code: str,
    *,
    fragment_count: int = 0,
) -> LigandPreparationResult:
    return LigandPreparationResult(
        ligand_id=ligand_id,
        original_smiles=smiles,
        status="failed",
        fragment_count=fragment_count,
        error_code=error_code,
    )


def _embed(molecule: Chem.Mol) -> str | None:
    molecule_with_hydrogens = Chem.AddHs(molecule)
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = 20_260_812
    parameters.useRandomCoords = True
    if AllChem.EmbedMolecule(molecule_with_hydrogens, parameters) != 0:
        return None
    if AllChem.MMFFHasAllMoleculeParams(molecule_with_hydrogens):
        AllChem.MMFFOptimizeMolecule(molecule_with_hydrogens, maxIters=200)
    else:
        AllChem.UFFOptimizeMolecule(molecule_with_hydrogens, maxIters=200)
    return str(Chem.MolToMolBlock(molecule_with_hydrogens))


def prepare_ligand(smiles: str, *, profile: Profile) -> LigandPreparationResult:
    """Prepare pH 7.4 ± 1.0 protomer/tautomer states with explicit failure codes."""
    ligand_id = hashlib.sha256(smiles.encode()).hexdigest()
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return _failure(smiles, ligand_id, "invalid_smiles")

    fragment_count = len(Chem.GetMolFrags(molecule))
    parent = rdMolStandardize.FragmentParent(rdMolStandardize.Cleanup(molecule))
    Chem.AssignStereochemistry(parent, cleanIt=True, force=True)
    unassigned = [
        center
        for center, label in Chem.FindMolChiralCenters(
            parent, includeUnassigned=True, useLegacyImplementation=False
        )
        if label == "?"
    ]
    uncertainty_flags: list[str] = []
    if unassigned:
        if profile == "production":
            return _failure(
                smiles,
                ligand_id,
                "undefined_stereochemistry",
                fragment_count=fragment_count,
            )
        uncertainty_flags.append("undefined_stereochemistry")

    molecular_weight = Descriptors.MolWt(parent)
    if profile == "production" and not 100 <= molecular_weight <= 900:
        return _failure(
            smiles,
            ligand_id,
            "molecular_weight_out_of_range",
            fragment_count=fragment_count,
        )

    parent_smiles = Chem.MolToSmiles(parent, canonical=True, isomericSmiles=True)
    try:
        protomers = protonate_smiles(
            parent_smiles,
            ph_min=6.4,
            ph_max=8.4,
            precision=1.0,
            max_variants=16,
        )
    except Exception:
        return _failure(
            smiles, ligand_id, "protonation_failed", fragment_count=fragment_count
        )
    if not protomers:
        return _failure(
            smiles, ligand_id, "protonation_failed", fragment_count=fragment_count
        )

    enumerator = rdMolStandardize.TautomerEnumerator()
    enumerator.SetMaxTautomers(16)
    unique: dict[str, Chem.Mol] = {}
    for protomer_smiles in protomers:
        protomer = Chem.MolFromSmiles(protomer_smiles)
        if protomer is None:
            continue
        for tautomer in enumerator.Enumerate(protomer):
            canonical = Chem.MolToSmiles(tautomer, canonical=True, isomericSmiles=True)
            unique.setdefault(canonical, tautomer)
            if len(unique) == 16:
                break
        if len(unique) == 16:
            break

    if not unique:
        return _failure(
            smiles, ligand_id, "protonation_failed", fragment_count=fragment_count
        )
    max_atoms = max(molecule_state.GetNumAtoms() for molecule_state in unique.values())
    if max_atoms > 128:
        return _failure(
            smiles, ligand_id, "ligand_too_large", fragment_count=fragment_count
        )
    if max_atoms > 56:
        uncertainty_flags.append("boltz_high_atom_count")

    states: list[LigandState] = []
    for canonical, molecule_state in sorted(unique.items()):
        mol_block = _embed(molecule_state)
        if mol_block is None:
            return _failure(
                smiles,
                ligand_id,
                "conformer_generation_failed",
                fragment_count=fragment_count,
            )
        state_id = hashlib.sha256(f"{canonical}|6.4|8.4".encode()).hexdigest()
        states.append(
            LigandState(
                ligand_state_id=state_id,
                canonical_smiles=canonical,
                ph_min=6.4,
                ph_max=8.4,
                atom_count=molecule_state.GetNumAtoms(),
                formal_charge=Chem.GetFormalCharge(molecule_state),
                mol_block=mol_block,
            )
        )

    return LigandPreparationResult(
        ligand_id=ligand_id,
        original_smiles=smiles,
        status="succeeded",
        fragment_count=fragment_count,
        states=states,
        uncertainty_flags=sorted(set(uncertainty_flags)),
    )
