from pathlib import Path

import pytest
from pydantic import ValidationError

from airti_tf.contracts import (
    ArtifactRecord,
    BoltzRecord,
    DockingRecord,
    LigandRecord,
    MDRecord,
    PocketRecord,
    RankedTarget,
    StageStatus,
    TargetRecord,
    TargetStatus,
)


BASE = {
    "schema_version": "1.0",
    "input_sha256": "a" * 64,
    "tool_version": "test-1.0",
}


def test_unsupported_target_cannot_have_numeric_score() -> None:
    with pytest.raises(ValidationError, match="cannot receive a numeric score"):
        TargetRecord(
            **BASE,
            record_id="target:P00001",
            uniprot_id="P00001",
            sequence="MPEPTIDE",
            status="unsupported",
            unsupported_reason="no_structure",
            calibrated_score=0.0,
        )


def test_unsupported_target_requires_reason() -> None:
    with pytest.raises(ValidationError, match="requires unsupported_reason"):
        TargetRecord(
            **BASE,
            record_id="target:P00001",
            uniprot_id="P00001",
            sequence="MPEPTIDE",
            status=TargetStatus.UNSUPPORTED,
        )


def test_core_records_share_provenance_fields() -> None:
    records = [
        LigandRecord(
            **BASE,
            record_id="ligand:L1",
            ligand_id="L1",
            canonical_smiles="CCO",
            status=StageStatus.SUCCEEDED,
        ),
        PocketRecord(
            **BASE,
            record_id="pocket:P1",
            pocket_id="P1",
            target_id="P00533",
            status=StageStatus.SUCCEEDED,
        ),
        DockingRecord(
            **BASE,
            record_id="dock:D1",
            docking_id="D1",
            ligand_state_id="L1:S1",
            pocket_id="P1",
            seed=11,
            status=StageStatus.SUCCEEDED,
            affinity_kcal_mol=-8.2,
        ),
        BoltzRecord(
            **BASE,
            record_id="boltz:B1",
            boltz_id="B1",
            ligand_state_id="L1:S1",
            target_id="P00533",
            seed=11,
            status=StageStatus.SUCCEEDED,
        ),
        MDRecord(
            **BASE,
            record_id="md:M1",
            md_id="M1",
            ligand_state_id="L1:S1",
            target_id="P00533",
            replica=1,
            status=StageStatus.SUCCEEDED,
            completed_ns=100.0,
        ),
        RankedTarget(
            **BASE,
            record_id="rank:R1",
            target_id="P00533",
            ligand_id="L1",
            rank=1,
            status=StageStatus.SUCCEEDED,
            candidate_priority=0.91,
        ),
        ArtifactRecord(
            **BASE,
            record_id="artifact:A1",
            artifact_id="A1",
            task_id="T1",
            path=Path("results/A1.json"),
            sha256="b" * 64,
            status=StageStatus.SUCCEEDED,
        ),
    ]

    assert all(record.schema_version == "1.0" for record in records)
    assert all(len(record.input_sha256) == 64 for record in records)
    assert all(record.tool_version == "test-1.0" for record in records)


def test_hash_fields_reject_non_sha256_values() -> None:
    with pytest.raises(ValidationError):
        LigandRecord(
            schema_version="1.0",
            record_id="ligand:L1",
            input_sha256="short",
            tool_version="test-1.0",
            ligand_id="L1",
            canonical_smiles="CCO",
            status=StageStatus.SUCCEEDED,
        )
