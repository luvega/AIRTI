import hashlib

import pytest
from pydantic import ValidationError

from airti_tf.stages import CofactorRecord, TargetPocketRow


def _ready_payload() -> dict[str, object]:
    sequence = "MPEPTIDE"
    return {
        "schema_version": "1.1",
        "target_id": "P08684",
        "family": "cytochrome_p450",
        "status": "ready",
        "environment": "soluble_construct",
        "orientation_source": None,
        "cofactors": [
            {
                "ccd_id": "HEM",
                "role": "essential",
                "parameter_id": "p450-ferric-thiolate-v1",
            }
        ],
        "sequence": sequence,
        "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
        "model_sequence": sequence,
        "model_sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
        "model_sequence_start": 1,
        "model_sequence_end": len(sequence),
        "structure_quality": 0.9,
        "structure_id": "3NXU",
        "structure_source": "pdb",
        "structure_path": "targets/P08684/3NXU.pdb",
        "calibration_path": "targets/P08684/calibration.json",
        "pocket_id": "P08684:p1",
        "receptor_pdbqt_path": "targets/P08684/receptor.pdbqt",
        "box": {"center": [0, 0, 0], "size": [20, 20, 20]},
        "background_affinities": [-5.0 - index / 100 for index in range(100)],
        "msa_path": "targets/P08684/P08684.a3m",
        "msa_database_version": "single-sequence-2026_02",
        "pocket_residues": [1, 2],
        "model_pocket_residues": [1, 2],
    }


def test_v11_target_contract_preserves_heme_and_environment() -> None:
    row = TargetPocketRow.model_validate(_ready_payload())

    assert row.schema_version == "1.1"
    assert row.environment == "soluble_construct"
    assert row.cofactors == [
        CofactorRecord(
            ccd_id="HEM",
            role="essential",
            parameter_id="p450-ferric-thiolate-v1",
        )
    ]


def test_membrane_target_requires_orientation_provenance() -> None:
    payload = _ready_payload()
    payload["environment"] = "membrane"

    with pytest.raises(ValidationError, match="orientation"):
        TargetPocketRow.model_validate(payload)


def test_v10_target_contract_defaults_to_soluble_without_cofactors() -> None:
    payload = _ready_payload()
    payload["schema_version"] = "1.0"
    payload.pop("environment")
    payload.pop("orientation_source")
    payload.pop("cofactors")

    row = TargetPocketRow.model_validate(payload)

    assert row.environment == "soluble"
    assert row.cofactors == []
