import json
from pathlib import Path

from airti_tf.targets.structures import StructureCandidate, choose_structure


def candidate(source: str, **overrides: object) -> StructureCandidate:
    payload: dict[str, object] = {
        "structure_id": f"{source}-1",
        "source": source,
        "coverage": 0.95,
        "sequence_identity": 1.0,
        "mainchain_missing_fraction": 0.0,
        "unsupported_chemistry": False,
        "has_ligand": False,
        "experimental_method": "x-ray" if source == "pdb" else None,
        "pae_supported": True,
    }
    payload.update(overrides)
    return StructureCandidate.model_validate(payload)


def test_prefers_usable_experimental_structure() -> None:
    chosen = choose_structure(
        [
            candidate("alphafold", coverage=0.99, confidence=0.88),
            candidate("pdb", coverage=0.85, resolution=2.1, has_ligand=True),
        ]
    )

    assert chosen.status == "ready"
    assert chosen.source == "pdb"
    assert chosen.structure_id == "pdb-1"


def test_fixture_selection_is_deterministic() -> None:
    payload = json.loads(Path("tests/fixtures/structure_candidates.json").read_text())
    candidates = [StructureCandidate.model_validate(item) for item in payload]

    first = choose_structure(candidates)
    second = choose_structure(list(reversed(candidates)))

    assert first == second
    assert first.structure_id == "1M17"


def test_apo_structure_must_meet_stricter_resolution_cutoff() -> None:
    result = choose_structure(
        [
            candidate("pdb", structure_id="apo-poor", resolution=2.9),
            candidate("alphafold", structure_id="af-good", confidence=0.91),
        ]
    )

    assert result.structure_id == "af-good"


def test_ligand_bound_cryo_em_structure_uses_membrane_aware_cutoff() -> None:
    chosen = choose_structure(
        [
            candidate(
                "pdb",
                structure_id="7A6F",
                experimental_method="electron_microscopy",
                resolution=3.5,
                has_ligand=True,
            ),
            candidate("alphafold", structure_id="af-abcb1", confidence=0.92),
        ]
    )

    assert chosen.structure_id == "7A6F"


def test_low_resolution_cryo_em_falls_back_to_alphafold() -> None:
    chosen = choose_structure(
        [
            candidate(
                "pdb",
                structure_id="cryo-poor",
                experimental_method="electron_microscopy",
                resolution=3.6,
                has_ligand=True,
            ),
            candidate("alphafold", structure_id="af-good", confidence=0.92),
        ]
    )

    assert chosen.structure_id == "af-good"


def test_no_usable_structure_is_preserved_as_unsupported() -> None:
    result = choose_structure([])

    assert result.status == "unsupported"
    assert result.unsupported_reason == "no_structure"
    assert result.score is None


def test_low_coverage_and_sequence_mismatch_have_explicit_reasons() -> None:
    low_coverage = choose_structure([candidate("pdb", coverage=0.69, resolution=2.0)])
    mismatch = choose_structure(
        [candidate("alphafold", sequence_identity=0.8, confidence=0.95)]
    )

    assert low_coverage.unsupported_reason == "low_coverage"
    assert mismatch.unsupported_reason == "sequence_mismatch"


def test_low_confidence_alphafold_is_not_scored_as_zero() -> None:
    result = choose_structure([candidate("alphafold", confidence=0.69)])

    assert result.status == "unsupported"
    assert result.unsupported_reason == "low_confidence"
    assert result.score is None


def test_unsupported_chemistry_fails_closed() -> None:
    result = choose_structure(
        [candidate("pdb", resolution=1.8, has_ligand=True, unsupported_chemistry=True)]
    )

    assert result.unsupported_reason == "unsupported_chemistry"
