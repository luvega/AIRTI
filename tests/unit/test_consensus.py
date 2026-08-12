import pytest

from airti_tf.ranking.consensus import (
    FROZEN_WEIGHTS,
    FrozenWeightError,
    TargetEvidence,
    rank_targets,
)


def target(**overrides: object) -> TargetEvidence:
    payload: dict[str, object] = {
        "target_id": "P00533",
        "ligand_id": "query-1",
        "status": "ready",
        "vina_score": 0.8,
        "docking_consistency": 0.9,
        "structure_quality": 0.85,
        "boltz_score": 0.75,
        "md_score": 0.7,
        "md_status": "stable",
        "successful_seeds": 3,
        "boltz_seed_spread": 0.05,
        "heavy_atom_count": 32,
    }
    payload.update(overrides)
    return TargetEvidence.model_validate(payload)


def test_missing_md_is_not_converted_to_zero() -> None:
    result = rank_targets(
        [target(md_score=None, md_status="failed")], stage="final"
    )

    assert result.ranked[0].md_score is None
    assert "missing_md" in result.ranked[0].uncertainty_flags
    assert result.ranked[0].evidence_tier == "partial_computational"


def test_unsupported_target_is_excluded_but_counted() -> None:
    result = rank_targets(
        [
            target(),
            target(
                target_id="Q99999",
                status="unsupported",
                vina_score=None,
                docking_consistency=None,
                structure_quality=None,
                boltz_score=None,
                md_score=None,
                md_status=None,
            ),
        ],
        stage="screen",
    )

    assert len(result.ranked) == 1
    assert result.coverage.unsupported == 1


def test_frozen_weights_reject_runtime_change() -> None:
    with pytest.raises(FrozenWeightError):
        FROZEN_WEIGHTS.final.vina = 0.9


def test_risk_penalty_is_capped_and_explained() -> None:
    result = rank_targets(
        [
            target(
                severe_clash=True,
                structure_low_confidence=True,
                boltz_seed_spread=0.8,
                heavy_atom_count=60,
                md_score=None,
                md_status="failed",
            )
        ],
        stage="final",
    )
    ranked = result.ranked[0]

    assert ranked.risk_penalty == pytest.approx(0.15)
    assert set(ranked.uncertainty_flags) >= {
        "severe_clash",
        "structure_low_confidence",
        "boltz_seed_disagreement",
        "large_ligand",
        "missing_md",
    }


def test_ties_are_broken_by_accession_after_evidence_and_quality() -> None:
    result = rank_targets(
        [target(target_id="Q99999"), target(target_id="P00533")], stage="final"
    )

    assert [item.target_id for item in result.ranked] == ["P00533", "Q99999"]
    assert all(item.weight_version == "retrospective-v1" for item in result.ranked)

