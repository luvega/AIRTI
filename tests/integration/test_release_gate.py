import json
from pathlib import Path

import pytest

from airti_tf.benchmark import (
    DataLeakageError,
    ReleaseMetrics,
    evaluate_release,
    freeze_weights,
    write_release_decision,
)


def metrics(**overrides: object) -> ReleaseMetrics:
    payload: dict[str, object] = {
        "blind_success_at_100": 0.30,
        "technical_success_rate": 0.95,
        "boltz_success_at_k": 0.40,
        "vina_success_at_k": 0.35,
        "successful_target_families": 3,
        "median_top20_jaccard": 0.70,
        "failures_with_error_codes": True,
        "report_integrity_passed": True,
    }
    payload.update(overrides)
    return ReleaseMetrics.model_validate(payload)


def test_release_gate_accepts_exact_boundaries() -> None:
    decision = evaluate_release(metrics())

    assert decision.status == "pass"
    assert all(item.passed for item in decision.criteria)


def test_release_gate_rejects_02999_success_at_100() -> None:
    decision = evaluate_release(metrics(blind_success_at_100=0.2999))

    assert decision.status == "fail"
    failed = {item.criterion for item in decision.criteria if not item.passed}
    assert failed == {"blind_success_at_100"}


def test_boltz_cannot_degrade_vina_baseline() -> None:
    decision = evaluate_release(
        metrics(boltz_success_at_k=0.34, vina_success_at_k=0.35)
    )

    assert decision.status == "fail"
    assert any(
        item.criterion == "boltz_not_worse_than_vina" and not item.passed
        for item in decision.criteria
    )


def test_blind_evaluation_requires_frozen_weights(tmp_path: Path) -> None:
    with pytest.raises(DataLeakageError, match="frozen"):
        evaluate_release(metrics(), dataset_role="blind", frozen_weights_path=tmp_path / "missing.yaml")


def test_release_decision_is_bound_to_frozen_weight_hash(tmp_path: Path) -> None:
    frozen = freeze_weights(
        {"final": {"vina": 0.25, "boltz": 0.30, "md": 0.30, "structure_quality": 0.15}},
        tmp_path / "frozen_weights.yaml",
    )
    decision = evaluate_release(
        metrics(), dataset_role="blind", frozen_weights_path=frozen.path
    )
    output = write_release_decision(decision, tmp_path / "release_decision.json")
    payload = json.loads(output.read_text(encoding="utf-8"))

    assert payload["status"] == "pass"
    assert payload["frozen_weights_sha256"] == frozen.sha256
    assert payload["decision_sha256"]

