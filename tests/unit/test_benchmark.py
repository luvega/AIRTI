import pytest

from airti_tf.benchmark import (
    DataLeakageError,
    bootstrap_success_at_k,
    fit_weights,
    reciprocal_rank,
    recall_at_k,
    success_at_k,
    top_k_jaccard,
)


def test_success_at_k_uses_any_known_human_target() -> None:
    truth = {"drug1": {"P00533", "P04637"}}
    ranked = {"drug1": ["P11111", "P00533", "P22222"]}

    assert success_at_k(truth, ranked, k=2) == 1.0
    assert reciprocal_rank(truth, ranked) == pytest.approx(0.5)
    assert recall_at_k(truth, ranked, k=2) == pytest.approx(0.5)


def test_blind_set_cannot_be_used_to_fit_weights() -> None:
    with pytest.raises(DataLeakageError):
        fit_weights(dataset_role="blind")


def test_validation_set_cannot_be_used_to_fit_weights() -> None:
    with pytest.raises(DataLeakageError):
        fit_weights(dataset_role="retrospective_validation")


def test_top_k_jaccard_is_deterministic() -> None:
    first = ["A", "B", "C", "D"]
    second = ["B", "C", "E", "F"]

    assert top_k_jaccard(first, second, k=4) == pytest.approx(2 / 6)


def test_bootstrap_interval_is_seeded_and_contains_point_estimate() -> None:
    truth = {
        "d1": {"A"},
        "d2": {"B"},
        "d3": {"C"},
        "d4": {"D"},
    }
    ranked = {
        "d1": ["A"],
        "d2": ["X"],
        "d3": ["C"],
        "d4": ["Y"],
    }

    first = bootstrap_success_at_k(truth, ranked, k=1, iterations=500, seed=20260812)
    second = bootstrap_success_at_k(truth, ranked, k=1, iterations=500, seed=20260812)

    assert first == second
    assert first.estimate == pytest.approx(0.5)
    assert first.lower <= first.estimate <= first.upper


def test_empty_truth_is_rejected() -> None:
    with pytest.raises(ValueError, match="truth"):
        success_at_k({}, {}, k=10)

