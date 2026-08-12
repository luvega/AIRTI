from pathlib import Path

import pytest

from airti_tf.state import InvalidTransitionError, StateStore


@pytest.fixture
def state_store(tmp_path: Path) -> StateStore:
    return StateStore(tmp_path / "job_status.sqlite")


def test_completed_task_is_not_claimed_twice(state_store: StateStore) -> None:
    task_id = state_store.register("dock", input_hash="abc")
    state_store.transition(task_id, "running")
    state_store.transition(task_id, "succeeded", output_hash="def")

    assert state_store.claim(task_id) is False


def test_failed_task_can_be_retried_with_audit_entry(state_store: StateStore) -> None:
    task_id = state_store.register("boltz", input_hash="abc")
    state_store.transition(task_id, "running")
    state_store.transition(task_id, "failed", error_code="OOM")

    assert state_store.retry(task_id, max_attempts=2) is True
    assert state_store.history(task_id)[-1].status == "pending"
    assert state_store.history(task_id)[-1].error_code == "OOM"


def test_retry_stops_at_attempt_limit(state_store: StateStore) -> None:
    task_id = state_store.register("boltz", input_hash="abc")
    assert state_store.claim(task_id)
    state_store.transition(task_id, "failed", error_code="OOM")
    assert state_store.retry(task_id, max_attempts=1) is False


def test_invalid_transition_is_rejected(state_store: StateStore) -> None:
    task_id = state_store.register("dock", input_hash="abc")

    with pytest.raises(InvalidTransitionError, match="pending -> succeeded"):
        state_store.transition(task_id, "succeeded", output_hash="def")


def test_registration_is_idempotent_and_artifacts_use_hash_identity(
    state_store: StateStore,
) -> None:
    first = state_store.register("dock", input_hash="abc")
    second = state_store.register("dock", input_hash="abc")
    state_store.register_artifact(first, sha256="f" * 64, path="one/result.json")
    state_store.register_artifact(first, sha256="f" * 64, path="renamed/result.json")

    assert first == second
    assert state_store.task_count() == 1
    assert state_store.artifact_count() == 1


def test_schema_contains_required_audit_tables(state_store: StateStore) -> None:
    assert state_store.table_names() >= {"runs", "tasks", "artifacts", "transitions"}
