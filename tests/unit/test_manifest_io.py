import json
from pathlib import Path

from airti_tf.manifest_io import (
    content_sha256,
    read_jsonl,
    write_artifact,
    write_jsonl_atomic,
)


def test_logically_identical_payloads_have_same_hash() -> None:
    assert content_sha256({"b": 2, "a": 1}) == content_sha256({"a": 1, "b": 2})


def test_jsonl_round_trip_uses_canonical_order(tmp_path: Path) -> None:
    path = tmp_path / "manifest.jsonl"
    records = [{"target": "P04637", "rank": 2}, {"rank": 1, "target": "P00533"}]

    digest = write_jsonl_atomic(path, records)

    assert read_jsonl(path) == records
    assert digest == content_sha256(path.read_bytes())
    first_line = path.read_text(encoding="utf-8").splitlines()[0]
    assert first_line == '{"rank":2,"target":"P04637"}'
    assert list(tmp_path.glob(".manifest.jsonl.*.tmp")) == []


def test_artifact_identity_is_content_hash_not_filename(tmp_path: Path) -> None:
    first = write_artifact(tmp_path / "first.json", {"answer": 42})
    second = write_artifact(tmp_path / "renamed.json", {"answer": 42})

    assert first.artifact_id == second.artifact_id
    assert first.artifact_id == first.sha256
    assert first.path != second.path
    assert json.loads(first.path.read_text(encoding="utf-8")) == {"answer": 42}

