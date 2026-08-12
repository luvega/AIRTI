"""Canonical hashing and crash-safe manifest writes."""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable


def canonical_json_bytes(payload: Any) -> bytes:
    """Serialize JSON-compatible data deterministically."""
    return json.dumps(
        payload,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")


def content_sha256(payload: Any) -> str:
    """Return the SHA-256 of raw bytes or canonical JSON data."""
    raw = payload if isinstance(payload, bytes) else canonical_json_bytes(payload)
    return hashlib.sha256(raw).hexdigest()


def _atomic_write(path: Path, raw: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    temporary_path = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(raw)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_path, path)
        directory_descriptor = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(directory_descriptor)
        finally:
            os.close(directory_descriptor)
    except BaseException:
        temporary_path.unlink(missing_ok=True)
        raise


def write_bytes_atomic(path: Path, raw: bytes) -> str:
    """Write arbitrary bytes atomically and return their content digest."""
    _atomic_write(path, raw)
    return content_sha256(raw)


def write_jsonl_atomic(path: Path, records: Iterable[dict[str, Any]]) -> str:
    """Write a canonical JSONL manifest atomically and return its digest."""
    raw = b"".join(canonical_json_bytes(record) + b"\n" for record in records)
    _atomic_write(path, raw)
    return content_sha256(raw)


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    """Read a JSONL manifest."""
    records: list[dict[str, Any]] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                decoded = json.loads(line)
                if not isinstance(decoded, dict):
                    raise ValueError(f"manifest row is not an object: {path}")
                records.append(decoded)
    return records


@dataclass(frozen=True)
class WrittenArtifact:
    artifact_id: str
    sha256: str
    path: Path


def write_artifact(path: Path, payload: Any) -> WrittenArtifact:
    """Write canonical JSON and identify the artifact by content, not path."""
    raw = canonical_json_bytes(payload) + b"\n"
    digest = content_sha256(raw)
    _atomic_write(path, raw)
    return WrittenArtifact(artifact_id=digest, sha256=digest, path=path)
