"""SQLite-backed task state with explicit, auditable transitions."""

from __future__ import annotations

import hashlib
import re
import sqlite3
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Literal, cast

TaskStatus = Literal["pending", "running", "succeeded", "failed", "skipped"]

ALLOWED_TRANSITIONS: set[tuple[TaskStatus, TaskStatus]] = {
    ("pending", "running"),
    ("running", "succeeded"),
    ("running", "failed"),
    ("failed", "pending"),
    ("pending", "skipped"),
}


class InvalidTransitionError(RuntimeError):
    """Raised for a transition outside the documented state machine."""


@dataclass(frozen=True)
class TransitionEvent:
    status: TaskStatus
    previous_status: TaskStatus
    error_code: str | None
    output_hash: str | None
    created_at: str


def _now() -> str:
    return datetime.now(UTC).isoformat()


class StateStore:
    """Durable store for runs, tasks, artifacts, and transition history."""

    def __init__(self, path: Path) -> None:
        self.path = path
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._initialize()

    def _connect(self) -> sqlite3.Connection:
        connection = sqlite3.connect(self.path, timeout=30)
        connection.row_factory = sqlite3.Row
        connection.execute("PRAGMA foreign_keys = ON")
        connection.execute("PRAGMA journal_mode = WAL")
        return connection

    def _initialize(self) -> None:
        with self._connect() as connection:
            connection.executescript(
                """
                CREATE TABLE IF NOT EXISTS runs (
                    run_id TEXT PRIMARY KEY,
                    created_at TEXT NOT NULL
                );
                CREATE TABLE IF NOT EXISTS tasks (
                    task_id TEXT PRIMARY KEY,
                    run_id TEXT NOT NULL REFERENCES runs(run_id),
                    stage TEXT NOT NULL,
                    input_hash TEXT NOT NULL,
                    status TEXT NOT NULL,
                    attempt_count INTEGER NOT NULL DEFAULT 0,
                    output_hash TEXT,
                    error_code TEXT,
                    created_at TEXT NOT NULL,
                    updated_at TEXT NOT NULL,
                    UNIQUE(run_id, stage, input_hash)
                );
                CREATE TABLE IF NOT EXISTS artifacts (
                    sha256 TEXT PRIMARY KEY,
                    task_id TEXT NOT NULL REFERENCES tasks(task_id),
                    path TEXT NOT NULL,
                    created_at TEXT NOT NULL
                );
                CREATE TABLE IF NOT EXISTS transitions (
                    transition_id INTEGER PRIMARY KEY AUTOINCREMENT,
                    task_id TEXT NOT NULL REFERENCES tasks(task_id),
                    from_status TEXT NOT NULL,
                    to_status TEXT NOT NULL,
                    error_code TEXT,
                    output_hash TEXT,
                    created_at TEXT NOT NULL
                );
                """
            )

    def register(self, stage: str, input_hash: str, run_id: str = "default") -> str:
        """Register a deterministic task and return its id."""
        task_id = hashlib.sha256(f"{run_id}\0{stage}\0{input_hash}".encode()).hexdigest()
        timestamp = _now()
        with self._connect() as connection:
            connection.execute(
                "INSERT OR IGNORE INTO runs(run_id, created_at) VALUES (?, ?)",
                (run_id, timestamp),
            )
            connection.execute(
                """
                INSERT OR IGNORE INTO tasks(
                    task_id, run_id, stage, input_hash, status, created_at, updated_at
                ) VALUES (?, ?, ?, ?, 'pending', ?, ?)
                """,
                (task_id, run_id, stage, input_hash, timestamp, timestamp),
            )
        return task_id

    def transition(
        self,
        task_id: str,
        status: TaskStatus,
        *,
        output_hash: str | None = None,
        error_code: str | None = None,
    ) -> None:
        """Apply one legal state change and append its audit entry."""
        with self._connect() as connection:
            connection.execute("BEGIN IMMEDIATE")
            row = connection.execute(
                "SELECT status, attempt_count, error_code FROM tasks WHERE task_id = ?",
                (task_id,),
            ).fetchone()
            if row is None:
                raise KeyError(task_id)
            previous = cast(TaskStatus, row["status"])
            if (previous, status) not in ALLOWED_TRANSITIONS:
                raise InvalidTransitionError(f"invalid transition: {previous} -> {status}")
            event_error = error_code
            if previous == "failed" and status == "pending" and event_error is None:
                event_error = cast(str | None, row["error_code"])
            attempts = int(row["attempt_count"]) + (status == "running")
            timestamp = _now()
            connection.execute(
                """
                UPDATE tasks
                SET status = ?, attempt_count = ?, output_hash = ?, error_code = ?, updated_at = ?
                WHERE task_id = ?
                """,
                (
                    status,
                    attempts,
                    output_hash,
                    None if status == "pending" else error_code,
                    timestamp,
                    task_id,
                ),
            )
            connection.execute(
                """
                INSERT INTO transitions(
                    task_id, from_status, to_status, error_code, output_hash, created_at
                ) VALUES (?, ?, ?, ?, ?, ?)
                """,
                (task_id, previous, status, event_error, output_hash, timestamp),
            )

    def claim(self, task_id: str) -> bool:
        """Atomically claim a pending task; completed tasks return false."""
        with self._connect() as connection:
            connection.execute("BEGIN IMMEDIATE")
            row = connection.execute(
                "SELECT status, attempt_count FROM tasks WHERE task_id = ?", (task_id,)
            ).fetchone()
            if row is None:
                raise KeyError(task_id)
            if row["status"] != "pending":
                return False
            timestamp = _now()
            connection.execute(
                """
                UPDATE tasks
                SET status = 'running', attempt_count = ?, updated_at = ?
                WHERE task_id = ?
                """,
                (int(row["attempt_count"]) + 1, timestamp, task_id),
            )
            connection.execute(
                """
                INSERT INTO transitions(task_id, from_status, to_status, created_at)
                VALUES (?, 'pending', 'running', ?)
                """,
                (task_id, timestamp),
            )
            return True

    def retry(self, task_id: str, *, max_attempts: int) -> bool:
        """Return a failed task to pending while attempts remain."""
        with self._connect() as connection:
            row = connection.execute(
                "SELECT status, attempt_count FROM tasks WHERE task_id = ?", (task_id,)
            ).fetchone()
        if row is None:
            raise KeyError(task_id)
        if row["status"] != "failed" or int(row["attempt_count"]) >= max_attempts:
            return False
        self.transition(task_id, "pending")
        return True

    def history(self, task_id: str) -> list[TransitionEvent]:
        """Return transition history in insertion order."""
        with self._connect() as connection:
            rows = connection.execute(
                """
                SELECT from_status, to_status, error_code, output_hash, created_at
                FROM transitions WHERE task_id = ? ORDER BY transition_id
                """,
                (task_id,),
            ).fetchall()
        return [
            TransitionEvent(
                status=cast(TaskStatus, row["to_status"]),
                previous_status=cast(TaskStatus, row["from_status"]),
                error_code=cast(str | None, row["error_code"]),
                output_hash=cast(str | None, row["output_hash"]),
                created_at=cast(str, row["created_at"]),
            )
            for row in rows
        ]

    def register_artifact(self, task_id: str, *, sha256: str, path: str) -> None:
        """Register an artifact using its content digest as primary identity."""
        if re.fullmatch(r"[0-9a-f]{64}", sha256) is None:
            raise ValueError("artifact sha256 must be 64 lowercase hexadecimal characters")
        with self._connect() as connection:
            connection.execute(
                """
                INSERT OR IGNORE INTO artifacts(sha256, task_id, path, created_at)
                VALUES (?, ?, ?, ?)
                """,
                (sha256, task_id, path, _now()),
            )

    def task_count(self) -> int:
        with self._connect() as connection:
            row = connection.execute("SELECT COUNT(*) AS count FROM tasks").fetchone()
        return int(row["count"])

    def artifact_count(self) -> int:
        with self._connect() as connection:
            row = connection.execute("SELECT COUNT(*) AS count FROM artifacts").fetchone()
        return int(row["count"])

    def table_names(self) -> set[str]:
        with self._connect() as connection:
            rows = connection.execute(
                "SELECT name FROM sqlite_master WHERE type = 'table'"
            ).fetchall()
        return {cast(str, row["name"]) for row in rows}
