"""PocketVina-compatible QuickVina2 and AutoDock Vina adapters."""

from __future__ import annotations

import re
import statistics
import subprocess
from collections.abc import Callable
from pathlib import Path
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator

from airti_tf.pockets.receptor import DockingBox

DockingStatus = Literal["succeeded", "failed"]
Backend = Literal["quickvina2", "vina"]


class InsufficientSeedsError(RuntimeError):
    """Raised when fewer than two of three fixed seeds complete."""


class DockingJob(BaseModel):
    model_config = ConfigDict(extra="forbid")

    job_id: str
    receptor_pdbqt: Path
    ligand_pdbqt: Path
    box: DockingBox
    output_dir: Path
    exhaustiveness: int = Field(default=8, gt=0)
    num_modes: int = Field(default=9, gt=0)
    cpu: int = Field(default=0, ge=0)


class DockingSeedResult(BaseModel):
    model_config = ConfigDict(extra="forbid")

    seed: int
    status: DockingStatus
    affinity_kcal_mol: float | None = None
    pose_count: int = Field(default=0, ge=0)
    error_code: str | None = None

    @model_validator(mode="after")
    def validate_state(self) -> "DockingSeedResult":
        if self.status == "succeeded" and self.affinity_kcal_mol is None:
            raise ValueError("successful docking seed requires affinity")
        if self.status == "failed" and not self.error_code:
            raise ValueError("failed docking seed requires error_code")
        return self


class DockingSummary(BaseModel):
    model_config = ConfigDict(extra="forbid")

    affinity_median: float
    affinity_range: float
    seed_success_count: int = Field(ge=2, le=3)
    successful_seeds: list[int]


class ExecutionResult(BaseModel):
    model_config = ConfigDict(extra="forbid")

    return_code: int
    stdout: str
    stderr: str
    timed_out: bool


Executor = Callable[[list[str], int], ExecutionResult]


def _build_command(executable: str, job: DockingJob, seed: int) -> list[str]:
    output = job.output_dir / f"{job.job_id}.seed{seed}.poses.pdbqt"
    log = job.output_dir / f"{job.job_id}.seed{seed}.log"
    command = [
        executable,
        "--receptor",
        str(job.receptor_pdbqt),
        "--ligand",
        str(job.ligand_pdbqt),
        "--center_x",
        f"{job.box.center[0]:.3f}",
        "--center_y",
        f"{job.box.center[1]:.3f}",
        "--center_z",
        f"{job.box.center[2]:.3f}",
        "--size_x",
        f"{job.box.size[0]:.3f}",
        "--size_y",
        f"{job.box.size[1]:.3f}",
        "--size_z",
        f"{job.box.size[2]:.3f}",
        "--exhaustiveness",
        str(job.exhaustiveness),
        "--num_modes",
        str(job.num_modes),
        "--seed",
        str(seed),
    ]
    if job.cpu:
        command.extend(["--cpu", str(job.cpu)])
    command.extend([
        "--out",
        str(output),
        "--log",
        str(log),
    ])
    return command


def build_quickvina_command(job: DockingJob, *, seed: int) -> list[str]:
    return _build_command("qvina2", job, seed)


def build_vina_command(job: DockingJob, *, seed: int) -> list[str]:
    return _build_command("vina", job, seed)


def parse_vina_log(path: Path, *, seed: int) -> DockingSeedResult:
    """Parse the standard Vina result table from a completed log."""
    rows = re.findall(
        r"^\s*\d+\s+(-?\d+(?:\.\d+)?)\s+[-+0-9.eE]+\s+[-+0-9.eE]+\s*$",
        path.read_text(encoding="utf-8"),
        flags=re.MULTILINE,
    )
    if not rows:
        return DockingSeedResult(seed=seed, status="failed", error_code="unparseable_output")
    affinities = [float(value) for value in rows]
    return DockingSeedResult(
        seed=seed,
        status="succeeded",
        affinity_kcal_mol=min(affinities),
        pose_count=len(affinities),
    )


def summarize_seeds(records: list[DockingSeedResult]) -> DockingSummary:
    """Summarize successful fixed seeds using the robust median."""
    successful = [record for record in records if record.status == "succeeded"]
    if len(successful) < 2:
        raise InsufficientSeedsError(f"only {len(successful)} of {len(records)} seeds succeeded")
    affinities = [
        record.affinity_kcal_mol
        for record in successful
        if record.affinity_kcal_mol is not None
    ]
    return DockingSummary(
        affinity_median=statistics.median(affinities),
        affinity_range=max(affinities) - min(affinities),
        seed_success_count=len(successful),
        successful_seeds=sorted(record.seed for record in successful),
    )


def classify_failure(return_code: int, stderr: str, *, timed_out: bool) -> str:
    """Map execution symptoms to stable technical error codes."""
    if timed_out:
        return "timeout"
    lowered = stderr.lower()
    if "pdbqt" in lowered and ("parse" in lowered or "parsing" in lowered):
        return "pdbqt_error"
    if "grid" in lowered or "box" in lowered:
        return "grid_error"
    if return_code in {137, 143}:
        return "resource_exhausted"
    return "nonzero_exit"


def _to_text(value: str | bytes | None) -> str:
    if value is None:
        return ""
    return value.decode(errors="replace") if isinstance(value, bytes) else value


def _subprocess_executor(command: list[str], timeout_seconds: int) -> ExecutionResult:
    try:
        result = subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
        )
    except subprocess.TimeoutExpired as error:
        return ExecutionResult(
            return_code=124,
            stdout=_to_text(error.stdout),
            stderr=_to_text(error.stderr),
            timed_out=True,
        )
    return ExecutionResult(
        return_code=result.returncode,
        stdout=result.stdout,
        stderr=result.stderr,
        timed_out=False,
    )


def run_docking_seed(
    job: DockingJob,
    *,
    seed: int,
    backend: Backend = "quickvina2",
    executor: Executor = _subprocess_executor,
    timeout_seconds: int = 600,
) -> DockingSeedResult:
    """Execute one seed through the backend-independent PocketVina contract."""
    job.output_dir.mkdir(parents=True, exist_ok=True)
    command = (
        build_quickvina_command(job, seed=seed)
        if backend == "quickvina2"
        else build_vina_command(job, seed=seed)
    )
    result = executor(command, timeout_seconds)
    if result.return_code != 0 or result.timed_out:
        return DockingSeedResult(
            seed=seed,
            status="failed",
            error_code=classify_failure(
                result.return_code, result.stderr, timed_out=result.timed_out
            ),
        )
    log_path = Path(command[command.index("--log") + 1])
    pose_path = Path(command[command.index("--out") + 1])
    if not log_path.is_file() or not pose_path.is_file():
        return DockingSeedResult(seed=seed, status="failed", error_code="output_missing")
    return parse_vina_log(log_path, seed=seed)
