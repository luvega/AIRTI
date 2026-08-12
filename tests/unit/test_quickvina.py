from pathlib import Path

import pytest

from airti_tf.pockets.receptor import DockingBox
from airti_tf.screening.quickvina import (
    DockingJob,
    DockingSeedResult,
    ExecutionResult,
    InsufficientSeedsError,
    build_quickvina_command,
    build_vina_command,
    classify_failure,
    parse_vina_log,
    run_docking_seed,
    summarize_seeds,
)


@pytest.fixture
def job(tmp_path: Path) -> DockingJob:
    return DockingJob(
        job_id="D1",
        receptor_pdbqt=Path("receptor.pdbqt"),
        ligand_pdbqt=Path("ligand.pdbqt"),
        box=DockingBox(center=(1, 2, 3), size=(18, 20, 22)),
        output_dir=tmp_path,
        exhaustiveness=8,
        num_modes=9,
    )


def docking(seed: int, affinity: float) -> DockingSeedResult:
    return DockingSeedResult(
        seed=seed,
        status="succeeded",
        affinity_kcal_mol=affinity,
        pose_count=9,
    )


def test_command_contains_locked_seed_and_box(job: DockingJob) -> None:
    command = build_quickvina_command(job, seed=29)

    assert command[0] == "qvina2"
    assert ["--seed", "29"] == command[command.index("--seed") : command.index("--seed") + 2]
    assert command[command.index("--center_x") + 1] == "1.000"
    assert command[command.index("--size_z") + 1] == "22.000"
    assert command[command.index("--exhaustiveness") + 1] == "8"


def test_vina_control_uses_same_job_contract(job: DockingJob) -> None:
    quickvina = build_quickvina_command(job, seed=11)
    vina = build_vina_command(job, seed=11)

    assert quickvina[1:] == vina[1:]
    assert vina[0] == "vina"


def test_parser_extracts_best_affinity_and_pose_count() -> None:
    result = parse_vina_log(Path("tests/fixtures/quickvina_output.txt"), seed=11)

    assert result.status == "succeeded"
    assert result.affinity_kcal_mol == pytest.approx(-8.2)
    assert result.pose_count == 3


def test_three_seed_summary_uses_median() -> None:
    result = summarize_seeds(
        [docking(11, -8.0), docking(29, -7.0), docking(47, -9.0)]
    )

    assert result.affinity_median == -8.0
    assert result.seed_success_count == 3
    assert result.affinity_range == 2.0


def test_at_least_two_seeds_are_required() -> None:
    records = [
        docking(11, -8.0),
        DockingSeedResult(seed=29, status="failed", error_code="timeout"),
        DockingSeedResult(seed=47, status="failed", error_code="nonzero_exit"),
    ]

    with pytest.raises(InsufficientSeedsError, match="1 of 3"):
        summarize_seeds(records)


def test_adapter_executes_and_parses_one_seed(job: DockingJob) -> None:
    fixture = Path("tests/fixtures/quickvina_output.txt").read_text(encoding="utf-8")

    def executor(command: list[str], timeout_seconds: int) -> ExecutionResult:
        assert timeout_seconds == 600
        Path(command[command.index("--log") + 1]).write_text(fixture, encoding="utf-8")
        Path(command[command.index("--out") + 1]).write_text("POSE", encoding="utf-8")
        return ExecutionResult(return_code=0, stdout="", stderr="", timed_out=False)

    result = run_docking_seed(job, seed=47, executor=executor)

    assert result.status == "succeeded"
    assert result.affinity_kcal_mol == pytest.approx(-8.2)


def test_adapter_classifies_failed_execution(job: DockingJob) -> None:
    def executor(command: list[str], timeout_seconds: int) -> ExecutionResult:
        del command, timeout_seconds
        return ExecutionResult(
            return_code=1, stdout="", stderr="PDBQT parsing error", timed_out=False
        )

    result = run_docking_seed(job, seed=11, executor=executor)

    assert result.status == "failed"
    assert result.error_code == "pdbqt_error"


@pytest.mark.parametrize(
    ("return_code", "stderr", "timed_out", "expected"),
    [
        (0, "", True, "timeout"),
        (1, "PDBQT parsing error", False, "pdbqt_error"),
        (1, "grid dimensions invalid", False, "grid_error"),
        (137, "", False, "resource_exhausted"),
        (1, "unexpected", False, "nonzero_exit"),
    ],
)
def test_failure_codes_are_specific(
    return_code: int, stderr: str, timed_out: bool, expected: str
) -> None:
    assert classify_failure(return_code, stderr, timed_out=timed_out) == expected
