import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

from airti_tf.cli import app
from airti_tf.runtime import GPUFact, HostFacts, ImageFact, run_preflight


@pytest.fixture
def fake_host() -> HostFacts:
    return HostFacts(
        data_root=Path("/data/airti-target-fishing"),
        cache_root=Path("/mnt/ssd4t/airti-target-fishing"),
        data_free_bytes=12 * 1024**4,
        cache_free_bytes=2 * 1024**4,
        commands={"docker": True, "docker_compose": True, "nextflow": False},
        versions={"docker_compose": "5.1.4"},
        nvidia_runtime=True,
        gpus=[GPUFact(memory_mib=49140, container_visible=True)],
        images={},
        reference_manifest_status="missing",
        open_file_limit=1_048_576,
    )


def test_production_preflight_reports_only_real_missing_layers(fake_host: HostFacts) -> None:
    result = run_preflight(profile="production", host=fake_host)

    assert result.ok is False
    assert set(result.blockers) >= {
        "nextflow_missing",
        "airti_images_missing",
    }
    assert "gpu_unavailable" not in result.blockers
    assert "data_space_below_1_tib" not in result.blockers
    assert "cache_space_below_500_gib" not in result.blockers
    assert result.exit_code == 2


def test_local_profile_allows_mock_execution(fake_host: HostFacts) -> None:
    result = run_preflight(profile="local", host=fake_host)

    assert result.ok is True
    assert result.exit_code == 0


def test_ready_production_host_passes(fake_host: HostFacts) -> None:
    fake_host.commands["nextflow"] = True
    fake_host.versions["nextflow"] = "24.10.4"
    fake_host.reference_manifest_status = "valid"
    fake_host.images = {
        name: ImageFact(command_ok=True, digest=f"sha256:{index:064x}")
        for index, name in enumerate(
            [
                "airti-targetlib-cpu",
                "airti-screening-cpu",
                "airti-boltz2-gpu",
                "airti-gromacs-gpu",
                "airti-orchestrator",
            ],
            start=1,
        )
    }

    result = run_preflight(profile="production", host=fake_host)

    assert result.ok is True
    assert result.exit_code == 0


def test_invalid_reference_manifest_is_configuration_error(fake_host: HostFacts) -> None:
    fake_host.reference_manifest_status = "invalid"

    result = run_preflight(profile="production", host=fake_host)

    assert result.ok is False
    assert result.exit_code == 3
    assert result.invalid == ["reference_manifest_invalid"]


def test_local_cli_writes_machine_readable_result(tmp_path: Path) -> None:
    output = tmp_path / "preflight.json"

    result = CliRunner().invoke(
        app,
        ["preflight", "--profile", "local", "--output", str(output)],
    )

    assert result.exit_code == 0
    payload = json.loads(output.read_text(encoding="utf-8"))
    assert payload["profile"] == "local"
    assert payload["ok"] is True
