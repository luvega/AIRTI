import json
import shutil
import subprocess
from pathlib import Path

import pytest
import yaml

PRODUCTION_IMAGE = "airti-tf:0.1.0-gpu"
EXPECTED_COMMANDS = [
    ["airti-tf", "version"],
    ["fpocket", "--help"],
    ["qvina2", "--help"],
    ["vina", "--help"],
    ["mk_prepare_ligand.py", "--help"],
    ["mk_prepare_receptor.py", "--help"],
    ["boltz", "--help"],
    ["gmx", "--version"],
    ["tleap", "-h"],
    ["antechamber", "-h"],
    ["parmed", "--version"],
]


def docker_available() -> bool:
    if shutil.which("docker") is None:
        return False
    result = subprocess.run(
        ["docker", "info"], capture_output=True, text=True, timeout=20
    )
    return result.returncode == 0


def image_exists(image: str) -> bool:
    result = subprocess.run(
        ["docker", "image", "inspect", image],
        capture_output=True,
        text=True,
        timeout=20,
    )
    return result.returncode == 0


def test_unified_dockerfile_pins_every_base_stage_by_digest() -> None:
    dockerfile = Path("containers/airti.Dockerfile")

    assert dockerfile.exists()
    from_lines = [
        line
        for line in dockerfile.read_text(encoding="utf-8").splitlines()
        if line.startswith("FROM ")
    ]
    assert from_lines
    assert all("@sha256:" in line for line in from_lines), from_lines


def test_cuda_conda_solver_override_matches_locked_cuda_runtime() -> None:
    dockerfile = Path("containers/airti.Dockerfile").read_text(encoding="utf-8")

    assert "CONDA_OVERRIDE_CUDA=12.8" in dockerfile


def test_boltz_cuda_dependencies_cannot_upgrade_locked_torch() -> None:
    dockerfile = Path("containers/airti.Dockerfile").read_text(encoding="utf-8")
    constraints = Path("containers/boltz.constraints.txt")

    assert constraints.exists()
    locked = set(constraints.read_text(encoding="utf-8").splitlines())
    assert {
        "torch==2.7.1",
        "cuequivariance==0.5.1",
        "cuequivariance-torch==0.5.1",
        "cuequivariance-ops-cu12==0.5.1",
        "cuequivariance-ops-torch-cu12==0.5.1",
    } <= locked
    assert "--constraint /tmp/boltz.constraints.txt" in dockerfile


def test_runtime_bypasses_broken_gromacs_activation_hook() -> None:
    dockerfile = Path("containers/airti.Dockerfile").read_text(encoding="utf-8")

    assert "ENV PATH=/opt/conda/bin:${PATH}" in dockerfile
    assert "MAMBA_DOCKERFILE_ACTIVATE=1" not in dockerfile
    assert "ENTRYPOINT []" in dockerfile


def test_non_root_scientific_caches_are_redirected_to_writable_runtime_paths() -> None:
    dockerfile = Path("containers/airti.Dockerfile").read_text(encoding="utf-8")

    assert "NUMBA_CACHE_DIR=/tmp/airti-cache/numba" in dockerfile
    assert "MPLCONFIGDIR=/tmp/airti-cache/matplotlib" in dockerfile
    assert "XDG_CACHE_HOME=/tmp/airti-cache/xdg" in dockerfile
    assert "TRITON_CACHE_DIR=/tmp/airti-cache/triton" in dockerfile
    assert "mkdir -p /tmp/airti-cache/numba" in dockerfile
    assert "chmod 1777 /tmp/airti-cache" in dockerfile
    assert "USER $MAMBA_USER" in dockerfile


def test_lockfile_declares_one_image_and_all_scientific_versions() -> None:
    lock = yaml.safe_load(Path("containers/images.lock.yaml").read_text())

    assert lock["tools"] == {
        "ambertools": "24.8",
        "autodock_vina": "1.2.7",
        "boltz": "2.2.1",
        "fpocket": "4.2.3",
        "gromacs": "2025.4",
        "meeko": "0.7.1",
        "parmed": "4.3.1",
        "qvina": "2.1.0",
    }
    assert set(lock["images"]) == {"production"}
    assert lock["images"]["production"]["tag"] == PRODUCTION_IMAGE
    digest = lock["images"]["production"]["digest"]
    assert digest.startswith("sha256:")
    assert Path(lock["images"]["production"]["sbom"]).is_file()
    assert lock["images"]["production"]["smoke_status"] == "passed"


def test_real_hardware_smoke_inputs_are_versioned() -> None:
    required = {
        "tests/fixtures/smoke/boltz/input.yaml",
        "tests/fixtures/smoke/docking/ligand.pdbqt",
        "tests/fixtures/smoke/docking/receptor.pdbqt",
        "tests/fixtures/smoke/gromacs/md.mdp",
        "tests/fixtures/smoke/gromacs/topol.top",
    }

    assert {path for path in required if not Path(path).is_file()} == set()


def test_hardware_smoke_script_runs_all_three_real_engines() -> None:
    script = Path("scripts/run_hardware_smoke.sh")

    assert script.is_file()
    source = script.read_text(encoding="utf-8")
    assert "set -euo pipefail" in source
    assert 'IMAGE="${AIRTI_IMAGE:-airti-tf:0.1.0-gpu}"' in source
    assert "qvina2 --receptor" in source
    assert "boltz predict" in source
    assert "gmx mdrun" in source
    assert "GPU available: True" in source
    assert "Finished mdrun" in source
    assert "sha256sum" in source


def test_boltz_model_artifacts_are_hash_locked() -> None:
    lock = yaml.safe_load(Path("containers/models.lock.yaml").read_text())

    assert lock["boltz2"]["version"] == "2.2.1"
    for artifact in ("confidence", "affinity"):
        assert lock["boltz2"][artifact]["sha256"].startswith("sha256:")


def test_production_profile_uses_one_unified_image_parameter() -> None:
    config = Path("workflow/nextflow.config").read_text(encoding="utf-8")

    assert config.count("params.production_image") >= 6
    assert "base_image" not in config
    assert "screening_image" not in config
    assert "boltz2_image" not in config
    assert "gromacs_image" not in config
    assert "latest" not in config


@pytest.mark.skipif(not docker_available(), reason="Docker daemon unavailable")
def test_unified_image_exposes_every_required_command() -> None:
    if not image_exists(PRODUCTION_IMAGE):
        pytest.skip(f"image not built: {PRODUCTION_IMAGE}")

    for command in EXPECTED_COMMANDS:
        result = subprocess.run(
            ["docker", "run", "--rm", PRODUCTION_IMAGE, *command],
            capture_output=True,
            text=True,
            timeout=120,
        )
        assert result.returncode == 0, json.dumps(
            {
                "image": PRODUCTION_IMAGE,
                "command": command,
                "stdout": result.stdout[-2000:],
                "stderr": result.stderr[-2000:],
            },
            indent=2,
        )


@pytest.mark.skipif(not docker_available(), reason="Docker daemon unavailable")
def test_unified_image_sees_rtx4090_and_cuda_torch() -> None:
    if not image_exists(PRODUCTION_IMAGE):
        pytest.skip(f"image not built: {PRODUCTION_IMAGE}")
    result = subprocess.run(
        [
            "docker",
            "run",
            "--rm",
            "--gpus",
            "all",
            PRODUCTION_IMAGE,
            "python",
            "-c",
            (
                "import torch; "
                "assert torch.cuda.is_available(); "
                "print(torch.cuda.get_device_name(0)); "
                "print(torch.version.cuda)"
            ),
        ],
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert result.returncode == 0, result.stderr
    assert "RTX 4090" in result.stdout
