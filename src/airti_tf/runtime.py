"""Runtime discovery and fail-closed production readiness checks."""

from __future__ import annotations

import hashlib
import json
import re
import resource
import shutil
import subprocess
from pathlib import Path
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field

Profile = Literal["local", "production"]

REQUIRED_IMAGES = ("airti-tf:0.1.0-gpu",)


class MutableModel(BaseModel):
    model_config = ConfigDict(validate_assignment=True)


class GPUFact(MutableModel):
    memory_mib: int = Field(ge=0)
    container_visible: bool = False


class ImageFact(MutableModel):
    command_ok: bool = False
    digest: str | None = None


class HostFacts(MutableModel):
    data_root: Path = Path("/data/airti-target-fishing")
    cache_root: Path = Path("/mnt/ssd4t/airti-target-fishing")
    data_free_bytes: int = 0
    cache_free_bytes: int = 0
    commands: dict[str, bool] = Field(default_factory=dict)
    versions: dict[str, str] = Field(default_factory=dict)
    nvidia_runtime: bool = False
    gpus: list[GPUFact] = Field(default_factory=list)
    images: dict[str, ImageFact] = Field(default_factory=dict)
    reference_manifest_status: Literal["valid", "missing", "invalid"] = "missing"
    open_file_limit: int = 0


class PreflightResult(BaseModel):
    profile: Profile
    ok: bool
    exit_code: Literal[0, 2, 3]
    blockers: list[str]
    invalid: list[str]
    checks: dict[str, bool]


def _version_tuple(value: str) -> tuple[int, ...]:
    match = re.search(r"(\d+)\.(\d+)(?:\.(\d+))?", value)
    if match is None:
        return ()
    return tuple(int(part) for part in match.groups(default="0"))


def _under(path: Path, root: Path) -> bool:
    try:
        path.resolve(strict=False).relative_to(root)
    except ValueError:
        return False
    return True


def _image_digest_locked(fact: ImageFact | None) -> bool:
    return fact is not None and fact.digest is not None and fact.digest.startswith("sha256:")


def run_preflight(*, profile: Profile, host: HostFacts) -> PreflightResult:
    """Evaluate host facts without mutating the host."""
    if profile == "local":
        return PreflightResult(
            profile=profile,
            ok=True,
            exit_code=0,
            blockers=[],
            invalid=[],
            checks={"mock_execution": True},
        )

    blockers: list[str] = []
    invalid: list[str] = []
    checks: dict[str, bool] = {}

    checks["data_root"] = _under(host.data_root, Path("/data"))
    checks["data_space"] = host.data_free_bytes >= 1024**4
    checks["cache_root"] = _under(host.cache_root, Path("/mnt/ssd4t"))
    checks["cache_space"] = host.cache_free_bytes >= 500 * 1024**3
    checks["nextflow"] = host.commands.get("nextflow", False)
    checks["nextflow_version"] = _version_tuple(host.versions.get("nextflow", "")) >= (
        24,
        10,
        0,
    )
    checks["docker"] = host.commands.get("docker", False)
    checks["docker_compose"] = host.commands.get("docker_compose", False)
    checks["nvidia_runtime"] = host.nvidia_runtime
    checks["gpu"] = any(
        gpu.memory_mib >= 40 * 1024 and gpu.container_visible for gpu in host.gpus
    )
    missing_images = [
        name
        for name in REQUIRED_IMAGES
        if name not in host.images or not host.images[name].command_ok
    ]
    checks["airti_images"] = not missing_images
    checks["image_digests"] = bool(host.images) and all(
        _image_digest_locked(host.images.get(name)) for name in REQUIRED_IMAGES
    )
    checks["reference_manifest"] = host.reference_manifest_status == "valid"
    checks["open_file_limit"] = host.open_file_limit >= 65_536

    blocker_by_check = {
        "data_root": "artifact_root_outside_data",
        "data_space": "data_space_below_1_tib",
        "cache_root": "cache_root_outside_ssd",
        "cache_space": "cache_space_below_500_gib",
        "nextflow": "nextflow_missing",
        "nextflow_version": "nextflow_version_unsupported",
        "docker": "docker_missing",
        "docker_compose": "docker_compose_missing",
        "nvidia_runtime": "nvidia_runtime_missing",
        "gpu": "gpu_unavailable",
        "airti_images": "airti_images_missing",
        "image_digests": "image_digests_unlocked",
        "reference_manifest": "reference_manifest_missing",
        "open_file_limit": "open_file_limit_below_65536",
    }
    for check, blocker in blocker_by_check.items():
        if not checks[check]:
            blockers.append(blocker)

    if not checks["nextflow"] and "nextflow_version_unsupported" in blockers:
        blockers.remove("nextflow_version_unsupported")
    if not checks["airti_images"] and "image_digests_unlocked" in blockers:
        blockers.remove("image_digests_unlocked")
    if host.reference_manifest_status == "invalid":
        blockers = [item for item in blockers if item != "reference_manifest_missing"]
        invalid.append("reference_manifest_invalid")

    exit_code: Literal[0, 2, 3] = 3 if invalid else (2 if blockers else 0)
    return PreflightResult(
        profile=profile,
        ok=exit_code == 0,
        exit_code=exit_code,
        blockers=blockers,
        invalid=invalid,
        checks=checks,
    )


def _run(command: list[str], *, timeout: int = 20) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        check=False,
        capture_output=True,
        text=True,
        timeout=timeout,
    )


def _free_bytes(path: Path) -> int:
    candidate = path
    while not candidate.exists() and candidate != candidate.parent:
        candidate = candidate.parent
    return shutil.disk_usage(candidate).free if candidate.exists() else 0


def _reference_manifest_status(path: Path) -> Literal["valid", "missing", "invalid"]:
    if not path.is_file():
        return "missing"
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return "invalid"
    if not isinstance(payload, dict) or not isinstance(payload.get("artifacts"), list):
        return "invalid"
    if payload.get("taxonomy_id") != 9606 or payload.get("proteome_id") != "UP000005640":
        return "invalid"
    for artifact in payload["artifacts"]:
        if not isinstance(artifact, dict):
            return "invalid"
        relative = artifact.get("path")
        expected = artifact.get("sha256")
        if not isinstance(relative, str) or not isinstance(expected, str):
            return "invalid"
        artifact_path = (path.parent / relative).resolve()
        try:
            artifact_path.relative_to(path.parent.resolve())
        except ValueError:
            return "invalid"
        if not artifact_path.is_file():
            return "invalid"
        observed = hashlib.sha256(artifact_path.read_bytes()).hexdigest()
        if observed != expected:
            return "invalid"
    return "valid"


def collect_host_facts(*, deep: bool) -> HostFacts:
    """Collect local facts; deep mode may start one disposable GPU container."""
    data_root = Path("/data/airti-target-fishing")
    cache_root = Path("/mnt/ssd4t/airti-target-fishing")
    docker_available = shutil.which("docker") is not None
    nextflow_available = shutil.which("nextflow") is not None
    commands = {
        "docker": docker_available,
        "docker_compose": False,
        "nextflow": nextflow_available,
    }
    versions: dict[str, str] = {}
    nvidia_runtime = False
    gpus: list[GPUFact] = []
    images: dict[str, ImageFact] = {}

    if nextflow_available:
        process = _run(["nextflow", "-version"])
        versions["nextflow"] = process.stdout + process.stderr
    if docker_available:
        compose = _run(["docker", "compose", "version"])
        commands["docker_compose"] = compose.returncode == 0
        versions["docker_compose"] = compose.stdout.strip()
        info = _run(["docker", "info", "--format", "{{json .Runtimes}}"])
        nvidia_runtime = info.returncode == 0 and '"nvidia"' in info.stdout
        if deep:
            for image in REQUIRED_IMAGES:
                inspect = _run(
                    ["docker", "image", "inspect", image, "--format", "{{.Id}}"]
                )
                if inspect.returncode != 0:
                    continue
                digest = inspect.stdout.strip()
                smoke = _run(
                    ["docker", "run", "--rm", image, "airti-tf", "version"],
                    timeout=120,
                )
                images[image] = ImageFact(
                    command_ok=smoke.returncode == 0,
                    digest=digest if digest.startswith("sha256:") else None,
                )

    nvidia_smi = shutil.which("nvidia-smi")
    if nvidia_smi is not None:
        query = _run(
            [
                nvidia_smi,
                "--query-gpu=memory.total",
                "--format=csv,noheader,nounits",
            ]
        )
        if query.returncode == 0:
            gpus = [
                GPUFact(memory_mib=int(line.strip()), container_visible=False)
                for line in query.stdout.splitlines()
                if line.strip().isdigit()
            ]

    if deep and docker_available and nvidia_runtime and gpus:
        smoke = _run(
            [
                "docker",
                "run",
                "--rm",
                "--gpus",
                "all",
                REQUIRED_IMAGES[0],
                "python",
                "-c",
                (
                    "import torch; assert torch.cuda.is_available(); "
                    "print(torch.cuda.get_device_name(0))"
                ),
            ],
            timeout=120,
        )
        if smoke.returncode == 0:
            gpus = [gpu.model_copy(update={"container_visible": True}) for gpu in gpus]

    return HostFacts(
        data_root=data_root,
        cache_root=cache_root,
        data_free_bytes=_free_bytes(data_root),
        cache_free_bytes=_free_bytes(cache_root),
        commands=commands,
        versions=versions,
        nvidia_runtime=nvidia_runtime,
        gpus=gpus,
        images=images,
        reference_manifest_status=_reference_manifest_status(
            data_root / "reference" / "reference_manifest.json"
        ),
        open_file_limit=resource.getrlimit(resource.RLIMIT_NOFILE)[0],
    )
