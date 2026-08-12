"""Typed, layered configuration with audited production overrides."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Any, Literal

import yaml
from pydantic import BaseModel, ConfigDict, Field

CONFIG_ROOT = Path(__file__).resolve().parents[2] / "configs"

LOCKED_PATHS = frozenset(
    {
        "screen.backend",
        "screen.seeds",
        "screen.exhaustiveness",
        "screen.background_probe_count",
        "routing.screen_top_n",
        "routing.boltz_top_n",
        "routing.md_top_n",
        "routing.md_replica_top_n",
        "boltz.diffusion_samples",
        "boltz.affinity_samples",
        "md.duration_ns",
        "md.timestep_fs",
        "md.forcefield",
        "md.ligand_forcefield",
        "md.water_model",
        "report.final_top_n",
    }
)


class LockedConfigurationError(ValueError):
    """Raised when production scientific defaults are changed silently."""


class StrictModel(BaseModel):
    """Reject unknown configuration keys."""

    model_config = ConfigDict(extra="forbid")


class ScreenSettings(StrictModel):
    backend: Literal["quickvina2"]
    seeds: list[int] = Field(min_length=3, max_length=3)
    exhaustiveness: int = Field(gt=0)
    background_probe_count: int = Field(gt=0)


class RoutingSettings(StrictModel):
    screen_top_n: int = Field(gt=0)
    boltz_top_n: int = Field(gt=0)
    md_top_n: int = Field(gt=0)
    md_replica_top_n: int = Field(gt=0)


class BoltzSettings(StrictModel):
    diffusion_samples: int = Field(gt=0)
    affinity_samples: int = Field(gt=0)


class MDSettings(StrictModel):
    duration_ns: int = Field(gt=0)
    timestep_fs: int = Field(gt=0)
    forcefield: str
    ligand_forcefield: str
    water_model: str


class ReportSettings(StrictModel):
    final_top_n: int = Field(gt=0)


class RuntimeSettings(StrictModel):
    executor: Literal["local", "slurm"]
    artifact_root: Path
    cache_root: Path


class Settings(StrictModel):
    schema_version: str
    screen: ScreenSettings
    routing: RoutingSettings
    boltz: BoltzSettings
    md: MDSettings
    report: ReportSettings
    runtime: RuntimeSettings
    override_audit: dict[str, dict[str, Any]] = Field(default_factory=dict)


def _read_yaml(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        payload = yaml.safe_load(handle) or {}
    if not isinstance(payload, dict):
        raise ValueError(f"configuration root must be a mapping: {path}")
    return payload


def _deep_merge(base: dict[str, Any], update: dict[str, Any]) -> dict[str, Any]:
    merged = deepcopy(base)
    for key, value in update.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _deep_merge(merged[key], value)
        else:
            merged[key] = deepcopy(value)
    return merged


def _flatten(payload: dict[str, Any], prefix: str = "") -> dict[str, Any]:
    flattened: dict[str, Any] = {}
    for key, value in payload.items():
        path = f"{prefix}.{key}" if prefix else key
        if isinstance(value, dict):
            flattened.update(_flatten(value, path))
        else:
            flattened[path] = value
    return flattened


def load_settings(
    profile: Literal["local", "production"],
    *,
    overrides: dict[str, Any] | None = None,
    allow_locked_override: bool = False,
    config_root: Path = CONFIG_ROOT,
) -> Settings:
    """Load defaults, a named profile, and explicit runtime overrides."""

    defaults = _read_yaml(config_root / "defaults.yaml")
    profile_values = _read_yaml(config_root / f"{profile}.yaml")
    values = _deep_merge(defaults, profile_values)
    audit: dict[str, dict[str, Any]] = {}

    if overrides:
        before = _flatten(values)
        requested = _flatten(overrides)
        for path, new_value in requested.items():
            if path in LOCKED_PATHS and before.get(path) != new_value:
                if not allow_locked_override:
                    raise LockedConfigurationError(
                        f"locked setting {path} requires allow_locked_override"
                    )
                audit[path] = {"from": before.get(path), "to": new_value}
        values = _deep_merge(values, overrides)

    values["override_audit"] = audit
    return Settings.model_validate(values)

