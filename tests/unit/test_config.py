from pathlib import Path

import pytest

from airti_tf.config import LockedConfigurationError, load_settings


def test_production_defaults_are_locked() -> None:
    cfg = load_settings("production")

    assert cfg.screen.seeds == [11, 29, 47]
    assert cfg.screen.background_probe_count == 100
    assert cfg.routing.screen_top_n == 300
    assert cfg.routing.boltz_top_n == 30
    assert cfg.routing.md_top_n == 10
    assert cfg.routing.md_replica_top_n == 3
    assert cfg.md.duration_ns == 100


def test_profile_overrides_unlocked_runtime_paths() -> None:
    cfg = load_settings("production")

    assert cfg.runtime.artifact_root == Path("/data/airti-target-fishing")
    assert cfg.runtime.cache_root == Path("/mnt/ssd4t/airti-target-fishing")
    assert cfg.runtime.executor == "local"


def test_locked_override_requires_explicit_permission() -> None:
    with pytest.raises(LockedConfigurationError, match="screen.seeds"):
        load_settings("production", overrides={"screen": {"seeds": [1, 2, 3]}})


def test_locked_override_is_auditable_when_allowed() -> None:
    cfg = load_settings(
        "production",
        overrides={"screen": {"seeds": [1, 2, 3]}},
        allow_locked_override=True,
    )

    assert cfg.screen.seeds == [1, 2, 3]
    assert cfg.override_audit == {"screen.seeds": {"from": [11, 29, 47], "to": [1, 2, 3]}}

