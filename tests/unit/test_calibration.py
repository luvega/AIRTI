import pytest

from airti_tf.screening.calibration import (
    BackgroundDistribution,
    InsufficientBackgroundError,
    ScreenHit,
    empirical_percentile,
    route_screen_candidates,
)


def hit(
    target: str,
    score: float,
    *,
    pocket: str = "P1",
    state: str = "S1",
    family: str = "kinase",
) -> ScreenHit:
    return ScreenHit(
        target_id=target,
        family=family,
        pocket_id=pocket,
        ligand_state_id=state,
        affinity_median=-5 - score,
        calibrated_score=score,
        seed_range=0.5,
        pose_consistency=0.8,
    )


def test_empirical_percentile_is_directionally_correct() -> None:
    background = [-10.0, -8.0, -7.0, -6.0, -5.0]

    assert empirical_percentile(query=-9.0, background=background) == pytest.approx(5 / 6)


def test_background_requires_95_valid_probe_values() -> None:
    with pytest.raises(InsufficientBackgroundError, match="94"):
        BackgroundDistribution(pocket_id="P1", affinities=[-6.0] * 94)


def test_top300_is_union_across_ligand_states_and_pockets() -> None:
    records = [
        hit("P00533", 0.95, pocket="P1", state="S1"),
        hit("P00533", 0.92, pocket="P2", state="S2"),
        hit("P04637", 0.90, pocket="P3", state="S1", family="transcription_factor"),
    ]

    routed = route_screen_candidates(records, top_n=300)

    assert len(routed) == 2
    assert len({record.target_id for record in routed}) == len(routed)
    assert routed[0].target_id == "P00533"
    assert routed[0].best_pocket_id == "P1"
    assert routed[0].best_state_id == "S1"
    assert routed[0].second_best_score == pytest.approx(0.92)


def test_family_cap_prevents_one_family_from_consuming_route() -> None:
    records = [hit(f"K{index:03}", 1 - index / 1000) for index in range(30)]
    records += [
        hit(f"G{index:03}", 0.8 - index / 1000, family=f"family-{index}")
        for index in range(30)
    ]

    routed = route_screen_candidates(
        records, top_n=40, primary_n=30, family_cap=15
    )

    assert len(routed) == 40
    assert sum(record.family == "kinase" for record in routed) <= 15
    assert len({record.family for record in routed}) >= 26


def test_ties_are_deterministic() -> None:
    routed = route_screen_candidates(
        [hit("P2", 0.8), hit("P1", 0.8)], top_n=2, primary_n=2
    )

    assert [record.target_id for record in routed] == ["P1", "P2"]

