from airti_tf.screening.calibration import ScreenHit, route_screen_candidates


def test_each_target_keeps_selection_reason() -> None:
    records = [
        ScreenHit(
            target_id=f"P{index:05}",
            family=f"family-{index % 20}",
            pocket_id=f"pocket-{index}",
            ligand_state_id="state-1",
            affinity_median=-8.0,
            calibrated_score=1 - index / 1000,
            seed_range=0.2,
            pose_consistency=0.9,
        )
        for index in range(350)
    ]

    routed = route_screen_candidates(records, top_n=300)

    assert len(routed) == 300
    assert all(record.selection_reason for record in routed)
