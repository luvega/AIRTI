from airti_tf.ranking.consensus import TargetEvidence, rank_targets


def evidence(index: int, *, complete: bool = True) -> TargetEvidence:
    score = 1.0 - index / 300.0
    return TargetEvidence(
        target_id=f"P{index:05}",
        ligand_id="query-1",
        status="ready",
        vina_score=score,
        docking_consistency=score,
        structure_quality=score,
        boltz_score=score if index <= 30 else None,
        md_score=score if index <= 10 and complete else None,
        md_status="stable" if index <= 10 and complete else None,
        successful_seeds=3,
        boltz_seed_spread=0.01,
        heavy_atom_count=30,
    )


def test_stages_route_300_to_30_to_10_and_final_top5() -> None:
    screen = rank_targets([evidence(i) for i in range(1, 301)], stage="screen")
    boltz_input = [item.source for item in screen.ranked[:30]]
    boltz = rank_targets(boltz_input, stage="boltz")
    final_input = [item.source for item in boltz.ranked[:10]]
    final = rank_targets(final_input, stage="final")

    assert len(screen.ranked) == 300
    assert len(boltz.ranked) == 30
    assert len(final.ranked) == 10
    assert [item.target_id for item in final.full_evidence[:5]] == [
        "P00001",
        "P00002",
        "P00003",
        "P00004",
        "P00005",
    ]


def test_partial_evidence_is_reported_separately_from_full_final_ranking() -> None:
    records = [evidence(1), evidence(2, complete=False)]
    result = rank_targets(records, stage="final")

    assert [item.target_id for item in result.full_evidence] == ["P00001"]
    assert [item.target_id for item in result.partial_evidence] == ["P00002"]
    assert result.partial_evidence[0].normalized_scores["md"] is None
