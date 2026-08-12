from pathlib import Path

from airti_tf.pockets.fpocket import (
    PocketCandidate,
    build_fpocket_command,
    parse_fpocket,
    qc_pocket,
    select_qualified_pockets,
)


def pocket(**overrides: object) -> PocketCandidate:
    payload: dict[str, object] = {
        "pocket_id": "P00533:abc:1",
        "target_id": "P00533",
        "rank": 1,
        "volume_a3": 180.0,
        "druggability": 0.5,
        "fpocket_score": 20.0,
        "residue_count": 8,
        "mean_plddt": 90.0,
        "exposed_fraction": 0.3,
        "severe_backbone_clash": False,
        "known_ligand_overlap": False,
    }
    payload.update(overrides)
    return PocketCandidate.model_validate(payload)


def test_pocket_ids_are_stable_across_parser_runs() -> None:
    fixture = Path("tests/fixtures/fpocket_output/sample_out")

    first = parse_fpocket(fixture, target_id="P00533")
    second = parse_fpocket(fixture, target_id="P00533")

    assert [item.pocket_id for item in first] == [item.pocket_id for item in second]
    assert [item.rank for item in first] == [1, 2]
    assert first[0].residue_count == 8


def test_rejects_buried_or_low_confidence_pocket() -> None:
    result = qc_pocket(
        pocket(volume_a3=55, druggability=0.02, mean_plddt=51, exposed_fraction=0.01)
    )

    assert result.status == "unsupported"
    assert result.unsupported_reason in {"volume_too_small", "low_confidence", "inaccessible"}
    assert result.score is None


def test_alphafold_pocket_requires_local_plddt_70() -> None:
    result = qc_pocket(pocket(mean_plddt=69), structure_source="alphafold")

    assert result.status == "unsupported"
    assert result.unsupported_reason == "low_confidence"


def test_selects_at_most_five_and_prioritizes_known_ligand() -> None:
    candidates = [pocket(pocket_id=f"p{i}", rank=i, druggability=0.9 - i / 20) for i in range(1, 8)]
    candidates[-1] = candidates[-1].model_copy(
        update={"known_ligand_overlap": True, "druggability": 0.1}
    )

    selected = select_qualified_pockets(candidates, limit=5)

    assert len(selected) == 5
    assert selected[0].pocket.pocket_id == "p7"
    assert all(result.status == "ready" for result in selected)


def test_fpocket_command_uses_documented_cli() -> None:
    assert build_fpocket_command(Path("target.pdb")) == ["fpocket", "-f", "target.pdb"]
