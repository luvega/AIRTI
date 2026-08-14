import json
from pathlib import Path

import yaml
from typer.testing import CliRunner

from airti_tf.cases.evaluation import evaluate_case_manifests, load_case_definition
from airti_tf.cli import app


def _write_jsonl(path: Path, rows: list[dict[str, object]]) -> None:
    path.write_text(
        "".join(json.dumps(row, sort_keys=True) + "\n" for row in rows),
        encoding="utf-8",
    )


def test_load_case_definition_rejects_query_target_leakage(tmp_path: Path) -> None:
    case_path = tmp_path / "case.yaml"
    case_path.write_text(
        yaml.safe_dump(
            {
                "schema_version": "1.0",
                "case_id": "nelfinavir-parent-v1",
                "query_id": "nelfinavir-parent",
                "query_inchikey": "QAGYKUNXZHXKMR-HKWSIXNMSA-N",
                "positive_anchors": [
                    {"target_id": "P08684", "gene_symbol": "CYP3A4", "tier": "gold"}
                ],
                "novelty_exclusions": ["P08684"],
                "upstream_visible_target_ids": ["P08684"],
            }
        ),
        encoding="utf-8",
    )

    try:
        load_case_definition(case_path)
    except ValueError as error:
        assert "anchor leakage" in str(error)
    else:
        raise AssertionError("positive anchors must not be upstream-visible")


def test_evaluate_case_reports_tiered_recall_and_stable_novel_hits(
    tmp_path: Path,
) -> None:
    case_path = tmp_path / "case.yaml"
    case_path.write_text(
        yaml.safe_dump(
            {
                "schema_version": "1.0",
                "case_id": "nelfinavir-parent-v1",
                "query_id": "nelfinavir-parent",
                "query_inchikey": "QAGYKUNXZHXKMR-HKWSIXNMSA-N",
                "positive_anchors": [
                    {"target_id": "P08684", "gene_symbol": "CYP3A4", "tier": "gold"},
                    {"target_id": "P20815", "gene_symbol": "CYP3A5", "tier": "gold"},
                    {"target_id": "P08183", "gene_symbol": "ABCB1", "tier": "silver"},
                ],
                "novelty_exclusions": ["P08684", "P20815", "P08183", "P08238"],
                "upstream_visible_target_ids": [],
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    screen = tmp_path / "screen.jsonl"
    boltz = tmp_path / "boltz.jsonl"
    md = tmp_path / "md.jsonl"
    _write_jsonl(
        screen,
        [
            {"query_id": "nelfinavir-parent", "target_id": "P08684", "screen_rank": 4},
            {"query_id": "nelfinavir-parent", "target_id": "P08183", "screen_rank": 20},
        ],
    )
    _write_jsonl(
        boltz,
        [
            {
                "query_id": "nelfinavir-parent",
                "target_id": "P08684",
                "boltz_status": "succeeded",
                "calibrated_score": 0.8,
                "pose_consistency": 0.8,
                "structure_quality": 0.9,
                "boltz_score": 0.9,
                "boltz_seed_success_count": 3,
                "ligand_atom_count": 40,
            },
            {
                "query_id": "nelfinavir-parent",
                "target_id": "Q99999",
                "boltz_status": "succeeded",
                "calibrated_score": 0.7,
                "pose_consistency": 0.7,
                "structure_quality": 0.8,
                "boltz_score": 0.7,
                "boltz_seed_success_count": 3,
                "ligand_atom_count": 40,
            },
        ],
    )
    _write_jsonl(
        md,
        [
            {
                "query_id": "nelfinavir-parent",
                "target_id": "Q99999",
                "md_status": "stable",
                "calibrated_score": 0.9,
                "pose_consistency": 0.8,
                "structure_quality": 0.9,
                "boltz_score": 0.9,
                "md_score": 0.95,
                "boltz_seed_success_count": 3,
                "ligand_atom_count": 40,
            },
            {
                "query_id": "nelfinavir-parent",
                "target_id": "P08238",
                "md_status": "stable",
                "calibrated_score": 0.8,
                "pose_consistency": 0.8,
                "structure_quality": 0.9,
                "boltz_score": 0.8,
                "md_score": 0.8,
                "boltz_seed_success_count": 3,
                "ligand_atom_count": 40,
            },
            {
                "query_id": "nelfinavir-parent",
                "target_id": "Q88888",
                "md_status": "unstable",
                "calibrated_score": 0.7,
                "pose_consistency": 0.8,
                "structure_quality": 0.9,
                "boltz_score": 0.7,
                "md_score": 0.7,
                "boltz_seed_success_count": 3,
                "ligand_atom_count": 40,
            },
        ],
    )

    result = evaluate_case_manifests(
        case_path=case_path,
        screen_manifest=screen,
        boltz_manifest=boltz,
        md_manifest=md,
    )

    assert result.case_id == "nelfinavir-parent-v1"
    assert result.stages["screen"].tiers["gold"].recall_at[10] == 0.5
    assert result.stages["screen"].tiers["silver"].recall_at[50] == 1.0
    assert result.stages["boltz"].anchor_ranks["P08684"] == 1
    assert result.exploratory_targets == ["Q99999"]


def test_evaluate_case_cli_writes_auditable_json(
    monkeypatch: object, tmp_path: Path
) -> None:
    case = tmp_path / "case.yaml"
    screen = tmp_path / "screen.jsonl"
    boltz = tmp_path / "boltz.jsonl"
    md = tmp_path / "md.jsonl"
    for path in (case, screen, boltz, md):
        path.write_text("{}\n", encoding="utf-8")
    output = tmp_path / "evaluation.json"

    class FakeResult:
        def model_dump(self, *, mode: str) -> dict[str, object]:
            assert mode == "json"
            return {"schema_version": "1.0", "case_id": "nelfinavir-parent-v1"}

        def model_dump_json(self, *, indent: int) -> str:
            assert indent == 2
            return '{"case_id":"nelfinavir-parent-v1"}'

    monkeypatch.setattr(  # type: ignore[attr-defined]
        "airti_tf.cases.evaluation.evaluate_case_manifests",
        lambda **_kwargs: FakeResult(),
    )
    result = CliRunner().invoke(
        app,
        [
            "evaluate-case",
            "--case",
            str(case),
            "--screen",
            str(screen),
            "--boltz",
            str(boltz),
            "--md",
            str(md),
            "--output",
            str(output),
        ],
    )

    assert result.exit_code == 0, result.stdout
    assert json.loads(output.read_text()) == {
        "case_id": "nelfinavir-parent-v1",
        "schema_version": "1.0",
    }
