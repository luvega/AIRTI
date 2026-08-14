import json
from pathlib import Path

import yaml

from airti_tf.cases.panel import build_case_panel
from airti_tf.manifest_io import write_jsonl_atomic


def test_case_panel_contains_anchors_family_controls_and_five_strata(
    tmp_path: Path,
) -> None:
    anchors = [f"A{index}" for index in range(6)]
    case = tmp_path / "case.yaml"
    case.write_text(
        yaml.safe_dump(
            {
                "schema_version": "1.0",
                "case_id": "synthetic",
                "query_id": "query",
                "query_inchikey": "A" * 14 + "-" + "B" * 10 + "-C",
                "positive_anchors": [
                    {
                        "target_id": target_id,
                        "gene_symbol": target_id,
                        "tier": "gold" if index < 2 else "silver",
                    }
                    for index, target_id in enumerate(anchors)
                ],
                "novelty_exclusions": anchors,
                "upstream_visible_target_ids": [],
            }
        ),
        encoding="utf-8",
    )
    spec = tmp_path / "panel.yaml"
    spec.write_text(
        yaml.safe_dump(
            {
                "schema_version": "1.0",
                "panel_id": "synthetic-panel",
                "panel_size": 64,
                "anchor_count": 6,
                "family_comparators_per_anchor": 3,
                "stratified_comparator_count": 40,
                "selection_seed": 20260814,
                "strata": {
                    "enzyme": 8,
                    "receptor": 8,
                    "transporter_or_channel": 8,
                    "chaperone_or_scaffold": 8,
                    "other": 8,
                },
                "comparison_proteins_are_true_negatives": False,
            }
        ),
        encoding="utf-8",
    )
    rows: list[dict[str, object]] = []
    for index, target_id in enumerate(anchors):
        family = f"anchor-family-{index}"
        rows.append({"target_id": target_id, "status": "ready", "family": family})
        rows.extend(
            {
                "target_id": f"F{index}{offset}",
                "status": "ready",
                "family": family,
            }
            for offset in range(3)
        )
    stratum_families = {
        "enzyme": "protein kinase family",
        "receptor": "G-protein coupled receptor family",
        "transporter_or_channel": "ABC transporter family",
        "chaperone_or_scaffold": "heat shock chaperone family",
        "other": "unclassified family",
    }
    for prefix, family in stratum_families.items():
        rows.extend(
            {
                "target_id": f"{prefix}-{index}",
                "status": "ready",
                "family": family,
            }
            for index in range(8)
        )
    reference = tmp_path / "targets.jsonl"
    output = tmp_path / "panel-targets.jsonl"
    audit = tmp_path / "panel-audit.json"
    write_jsonl_atomic(reference, rows)

    summary = build_case_panel(
        case_path=case,
        panel_spec_path=spec,
        target_manifest=reference,
        output_manifest=output,
        audit_output=audit,
    )

    selected = [json.loads(line) for line in output.read_text().splitlines()]
    audit_payload = json.loads(audit.read_text())
    assert summary.target_count == 64
    assert {row["target_id"] for row in selected} >= set(anchors)
    assert len(audit_payload["anchors"]) == 6
    assert sum(
        len(value) for value in audit_payload["family_comparators"].values()
    ) == 18
    assert {
        key: len(value) for key, value in audit_payload["strata"].items()
    } == {
        "enzyme": 8,
        "receptor": 8,
        "transporter_or_channel": 8,
        "chaperone_or_scaffold": 8,
        "other": 8,
    }


def test_case_panel_honors_frozen_family_comparator_ids(tmp_path: Path) -> None:
    anchors = [f"A{index}" for index in range(2)]
    case = tmp_path / "case.yaml"
    case.write_text(
        yaml.safe_dump(
            {
                "schema_version": "1.0",
                "case_id": "frozen-comparators",
                "query_id": "query",
                "query_inchikey": "A" * 14 + "-" + "B" * 10 + "-C",
                "positive_anchors": [
                    {
                        "target_id": target_id,
                        "gene_symbol": target_id,
                        "tier": "gold",
                    }
                    for target_id in anchors
                ],
                "novelty_exclusions": anchors,
                "upstream_visible_target_ids": [],
            }
        ),
        encoding="utf-8",
    )
    spec = tmp_path / "panel.yaml"
    spec.write_text(
        yaml.safe_dump(
            {
                "schema_version": "1.0",
                "panel_id": "frozen-panel",
                "panel_size": 8,
                "anchor_count": 2,
                "family_comparators_per_anchor": 2,
                "stratified_comparator_count": 2,
                "selection_seed": 7,
                "strata": {"other": 2},
                "family_comparators": {
                    "A0": ["X0", "X1"],
                    "A1": ["X2", "X3"],
                },
                "comparison_proteins_are_true_negatives": False,
            }
        ),
        encoding="utf-8",
    )
    rows = [
        {"uniprot_id": target_id, "reviewed": True, "protein_family": family}
        for target_id, family in [
            ("A0", "DDI1 family"),
            ("A1", "Peptidase M50A family"),
            ("X0", "DDI1 family"),
            ("X1", "Peptidase A2 family"),
            ("X2", "Peptidase S8 family"),
            ("X3", "Peptidase M48A family"),
            ("O0", "unclassified family"),
            ("O1", "unclassified family"),
        ]
    ]
    reference = tmp_path / "targets.jsonl"
    output = tmp_path / "panel-targets.jsonl"
    audit = tmp_path / "panel-audit.json"
    write_jsonl_atomic(reference, rows)

    summary = build_case_panel(
        case_path=case,
        panel_spec_path=spec,
        target_manifest=reference,
        output_manifest=output,
        audit_output=audit,
    )

    audit_payload = json.loads(audit.read_text())
    assert summary.target_count == 8
    assert audit_payload["family_comparators"] == {
        "A0": ["X0", "X1"],
        "A1": ["X2", "X3"],
    }
    assert audit_payload["family_comparator_selection"] == "frozen_ids"
