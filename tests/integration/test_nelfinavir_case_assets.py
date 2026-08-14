from pathlib import Path

import yaml

from airti_tf.cases.evaluation import load_case_definition


CASE_ROOT = Path("cases/nelfinavir")


def test_nelfinavir_parent_case_assets_are_frozen_and_blinded() -> None:
    definition = load_case_definition(CASE_ROOT / "case.yaml")
    query = (CASE_ROOT / "input" / "nelfinavir.smi").read_text().strip()
    query_line = next(
        line for line in query.splitlines() if not line.startswith("#")
    )
    query_fields = query_line.split()
    panel = yaml.safe_load((CASE_ROOT / "panel.yaml").read_text())

    assert definition.query_inchikey == "QAGYKUNXZHXKMR-HKWSIXNMSA-N"
    assert query_fields[1] == "nelfinavir-parent"
    assert query_fields[2].startswith("states=")
    assert len(query_fields[2][len("states=") :].split("|")) == 2
    assert "M8" not in query
    assert definition.upstream_visible_target_ids == []
    assert {anchor.target_id for anchor in definition.positive_anchors} == {
        "P08684",
        "P20815",
        "P08183",
        "Q5TDH0",
        "P07900",
        "O43462",
    }
    assert panel["panel_size"] == 64
    assert panel["family_comparators_per_anchor"] == 3
    assert panel["stratified_comparator_count"] == 40
    assert panel["selection_seed"] == 20260814


def test_nelfinavir_production_config_preserves_frozen_scientific_defaults() -> None:
    config = yaml.safe_load((CASE_ROOT / "production.yaml").read_text())

    assert config["ligand"] == {
        "state_policy": "curated_input",
        "formal_charges": [0, 1],
        "experimental_pka": [6.00, 11.06],
        "pka_reference_doi": "10.1016/S0378-4347(97)00193-X",
    }
    assert config["routing"] == {
        "screen_top_n": 300,
        "boltz_top_n": 30,
        "md_top_n": 10,
        "md_replica_top_n": 3,
    }
    assert config["md"]["duration_ns"] == 100
    assert config["md"]["membrane"]["lipids"] == ["POPC", "CHL1"]
    assert config["md"]["membrane"]["ratio"] == [4, 1]
    assert config["md"]["membrane"]["salt_molar"] == 0.15
