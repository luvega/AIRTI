from pathlib import Path

from airti_tf.screening.background import (
    BackgroundCandidate,
    build_chembl_candidate,
    select_background_panel,
    validate_panel_tsv,
    write_panel_tsv,
)


def candidate(index: int) -> BackgroundCandidate:
    molecular_weight = 100 + index * 4
    return BackgroundCandidate(
        source_id=f"CHEMBL{index:06}",
        canonical_smiles=f"C{index}",
        molecular_weight=molecular_weight,
        clogp=float(index % 7),
        hbd=index % 4,
        hba=index % 9,
        rotatable_bonds=index % 12,
        formal_charge=(index % 3) - 1,
        fingerprint=(index.to_bytes(4, "little") * 64),
    )


def test_panel_selects_ten_from_each_weight_stratum() -> None:
    candidates = [candidate(index) for index in range(250)]

    panel = select_background_panel(candidates, panel_size=100, strata=10, seed=20260812)

    assert len(panel) == 100
    assert len({item.source_id for item in panel}) == 100
    assert {item.mw_stratum for item in panel} == set(range(10))
    assert all(sum(item.mw_stratum == stratum for item in panel) == 10 for stratum in range(10))


def test_panel_selection_is_deterministic() -> None:
    candidates = [candidate(index) for index in range(250)]

    first = select_background_panel(candidates, panel_size=100, strata=10, seed=20260812)
    second = select_background_panel(list(reversed(candidates)), panel_size=100, strata=10, seed=20260812)

    assert [item.source_id for item in first] == [item.source_id for item in second]


def test_chemistry_filter_rejects_metals_reactivity_and_undefined_stereo() -> None:
    assert build_chembl_candidate("METAL", "[Na+].[Cl-]") is None
    assert build_chembl_candidate("REACTIVE", "O=C(Cl)c1ccccc1") is None
    assert build_chembl_candidate("STEREO", "CC(O)c1ccccc1") is None
    assert build_chembl_candidate("CHEMBL25", "CC(=O)Oc1ccccc1C(=O)O") is not None


def test_panel_tsv_round_trip_validates_exact_counts(tmp_path: Path) -> None:
    panel = select_background_panel(
        [candidate(index) for index in range(250)],
        panel_size=100,
        strata=10,
        seed=20260812,
    )
    output = tmp_path / "panel.smi"

    digest = write_panel_tsv(output, panel)
    summary = validate_panel_tsv(output)

    assert len(digest) == 64
    assert summary == {"probe_count": 100, "strata": 10, "per_stratum": 10}
