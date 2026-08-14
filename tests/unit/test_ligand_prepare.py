import pytest

from airti_tf.ligands.prepare import prepare_ligand, validate_query_batch


def test_undefined_stereochemistry_blocks_production() -> None:
    result = prepare_ligand("CC(O)C(=O)O", profile="production")

    assert result.status == "failed"
    assert result.error_code == "undefined_stereochemistry"
    assert result.states == []


def test_query_count_is_limited_to_five() -> None:
    with pytest.raises(ValueError, match="1 to 5"):
        validate_query_batch([f"C{'C' * index}" for index in range(6)])


def test_query_batch_rejects_duplicates() -> None:
    with pytest.raises(ValueError, match="unique"):
        validate_query_batch(["CCO", "CCO"])


def test_explicit_stereo_generates_deterministic_ph74_states() -> None:
    first = prepare_ligand("CC[C@H](O)C(=O)O", profile="production")
    second = prepare_ligand("CC[C@H](O)C(=O)O", profile="production")

    assert first.status == "succeeded"
    assert first.states
    assert [state.ligand_state_id for state in first.states] == [
        state.ligand_state_id for state in second.states
    ]
    assert all(state.ph_min == 6.4 and state.ph_max == 8.4 for state in first.states)
    assert all(state.mol_block for state in first.states)


def test_salt_is_removed_but_mapping_to_original_is_retained() -> None:
    result = prepare_ligand("CCN.Cl", profile="local")

    assert result.status == "succeeded"
    assert result.original_smiles == "CCN.Cl"
    assert result.fragment_count == 2
    assert all("Cl" not in state.canonical_smiles for state in result.states)


def test_boltz_atom_limits_are_explicit() -> None:
    result = prepare_ligand("C" * 60, profile="production")

    assert result.status == "succeeded"
    assert "boltz_high_atom_count" in result.uncertainty_flags
    assert all(state.atom_count <= 128 for state in result.states)


def test_curated_nelfinavir_states_override_generic_pka_enumeration() -> None:
    neutral = (
        "CC1=C(C=CC=C1O)C(=O)N[C@@H](CSC2=CC=CC=C2)"
        "[C@@H](CN3C[C@H]4CCCC[C@H]4C[C@H]3C(=O)NC(C)(C)C)O"
    )
    cation = neutral.replace("CN3C", "C[NH+]3C")

    result = prepare_ligand(
        neutral,
        profile="production",
        explicit_state_smiles=[neutral, cation],
    )

    assert result.status == "succeeded"
    assert [state.formal_charge for state in result.states] == [0, 1]
    assert all("[O-]" not in state.canonical_smiles for state in result.states)
    assert "curated_protonation_states" in result.uncertainty_flags


def test_curated_state_must_preserve_parent_connectivity() -> None:
    result = prepare_ligand(
        "CCO",
        profile="local",
        explicit_state_smiles=["CCN"],
    )

    assert result.status == "failed"
    assert result.error_code == "explicit_state_connectivity_mismatch"
