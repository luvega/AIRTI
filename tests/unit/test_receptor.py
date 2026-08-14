from pathlib import Path

import pytest

from airti_tf.pockets.receptor import (
    AtomCoordinate,
    BoxTooLargeError,
    build_box,
    build_meeko_command,
)


def test_receptor_box_contains_all_pocket_atoms() -> None:
    atoms = [
        AtomCoordinate(coord=(-2.0, 1.0, 3.0)),
        AtomCoordinate(coord=(4.0, 6.0, 8.0)),
    ]

    box = build_box(atoms, padding_a=5.0, min_size_a=18.0, max_size_a=32.0)

    assert all(box.contains(atom.coord) for atom in atoms)
    assert box.size == (18.0, 18.0, 18.0)


def test_oversized_box_is_rejected_not_silently_clipped() -> None:
    atoms = [AtomCoordinate(coord=(0, 0, 0)), AtomCoordinate(coord=(30, 0, 0))]

    with pytest.raises(BoxTooLargeError, match="40.0"):
        build_box(atoms, padding_a=5.0, min_size_a=18.0, max_size_a=32.0)


def test_meeko_command_writes_pdbqt_and_exact_box() -> None:
    box = build_box(
        [AtomCoordinate(coord=(0, 0, 0)), AtomCoordinate(coord=(2, 2, 2))],
        padding_a=5.0,
        min_size_a=18.0,
        max_size_a=32.0,
    )

    command = build_meeko_command(
        input_pdb=Path("target.pdb"), output_prefix=Path("prepared/target"), box=box
    )

    assert command[:3] == ["mk_prepare_receptor.py", "-i", "target.pdb"]
    assert "--write_pdbqt" in command
    assert command[command.index("--box_center") + 1 : command.index("--box_center") + 4] == [
        "1.000",
        "1.000",
        "1.000",
    ]


def test_meeko_command_accepts_named_sdf_and_complete_json_templates() -> None:
    box = build_box(
        [AtomCoordinate(coord=(0, 0, 0))],
        padding_a=5.0,
        min_size_a=18.0,
        max_size_a=32.0,
    )

    command = build_meeko_command(
        input_pdb=Path("target.pdb"),
        output_prefix=Path("prepared/target"),
        box=box,
        residue_templates={
            "HEM": Path("HEM-template.json"),
            "FAD": Path("FAD.sdf"),
        },
    )

    additions = [
        command[index + 1]
        for index, value in enumerate(command)
        if value == "--add_templates"
    ]
    assert additions == ["FAD:FAD.sdf", "HEM-template.json"]
