"""Pocket-derived Vina boxes and Meeko receptor commands."""

from __future__ import annotations

from pathlib import Path

from pydantic import BaseModel, ConfigDict


class BoxTooLargeError(ValueError):
    """Raised when a pocket cannot be represented by one allowed box."""


class AtomCoordinate(BaseModel):
    model_config = ConfigDict(extra="forbid")

    coord: tuple[float, float, float]


class DockingBox(BaseModel):
    model_config = ConfigDict(extra="forbid")

    center: tuple[float, float, float]
    size: tuple[float, float, float]

    def contains(self, coord: tuple[float, float, float]) -> bool:
        return all(
            abs(value - center) <= length / 2 + 1e-9
            for value, center, length in zip(coord, self.center, self.size, strict=True)
        )


def build_box(
    pocket_atoms: list[AtomCoordinate],
    *,
    padding_a: float,
    min_size_a: float,
    max_size_a: float,
) -> DockingBox:
    """Build an enclosing box; never silently clip oversized pockets."""
    if not pocket_atoms:
        raise ValueError("at least one pocket atom is required")
    axes = list(zip(*(atom.coord for atom in pocket_atoms), strict=True))
    minima = tuple(min(axis) for axis in axes)
    maxima = tuple(max(axis) for axis in axes)
    center = tuple((low + high) / 2 for low, high in zip(minima, maxima, strict=True))
    raw_size = tuple(
        high - low + 2 * padding_a for low, high in zip(minima, maxima, strict=True)
    )
    if any(length > max_size_a for length in raw_size):
        largest = max(raw_size)
        raise BoxTooLargeError(
            f"required box dimension {largest:.1f} A exceeds maximum {max_size_a:.1f} A"
        )
    size = tuple(max(length, min_size_a) for length in raw_size)
    return DockingBox(center=center, size=size)


def build_meeko_command(
    *, input_pdb: Path, output_prefix: Path, box: DockingBox
) -> list[str]:
    """Build a Meeko receptor-preparation command with explicit output paths."""
    center = [f"{value:.3f}" for value in box.center]
    size = [f"{value:.3f}" for value in box.size]
    return [
        "mk_prepare_receptor.py",
        "-i",
        str(input_pdb),
        "--write_pdbqt",
        str(output_prefix.with_suffix(".pdbqt")),
        "--write_vina_box",
        str(output_prefix.with_suffix(".box.txt")),
        "--box_center",
        *center,
        "--box_size",
        *size,
    ]

