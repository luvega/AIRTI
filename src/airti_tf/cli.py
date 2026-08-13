"""Command-line entry point for AIRTI Target Fishing."""

from enum import StrEnum
from pathlib import Path

import typer

from airti_tf import __version__
from airti_tf.manifest_io import write_artifact
from airti_tf.runtime import collect_host_facts, run_preflight

app = typer.Typer(no_args_is_help=True)
targets_app = typer.Typer(no_args_is_help=True)
app.add_typer(targets_app, name="targets")


class Profile(StrEnum):
    LOCAL = "local"
    PRODUCTION = "production"


class StructureSource(StrEnum):
    PDB = "pdb"
    ALPHAFOLD = "alphafold"


@app.callback()
def main() -> None:
    """Run the AIRTI human-proteome target-fishing pipeline."""


@app.command()
def version() -> None:
    """Print the package version."""
    typer.echo(f"airti-tf {__version__}")


@app.command()
def preflight(
    profile: Profile = typer.Option(Profile.LOCAL),
    output: Path = typer.Option(Path("preflight.json")),
) -> None:
    """Check whether a local or production run may start."""
    host = collect_host_facts(deep=profile == Profile.PRODUCTION)
    result = run_preflight(profile=profile.value, host=host)
    write_artifact(output, result.model_dump(mode="json"))
    typer.echo(result.model_dump_json(indent=2))
    if result.exit_code:
        raise typer.Exit(result.exit_code)


@targets_app.command("fetch-uniprot")
def fetch_uniprot(
    proteome: str = typer.Option("UP000005640"),
    release: str = typer.Option(...),
    output: Path = typer.Option(..., "--out"),
) -> None:
    """Fetch a release-pinned human canonical proteome snapshot."""
    from airti_tf.sources.uniprot import fetch_uniprot_snapshot

    snapshot = fetch_uniprot_snapshot(
        proteome=proteome,
        release=release,
        output_dir=output,
    )
    typer.echo(snapshot.manifest_path)


@targets_app.command("compile-reference")
def compile_reference(
    root: Path = typer.Option(..., "--root", exists=True, file_okay=False),
    proteome: Path = typer.Option(..., "--proteome", exists=True, readable=True),
    targets: Path = typer.Option(..., "--targets", exists=True, readable=True),
    output: Path = typer.Option(..., "--output"),
) -> None:
    """Validate and hash-lock a human whole-proteome target library."""
    from airti_tf.targets.reference import compile_reference_bundle

    summary = compile_reference_bundle(
        root=root,
        proteome_manifest=proteome,
        target_manifest=targets,
        output_manifest=output,
    )
    typer.echo(summary.model_dump_json(indent=2))


@targets_app.command("calibrate-pocket")
def calibrate_pocket(
    panel: Path = typer.Option(..., "--panel", exists=True, readable=True),
    receptor: Path = typer.Option(..., "--receptor", exists=True, readable=True),
    pocket_id: str = typer.Option(..., "--pocket-id"),
    output: Path = typer.Option(..., "--output"),
    asset_dir: Path = typer.Option(..., "--asset-dir"),
    center_x: float = typer.Option(..., "--center-x"),
    center_y: float = typer.Option(..., "--center-y"),
    center_z: float = typer.Option(..., "--center-z"),
    size_x: float = typer.Option(..., "--size-x", min=1.0),
    size_y: float = typer.Option(..., "--size-y", min=1.0),
    size_z: float = typer.Option(..., "--size-z", min=1.0),
    expected_probes: int = typer.Option(100, "--expected-probes", min=1),
    minimum_successful: int = typer.Option(95, "--minimum-successful", min=1),
    workers: int = typer.Option(8, "--workers", min=1),
) -> None:
    """Build a fixed-seed empirical docking background for one pocket."""
    from airti_tf.pockets.receptor import DockingBox
    from airti_tf.screening.calibration_build import calibrate_pocket_background

    summary = calibrate_pocket_background(
        panel,
        receptor_pdbqt=receptor,
        box=DockingBox(
            center=(center_x, center_y, center_z),
            size=(size_x, size_y, size_z),
        ),
        pocket_id=pocket_id,
        output=output,
        asset_dir=asset_dir,
        expected_probe_count=expected_probes,
        minimum_successful_probes=minimum_successful,
        workers=workers,
    )
    typer.echo(summary.model_dump_json(indent=2))


@targets_app.command("build-pilot")
def build_pilot(
    root: Path = typer.Option(..., "--root", exists=True, file_okay=False),
    proteome: Path = typer.Option(..., "--proteome", exists=True, readable=True),
    target_id: str = typer.Option(..., "--target-id"),
    family: str = typer.Option(..., "--family"),
    structure_id: str = typer.Option(..., "--structure-id"),
    structure_source: StructureSource = typer.Option(..., "--structure-source"),
    structure: Path = typer.Option(..., "--structure", exists=True, readable=True),
    ligand_resname: str = typer.Option(..., "--ligand-resname"),
    receptor: Path = typer.Option(..., "--receptor", exists=True, readable=True),
    calibration: Path = typer.Option(..., "--calibration", exists=True, readable=True),
    center_x: float = typer.Option(..., "--center-x"),
    center_y: float = typer.Option(..., "--center-y"),
    center_z: float = typer.Option(..., "--center-z"),
    size_x: float = typer.Option(..., "--size-x", min=1.0),
    size_y: float = typer.Option(..., "--size-y", min=1.0),
    size_z: float = typer.Option(..., "--size-z", min=1.0),
    structure_quality: float = typer.Option(
        ..., "--structure-quality", min=0.0, max=1.0
    ),
    msa_version: str = typer.Option(..., "--msa-version"),
    output: Path = typer.Option(..., "--output"),
) -> None:
    """Create a full-coverage manifest with one real pilot target ready."""
    from airti_tf.pockets.receptor import DockingBox
    from airti_tf.targets.reference import build_pilot_target_manifest

    summary = build_pilot_target_manifest(
        root=root,
        proteome_manifest=proteome,
        target_id=target_id,
        family=family,
        structure_id=structure_id,
        structure_source=structure_source.value,
        structure_pdb=structure,
        ligand_resname=ligand_resname,
        receptor_pdbqt=receptor,
        calibration_json=calibration,
        box=DockingBox(
            center=(center_x, center_y, center_z),
            size=(size_x, size_y, size_z),
        ),
        structure_quality=structure_quality,
        msa_database_version=msa_version,
        output_manifest=output,
    )
    typer.echo(summary.model_dump_json(indent=2))


@app.command("prepare-ligands")
def prepare_ligands(
    input_path: Path = typer.Argument(..., exists=True, readable=True),
    output: Path = typer.Option(..., "--output"),
    asset_dir: Path = typer.Option(Path("prepared_ligands"), "--asset-dir"),
    profile: Profile = typer.Option(Profile.PRODUCTION),
    max_molecules: int = typer.Option(5, min=1, max=5),
) -> None:
    """Prepare one to five query ligands and their docking assets."""
    from airti_tf.stages import prepare_ligand_bundle

    summary = prepare_ligand_bundle(
        input_path,
        output_manifest=output,
        asset_dir=asset_dir,
        profile=profile.value,
        max_molecules=max_molecules,
    )
    typer.echo(summary.model_dump_json(indent=2))
    if summary.failed_query_count:
        raise typer.Exit(3)


@app.command()
def screen(
    ligands: Path = typer.Option(..., "--ligands", exists=True, readable=True),
    targets: Path = typer.Option(..., "--targets", exists=True, readable=True),
    output: Path = typer.Option(..., "--output"),
    asset_dir: Path = typer.Option(Path("screen_assets"), "--asset-dir"),
    top_n: int = typer.Option(300, "--top-n", min=1, max=500),
) -> None:
    """Run calibrated whole-library docking."""
    from airti_tf.stages import screen_ligand_bundle

    summary = screen_ligand_bundle(
        ligands,
        targets,
        output_manifest=output,
        asset_dir=asset_dir,
        top_n=top_n,
    )
    typer.echo(summary.model_dump_json(indent=2))
    if summary.candidate_count == 0:
        raise typer.Exit(2)


@app.command("refine-boltz")
def refine_boltz(
    candidates: Path = typer.Option(..., "--candidates", exists=True, readable=True),
    output: Path = typer.Option(..., "--output"),
    asset_dir: Path = typer.Option(Path("boltz_assets"), "--asset-dir"),
    profile: Profile = typer.Option(Profile.PRODUCTION),
    top_n: int = typer.Option(30, "--top-n", min=1, max=300),
    cache: Path = typer.Option(Path("/models/boltz"), "--cache"),
) -> None:
    """Run multi-seed Boltz-2 refinement."""
    from airti_tf.stages import refine_boltz_bundle

    summary = refine_boltz_bundle(
        candidates,
        output_manifest=output,
        asset_dir=asset_dir,
        profile=profile.value,
        top_n=top_n,
        cache_path=cache,
    )
    typer.echo(summary.model_dump_json(indent=2))
    if summary.succeeded_candidate_count == 0:
        raise typer.Exit(2)


@app.command("run-md")
def run_md(
    candidates: Path = typer.Option(..., "--candidates", exists=True, readable=True),
    output: Path = typer.Option(..., "--output"),
    asset_dir: Path = typer.Option(Path("md_assets"), "--asset-dir"),
    top_n: int = typer.Option(10, "--top-n", min=1, max=10),
) -> None:
    """Run checkpoint-aware molecular dynamics."""
    from airti_tf.stages import run_md_bundle

    summary = run_md_bundle(
        candidates,
        output_manifest=output,
        asset_dir=asset_dir,
        top_n=top_n,
    )
    typer.echo(summary.model_dump_json(indent=2))
    if summary.succeeded_candidate_count == 0:
        raise typer.Exit(2)


@app.command("render-report")
def render_report_command(
    candidates: Path = typer.Option(..., "--candidates", exists=True, readable=True),
    output: Path = typer.Option(..., "--output"),
    state_db: Path = typer.Option(Path("job_status.sqlite"), "--state-db"),
    project_id: str | None = typer.Option(None, "--project-id"),
) -> None:
    """Render the traceable computation-only report."""
    from airti_tf.stages import render_report_bundle

    summary = render_report_bundle(
        candidates,
        output_dir=output,
        state_db=state_db,
        project_id=project_id,
    )
    typer.echo(summary.model_dump_json(indent=2))


if __name__ == "__main__":
    app()
