"""Command-line entry point for AIRTI Target Fishing."""

from enum import StrEnum
from pathlib import Path

import typer

from airti_tf import __version__
from airti_tf.manifest_io import write_artifact
from airti_tf.runtime import collect_host_facts, run_preflight

app = typer.Typer(no_args_is_help=True)


class Profile(StrEnum):
    LOCAL = "local"
    PRODUCTION = "production"


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


if __name__ == "__main__":
    app()
