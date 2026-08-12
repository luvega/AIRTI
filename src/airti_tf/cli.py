"""Command-line entry point for AIRTI Target Fishing."""

import typer

from airti_tf import __version__

app = typer.Typer(no_args_is_help=True)


@app.callback()
def main() -> None:
    """Run the AIRTI human-proteome target-fishing pipeline."""


@app.command()
def version() -> None:
    """Print the package version."""
    typer.echo(f"airti-tf {__version__}")


if __name__ == "__main__":
    app()
