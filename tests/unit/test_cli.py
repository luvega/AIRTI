from typer.testing import CliRunner

from airti_tf.cli import app


def test_version_command_reports_package_version() -> None:
    result = CliRunner().invoke(app, ["version"])

    assert result.exit_code == 0
    assert result.stdout.strip() == "airti-tf 0.1.0"


def test_target_library_build_commands_are_exposed() -> None:
    result = CliRunner().invoke(app, ["targets", "--help"])

    assert result.exit_code == 0
    assert "calibrate-pocket" in result.stdout
    assert "build-pilot" in result.stdout
