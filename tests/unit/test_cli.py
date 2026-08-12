from typer.testing import CliRunner

from airti_tf.cli import app


def test_version_command_reports_package_version() -> None:
    result = CliRunner().invoke(app, ["version"])

    assert result.exit_code == 0
    assert result.stdout.strip() == "airti-tf 0.1.0"
