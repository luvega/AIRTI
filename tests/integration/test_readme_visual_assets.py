import struct
from pathlib import Path


PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"


def png_dimensions(path: Path) -> tuple[int, int]:
    payload = path.read_bytes()
    assert payload[:8] == PNG_SIGNATURE
    assert payload[12:16] == b"IHDR"
    return struct.unpack(">II", payload[16:24])


def test_github_visual_assets_are_valid_pngs_with_expected_aspect_ratios() -> None:
    icon = Path("docs/assets/airti-icon.png")
    workflow = Path("docs/assets/airti-workflow.png")

    assert icon.is_file() and icon.stat().st_size > 100_000
    assert workflow.is_file() and workflow.stat().st_size > 100_000
    icon_width, icon_height = png_dimensions(icon)
    workflow_width, workflow_height = png_dimensions(workflow)
    assert icon_width == icon_height
    assert icon_width >= 1024
    assert workflow_width > workflow_height
    assert workflow_width >= 1024
