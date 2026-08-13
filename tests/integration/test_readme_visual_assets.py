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


def test_readme_references_visual_assets_with_alt_text_and_captions() -> None:
    readme = Path("README.md").read_text(encoding="utf-8")

    assert 'src="docs/assets/airti-icon.png"' in readme
    assert 'alt="AIRTI 反向钓靶项目图标"' in readme
    assert (
        "![AIRTI 全人蛋白组反向钓靶流程示意图]"
        "(docs/assets/airti-workflow.png)"
    ) in readme
    assert readme.count("图注：") >= 2
    assert (
        "不表示全人蛋白组 ready 覆盖、100 ns MD 或湿实验靶点确认已经完成"
        in readme
    )
    assert "中央小分子进入开放蛋白结合口袋" in readme
    assert "外围节点网络和弧形轨迹" in readme
    assert "中央发光小分子" not in readme
