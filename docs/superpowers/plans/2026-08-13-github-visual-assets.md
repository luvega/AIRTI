# AIRTI GitHub Visual Assets Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Generate a text-free AIRTI project icon and full workflow illustration, document both with precise Chinese captions in the README, and publish the validated assets to GitHub.

**Architecture:** Two independently generated PNG assets live under `docs/assets/` and are consumed only by `README.md`. A standard-library integration test verifies PNG headers, dimensions, aspect ratios, README references, alt text, and caption presence; human visual inspection verifies the semantic sequence and absence of embedded text.

**Tech Stack:** Built-in ImageGen, PNG, GitHub-flavored Markdown/HTML, Python standard library, pytest, Git, GitHub CLI.

---

## File map

- Create `docs/assets/airti-icon.png`: square repository icon with the reverse-target-fishing visual metaphor.
- Create `docs/assets/airti-workflow.png`: wide, left-to-right, text-free workflow illustration.
- Create `tests/integration/test_readme_visual_assets.py`: deterministic file and README contract checks without adding an image-library dependency.
- Modify `README.md`: display the icon and workflow image with accessible alternative text and Chinese captions.
- Retain `docs/superpowers/specs/2026-08-13-github-visual-assets-design.md`: approved visual specification and evidence-boundary source.

### Task 1: Add the visual-asset file contract

**Files:**
- Create: `tests/integration/test_readme_visual_assets.py`

- [ ] **Step 1: Write the failing PNG contract test**

```python
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
```

- [ ] **Step 2: Run the test and verify that the missing assets fail closed**

Run:

```bash
.venv/bin/python -m pytest tests/integration/test_readme_visual_assets.py -q
```

Expected: FAIL because `docs/assets/airti-icon.png` and `docs/assets/airti-workflow.png` do not yet exist.

### Task 2: Generate and validate both text-free images

**Files:**
- Create: `docs/assets/airti-icon.png`
- Create: `docs/assets/airti-workflow.png`

- [ ] **Step 1: Generate the square project icon with built-in ImageGen**

Use this production prompt as a separate image-generation request:

```text
Use case: logo-brand
Asset type: square GitHub repository icon and README emblem
Primary request: Create a sophisticated scientific emblem for a human-proteome reverse target fishing platform. At the center, a small luminous cyan molecular structure is entering and fitting into an open abstract protein binding pocket made of smooth indigo and violet molecular surfaces. Around it, sparse neural-network nodes and one elegant orbital particle trajectory suggest AI assessment and molecular dynamics.
Scene/backdrop: complete deep navy-to-near-black background, no transparency
Style/medium: premium three-dimensional biomedical visualization, clean geometric silhouette, high scientific credibility, restrained cinematic glow
Composition/framing: perfectly square, centered subject, generous safe padding, readable at 64 to 128 pixels, no border
Lighting/mood: cyan, blue and violet bioluminescent accents; precise, calm, advanced research mood
Constraints: absolutely no text, no letters, no numbers, no glyphs, no labels, no logo wordmark, no watermark; do not create a closed bullseye, medal, check mark, trophy, or any symbol implying confirmed success
Avoid: clutter, human figures, laboratory branding, UI panels, charts, typographic shapes, random pseudo-writing
```

Save the selected generated PNG as `docs/assets/airti-icon.png`.

- [ ] **Step 2: Generate the wide workflow illustration with built-in ImageGen**

Use this production prompt as a second, independent image-generation request:

```text
Use case: scientific-educational
Asset type: wide GitHub README workflow illustration without labels
Primary request: Show a complete human-proteome reverse target fishing workflow as one continuous left-to-right scientific scene. On the far left, one to five distinct small molecules enter a vast network of diverse human protein structures and binding pockets. Next, many molecular poses stream into illuminated pockets for large-scale docking. The candidates then narrow into a neural-network-like computational core and emerge as a few high-quality protein-ligand complexes. Those complexes enter circular dynamic trajectories surrounded by subtle solvent particles and conformational motion. On the far right, only a few ranked candidate complexes remain, separated by visible space from a restrained laboratory assay vessel and detection ripples that represent a future, independent experimental validation step.
Scene/backdrop: seamless deep navy scientific space with subtle depth; each stage distinct through spatial grouping but connected by cyan-violet light flow
Style/medium: premium three-dimensional biomedical data visualization matching the AIRTI icon, scientifically plausible proteins and molecules, polished GitHub hero illustration
Composition/framing: panoramic landscape, strong left-to-right direction, progressively narrowing visual funnel, balanced negative space, no boxed panels
Lighting/mood: cyan, electric blue and violet glow on a near-black background, precise and credible rather than fantastical
Constraints: absolutely no text, no letters, no numbers, no labels, no legends, no axis marks, no glyphs, no watermark; direction must be expressed only through position, narrowing streams and elegant light-flow arrows; keep the laboratory validation motif visually separate from the computational output
Avoid: readable or pseudo-readable writing, software logos, dashboards, user interfaces, closed bullseyes, check marks, awards, definitive-success symbolism, excessive fantasy effects
```

Save the selected generated PNG as `docs/assets/airti-workflow.png`.

- [ ] **Step 3: Inspect both project files visually**

Load both files with the local image viewer at original detail. Confirm:

- neither image contains readable or pseudo-readable characters;
- the icon remains legible at thumbnail scale and has complete edges;
- the workflow reads left to right as input, proteome, docking, AI refinement, dynamics, candidate output, and separate experimental follow-up;
- no check mark, bullseye, medal, trophy, or other confirmation symbol appears;
- both images share the same dark navy, cyan, blue, and violet visual language.

If one image fails, perform one targeted ImageGen edit or regenerate only that image, preserving the other approved asset.

- [ ] **Step 4: Run the PNG contract test**

Run:

```bash
.venv/bin/python -m pytest tests/integration/test_readme_visual_assets.py -q
```

Expected: `1 passed`.

- [ ] **Step 5: Commit the validated assets and contract**

```bash
git add docs/assets/airti-icon.png docs/assets/airti-workflow.png tests/integration/test_readme_visual_assets.py
git commit -m "docs: add AIRTI visual assets"
```

### Task 3: Integrate images and captions into the README

**Files:**
- Modify: `tests/integration/test_readme_visual_assets.py`
- Modify: `README.md`

- [ ] **Step 1: Extend the test with the README consumption contract**

Append:

```python
def test_readme_references_visual_assets_with_alt_text_and_captions() -> None:
    readme = Path("README.md").read_text(encoding="utf-8")

    assert 'src="docs/assets/airti-icon.png"' in readme
    assert 'alt="AIRTI 反向钓靶项目图标"' in readme
    assert "![AIRTI 全人蛋白组反向钓靶流程示意图](docs/assets/airti-workflow.png)" in readme
    assert readme.count("图注：") >= 2
    assert "不表示全人蛋白组 ready 覆盖、100 ns MD 或湿实验靶点确认已经完成" in readme
```

- [ ] **Step 2: Run the new README test and verify that it fails**

Run:

```bash
.venv/bin/python -m pytest tests/integration/test_readme_visual_assets.py::test_readme_references_visual_assets_with_alt_text_and_captions -q
```

Expected: FAIL because the README does not yet reference the images.

- [ ] **Step 3: Add the project icon directly below the README title**

Insert below `# AIRTI Target Fishing`:

```html
<p align="center">
  <img src="docs/assets/airti-icon.png" width="220" alt="AIRTI 反向钓靶项目图标">
</p>

<p align="center"><em>图注：中央发光小分子进入抽象蛋白结合口袋，外围神经网络节点和环形粒子轨迹分别表示 AI 复合物精评与分子动力学复核；图形表达计算候选生成，不表示实验靶点确认。</em></p>
```

- [ ] **Step 4: Add the workflow image after the opening evidence-boundary paragraph**

Insert after `项目只输出候选优先级与可追溯计算证据，不把纯计算结果表述为实验确认的直接靶点。`:

```markdown
![AIRTI 全人蛋白组反向钓靶流程示意图](docs/assets/airti-workflow.png)

> 图注：图中从左到右表示 1–5 个查询小分子进入人源 canonical 蛋白结构与口袋库，依次经过背景校准的多种子批量对接、Boltz-2 多种子复合物精评和 GROMACS 分子动力学复核，最终形成少量可追溯候选，并转交独立湿实验验证。该图是目标生产流程示意，不表示全人蛋白组 ready 覆盖、100 ns MD 或湿实验靶点确认已经完成；当前通过范围以“当前阶段”和验证报告为准。
```

- [ ] **Step 5: Run README and Markdown integrity checks**

Run:

```bash
.venv/bin/python -m pytest tests/integration/test_readme_visual_assets.py -q
git diff --check
```

Expected: `2 passed`; `git diff --check` emits no output.

- [ ] **Step 6: Commit the README integration**

```bash
git add README.md tests/integration/test_readme_visual_assets.py
git commit -m "docs: present AIRTI workflow in README"
```

### Task 4: Verify and publish the GitHub update

**Files:**
- Verify: `README.md`
- Verify: `docs/assets/airti-icon.png`
- Verify: `docs/assets/airti-workflow.png`
- Verify: `tests/integration/test_readme_visual_assets.py`

- [ ] **Step 1: Run focused and full repository validation**

Run:

```bash
.venv/bin/python -m pytest tests/integration/test_readme_visual_assets.py -q
.venv/bin/python -m pytest -q
.venv/bin/ruff check src tests
.venv/bin/python -m mypy src/airti_tf
git diff --check
```

Expected: visual tests pass; the existing full suite remains green with only previously documented skips/warnings; Ruff and mypy pass; no whitespace errors.

- [ ] **Step 2: Confirm that only intended files and commits are present**

Run:

```bash
git status -sb
git log -5 --oneline
git diff origin/codex/airti-production-readiness...HEAD --stat
```

Expected: no unstaged changes; the branch contains the approved design, implementation plan, two visual assets, README integration, and tests.

- [ ] **Step 3: Push the existing feature branch**

```bash
git push -u origin codex/airti-production-readiness
```

Expected: the remote branch advances to the local HEAD.

- [ ] **Step 4: Update the existing pull request or open one draft pull request**

Run:

```bash
gh pr list --head codex/airti-production-readiness --json number,url,state,title
```

If an open PR exists, retain it and report its URL. If none exists, create one draft PR against the repository default branch with:

```bash
gh pr create --draft \
  --head codex/airti-production-readiness \
  --title "Add production-ready AIRTI target-fishing pipeline and visuals" \
  --body $'## Summary\n\n- connect the production target-fishing stage adapters and lock the unified GPU runtime\n- record the bounded EGFR adapter pilot evidence\n- add text-free AIRTI icon and full workflow visuals with explicit README captions\n\n## Evidence boundary\n\nThe EGFR run is a single-ready-target adapter pilot. It is not a full-human-proteome retrieval validation, a completed 100 ns MD result, or experimental target confirmation.\n\n## Validation\n\n- `.venv/bin/python -m pytest -q`\n- `.venv/bin/ruff check src tests`\n- `.venv/bin/python -m mypy src/airti_tf`\n- unified-image QuickVina2, Boltz-2, and GROMACS GPU hardware smoke'
```

The PR body must summarize the production adapters, locked unified GPU image, EGFR pilot evidence, text-free README visuals, evidence boundaries, and validation commands. Expected: exactly one open draft PR represents the branch.
