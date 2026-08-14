# AIRTI Academic Visual Refresh Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace AIRTI's dark cinematic README images with a coherent, text-free, journal-style scientific graphical abstract system and publish the validated update to GitHub.

**Architecture:** Keep the stable README asset paths and replace only the two referenced PNG files. Use the built-in ImageGen tool independently for the square icon and wide workflow, validate each visually and structurally, update the stale icon caption through a focused contract test, then run the repository-wide quality gates before publishing through a pull request.

**Tech Stack:** Built-in ImageGen, PNG assets, Markdown, Python standard-library integration tests, pytest, Ruff, mypy, Git, GitHub CLI.

---

## File map

- Modify `tests/integration/test_readme_visual_assets.py`: encode the new academic-caption contract while preserving path, dimensions, alt text, captions, and evidence-boundary checks.
- Modify `README.md`: replace the old "发光" icon caption with wording that matches the flat academic illustration.
- Replace `docs/assets/airti-icon.png`: square, white-background academic project emblem.
- Replace `docs/assets/airti-workflow.png`: wide, white-background scientific graphical abstract.
- Preserve `docs/superpowers/specs/2026-08-14-academic-visual-refresh-design.md`: approved source of truth for composition, palette, and evidence boundaries.

### Task 1: Update the README visual contract

**Files:**
- Modify: `tests/integration/test_readme_visual_assets.py`
- Modify: `README.md:1-16`

- [ ] **Step 1: Add the new caption assertions**

Use `apply_patch` to add these assertions at the end of `test_readme_references_visual_assets_with_alt_text_and_captions`:

```python
    assert "中央小分子进入开放蛋白结合口袋" in readme
    assert "外围节点网络和弧形轨迹" in readme
    assert "中央发光小分子" not in readme
```

- [ ] **Step 2: Run the focused test and verify the new contract fails**

Run:

```bash
.venv/bin/pytest tests/integration/test_readme_visual_assets.py -q
```

Expected: one test fails because README still contains `中央发光小分子` and does not yet contain the new academic caption phrases.

- [ ] **Step 3: Replace only the stale icon caption**

Use `apply_patch` to replace the icon caption with:

```html
<p align="center"><em>图注：中央小分子进入开放蛋白结合口袋，外围节点网络和弧形轨迹分别表示 AI 复合物精评与分子动力学复核；图形表达计算候选生成，不表示实验靶点确认。</em></p>
```

Do not alter the workflow caption or its evidence-boundary sentence.

- [ ] **Step 4: Run the focused test and verify it passes**

Run:

```bash
.venv/bin/pytest tests/integration/test_readme_visual_assets.py -q
```

Expected: `2 passed`.

- [ ] **Step 5: Commit the caption contract**

```bash
git add README.md tests/integration/test_readme_visual_assets.py
git commit -m "docs: align captions with academic visuals"
```

### Task 2: Regenerate the academic project icon

**Files:**
- Replace: `docs/assets/airti-icon.png`

- [ ] **Step 1: Generate the icon in built-in mode**

Invoke the built-in ImageGen tool once with no referenced image and this prompt:

```text
Use case: logo-brand
Asset type: square GitHub README emblem for a scientific software project
Primary request: Create a restrained journal-style scientific emblem for a human-proteome reverse target-fishing workflow. Show one small ball-and-stick ligand entering an open protein binding pocket. Add only a few connected nodes behind the pocket to suggest AI complex assessment and one thin curved trajectory around part of the complex to suggest molecular-dynamics review.
Scene/backdrop: pure clean white background with generous empty space
Style/medium: flat vector-like scientific graphical abstract suitable for a Nature or Cell methods figure; crisp two-dimensional shapes, fine consistent outlines, subtle depth only where needed for scientific legibility
Composition/framing: square, centered, approximately twelve percent safe margin, one dominant ligand-pocket relationship readable at 64 to 128 pixels
Color palette: low-saturation academic blue #3B6F9C, blue-green #4C9A91, dark blue-gray #334E5C, pale blue-gray #DCE8EC, with at most one small muted warm accent atom
Constraints: open non-circular pocket; scientifically plausible molecule and pocket; absolutely no text, letters, numbers, labels, legends, glyphs, wordmarks, watermark, pseudo-writing, check marks, bullseyes, medals, trophies, or symbols implying confirmed success
Avoid: dark background, neon glow, cinematic lighting, glossy 3D rendering, metallic surfaces, lens flare, particle clouds, dashboards, interface panels, commercial advertising aesthetics, excessive detail
```

- [ ] **Step 2: Inspect the generated icon before copying it**

Open the generated file with `view_image` at original detail. Confirm all of the following:

- white background and flat journal-figure treatment;
- no text, pseudo-text, watermark, bullseye, or success symbol;
- ligand clearly enters an open protein pocket;
- nodes and one curved trajectory remain subordinate to the main subject;
- no dark cinematic background, neon glow, or glossy 3D material.

If exactly one criterion fails, issue one targeted ImageGen edit that preserves all passing elements and changes only the failed feature. Inspect the edited result again.

- [ ] **Step 3: Replace the stable project asset**

Copy the selected built-in ImageGen output from its `$CODEX_HOME/generated_images/...` location to:

```text
docs/assets/airti-icon.png
```

The user explicitly authorized replacement, so keep the stable filename rather than creating a versioned sibling.

- [ ] **Step 4: Validate the icon contract**

Run:

```bash
file docs/assets/airti-icon.png
.venv/bin/pytest tests/integration/test_readme_visual_assets.py::test_github_visual_assets_are_valid_pngs_with_expected_aspect_ratios -q
```

Expected: `file` identifies a PNG image; the focused pytest test passes; the icon is square and at least 1024 px per side.

- [ ] **Step 5: Commit the new icon**

```bash
git add docs/assets/airti-icon.png
git commit -m "docs: replace AIRTI icon with academic artwork"
```

### Task 3: Regenerate the academic workflow illustration

**Files:**
- Replace: `docs/assets/airti-workflow.png`

- [ ] **Step 1: Generate the workflow in built-in mode**

Invoke the built-in ImageGen tool once with no referenced image and this prompt:

```text
Use case: scientific-educational
Asset type: wide text-free graphical abstract for a GitHub README and research presentation
Primary request: Create a coherent left-to-right scientific graphical abstract for a human-proteome reverse target-fishing workflow. Arrange six visually distinct but unboxed stages: one to five different small-molecule inputs; a diverse human protein and binding-pocket library; many protein-ligand docking poses; a compact node network that narrows to a few refined complexes; one complex surrounded by sparse water dots and thin conformational trajectories for molecular-dynamics review; then a few ranked candidate complexes. Place a restrained assay plate or laboratory vessel farther to the right as an independent future validation entry, separated from computational candidates by visible white space and one thin muted-orange dashed connector.
Scene/backdrop: pure white background with ample negative space
Style/medium: flat vector-like scientific graphical abstract suitable for a Nature or Cell methods figure; precise two-dimensional biomedical illustration, clean consistent outlines, restrained detail, subtle shallow layering only for clarity
Composition/framing: panoramic approximately two-to-one landscape; strong left-to-right flow; visual density progressively narrows from the proteome library to a few candidates; consistent thin arrows point only right; no boxed panels
Color palette: computational stages use low-saturation academic blue #3B6F9C, blue-green #4C9A91, dark blue-gray #334E5C and pale blue-gray #DCE8EC; reserve muted orange #D58A4A only for the separated experimental-validation motif and dashed handoff
Constraints: scientifically plausible proteins, pockets, ball-and-stick ligands, network and solvent motifs; communicate input, screening, refinement, dynamics, ranked candidates and independent validation through objects and spacing alone; absolutely no text, letters, numbers, labels, legends, axes, ticks, glyphs, wordmarks, watermark or pseudo-writing; no check mark, bullseye, positive-result burst, medal, trophy or symbol implying target confirmation
Avoid: dark background, neon or volumetric glow, cinematic 3D rendering, glossy surfaces, lens flare, decorative particle rain, software interface, dashboard, cards, charts, commercial marketing style, reverse arrows, crowded composition
```

- [ ] **Step 2: Inspect the generated workflow before copying it**

Open the generated file with `view_image` at original detail. Confirm:

- white background and coherent flat academic style matching the icon;
- no text, pseudo-text, labels, watermark, boxes, UI, or chart fragments;
- all arrows and stage relationships run left to right;
- visual density narrows from proteome-scale screening to few candidates;
- molecular dynamics is expressed by solvent points and trajectories;
- orange is restricted to a visibly separated experimental-validation entry;
- no symbol implies wet-lab confirmation has already occurred.

If one criterion fails, issue one targeted ImageGen edit preserving the rest, then inspect the edited result again.

- [ ] **Step 3: Replace the stable project asset**

Copy the selected output to:

```text
docs/assets/airti-workflow.png
```

- [ ] **Step 4: Validate the workflow contract**

Run:

```bash
file docs/assets/airti-workflow.png
.venv/bin/pytest tests/integration/test_readme_visual_assets.py -q
```

Expected: `file` identifies a PNG; `2 passed`; the workflow is wider than high and at least 1024 px wide. Visually confirm the design target of at least 1536 px width and approximately 2:1 aspect ratio.

- [ ] **Step 5: Commit the new workflow**

```bash
git add docs/assets/airti-workflow.png
git commit -m "docs: replace workflow with academic graphical abstract"
```

### Task 4: Run repository-wide validation

**Files:**
- Verify: all changed files on `codex/academic-visual-refresh`

- [ ] **Step 1: Check the exact change scope**

Run:

```bash
git status --short --branch
git diff --stat origin/main...HEAD
git diff --check origin/main...HEAD
```

Expected: only the approved design/plan documents, README, the focused integration test, and two PNG assets differ; `git diff --check` prints no errors.

- [ ] **Step 2: Run the full Python test suite**

```bash
.venv/bin/pytest -q
```

Expected: all tests pass; the prior environment may continue to report three documented skips and the existing non-fatal Biopython/PyRosetta warning.

- [ ] **Step 3: Run lint and type checks**

```bash
.venv/bin/ruff check src tests
.venv/bin/python -m mypy src/airti_tf
```

Expected: Ruff reports all checks passed; mypy reports no issues.

- [ ] **Step 4: Verify the final images at original resolution**

Open both repository files with `view_image` at original detail and compare them against Sections 2–4 and 6 of the approved design specification. Do not publish if either image contains text, has a dark cinematic background, or blurs the computational-versus-experimental boundary.

### Task 5: Publish and merge the refresh

**Files:**
- Publish: branch `codex/academic-visual-refresh`

- [ ] **Step 1: Push the validated branch**

```bash
git push -u origin codex/academic-visual-refresh
```

Expected: the remote branch advances to the local HEAD.

- [ ] **Step 2: Create a pull request targeting `main`**

Use GitHub CLI to create a draft PR with this summary:

```text
- replace the cinematic README icon and workflow with text-free journal-style scientific artwork
- update the icon caption to match the flat graphical-abstract composition
- preserve the computational-candidate versus wet-lab-confirmation evidence boundary
- validate PNG contracts and the full repository quality gates
```

- [ ] **Step 3: Mark ready and merge after confirming the remote diff**

Confirm the PR contains only the intended files and is mergeable. Then mark it ready and merge with a merge commit, preserving the design, caption, and asset commits.

- [ ] **Step 4: Verify GitHub `main`**

Run:

```bash
git fetch origin main
gh pr view --json number,url,state,mergedAt,mergeCommit
gh api 'repos/luvega/AIRTI/contents/docs/assets/airti-icon.png?ref=main' --jq '[.name,.size,.sha] | @tsv'
gh api 'repos/luvega/AIRTI/contents/docs/assets/airti-workflow.png?ref=main' --jq '[.name,.size,.sha] | @tsv'
```

Expected: PR state is `MERGED`; `origin/main` contains the merge commit; both PNG assets exist under their stable paths on `main`.
