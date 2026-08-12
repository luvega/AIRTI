# AIRTI Target Fishing

AIRTI Target Fishing 是面向 1–5 个小分子的全人蛋白组反向钓靶内部计算流水线。首版主链为结构与口袋质控、QuickVina2 背景校准初筛、Boltz-2 精评、GROMACS 100 ns 分子动力学和分阶段共识排序。

项目只输出候选优先级与可追溯计算证据，不把纯计算结果表述为实验确认的直接靶点。

## 当前阶段

- 全人 canonical 蛋白清单、结构选择、口袋质控和配体准备已实现；
- QuickVina2 三种子初筛、背景校准、Boltz-2 精评和 GROMACS 100 ns 复核适配器已实现；
- 分阶段共识排序、工件追溯、报告措辞门控和 Nextflow DSL2 主链已实现；
- 现有 Docker/NVIDIA 蛋白质与多肽设计平台已审计并作为运行基础复用；
- 大数据与缓存计划分别使用 `/data/airti-target-fishing` 和 `/mnt/ssd4t/airti-target-fishing`。

研究依据、课题设计和实施规格分别见：

- [`ai-reverse-target-fishing/AI反向钓靶项目设计与课题调研报告.md`](ai-reverse-target-fishing/AI反向钓靶项目设计与课题调研报告.md)
- [`docs/superpowers/specs/2026-08-10-human-proteome-reverse-target-fishing-design.md`](docs/superpowers/specs/2026-08-10-human-proteome-reverse-target-fishing-design.md)
- [`docs/environment/2026-08-12-existing-platform-audit.md`](docs/environment/2026-08-12-existing-platform-audit.md)

## 开发入口

```bash
python -m venv --system-site-packages .venv
.venv/bin/pip install --no-deps -e .
.venv/bin/pytest -q
.venv/bin/ruff check src tests
.venv/bin/python -m mypy src/airti_tf
.venv/bin/airti-tf version
```

Nextflow 的 `test` profile 使用模拟工具验证 DAG、交付目录与 `-resume`，不消耗 GPU：

```bash
nextflow run workflow/main.nf \
  -profile test \
  --queries tests/fixtures/ligands.smi \
  --outdir results/mock
```

生产运行必须通过 `airti-tf preflight --profile production`，且需要 AIRTI 专用工具镜像。
真实计算尚需固定并构建 screening、Boltz-2 与 GROMACS 镜像；当前版本不应被表述为已完成生产级全蛋白组验证。
