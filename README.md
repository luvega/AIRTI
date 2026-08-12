# AIRTI Target Fishing

AIRTI Target Fishing 是面向 1–5 个小分子的全人蛋白组反向钓靶内部计算流水线。首版主链为结构与口袋质控、QuickVina2 背景校准初筛、Boltz-2 精评、GROMACS 100 ns 分子动力学和分阶段共识排序。

项目只输出候选优先级与可追溯计算证据，不把纯计算结果表述为实验确认的直接靶点。

## 当前阶段

- 设计规格已确认；
- 现有 Docker/NVIDIA 蛋白设计平台已审计；
- 工程实现按 `docs/superpowers/plans/` 中的计划推进；
- 大数据与缓存计划分别使用 `/data/airti-target-fishing` 和 `/mnt/ssd4t/airti-target-fishing`。

## 开发入口

```bash
python -m venv --system-site-packages .venv
.venv/bin/pip install --no-deps -e .
.venv/bin/pytest -q
.venv/bin/airti-tf version
```

生产运行必须通过 `airti-tf preflight --profile production`，且需要 AIRTI 专用工具镜像。

