# AIRTI Target Fishing

AIRTI Target Fishing 是面向 1–5 个小分子的全人蛋白组反向钓靶内部计算流水线。首版主链为结构与口袋质控、QuickVina2 背景校准初筛、Boltz-2 精评、GROMACS 100 ns 分子动力学和分阶段共识排序。

项目只输出候选优先级与可追溯计算证据，不把纯计算结果表述为实验确认的直接靶点。

## 当前阶段

- 全人 canonical 蛋白清单、结构选择、口袋质控和配体准备的 Python 组件已实现；
- QuickVina2 三种子初筛、背景校准、Boltz-2 精评和 GROMACS 100 ns 复核适配器已实现；
- 分阶段共识排序、工件追溯和报告措辞门控已实现；Nextflow DSL2 已形成可恢复的编排骨架和模拟工具测试链；
- 现有 Docker/NVIDIA 蛋白质与多肽设计平台已审计并作为运行基础复用；
- QuickVina2、Boltz-2、CUDA GROMACS、AmberTools、fpocket 与 Meeko 已封装在同一固定版本镜像 `airti-tf:0.1.0-gpu`；
- 统一镜像已在当前 RTX 4090 节点通过真实 QuickVina2 搜索、Boltz-2 结构/亲和力推理和 1000 步 GROMACS GPU smoke；
- 生产 Nextflow 所调用的 `prepare-ligands`、`screen`、`refine-boltz`、`run-md` 和 `render-report` 尚未全部接入 CLI；在该桥接完成前，不把模拟编排通过解释为真实端到端通过；
- 大数据与缓存计划分别使用 `/data/airti-target-fishing` 和 `/mnt/ssd4t/airti-target-fishing`。

研究依据、课题设计和实施规格分别见：

- [`ai-reverse-target-fishing/AI反向钓靶项目设计与课题调研报告.md`](ai-reverse-target-fishing/AI反向钓靶项目设计与课题调研报告.md)
- [`docs/superpowers/specs/2026-08-10-human-proteome-reverse-target-fishing-design.md`](docs/superpowers/specs/2026-08-10-human-proteome-reverse-target-fishing-design.md)
- [`docs/environment/2026-08-12-existing-platform-audit.md`](docs/environment/2026-08-12-existing-platform-audit.md)
- [`docs/validation/v0.1.0-hardware-smoke-report.md`](docs/validation/v0.1.0-hardware-smoke-report.md)
- [`docs/validation/v0.1.0-orchestration-smoke-report.md`](docs/validation/v0.1.0-orchestration-smoke-report.md)
- [`docs/operations/runbook.md`](docs/operations/runbook.md)

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
./scripts/run_orchestration_smoke.sh
```

该脚本将 10 个查询分成两个 5 分子批次，验证六阶段身份传递、报告哈希、unsupported 缺失值语义和无重算恢复。输出带有 `orchestration_mock_only` 标志，不产生真实靶点排序。

## 统一镜像与硬件 smoke

```bash
docker build -f containers/airti.Dockerfile -t airti-tf:0.1.0-gpu .
./scripts/run_hardware_smoke.sh
```

镜像内容标识、工具版本、SBOM 和模型检查点哈希见 `containers/images.lock.yaml`、`containers/models.lock.yaml` 与 `docs/sbom/`。模型和 CCD 数据不写入镜像，默认缓存在 `/mnt/ssd4t/airti-target-fishing/boltz`。

生产运行必须先执行：

```bash
airti-tf preflight --profile production --output preflight.json
```

当前通过的是工具链、GPU 执行路径和模拟编排验证，不等同于全人蛋白组端到端验证，也不构成靶点确认。启动首批真实分子前仍需完成生产 CLI 桥接，并冻结人源 canonical 蛋白结构/口袋库、背景探针对接分布和 10 例真实检索基准。单次生产任务仅支持 1–5 个查询分子。
