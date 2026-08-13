# AIRTI 单节点生产运行手册

## 1. 适用范围

本手册适用于人源 canonical 蛋白组反向钓靶流程的单节点 Docker/NVIDIA 部署。单个任务接收 1–5 个小分子，CPU 阶段可并行，Boltz-2 与 GROMACS 在同一 GPU 上串行运行。首版不包含湿实验、客户 Web UI、多物种参考库或单任务超过 5 个分子的生产支持。

计算输出是候选靶点优先级及其可追溯证据，不等同于直接结合或作用机制已经实验确认。

## 2. 固定资产

- 代码：与运行记录中的 Git commit 一致；
- 镜像：`airti-tf:0.1.0-gpu`，内容标识见 `containers/images.lock.yaml`；
- 模型：`/mnt/ssd4t/airti-target-fishing/boltz`，哈希见 `containers/models.lock.yaml`；
- 参考库：`/data/airti-target-fishing/reference`；
- 运行工件：`/data/airti-target-fishing/runs/<project_id>`；
- 临时缓存：`/mnt/ssd4t/airti-target-fishing/cache`。

镜像包含软件，不包含 Boltz 权重、全人蛋白结构库、口袋库或任务结果。升级任一工具、模型、参考库或排序权重均视为新版本，必须重新执行 benchmark 和 smoke。

## 3. 初次部署

```bash
docker build -f containers/airti.Dockerfile -t airti-tf:0.1.0-gpu .
./scripts/run_hardware_smoke.sh
```

硬件 smoke 必须同时满足：QuickVina2 写出构象；GROMACS 日志出现 `1 GPU selected for this run.` 和 `Finished mdrun`；Boltz 日志出现 `GPU available: True` 且结构与亲和力两个阶段的失败数均为 0；模型 SHA-256 校验通过。

在不消耗 GPU 的条件下，可运行编排 smoke：

```bash
./scripts/run_orchestration_smoke.sh
```

该脚本按服务上限把 10 个查询拆成两个 5 分子批次，并验证首次执行、`-resume` 全缓存、交付清单和工件哈希。其输出均标记为 `orchestration_mock_only`，仅能证明编排行为，不能用于判断靶点检索性能。

## 4. 生产前检查

```bash
airti-tf preflight \
  --profile production \
  --output /data/airti-target-fishing/runs/<project_id>/preflight.json
```

预检采用 fail-closed 策略。GPU、Docker/NVIDIA runtime、Nextflow、统一镜像、镜像内容标识、参考清单、磁盘空间或文件句柄任一必要条件不满足时，不得启动生产流程。

输入需保存原始文件及 SHA-256。SMILES/SDF 中应含稳定分子标识；未定义立体化学、超出分子量范围、无法质子化、构象生成失败或原子数超限的分子按明确错误码保留，不得静默删除或替换。

## 5. 启动与恢复

正式运行前冻结参数文件，并记录其 SHA-256。GPU 并发固定为 1。

```bash
nextflow run workflow/main.nf \
  -profile production \
  -params-file configs/production-frozen.yaml \
  --queries /data/airti-target-fishing/runs/<project_id>/input/queries.sdf \
  --outdir /data/airti-target-fishing/runs/<project_id>/delivery \
  -with-report /data/airti-target-fishing/runs/<project_id>/audit/nextflow-report.html \
  -with-trace /data/airti-target-fishing/runs/<project_id>/audit/nextflow-trace.tsv
```

中断后使用完全相同的代码、镜像、参数和输入执行 `-resume`。不得在恢复时修改随机种子、Boltz 采样数、对接 exhaustiveness、MD 步数或排序权重。GROMACS 仅从匹配的 `md.cpt` 继续，并保留原日志。

## 6. 结果复核

交付前至少核对：

1. 参考库版本、可计算覆盖率和 unsupported 原因统计；
2. QuickVina2 三种子至少 2 个成功，背景分布每口袋至少 95 个有效探针；
3. Boltz-2 至少 2/3 种子成功、无严重碰撞，并满足口袋约束门；
4. Top 10 MD 完成时长、checkpoint、稳定性指标和异常轨迹；
5. 每个报告数值均能定位到工件 ID 与 SHA-256；
6. 报告不使用“确认靶点”“已证实直接结合”等超出计算证据的表述。

故障分类与恢复规则见 `failure-catalog.md`，保留策略见 `data-retention.md`。

## 7. 当前生产边界

截至 2026-08-13，统一镜像、三引擎硬件 smoke、任意非 root UID 检查和 10 查询模拟编排 smoke 已通过。生产 DAG 的五个 CLI 合同均已桥接；EGFR/厄洛替尼试点完成了真实配体准备、背景校准对接、三种子初筛、三种子 Boltz 精评，以及最小化、100 ps NVT 和 500 ps NPT 的体系构建恢复 smoke。

当前人源清单含 20,416 个条目，但试点参考快照仅有 EGFR 一个 `ready` 靶点，其余 20,415 个均明确标记为 `unsupported`。此外，`run-md` 默认生产 runner、100 ns 轨迹分析和 10 例真实端到端检索尚未形成通过记录。因此，现阶段可以继续扩建参考库并开展受控适配器实测，但不得启动无人值守生产任务，也不得把 EGFR 单靶试点表述为全人蛋白组检索验证或实验靶点确认。详细证据见 `docs/validation/2026-08-13-egfr-adapter-pilot.md`。
