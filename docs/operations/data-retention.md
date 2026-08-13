# AIRTI 数据保留与归档规则

## 永久或项目期保留

- 原始查询文件、知情授权范围内的项目元数据及其 SHA-256；
- 冻结参数、代码 commit、镜像内容标识、模型与参考库清单；
- 全阶段 JSONL/Parquet 结果、失败与 unsupported 记录；
- Boltz Top 30 结构、Top 10 MD 输入、checkpoint、日志和最终轨迹；
- 最终报告、报告清单、工件清单、Nextflow trace/report 和发布门判定。

## 可重建缓存

- Boltz CCD 与模型缓存：跨项目共享，只读使用，按锁文件校验；
- MSA：按 UniProt ID、序列哈希和数据库版本寻址；
- 受体 PDBQT、口袋网格与背景分布：跟随参考库版本保留；
- Nextflow work：在报告发布并验证 `-resume` 后方可清理。

## 清理规则

临时下载、失败的未完成中间文件和可重建 scratch 可在审计完成后清理，但不得删除对应错误日志和状态历史。MD checkpoint 在轨迹归档并通过完整性校验前不得删除。任何清理操作都以明确项目目录为边界，不对 `/data`、`/mnt/ssd4t` 或工作区根目录执行递归删除。

归档包至少包含 `manifest.json` 和 `artifacts.sha256`。恢复归档时先验证所有哈希，再打开报告或继续计算；哈希不一致的归档视为损坏，不进入结果发布。
