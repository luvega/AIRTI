# Nelfinavir 人源反向钓靶实测与生产就绪记录

日期：2026-08-14  
案例：`nelfinavir-parent-v1`  
结论口径：本记录验证工程链路与输入冻结，不把计算排名表述为直接结合证据。

## 冻结输入

| 项目 | 冻结值 |
| --- | --- |
| 化合物 | Nelfinavir 母体，不含 M8 |
| PubChem / ChEMBL | CID 64143 / CHEMBL584 |
| InChIKey | `QAGYKUNXZHXKMR-HKWSIXNMSA-N` |
| 配体状态 | 人工审计的中性态、叔胺阳离子态；形式电荷 0、+1 |
| 配体 manifest | `/data/airti-target-fishing/cases/nelfinavir/ligand-prep-v2/prepared_ligands.jsonl` |
| 配体 manifest SHA-256 | `39c2bc6c504e74eab24a4115de93f29d9c756dae6084d345c88467cf808d36aa` |

通用质子化枚举曾产生不合理的酚盐和酰胺负离子，已经废弃且不得用于本案例。
案例输入依据实验 pKa 6.00、11.06 冻结两个可审计状态；每个状态仍独立保留三维
构象、PDBQT、SDF 和哈希。

## 人源参考库门禁

| 资产 | 条目数 | SHA-256 |
| --- | ---: | --- |
| reviewed human canonical manifest | 20,416 | `15d64442553a850d6212856b06cbcabbfc608fe66bc56baee9e40c07da2a10b4` |
| reference summary | 20,416 | `ad69e21502a4efcc0e34f51b1cb13f161a2efa96d4d38595355a2781f96e6e59` |
| 64 蛋白工程面板 | 64 | `5b9311976c7c69b0005c3fb3317ee4ee0172c6e0771b2c3c76d6bb84e13250a7` |
| 64 面板审计 | 64 | `924e2979f21435e7c60d8fc343d73cc3eb79104abc6c852932ffb1834009237b` |
| CYP3A4/HSP90AA1 代表集 | 2 | `90572a4ef9e8bfa80c150268786895fde9ab807f6ce6fe509fa7641aefe41ff7` |

20,416 条目是 reviewed、Homo sapiens、canonical 的冻结序列范围。正式发现阶段只有
在每个条目均被标记为 `ready` 或带明确原因的 `unsupported` 后才能通过门禁；64 面板
和双靶点实测均不能替代全库门禁。

## 真实引擎实测

### 代表性受体构建与背景校准

使用单一 `airti-tf:0.2.0-gpu` 镜像，对 CYP3A4（P08684）和 HSP90AA1
（P07900）执行真实结构获取、fpocket、Meeko 受体准备，以及每口袋 100 个背景探针、
3 个固定种子的 QuickVina2 校准。结果目录为：

`/data/airti-target-fishing/cases/nelfinavir/target-build-smoke-v8`

此处将在任务完成后记录门禁结果、校准成功数、结构/口袋选择和产物哈希。

### Nelfinavir 代表性筛选

此处将在受体校准通过门禁后记录两个冻结配体状态对两个代表性受体的真实
QuickVina2 输出。该小规模运行只验证数据契约、HEM 受体链路、经验背景归一化和
断点续跑，不用于报告全蛋白组候选排名。

## 已验证与尚未执行的边界

已验证：20,416 条人源 reviewed canonical 输入构建、案例和配体状态冻结、生产受体
构建实现、HEM 对接模板生成、单镜像 Nextflow 编排、可恢复检查点以及真实配体准备。

尚未执行：64 蛋白完整对接、20,416 蛋白完整受体资产构建与 Nelfinavir 全库筛选、
Top 30 Boltz-2、Top 10/Top 3 的 100 ns MD。CYP 的 MD 还必须等待经审计的
`p450-ferric-thiolate-v1` 参数适配器；缺失时流程按设计失败关闭。

