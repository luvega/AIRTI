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

### 单镜像与硬件门禁

| 项目 | 实测值 |
| --- | --- |
| 镜像 | `airti-tf:0.2.0-gpu` |
| 本地内容标识 | `sha256:e4963bfd1fb054e657c40b169c9b64de485cb75a5ee46ab1c58413a2557df147` |
| 源码修订 | `eb8591885d6351cd19864a966f3b1fb5770c44a7` |
| 镜像大小 | 16,986,092,598 bytes |
| SPDX 2.3 包数 | 672 |
| SBOM SHA-256 | `dc823ad1762ec9f35d289fe3e7f4016d9ff80528cdc58c0b7dcabcaf38a744eb` |
| 硬件 | NVIDIA GeForce RTX 4090 |
| 硬件冒烟 | QuickVina2、CUDA GROMACS、Boltz-2 全部通过 |

硬件记录位于
`/data/airti-target-fishing/validation/v0.2.0-hardware-smoke-eb85918`；
其 `artifacts.sha256` 文件哈希为
`5f9e62eb09147e708288a0f11cee4e7d21469520ce47fec3dbac1e5951c34d6c`。
Boltz-2 结构与亲和力模型检查点均通过冻结哈希核验。

### 代表性受体构建与背景校准

使用单一 `airti-tf:0.2.0-gpu` 镜像，对 CYP3A4（P08684）和 HSP90AA1
（P07900）执行真实结构获取、fpocket、Meeko 受体准备，以及每口袋 100 个背景探针、
3 个固定种子的 QuickVina2 校准。结果目录为：

`/data/airti-target-fishing/cases/nelfinavir/target-build-smoke-v8`

| 靶点 | 结构 / 口袋 | 环境 | 成功探针 | 状态数 / 种子任务 | 背景中位数 |
| --- | --- | --- | ---: | ---: | ---: |
| CYP3A4 | PDB 5VCC / fpocket 1 | membrane，保留 HEM | 95/100 | 318 / 954 | −8.4 kcal/mol |
| HSP90AA1 | PDB 7KRJ / fpocket 32 | soluble | 95/100 | 318 / 954 | −5.6 kcal/mol |

两个靶点均为 `ready`，独立门禁结果为 `passed=true`、2/2 覆盖、0 failed、
0 unsupported。靶点 manifest SHA-256 为
`77a836a1c63a6f2f52d4c09da1cc7e570cca1ab50bbd64c69637f2efa8a399bb`，
门禁 JSON SHA-256 为
`24a8ac103cb0888903181dd43260247a6e89b954660688fe3919efe090c3ff34`。
CYP3A4 与 HSP90AA1 校准 JSON SHA-256 分别为
`93e5788cd1579ef4e1fa13275c787a2eb8e88f70dc08f8266ac2a91deb494201`
和 `ddeb935bfc11a5b751192ead1b48843d9e56114370add4e13195646e11db60c8`。

本轮以未限制单个 QuickVina 进程 CPU 数的候选镜像启动，实测出现明显过度订阅。
最终镜像已将背景校准子进程固定为 1 CPU，并由 `calibration_workers` 统一控制并发；
构建指纹同步升级，旧检查点不会被新版本误复用。

### Nelfinavir 代表性筛选

两个冻结配体状态对两个代表性受体共形成 4 个对接组合、12 个固定种子任务；全部
成功，0 failed。筛选 manifest SHA-256 为
`65b6782677409d82648d0d87d92a535b65488919fc902ff0b798877f175fe441`。

| 代表集顺序 | 靶点 | 入选状态 | 三种子中位数 | 经验分位 | 种子范围 | 构象一致性 |
| ---: | --- | --- | ---: | ---: | ---: | ---: |
| 1 | CYP3A4 | 中性态 | −11.8 kcal/mol | 0.96875 | 0.1 | 0.653192 |
| 2 | HSP90AA1 | 中性态 | −7.1 kcal/mol | 0.93750 | 1.7 | 0.000000 |

两种状态均真实完成计算；表中“入选状态”是每个靶点跨状态路由后的最佳状态。CYP3A4
在这个只含两个已知参照的工程集内位居首位，说明正对照链路可工作。HSP90AA1 虽有
较高经验分位，但三种子构象不一致，必须降级解释。该小规模运行只验证数据契约、
HEM 受体链路与经验背景归一化，不用于报告全蛋白组候选排名或新靶点发现。

阳离子态也完成全部种子，其 CYP3A4 和 HSP90AA1 三种子中位数分别为 −11.4 和
−6.9 kcal/mol，均弱于各靶点入选的中性态。

### 代表性 Boltz-2 精评

对上述两个代表集候选以生产参数执行种子 11/29/47；每个种子包含 3 个扩散结构
样本与 3 个亲和力样本。精评 manifest SHA-256 为
`72ff85cff7f2abbd88b407a1a788a2ced4c48ae149375b470aff5e15fe596734`。

| 靶点 | 状态 | 成功种子 | 关键结果 |
| --- | --- | ---: | --- |
| CYP3A4 | succeeded | 3/3 | score 0.796325；confidence 0.913126；ligand ipTM 0.886346；口袋约束覆盖 0.756757；亲和力概率 0.629070 |
| HSP90AA1 | failed | 1/3 | `insufficient_successful_seeds`；种子 11/29/47 的口袋覆盖分别为 0.25、0.50、0.45 |

CYP3A4 的 Boltz 输入明确包含蛋白链 A、Nelfinavir 链 B 和 `CCD: HEM` 链 C，三个
种子均无严重原子碰撞。HSP90AA1 的种子 11 和 47 未达到预设的 0.50 口袋覆盖门槛，
与其对接阶段构象一致性为 0 的信号方向一致，因此按预注册规则失败关闭；这不是
程序崩溃，也不应通过放宽阈值补救。

## 已验证与尚未执行的边界

已验证：20,416 条人源 reviewed canonical 输入构建、案例和配体状态冻结、生产受体
构建实现、HEM 对接模板生成、单镜像 Nextflow 编排、可恢复检查点以及真实配体准备。

尚未执行：64 蛋白完整对接、20,416 蛋白完整受体资产构建与 Nelfinavir 全库筛选、
全库 Top 30 Boltz-2、Top 10/Top 3 的 100 ns MD。本记录只完成代表集 Top 2 的
Boltz-2 精评。CYP 的 MD 还必须等待经审计的
`p450-ferric-thiolate-v1` 参数适配器；缺失时流程按设计失败关闭。
