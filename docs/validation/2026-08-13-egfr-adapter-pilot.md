# AIRTI EGFR 适配器实测试点报告

**验证日期**：2026-08-13  
**查询分子**：厄洛替尼（erlotinib）  
**试点靶点**：EGFR（UniProt P00533）  
**结论**：配体准备、QuickVina2 背景校准初筛和 Boltz-2 三种子精评已在真实二进制与 GPU 环境中完成；AmberTools/ParmEd/GROMACS 体系经电荷残差修复后完成最小化、100 ps NVT 和 500 ps NPT。该结果证明单 ready 靶点适配链可运行，但不能证明全人蛋白组检索性能、100 ns 动力学稳定性或实验靶点关系。

## 1. 证据范围与冻结输入

本次试点使用 UniProt 2026_02 人源 canonical 清单。参考清单包含 20,416 个蛋白条目，其中 EGFR 为唯一 `ready` 靶点，20,415 个条目为 `unsupported`，失败条目为 0。EGFR 采用 1M17 晶体结构定义 AQ4 口袋，Boltz 输入为 P00533 第 695–1022 位的 328 aa 激酶结构域。

| 工件 | SHA-256 |
|---|---|
| `reference_manifest.json` | `aec6da544de45ddb2055a30882dfaee2e69332d1d9ae0bb2be508d51f5a65d47` |
| `human_canonical_targets.jsonl` | `5768036cdde912fdf621bd2c1062830e0d096402708b83855fabf16d5d753e10` |
| `background_calibration_v2.json` | `3e401c52a372e5589b6302868356327ef406ec9b264485e9f0d1003dcf0944e1` |
| `screened_candidates.jsonl` | `b5824b63b51b36a57d611ec448c99741bbe764ba1db699f21e7e55bd487a63b3` |
| `boltz_candidates.jsonl` | `b8716e813c4b9fbaea808d504c36e69f992c335d32f001b78f0f7742c5dc11e7` |

参考快照位于 `/data/airti-target-fishing/reference/pilots/2026_02-egfr`，运行工件位于 `/data/airti-target-fishing/runs/EGFR-PILOT-20260813`。上述绝对路径是当前节点的审计位置，不属于可移植输入合同。

## 2. QuickVina2 初筛

背景校准面板包含 100 个探针，其中 95 个获得有效结果，5 个明确记录为失败；有效探针的口袋亲和力中位数为 −7.4 kcal/mol。厄洛替尼经标准化后形成 11 个可计算状态，11/11 状态均完成对接，未出现静默丢弃。

最终候选的三种子亲和力中位数为 −7.6 kcal/mol，种子范围为 0.1 kcal/mol，背景校准分数为 0.5625。修正后的构象一致性为 0.963331；该值只比较每个种子的最佳 `MODEL`，不再错误地混合 PDBQT 中的多个输出模式。

由于本试点只有一个 ready 靶点，`screen_rank=1` 仅表示 EGFR 在这一单靶适配器集合中成功完成计算，不能解释为其在 20,416 个蛋白中的全局排名。

## 3. Boltz-2 精评

固定种子 11、29 和 47 均成功完成结构与亲和力推理，失败种子数为 0。三种子聚合结果如下：

| 指标 | 结果 |
|---|---:|
| 成功种子 | 3/3 |
| 置信度中位数 | 0.937473 |
| ligand iPTM 中位数 | 0.980440 |
| 亲和力概率中位数 | 0.804342 |
| affinity prediction 中位数 | −1.261638 |
| 口袋约束满足率中位数 | 1.000000 |
| 置信度种子范围 | 0.001340 |
| AIRTI Boltz 分数 | 0.930564 |

三个代表构象均无严重原子碰撞，最小蛋白–配体原子距离分别为 2.901、2.938 和 2.903 Å。代表结构来自种子 11。Boltz 输出属于模型预测证据；高置信度、口袋约束满足和预测亲和力均不能单独等同于直接结合或功能调控已经得到实验确认。

## 4. GROMACS 体系构建与恢复 smoke

体系构建从 Boltz 复合物提取 2,621 个蛋白原子和 29 个配体重原子。配体经补氢后包含 52 个原子，采用 AM1-BCC 电荷和 GAFF2 参数；蛋白采用 ff19SB，溶剂采用 TIP3P，并加入 0.15 M NaCl。

首次自动构建在 `gmx grompp` 最小化预处理处失败。原因是 Amber→GROMACS 转换后，中性配体残留约 +0.002 e 的数值电荷，PME 预处理按 fail-closed 策略拒绝继续。本次修复未使用 `-maxwarn` 绕过警告，而是在以下约束下进行确定性归一化：仅允许一个 `MOL` 配体；残差绝对值不得超过 0.05 e；允许范围内的残差均匀分配到配体原子，并写出 `charge_normalization.json`。本例总修正量为 −0.001999997 e，单原子修正量为 −3.84615×10⁻⁵ e。

修复后从审计检查点手工恢复并获得以下结果：

| 阶段 | 结果 |
|---|---|
| 最小化 | 1,234 步收敛；Fmax 913.705 kJ mol⁻¹ nm⁻¹；势能 −4.82883×10⁵ kJ/mol |
| NVT | 100 ps 完成；334.755 ns/day |
| NPT | 500 ps 完成；475.532 ns/day |
| NPT 温度 | 均值 299.977 K；RMSD 1.602 K |
| NPT 压力 | 均值 −7.273 bar；RMSD 174.993 bar |
| NPT 密度 | 均值 1031.14 kg/m³；RMSD 20.682 kg/m³ |

原 `system_build_status.json` 保留首次失败状态，没有被事后改写为成功；恢复过程由后续 `min.log`、`nvt.log`、`npt.log` 和 `charge_normalization.json` 共同证明。500 ps NPT 仅用于体系构建与执行路径 smoke，时间不足以支持“体系已充分平衡”或“配体稳定结合”的判断。默认 `run-md` 尚未接通 100 ns 生产执行与轨迹指标提取，因此本次不产生 MD 分数。

## 5. 统一镜像复验

针对 Boltz 的 Numba、Matplotlib、XDG 和 Triton 缓存权限问题，统一镜像将四个缓存目录固定为 `/tmp/airti-cache/*` 并设置为 `1777`。镜像使用任意 UID `1000:1000` 完成写入检查后，在 RTX 4090 节点通过 QuickVina2、GROMACS GPU 和 Boltz GPU 三引擎硬件 smoke。

| 项目 | 记录 |
|---|---|
| 镜像 | `airti-tf:0.1.0-gpu` |
| 镜像内容标识 | `sha256:c8085c3ff5bd1080aaa4b4d291a0112fb2612d09969f7d982ebf30eda74520a9` |
| AIRTI 源码提交 | `aa357690fc3db65aac26e5940246584b7804fd6b` |
| 镜像大小 | 16,815,059,528 bytes |
| SBOM SHA-256 | `f619f87fc643454dee4e7767aef924e4cdb12605f2b810ccffae769afd251948` |
| 硬件 smoke | `/data/airti-target-fishing/runs/IMAGE-SMOKE-20260813-AA35769`；通过 |

## 6. 当前判断与下一门禁

本试点支持以下工程判断：真实配体可以沿 CLI 适配链进入三种子对接与三种子 Boltz；计算失败能够按阶段保留；统一镜像可在目标 GPU 和外部非 root UID 下运行；GROMACS 体系构建的微小电荷残差可以在严格阈值内审计修复。

本试点不支持以下结论：EGFR 已从全人蛋白组中被正确检索；当前排序具有可量化的召回率或校准度；厄洛替尼–EGFR 结合获得了新的实验确认；500 ps 平衡可以替代 100 ns 生产 MD。

进入下一阶段需依次满足：

1. 扩展人源结构/口袋库的 `ready` 覆盖，并为每个新增口袋形成独立背景校准；
2. 接通默认 100 ns GROMACS runner、checkpoint 续跑和轨迹指标提取，先以缩短时长参数完成自动化集成测试，再运行冻结的 100 ns 协议；
3. 使用两个不超过 5 个分子的真实批次完成 10 例检索基准，报告 Top-k 召回、失败率、校准和分家族结果；
4. 将计算候选用于后续靶向结合、竞争或遗传学实验设计；实验结果与计算证据分开记录。
