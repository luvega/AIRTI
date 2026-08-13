# AIRTI 故障分类与恢复规则

| 阶段 | 稳定错误码或症状 | 自动处理 | 发布影响 |
|---|---|---|---|
| 配体准备 | `invalid_smiles`、`undefined_stereochemistry`、`molecular_weight_out_of_range` | 不重试；保留原输入与失败记录 | 该分子不进入筛选 |
| 配体准备 | `protonation_failed`、`conformer_generation_failed` | 可人工检查输入；不得自动改变化学结构 | 未解决前不评分 |
| 结构/口袋 | `no_structure`、`low_coverage`、`low_confidence`、`box_too_large` | 不用数值零分替代；标为 unsupported | 计入覆盖率分母 |
| QuickVina2 | `timeout`、`resource_exhausted`、非零退出 | 同参数最多重试 1 次；可增加资源，不改变科学参数 | 少于 2/3 种子成功则无筛选分数 |
| QuickVina2 | `pdbqt_error`、`grid_error` | 不盲目重试；回到受体/网格构建 | 修复并生成新工件哈希后重跑 |
| Boltz-2 | `cuda_oom` | 降低并发后最多重试 1 次；采样参数不变 | 再次失败时保留缺失证据 |
| Boltz-2 | `missing_msa`、`invalid_yaml`、`ligand_too_large` | 回到输入构建；不以单序列模式替代生产 MSA | 不产生 Boltz 数值分数 |
| Boltz-2 | `nan_output`、`constraint_violation`、严重碰撞 | 该种子失败；至少 2 个非碰撞种子方可形成共识 | 不满足时不进入完整证据层 |
| GROMACS | 退出码 137/143、节点中断 | 仅从匹配 checkpoint 使用 `-cpi -append` 恢复 | 未达到规定时长时标注 partial |
| GROMACS | 参数化或拓扑转换失败 | 不替换力场；人工复核 GAFF2/电荷/残基 | 无 MD 分数 |
| 报告 | 工件缺失或 SHA-256 不匹配 | 禁止发布；定位上游工件，不自动改清单 | 硬阻断 |
| 预检 | 镜像、GPU、参考清单、空间或版本不满足 | 不启动流程 | 硬阻断 |

所有重试均写入任务状态历史。资源重试不得改变模型、随机种子集合或评分规则。人工修复会产生新输入工件，必须使用新的哈希和任务 ID，不能覆盖原失败记录。
