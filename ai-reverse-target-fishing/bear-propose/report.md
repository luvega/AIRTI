---
skill: bear-propose
topic: "AI 反向钓靶服务课题立项前评估"
date: 2026-08-10
generated_at: 2026-08-10T12:26:00+08:00
query_count: 16
result_count: 160
empty_result_count: 0
output_files: {markdown: "report.md", html: "report.html", bibtex: "references.bib"}
source_policy: "All papers, authors, and claims come from this session's sci search results."
queries:
  - {id: q1-q12, label: "撞车检测", query: "见 bear-scoop/report.md", mode: "ultra_low/low", result_count: 120, useful_count: 8}
  - {id: q13, label: "安静区：前瞻与校准", query: "prospective target deconvolution benchmark uncertainty calibration top-k drug target prediction", mode: low, result_count: 10, useful_count: 3}
  - {id: q14, label: "安静区：实验反馈", query: "closed-loop active learning chemoproteomics target identification small molecule", mode: low, result_count: 10, useful_count: 2}
  - {id: q15, label: "挑战：泄漏与泛化", query: "data leakage scaffold split cold target generalization drug target interaction prediction benchmark", mode: low, result_count: 10, useful_count: 3}
  - {id: q16, label: "挑战：对接假阳性", query: "reverse docking target fishing limitations receptor conformation scoring false positives", mode: low, result_count: 10, useful_count: 4}
---

# AI 反向钓靶服务课题立项前评估

> 一句话结论：计算模块拥挤但可直接复用；服务课题的相对安静区是“真实客户情境 + 置信度/拒判 + 正交实验 + 可追溯交付”，最需预先控制的是域外泛化和把结构分数误当机制。
> 可信边界：所有结果来自本次 sci search；安静区是相对判断，不等于没有人在做。

## 1. 一眼结论

| 指标 | 值 |
|---|---|
| 直接撞车 | 1 篇高邻近原创工作 |
| 安静区有效支撑 | 5 篇 |
| 高威胁挑战 | 2 类 |
| 总体信号 | 适合服务立项；算法创新不是必要条件 |

撞车风险对“新算法论文”较高，对“服务体系”不构成阻断。安静区中的前瞻验证和置信度量化已有部分支撑，真正持续的实验反馈闭环证据较弱，需要把它作为逐步积累的运营能力，而非首期必达算法承诺。

## 2. 核心发现

### 综合判断

| 维度 | 结论 | 依据 |
|---|---|---|
| 撞车 | 高，但主要意味着方法成熟 | E1 MAI-TargetFisher 与总体计算路线高度接近 |
| 安静区 | 中等支撑 | E2/E3 支持校准与应用导向外部验证；E4/E5 支持实验去卷积 |
| 关键挑战 | 高 | E6/E7/E8 显示对接评分、结构构象和基准偏倚不可忽视 |
| 立项解释 | 服务可行，论文创新后置 | 以案例交付、实验闭环和证据等级衡量 |

### 安静区支撑

| 层级 | 文献 ID | 判断 |
|---|---|---|
| 直接支撑 | E2, E3 | 支持不确定性引导实验和应用导向外部验证 |
| 部分支撑 | E4, E5 | 支持蛋白组去卷积与结合位点识别，但不是主动学习闭环 |
| 间接相关 | E6 | 支持严格验证策略与适用域报告 |

### 潜在挑战

| 挑战 | 威胁 | 回应思路 |
|---|---|---|
| 冷启动、数据稀疏与热门靶点偏倚 | 高 | 时间/骨架/靶点切分，报告适用域与拒判 |
| 口袋、构象与打分误差 | 高 | 先情境过滤，再多构象/共识结构复核，最终实验确认 |
| 实验假阴性与路线不匹配 | 中高 | 按化合物性质选择无标记、探针或靶向实验，保留技术失败状态 |

## 3. 研究动作建议

| 场景 | 建议动作 | 依据 |
|---|---|---|
| 首期立项 | 建 S1/S2/S3 三档服务与 L0–L4 证据等级 | 避免计算排名越级解释 |
| MVP | 两个已知对照 + 一个盲样；盲样前冻结 Top-k | E3/E6 强调真实、非重叠验证 |
| 模型栈 | 复用成熟配体中心、多模态和结构方法 | E1 表明集成框架具备可实施性 |
| 实验 | 至少一条蛋白组/细胞占有路线 + 一条遗传因果路线 | E4/E5 及化学探针验证文献 |
| 论文支线 | 等内部案例库形成后研究风险—覆盖与实验选择 | 当前闭环直接证据较少 |

## 4. 证据表

| ID | 类型 | 强度 | 作者年份 | 标题 | DOI | 关系说明 | Query |
|---|---|---|---|---|---|---|---|
| E1 | 撞车 | 强 | Li et al. 2025 | MAI-TargetFisher | 10.1038/s41401-024-01444-z | 计算骨架高度重叠，可直接作成熟基线 | q12 |
| E2 | 直接支撑 | 强 | Zhao et al. 2025 | Evidential DTI prediction | 10.1038/s41467-025-62235-6 | 不确定性用于优先安排实验，并有前瞻结合验证 | q13 |
| E3 | 直接支撑 | 强 | Daina & Zoete 2024 | Testing predictive power of reverse screening | 10.1038/s42004-024-01179-2 | 支持应用导向、外部非重叠评测 | q10/q13 |
| E4 | 部分支撑 | 强 | Piazza et al. 2020 | LiP-Quant | 10.1038/s41467-020-18071-x | 支持机器学习化学蛋白组去卷积 | q14 |
| E5 | 部分支撑 | 强 | Vu et al. 2022 | Validating small molecule chemical probes | 10.1146/annurev-biochem-032620-105344 | 支持靶点占有、选择性与功能验证组合 | q14 |
| E6 | 方法批评 | 强 | Mathai et al. 2019 | Validation strategies for target prediction methods | 10.1093/bib/bbz026 | 指出骨架和靶点家族偏倚会扭曲性能 |
| E7 | 方法批评 | 强 | Luo et al. 2024 | Benchmarking reverse docking | 10.1002/pro.5167 | 显示全蛋白反向对接 Top-k 边界 | q16 |
| E8 | 方法批评 | 强 | Wang et al. 2016 | Evaluation of ten docking programs | 10.1039/C6CP01555G | 打分难以稳定排序全数据集亲和力 | q16 |

## 5. 详细证据

### E2. Zhao et al. 2025

**标题**：[Evidential deep learning-based drug-target interaction prediction](https://doi.org/10.1038/s41467-025-62235-6)
**类型**：安静区直接支撑；**强度**：强；**对应查询**：q13。
**关系说明**：支持用不确定性信息辅助安排候选实验优先级；其前瞻案例是可参考模板，但不等于已经解决所有客户域外问题。

### E3. Daina & Zoete 2024

**标题**：[Testing the predictive power of reverse screening](https://doi.org/10.1038/s42004-024-01179-2)
**类型**：安静区直接支撑；**强度**：强；**对应查询**：q10/q13。
**关系说明**：支持大规模外部非重叠数据与应用导向评估，适合转化为服务验收标准。

### E7. Luo et al. 2024

**标题**：[Benchmarking reverse docking through AlphaFold2 human proteome](https://doi.org/10.1002/pro.5167)
**类型**：最高威胁方法批评；**强度**：强；**对应查询**：q16。
**关系说明**：最佳 Top-100 命中率仍有限，要求对接只作富集层。

## 6. 检索透明度

| Query ID | 方向 | Query | Mode | 返回 | 采用 | 说明 |
|---|---|---|---|---:|---:|---|
| q1–q12 | 撞车 | 六角度及复跑 | mixed | 120 | 8 | 计算路线拥挤 |
| q13 | 前瞻与校准 | prospective ... uncertainty ... | low | 10 | 3 | 找到 EviDTI 与可靠性工作 |
| q14 | 实验反馈 | closed-loop active learning chemoproteomics ... | low | 10 | 2 | 找到化学蛋白组，未找到完全同构闭环 |
| q15 | 泛化挑战 | data leakage scaffold split cold target ... | low | 10 | 3 | 近期结果含预印本，核心报告改用已同行评议的验证综述 |
| q16 | 对接挑战 | reverse docking ... limitations ... | low | 10 | 4 | 找到反向对接与评分限制 |

### 空结果与缺口

| 方向 | Query | 说明 |
|---|---|---|
| 完整主动学习闭环 | q14 | 未检索到与本服务完全同构且具多轮客户实验反馈的强直接证据；只判为弱安静区信号 |

## 7. 可复用数据

```json
{"collision_risk":"high_for_algorithm_low_for_service","quiet_zone_support":"moderate","closed_loop_support":"weak_to_moderate","highest_threats":["domain shift","docking score overinterpretation"],"recommended_decision":"proceed as service MVP"}
```
