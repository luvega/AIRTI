---
skill: bear-scoop
topic: "面向表型活性小分子的多证据 AI 反向钓靶与正交验证服务"
date: 2026-08-10
generated_at: 2026-08-10T12:25:00+08:00
query_count: 12
result_count: 120
empty_result_count: 0
output_files:
  markdown: "report.md"
  html: "report.html"
  bibtex: "references.bib"
source_policy: "All papers, authors, and claims come from this session's sci search results."
queries:
  - {id: q1, label: "字面表述", query: '"AI target fishing" small molecule unknown target multimodal', mode: ultra_low, result_count: 10, useful_count: 3}
  - {id: q2, label: "方法中心", query: "multimodal machine learning drug target interaction reverse docking knowledge graph target identification", mode: ultra_low, result_count: 10, useful_count: 4}
  - {id: q3, label: "问题中心", query: "target deconvolution phenotypic screening unknown mechanism small molecule", mode: ultra_low, result_count: 10, useful_count: 4}
  - {id: q4, label: "结论中心", query: "prospective validation AI target prediction wet lab hit rate", mode: ultra_low, result_count: 10, useful_count: 2}
  - {id: q5, label: "相邻领域", query: "proteome-wide drug target identification chemoproteomics machine learning", mode: ultra_low, result_count: 10, useful_count: 4}
  - {id: q6, label: "可能标题", query: '"AI-guided" multimodal framework small molecule target deconvolution', mode: ultra_low, result_count: 10, useful_count: 2}
  - {id: q7, label: "字面复跑", query: "AI target fishing drug target prediction small molecule", mode: low, result_count: 10, useful_count: 5}
  - {id: q8, label: "方法复跑", query: "multimodal drug target prediction knowledge graph structure gene expression", mode: low, result_count: 10, useful_count: 5}
  - {id: q9, label: "问题复跑", query: "target deconvolution phenotypic screening chemical proteomics small molecule", mode: low, result_count: 10, useful_count: 5}
  - {id: q10, label: "结论复跑", query: "prospective validation target prediction reverse screening wet lab", mode: low, result_count: 10, useful_count: 3}
  - {id: q11, label: "相邻复跑", query: "machine learning chemoproteomics target identification LiP-MS", mode: low, result_count: 10, useful_count: 4}
  - {id: q12, label: "标题复跑", query: "proteome-wide drug target prediction AI physical modeling target fishing", mode: low, result_count: 10, useful_count: 4}
---

# AI 反向钓靶服务的邻近程度地图

> 一句话结论：计算钓靶的主要模块均已有成熟先例；最接近的是 MAI-TargetFisher，但对服务项目而言这代表可复用技术基础，不构成立项障碍。
> 可信边界：所有结果来自本次 sci search；未检索到的内容没有补写。

## 1. 一眼结论

| 指标 | 值 |
|---|---|
| 有效文献 | 8 |
| 检索方向 | 6 个正交角度 + 6 次 low 复跑 |
| 总体信号 | 计算方法拥挤，实验闭环服务仍可差异化 |
| 空结果 | 0 |

“多模态 AI + 结构反向筛选”不是空白；MAI-TargetFisher 已形成全蛋白组 AI/物理模型框架。MDTips、PLATO、SwissTargetPrediction 和 STarFish 覆盖了多模态、平台化、配体中心和天然产物路线。客户情境约束、证据追溯、正交实验和服务分级才是本项目应重点建设的部分。

## 2. 核心发现

### 撞车地图

```text
── 撞车地图 ───────── 多证据 AI 反向钓靶服务 ──

[!] 直接撞车/最近邻
    E1 Li et al. 2025 — MAI-TargetFisher
    邻近度  █████████░  7/8

    方法孪生
    E2 Xia et al. 2023 — MDTips                 5/8
    E3 Luo et al. 2024 — AlphaFold2 reverse docking  6/8

    问题孪生
    E4 Daina & Zoete 2024 — reverse screening  6/8
    E5 Ciriaco et al. 2022 — PLATO              5/8
    E6 Piazza et al. 2020 — LiP-Quant           6/8

    邻居
    E7 Cockroft et al. 2019 — STarFish          3/8
    E8 Daina et al. 2019 — SwissTargetPrediction 3/8
──────────────────────────────────────────────
```

### 四轴评分

| ID | 问题 | 方法 | 结论 | 场景 | 总分 | 分层 | 为什么在这一层 |
|---|---:|---:|---:|---:|---:|---|---|
| E1 | 2 | 2 | 2 | 1 | 7 | 直接撞车 | 同为小分子到蛋白组候选，并融合 AI 与物理模型 |
| E2 | 1 | 2 | 1 | 1 | 5 | 方法孪生 | 多模态 DTI 与候选靶点预测相同，缺少服务实验闭环 |
| E3 | 2 | 2 | 1 | 1 | 6 | 方法孪生 | 同属蛋白组反向结构筛选，主要是基准方法研究 |
| E4 | 2 | 1 | 2 | 1 | 6 | 问题孪生 | 同样从化合物反推靶点，但主要依赖配体相似性 |
| E5 | 2 | 1 | 1 | 1 | 5 | 问题孪生 | 提供同类 target-fishing 输出，方法较单一 |
| E6 | 2 | 1 | 1 | 2 | 6 | 问题孪生 | 解决同一靶点去卷积问题，但以实验蛋白组学为主 |
| E7 | 1 | 1 | 1 | 0 | 3 | 邻居 | 面向天然产物的集成模型，可作应用基线 |
| E8 | 1 | 1 | 1 | 0 | 3 | 邻居 | 成熟配体中心工具，可作候选源而非完整流程 |

## 3. 研究动作建议

| 场景 | 建议动作 | 依据 |
|---|---|---|
| 服务立项 | 直接复用/对接 E1、E2、E5、E8 的输出，不以重新训练网络为前置条件 | 方法模块已有成熟先例 |
| 产品定义 | 把差异化写成“情境过滤 + 证据卡 + L0–L4 验证等级 + 实验闭环” | 近邻工作多停留在候选排序或单一实验路线 |
| 客户报告 | 明确“候选靶点”与“确认靶点”的证据边界 | E3/E4 表明计算排名有适用域 |
| 后续论文 | 案例积累后研究校准、拒判和实验选择，而非先承诺算法创新 | 当前更安静的方向是服务数据闭环 |

## 4. 证据表

| ID | 类型 | 强度 | 作者年份 | 标题 | DOI | 关系说明 | Query |
|---|---|---|---|---|---|---|---|
| E1 | 直接撞车 | 强 | Li et al. 2025 | MAI-TargetFisher | 10.1038/s41401-024-01444-z | 最接近全蛋白组 AI 与物理模型集成 | q12 |
| E2 | 方法孪生 | 强 | Xia et al. 2023 | MDTips | 10.1093/bioinformatics/btad411 | 多模态数据融合生成候选靶点 | q2/q8 |
| E3 | 方法孪生 | 强 | Luo et al. 2024 | Benchmarking reverse docking through AlphaFold2 human proteome | 10.1002/pro.5167 | 蛋白组反向对接及性能边界 | q2 |
| E4 | 问题孪生 | 强 | Daina & Zoete 2024 | Testing the predictive power of reverse screening | 10.1038/s42004-024-01179-2 | 大规模外部集验证化合物到靶点反推 | q10 |
| E5 | 问题孪生 | 中 | Ciriaco et al. 2022 | PLATO | 10.3390/ijms23095245 | 平台化 target fishing 与生物活性画像 | q7 |
| E6 | 问题孪生 | 强 | Piazza et al. 2020 | LiP-Quant | 10.1038/s41467-020-18071-x | 机器学习化学蛋白组学靶点去卷积 | q5/q11 |
| E7 | 邻居 | 中 | Cockroft et al. 2019 | STarFish | 10.1021/acs.jcim.9b00489 | 天然产物外部集与模型集成 | q1 |
| E8 | 邻居 | 中 | Daina et al. 2019 | SwissTargetPrediction | 10.1093/nar/gkz382 | 可直接复用的配体中心候选源 | q12 |

## 5. 详细证据

### E1. Li et al. 2025

**标题**：[MAI-TargetFisher](https://doi.org/10.1038/s41401-024-01444-z)
**类型**：直接撞车；**强度**：强；**对应查询**：q12。
**关系说明**：覆盖 82% 蛋白编码基因的结构与口袋，并融合 AI 和物理模型；最适合作为本项目的现成计算骨架或外部基线。

### E2. Xia et al. 2023

**标题**：[MDTips](https://doi.org/10.1093/bioinformatics/btad411)
**类型**：方法孪生；**强度**：强；**对应查询**：q2/q8。
**关系说明**：融合知识图谱、表达和结构数据，支持多模态候选生成，但不替代客户样本的实验验证。

### E3. Luo et al. 2024

**标题**：[Benchmarking reverse docking through AlphaFold2 human proteome](https://doi.org/10.1002/pro.5167)
**类型**：方法孪生；**强度**：强；**对应查询**：q2。
**关系说明**：其最佳流程 Top-100 命中率为 27.8%，直接限定了结构筛选在服务中的角色。

### E4–E8

E4/E5/E8 共同说明配体中心反向筛选已有可复用平台；E6 提供独立实验路线；E7 显示天然产物存在域外性能下降，要求项目报告适用域和分歧。

## 6. 检索透明度

| Query ID | 方向 | Query | Mode | 返回 | 采用 | 说明 |
|---|---|---|---|---:|---:|---|
| q1/q7 | 字面表述 | AI target fishing ... | ultra_low + low | 20 | 3 | 找到平台与天然产物路线 |
| q2/q8 | 方法中心 | multimodal ... reverse docking ... | ultra_low + low | 20 | 3 | 找到 MDTips 与反向对接基准 |
| q3/q9 | 问题中心 | target deconvolution ... | ultra_low + low | 20 | 2 | 找到实验去卷积路线 |
| q4/q10 | 结论中心 | prospective validation ... | ultra_low + low | 20 | 2 | 找到外部验证工作 |
| q5/q11 | 相邻领域 | chemoproteomics ... LiP-MS | ultra_low + low | 20 | 2 | 找到 LiP-Quant |
| q6/q12 | 可能标题 | proteome-wide ... target fishing | ultra_low + low | 20 | 3 | 命中 MAI-TargetFisher |

### 空结果与缺口

| 方向 | Query | 说明 |
|---|---|---|
| 无完全空结果 | — | 未检索到不等于不存在；专利、中文文献和未公开服务可能遗漏 |

## 7. 可复用数据

```json
{"closest":"10.1038/s41401-024-01444-z","collision_risk_for_algorithm":"high","project_interpretation":"mature reusable methods; service project can proceed","service_gap":["context filtering","uncertainty and abstention","orthogonal validation","traceable delivery"]}
```
