# Nelfinavir 人源反向钓靶案例

本案例仅使用 Nelfinavir 母体游离碱，不包含活性代谢物 M8。输入结构固定为
PubChem CID 64143 / ChEMBL CHEMBL584 对应立体异构体；若实验来源为甲磺酸盐，
配体准备阶段仅去除盐片段，不改变母体立体化学。

冻结身份为 InChIKey `QAGYKUNXZHXKMR-HKWSIXNMSA-N`，原始输入见
`input/nelfinavir.smi`。每次运行必须保留输入、配置、靶点 manifest、镜像和模型
哈希；不允许在看到锚点召回结果后修改口袋阈值或排序权重。

配体状态不采用通用 pKa 模型的无约束输出。[Wu 等的实验报告](https://doi.org/10.1016/S0378-4347(97)00193-X)给出 Nelfinavir pKa
6.00 和 11.06；因此输入行用
`states=` 明确冻结中性态和叔胺阳离子态，酚羟基与两个酰胺保持质子化。该覆盖会
在配体 manifest 中标记 `curated_protonation_states`，便于复核并避免将不合理的
酚盐或酰胺负离子送入全库筛选。

评价采用双轨设计：CYP3A4、CYP3A5 为金标准参照；ABCB1、DDI2、HSP90AA1、
MBTPS2 为证据类型不等价的银标准参照。锚点文件只允许在全部计算排名写出后由
`airti-tf evaluate-case` 读取，不得参与结构选择、口袋选择、参数调整或排序。

第一阶段的 64 蛋白面板仅用于工程适配，不代表全蛋白组发现。第二阶段必须使用
20,416 个 reviewed human canonical 条目均有明确状态的冻结参考库。所有输出均为
计算候选，不能表述为已确认直接结合或已证实作用机制。

## 固定计算路线

1. QuickVina2 使用种子 11/29/47、exhaustiveness 8，在每个口袋的 100 探针经验背景上校准，保留 Top 300；
2. Boltz-2 对 Top 30 进行 3 次复合物/亲和力采样；
3. GROMACS 对 Top 10 各执行 100 ns，最终 Top 3 增加两个独立重复；
4. 膜蛋白使用 ff19SB/Lipid21/GAFF2/TIP3P、POPC:CHL1=4:1、0.15 M NaCl 和 5 ns 限制性平衡；
5. 输出 Top 5 计算候选，并用 `case.yaml` 生成金/银标准召回和新候选清单。

CYP 受体对接所需 HEM 模板由 RCSB CCD 经 Meeko chemical-component 生成器构建，原始 CCD、受控电荷字段变换和哈希全部留档。该模板不等于 MD 参数；缺少经审计的 `p450-ferric-thiolate-v1` 时，CYP 的 MD 必须以 `cofactor_parameter_adapter_missing` 失败关闭。

## 执行入口

64 面板生成、受体构建和全库门禁命令见 `docs/operations/runbook.md`。面板或正式全库通过门禁后，可使用同一个生产镜像启动：

```bash
nextflow run workflow/main.nf -profile production \
  --queries cases/nelfinavir/input/nelfinavir.smi \
  --target_manifest /data/airti-target-fishing/reference/nelfinavir-v1/targets.jsonl \
  --target_assets /data/airti-target-fishing/reference/nelfinavir-v1/targets \
  --case_definition cases/nelfinavir/case.yaml \
  --screen_top_n 300 --boltz_top_n 30 --md_top_n 10 \
  --md_protocol production \
  --outdir /data/airti-target-fishing/runs/NELFINAVIR-V1/delivery \
  -with-report /data/airti-target-fishing/runs/NELFINAVIR-V1/audit/nextflow-report.html \
  -with-trace /data/airti-target-fishing/runs/NELFINAVIR-V1/audit/nextflow-trace.tsv
```

锚点评价分支只读取已经写出的 screen、Boltz 和 MD manifest。若需单独重算评价，可执行：

```bash
airti-tf evaluate-case \
  --case cases/nelfinavir/case.yaml \
  --screen <screened_candidates.jsonl> \
  --boltz <boltz_candidates.jsonl> \
  --md <md_candidates.jsonl> \
  --output <case_evaluation.json>
```
