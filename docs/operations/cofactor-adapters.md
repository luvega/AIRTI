# 金属辅因子参数适配器

AIRTI 不会把 HEM 当作普通有机配体自动套用 GAFF2。细胞色素 P450 的铁—卟啉—轴向半胱氨酸体系涉及明确的氧化态、电子自旋态和金属配位键；未通过适配器门禁时，MD 必须以 `cofactor_parameter_adapter_missing` 失败关闭。

受体构建阶段所用的 `HEM_meeko_template.json` 仅服务于 QuickVina2 PDBQT。构建器同时保存原始 `HEM.cif`、`HEM_meeko_source.cif` 和 `HEM_meeko_template.provenance.json`；后者记录 RCSB CCD 哈希及为使 Meeko chemical-component 生成器处理显式酸性氢和未知 Fe 电荷所作的三处字段变换。该过程生成的是对接原子类型与部分电荷，不能替代下面的金属中心 MD 参数。

## P450 适配器目录

生产数据根目录下应建立：

```text
/data/airti-target-fishing/reference/cofactors/
└── p450-ferric-thiolate-v1/
    ├── adapter.json
    ├── cofactor.frcmod
    └── cofactor.lib
```

`adapter.json` 的最小合同如下：

```json
{
  "schema_version": "1.0",
  "parameter_id": "p450-ferric-thiolate-v1",
  "ccd_id": "HEM",
  "chemical_state": "ferric P450 heme, explicitly documented spin/protonation state",
  "source": "MCPB.py project DOI, parameter-generation record and local audit path",
  "frcmod_sha256": "64位小写十六进制SHA-256",
  "library_sha256": "64位小写十六进制SHA-256",
  "leap_lines": []
}
```

其中 `source` 必须指向实际采用的参数来源或本项目 MCPB.py/QM 生成记录；不得用空文件、通用血红蛋白参数或未经说明的网络文件代替 P450 化学态。AIRTI 会核对参数身份、来源字段和两个文件的 SHA-256。`leap_lines` 仅用于经审计的额外 LEaP 命令。

适配器通过后，膜体系使用同一统一镜像中的 AmberTools、PACKMOL-Memgen、ParmEd 和 CUDA GROMACS；体系固定为 ff19SB/Lipid21/GAFF2/TIP3P、POPC:CHL1=4:1、NaCl 0.15 M。参数文件本身属于版本化参考数据，不写入通用镜像。
