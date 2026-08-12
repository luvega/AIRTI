# Human Proteome Reverse Target Fishing Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 将已确认的全人蛋白组 AI 反向钓靶设计实现为可复现、可恢复、可审计的内部批处理流水线，稳定处理每项目 1–5 个小分子并交付 Top 5 候选靶点及证据链。

**Architecture:** 使用 Nextflow DSL2 编排六个阶段：人源 canonical 蛋白与结构/口袋库构建、配体标准化、QuickVina2 三种子初筛与背景探针校准、Boltz-2 三种子精评、GROMACS 100 ns 分子动力学、分阶段共识排序与 Markdown 报告。Python 包负责契约、适配器、状态、排序和报告；SQLite 与内容寻址文件保存追踪信息；首版复用既有单机 Docker/NVIDIA 平台，未来按需增加 Apptainer/HPC profile。

**Tech Stack:** Python 3.11（生产镜像）、Pydantic 2、Typer、pytest、RDKit、BioPython、pandas、DuckDB/SQLite、Nextflow DSL2、Docker Compose/NVIDIA runtime、fpocket、Meeko、QuickVina2、Boltz-2、GROMACS 2025、Jinja2；Apptainer 作为未来 HPC 扩展。

---

## 0. 实施约束与里程碑

本计划实现并验证计算服务，不在本轮直接启动全人蛋白组生产计算。2026-08-12 对既有蛋白质/多肽设计平台完成宿主级审计：

- 宿主有 72 logical CPUs、503 GiB 内存和一张 49140 MiB 的 RTX 4090；
- Docker NVIDIA runtime 可用，既有 AF3 smoke 已返回 `CudaDevice(id=0)`；
- `/data` 可用约 12 TiB，`/mnt/ssd4t` 可用约 2.5 TiB；
- `/data/protein-design` 已有 AF2/AF3、BindCraft、Foundry、RFpeptide、PepMimic、Rosetta 与 benchmark 平台；
- Nextflow 以及 AIRTI 所需的 fpocket、Meeko、QuickVina2、Boltz-2、GROMACS 尚未安装。

因此，`local` 配置执行单元测试和模拟适配器；`production` 复用现有 Docker/GPU/大盘基础，但必须先补齐五个 AIRTI 专用镜像、Nextflow 和参考库。完整证据见 `docs/environment/2026-08-12-existing-platform-audit.md`。

| 里程碑 | 范围 | 通过条件 |
|---|---|---|
| M0 | 工程骨架、契约、状态与运行预检 | CPU 单元测试通过；生产预检能明确阻断当前主机 |
| M1 | 全人蛋白组、结构与口袋库 | canonical 清单可复现；所有蛋白均有 `ready` 或 `unsupported` 状态 |
| M2 | 配体、QuickVina2、背景校准 | 三种子结果完整；每分子可形成校准后的 Top 300 并集 |
| M3 | Boltz-2、GROMACS 与分阶段排序 | Top 30/10/3 路由正确；断点续跑不重复已完成任务 |
| M4 | Nextflow、报告与审计 | 10 例 smoke 端到端完成，报告可追溯到每个工件 |
| M5 | 回溯、盲测与发布 | 约 100 例回溯、20–30 例盲测完成；Success@100 ≥ 30%，技术成功率 ≥ 95% |

## 1. 目标文件结构

~~~text
.
├── README.md
├── pyproject.toml
├── environment.yml
├── configs/
│   ├── defaults.yaml
│   ├── local.yaml
│   ├── production.yaml
│   └── benchmark.yaml
├── containers/
│   ├── base.Dockerfile
│   ├── screening.Dockerfile
│   ├── boltz2.Dockerfile
│   └── gromacs.Dockerfile
├── compose/
│   └── docker-compose.yml
├── data/
│   ├── reference/background_probes_v1.smi
│   └── forcefields/amberff19sb.ff/
├── scripts/
│   ├── build_background_panel.py
│   └── run_smoke.sh
├── src/airti_tf/
│   ├── __init__.py
│   ├── cli.py
│   ├── config.py
│   ├── contracts.py
│   ├── manifest_io.py
│   ├── state.py
│   ├── runtime.py
│   ├── sources/uniprot.py
│   ├── targets/structures.py
│   ├── pockets/fpocket.py
│   ├── pockets/receptor.py
│   ├── ligands/prepare.py
│   ├── screening/quickvina.py
│   ├── screening/calibration.py
│   ├── refinement/boltz2.py
│   ├── simulation/gromacs.py
│   ├── simulation/analysis.py
│   ├── ranking/consensus.py
│   └── reporting/render.py
├── templates/target_fishing_report.md.j2
├── workflow/
│   ├── main.nf
│   ├── nextflow.config
│   └── modules/
│       ├── target_library.nf
│       ├── ligand_prep.nf
│       ├── screen.nf
│       ├── refine.nf
│       ├── md.nf
│       └── report.nf
└── tests/
    ├── fixtures/
    ├── unit/
    └── integration/
~~~

## 2. 实施任务

### Task 1: 建立可测试的 Python 工程骨架

**Files:**
- Create: `pyproject.toml`
- Create: `environment.yml`
- Create: `README.md`
- Create: `src/airti_tf/__init__.py`
- Create: `src/airti_tf/cli.py`
- Create: `tests/unit/test_cli.py`

- [ ] **Step 1: 若目录尚无 Git 元数据，初始化仓库**

Run:

~~~bash
test -d .git || git init -b codex/airti-target-fishing
git branch --show-current
~~~

Expected: 新仓库位于 `codex/airti-target-fishing` 分支；已有仓库时不覆盖历史或切换已有分支。

- [ ] **Step 2: 先写失败的 CLI 测试**

~~~python
from typer.testing import CliRunner

from airti_tf.cli import app


def test_version_command_reports_package_version():
    result = CliRunner().invoke(app, ["version"])
    assert result.exit_code == 0
    assert result.stdout.strip() == "airti-tf 0.1.0"
~~~

Run: `pytest tests/unit/test_cli.py -q`
Expected: FAIL，提示 `airti_tf` 或 `airti_tf.cli` 不存在。

- [ ] **Step 3: 创建最小包配置和 CLI**

`pyproject.toml` 的核心内容：

~~~toml
[build-system]
requires = ["setuptools>=68"]
build-backend = "setuptools.build_meta"

[project]
name = "airti-tf"
version = "0.1.0"
requires-python = ">=3.11,<3.14"
dependencies = [
  "typer>=0.16,<0.22",
  "pydantic>=2.10,<3",
  "pydantic-settings>=2.10,<3",
  "pyyaml>=6.0,<7",
  "pandas>=2.3,<3",
  "jinja2>=3.1,<4",
]

[project.optional-dependencies]
chem = ["rdkit>=2025.3"]
dev = ["pytest>=8.4,<9", "pytest-cov>=6.2,<7", "ruff>=0.12,<0.13", "mypy>=1.17,<2"]

[project.scripts]
airti-tf = "airti_tf.cli:app"

[tool.pytest.ini_options]
testpaths = ["tests"]
~~~

`src/airti_tf/cli.py`：

~~~python
import typer

app = typer.Typer(no_args_is_help=True)


@app.command()
def version() -> None:
    typer.echo("airti-tf 0.1.0")
~~~

`environment.yml` 固定 Python 3.11，并仅放 Python 层依赖；重型计算工具由独立容器提供。

- [ ] **Step 4: 安装开发包并验证**

Run:

~~~bash
python -m venv --system-site-packages .venv
.venv/bin/pip install --no-deps -e .
.venv/bin/pytest tests/unit/test_cli.py -q
~~~

Expected: `1 passed`。

- [ ] **Step 5: 提交**

~~~bash
git add pyproject.toml environment.yml README.md src/airti_tf tests/unit/test_cli.py
git commit -m "build: scaffold reverse target fishing package"
~~~

### Task 2: 固化配置与跨阶段数据契约

**Files:**
- Create: `configs/defaults.yaml`
- Create: `configs/local.yaml`
- Create: `configs/production.yaml`
- Create: `src/airti_tf/config.py`
- Create: `src/airti_tf/contracts.py`
- Create: `tests/unit/test_config.py`
- Create: `tests/unit/test_contracts.py`

- [ ] **Step 1: 先写配置覆盖与契约校验测试**

测试必须覆盖：

~~~python
def test_production_defaults_are_locked(load_settings):
    cfg = load_settings("production")
    assert cfg.screen.seeds == [11, 29, 47]
    assert cfg.screen.background_probe_count == 100
    assert cfg.routing.screen_top_n == 300
    assert cfg.routing.boltz_top_n == 30
    assert cfg.routing.md_top_n == 10
    assert cfg.routing.md_replica_top_n == 3
    assert cfg.md.duration_ns == 100


def test_unsupported_target_cannot_have_numeric_score(TargetRecord):
    with pytest.raises(ValueError):
        TargetRecord(
            uniprot_id="P00001",
            sequence="MPEPTIDE",
            status="unsupported",
            unsupported_reason="no_structure",
            calibrated_score=0.0,
        )
~~~

Run: `pytest tests/unit/test_config.py tests/unit/test_contracts.py -q`
Expected: FAIL，因为加载器和模型尚不存在。

- [ ] **Step 2: 定义 Pydantic 模型**

至少定义以下模型及枚举：

~~~python
class TargetStatus(StrEnum):
    READY = "ready"
    UNSUPPORTED = "unsupported"
    FAILED = "failed"


class StageStatus(StrEnum):
    PENDING = "pending"
    RUNNING = "running"
    SUCCEEDED = "succeeded"
    FAILED = "failed"
    SKIPPED = "skipped"


class TargetRecord(BaseModel):
    uniprot_id: str
    sequence: str
    status: TargetStatus
    unsupported_reason: str | None = None
    calibrated_score: float | None = None

    @model_validator(mode="after")
    def preserve_missingness(self):
        if self.status != TargetStatus.READY and self.calibrated_score is not None:
            raise ValueError("unsupported or failed targets cannot receive a numeric score")
        return self
~~~

同时为 `LigandRecord`、`PocketRecord`、`DockingRecord`、`BoltzRecord`、`MDRecord`、`RankedTarget` 和 `ArtifactRecord` 定义 schema 版本字段、主键、输入哈希、工具版本和状态。

- [ ] **Step 3: 写入锁定配置**

`configs/defaults.yaml`：

~~~yaml
schema_version: "1.0"
screen:
  backend: quickvina2
  seeds: [11, 29, 47]
  exhaustiveness: 8
  background_probe_count: 100
routing:
  screen_top_n: 300
  boltz_top_n: 30
  md_top_n: 10
  md_replica_top_n: 3
boltz:
  diffusion_samples: 3
  affinity_samples: 3
md:
  duration_ns: 100
  timestep_fs: 2
  forcefield: amber_ff19sb
  ligand_forcefield: gaff2
  water_model: tip3p
report:
  final_top_n: 5
~~~

生产配置只允许通过显式 `--allow-locked-override` 修改上述锁定项，且必须记录到运行清单。

- [ ] **Step 4: 运行测试与静态检查**

Run:

~~~bash
.venv/bin/pytest tests/unit/test_config.py tests/unit/test_contracts.py -q
.venv/bin/ruff check src tests
.venv/bin/mypy src/airti_tf
~~~

Expected: 全部通过，无类型错误。

- [ ] **Step 5: 提交**

~~~bash
git add configs src/airti_tf/config.py src/airti_tf/contracts.py tests/unit/test_config.py tests/unit/test_contracts.py
git commit -m "feat: define pipeline configuration and contracts"
~~~

### Task 3: 实现内容寻址清单与 SQLite 状态机

**Files:**
- Create: `src/airti_tf/manifest_io.py`
- Create: `src/airti_tf/state.py`
- Create: `tests/unit/test_manifest_io.py`
- Create: `tests/unit/test_state.py`

- [ ] **Step 1: 写失败测试**

覆盖同一输入产生同一 SHA-256、原子写入、状态迁移合法性和重复任务幂等性：

~~~python
def test_completed_task_is_not_claimed_twice(state_store):
    task_id = state_store.register("dock", input_hash="abc")
    state_store.transition(task_id, "running")
    state_store.transition(task_id, "succeeded", output_hash="def")
    assert state_store.claim(task_id) is False


def test_failed_task_can_be_retried_with_audit_entry(state_store):
    task_id = state_store.register("boltz", input_hash="abc")
    state_store.transition(task_id, "running")
    state_store.transition(task_id, "failed", error_code="OOM")
    assert state_store.retry(task_id, max_attempts=2) is True
    assert state_store.history(task_id)[-1].status == "pending"
~~~

Run: `pytest tests/unit/test_manifest_io.py tests/unit/test_state.py -q`
Expected: FAIL。

- [ ] **Step 2: 实现 manifest JSONL 与状态表**

SQLite 至少包含 `runs`、`tasks`、`artifacts`、`transitions` 四张表。允许的迁移仅为：

~~~text
pending -> running
running -> succeeded
running -> failed
failed -> pending
pending -> skipped
~~~

所有文件先写同目录临时文件、`fsync` 后执行原子 `replace`；工件主键为内容 SHA-256，不得用文件名代替。

- [ ] **Step 3: 验证并提交**

Run: `.venv/bin/pytest tests/unit/test_manifest_io.py tests/unit/test_state.py -q`
Expected: 所有测试通过。

~~~bash
git add src/airti_tf/manifest_io.py src/airti_tf/state.py tests/unit/test_manifest_io.py tests/unit/test_state.py
git commit -m "feat: add auditable manifests and resumable state"
~~~

### Task 4: 建立本地与生产运行预检门

**Files:**
- Create: `src/airti_tf/runtime.py`
- Create: `tests/unit/test_runtime.py`
- Modify: `src/airti_tf/cli.py`

- [ ] **Step 1: 写失败的预检测试**

~~~python
def test_production_preflight_blocks_insufficient_host(fake_host):
    fake_host.data_free_bytes = 12 * 1024**4
    fake_host.cache_free_bytes = 2 * 1024**4
    fake_host.commands = {"docker": True, "nextflow": False}
    fake_host.gpus = [{"memory_mib": 49140, "container_visible": True}]
    result = run_preflight(profile="production", host=fake_host)
    assert result.ok is False
    assert set(result.blockers) >= {
        "nextflow_missing",
        "airti_images_missing",
    }


def test_local_profile_allows_mock_execution(fake_host):
    result = run_preflight(profile="local", host=fake_host)
    assert result.ok is True
~~~

Run: `pytest tests/unit/test_runtime.py -q`
Expected: FAIL。

- [ ] **Step 2: 实现 `airti-tf preflight`**

生产门控精确检查：

- 大工件根目录位于 `/data` 且可用空间 ≥ 1 TiB；
- 模型/高速缓存根目录位于 `/mnt/ssd4t` 且可用空间 ≥ 500 GiB；
- Nextflow ≥ 24.10；
- Docker Compose 与 NVIDIA runtime 可用；
- GPU 至少 1 张、显存 ≥ 40 GiB，并能在容器内运行计算 smoke；
- 五个 AIRTI 专用镜像可用且命令契约通过；
- 生产镜像 digest 均已锁定；
- 参考库 manifest 校验和全部通过；
- `ulimit -n` ≥ 65536。

命令返回码：通过为 0；可修复阻断为 2；配置或清单损坏为 3。输出同时生成 `preflight.json`。

- [ ] **Step 3: 测试当前主机的预期阻断**

Run:

~~~bash
.venv/bin/airti-tf preflight --profile local
.venv/bin/airti-tf preflight --profile production
~~~

Expected: local 返回 0；当前 production 返回 2，并报告 Nextflow 与 AIRTI 专用镜像阻断；GPU、Docker 和容量检查通过。

- [ ] **Step 4: 提交**

~~~bash
git add src/airti_tf/runtime.py src/airti_tf/cli.py tests/unit/test_runtime.py
git commit -m "feat: gate local and production runtimes"
~~~

### Task 5: 构建版本化的人源 canonical 蛋白清单

**Files:**
- Create: `src/airti_tf/sources/uniprot.py`
- Create: `tests/fixtures/uniprot_human_sample.tsv`
- Create: `tests/unit/test_uniprot.py`
- Modify: `src/airti_tf/cli.py`

- [ ] **Step 1: 写解析与去重失败测试**

测试夹具包含 reviewed/unreviewed、isoform 和重复 accession；断言只保留人源 reference proteome 的 canonical accession，一条 accession 一行，序列只含标准氨基酸字母。

~~~python
def test_builds_canonical_manifest(uniprot_fixture):
    records = parse_uniprot_tsv(uniprot_fixture, release="2026_03")
    assert [r.uniprot_id for r in records] == ["P00533", "P04637"]
    assert all(r.taxonomy_id == 9606 for r in records)
    assert all(r.sequence_sha256 for r in records)
~~~

Run: `pytest tests/unit/test_uniprot.py -q`
Expected: FAIL。

- [ ] **Step 2: 实现下载与快照固定**

命令：

~~~bash
airti-tf targets fetch-uniprot \
  --proteome UP000005640 \
  --release 2026_03 \
  --out work/reference/uniprot
~~~

下载必须保存请求 URL、HTTP ETag、返回日期、原始压缩文件 SHA-256、记录数和解析器版本。若指定 release 与响应元数据不一致，立即失败。网络测试使用录制夹具，CI 不访问在线 API。

- [ ] **Step 3: 生成清单并核对不变量**

输出 `human_canonical_proteome.jsonl` 与 `human_canonical_proteome.summary.json`。不变量：

- accession 唯一；
- taxonomy ID 全为 9606；
- 不丢弃 unreviewed canonical 蛋白；
- isoform 仅记录为别名，不作为独立主索引；
- 每条记录有序列哈希。

- [ ] **Step 4: 测试并提交**

Run: `.venv/bin/pytest tests/unit/test_uniprot.py -q`
Expected: 所有测试通过。

~~~bash
git add src/airti_tf/sources/uniprot.py src/airti_tf/cli.py tests/fixtures/uniprot_human_sample.tsv tests/unit/test_uniprot.py
git commit -m "feat: build versioned human canonical proteome manifest"
~~~

### Task 6: 选择结构并显式保留 unsupported 蛋白

**Files:**
- Create: `src/airti_tf/targets/structures.py`
- Create: `tests/fixtures/structure_candidates.json`
- Create: `tests/unit/test_structures.py`

- [ ] **Step 1: 写结构优先级与缺失状态测试**

~~~python
def test_prefers_usable_experimental_structure():
    chosen = choose_structure(
        [
            candidate("alphafold", coverage=0.99, confidence=0.88),
            candidate("pdb", coverage=0.85, resolution=2.1, has_ligand=True),
        ]
    )
    assert chosen.source == "pdb"


def test_no_usable_structure_is_preserved_as_unsupported():
    result = choose_structure([])
    assert result.status == "unsupported"
    assert result.unsupported_reason == "no_structure"
    assert result.score is None
~~~

Run: `pytest tests/unit/test_structures.py -q`
Expected: FAIL。

- [ ] **Step 2: 实现确定性结构选择规则**

候选结构先做硬门控：序列映射一致、目标域覆盖 ≥ 70%、无大范围主链缺失、无无法解释的共价修饰。通过后按以下顺序选择：

1. 含已知配体且分辨率 ≤ 3.0 Å 的实验结构；
2. 不含配体但分辨率 ≤ 2.8 Å、目标域完整的实验结构；
3. AlphaFold 结构，口袋邻域平均 pLDDT ≥ 70 且 PAE 支持域构象；
4. 均不满足时标为 `unsupported`，原因限定为 `no_structure`、`low_coverage`、`low_confidence`、`sequence_mismatch` 或 `unsupported_chemistry`。

结构标准化只移除与目标无关的晶体添加物；金属、辅因子、结构水的保留决策写入工件 manifest，禁止静默删除。

- [ ] **Step 3: 测试与提交**

Run: `.venv/bin/pytest tests/unit/test_structures.py -q`
Expected: 所有测试通过。

~~~bash
git add src/airti_tf/targets/structures.py tests/fixtures/structure_candidates.json tests/unit/test_structures.py
git commit -m "feat: select target structures with explicit coverage states"
~~~

### Task 7: 发现、质控并版本化蛋白口袋

**Files:**
- Create: `src/airti_tf/pockets/fpocket.py`
- Create: `tests/fixtures/fpocket_output/`
- Create: `tests/unit/test_fpocket.py`
- Modify: `src/airti_tf/cli.py`

- [ ] **Step 1: 写 fpocket 解析与口袋门控失败测试**

~~~python
def test_pocket_ids_are_stable_across_parser_runs(fpocket_fixture):
    first = parse_fpocket(fpocket_fixture, target_id="P00533")
    second = parse_fpocket(fpocket_fixture, target_id="P00533")
    assert [p.pocket_id for p in first] == [p.pocket_id for p in second]


def test_rejects_buried_or_low_confidence_pocket():
    result = qc_pocket(
        pocket(volume_a3=55, druggability=0.02, mean_plddt=51, exposed_fraction=0.01)
    )
    assert result.status == "unsupported"
    assert result.unsupported_reason in {"volume_too_small", "low_confidence", "inaccessible"}
~~~

Run: `pytest tests/unit/test_fpocket.py -q`
Expected: FAIL。

- [ ] **Step 2: 封装 fpocket CLI 与解析器**

适配器实际执行：

~~~bash
fpocket -f standardized_target.pdb
~~~

记录 fpocket 版本、镜像 digest、输入结构哈希、标准输出、标准错误和运行时间。`pocket_id` 使用 `target_id + structure_hash + ranked_fpocket_number` 生成，不依赖临时目录。

- [ ] **Step 3: 固化口袋 QC**

硬门控为：体积 ≥ 80 Å³、至少 6 个口袋邻域残基、AlphaFold 来源的邻域平均 pLDDT ≥ 70、网格不与主链严重冲突。保留最多 5 个合格口袋；已知配体口袋优先，其余按 fpocket druggability、体积、疏水/极性平衡和结构置信度的固定权重排序。无合格口袋的蛋白标为 `unsupported:no_qualified_pocket`。

- [ ] **Step 4: 测试与提交**

Run: `.venv/bin/pytest tests/unit/test_fpocket.py -q`
Expected: 所有测试通过。

~~~bash
git add src/airti_tf/pockets/fpocket.py src/airti_tf/cli.py tests/fixtures/fpocket_output tests/unit/test_fpocket.py
git commit -m "feat: detect and quality-control target pockets"
~~~

### Task 8: 生成受体网格并标准化 1–5 个查询配体

**Files:**
- Create: `src/airti_tf/pockets/receptor.py`
- Create: `src/airti_tf/ligands/prepare.py`
- Create: `tests/fixtures/ligands.smi`
- Create: `tests/unit/test_receptor.py`
- Create: `tests/unit/test_ligand_prepare.py`

- [ ] **Step 1: 写失败的配体和网格测试**

~~~python
def test_receptor_box_contains_all_pocket_atoms():
    box = build_box(pocket_atoms, padding_a=5.0, min_size_a=18.0, max_size_a=32.0)
    assert all(box.contains(atom.coord) for atom in pocket_atoms)


def test_undefined_stereochemistry_blocks_production():
    result = prepare_ligand("CC(O)C(=O)O", profile="production")
    assert result.status == "failed"
    assert result.error_code == "undefined_stereochemistry"


def test_query_count_is_limited_to_five():
    with pytest.raises(ValueError, match="1 to 5"):
        validate_query_batch([f"C{i}" for i in range(6)])
~~~

Run: `pytest tests/unit/test_receptor.py tests/unit/test_ligand_prepare.py -q`
Expected: FAIL。

- [ ] **Step 2: 实现受体准备**

使用 Meeko `mk_prepare_receptor.py` 生成 PDBQT 与 Vina box。网格中心取口袋原子几何中心，边界为口袋包围盒各向外扩 5 Å，单边限制在 18–32 Å；超过 32 Å 的口袋拆分或标为 `box_too_large`，不得静默截断。产物同时保存 `box.txt` 与 JSON 坐标。

- [ ] **Step 3: 实现配体标准化**

依次执行：解析 SMILES/SDF、去盐但保留主成分映射、标准电荷、互变异构体枚举、pH 7.4 ± 1.0 质子化、立体化学检查、ETKDGv3 三维构象和 MMFF94s 最小化。生产模式下：

- 输入必须是 1–5 个唯一分子；
- 分子量 100–900 Da；
- Boltz 路由前检查总原子数 ≤ 128，并对 >56 原子的结果加高不确定性标记；
- 未定义立体中心必须由项目输入明确指定，不能自动任选；
- 每个微观状态有确定性 `ligand_state_id`，最终报告能映射回原始输入。

- [ ] **Step 4: 测试与提交**

Run: `.venv/bin/pytest tests/unit/test_receptor.py tests/unit/test_ligand_prepare.py -q`
Expected: 所有测试通过。

~~~bash
git add src/airti_tf/pockets/receptor.py src/airti_tf/ligands/prepare.py tests/fixtures/ligands.smi tests/unit/test_receptor.py tests/unit/test_ligand_prepare.py
git commit -m "feat: prepare receptor grids and query ligands"
~~~

### Task 9: 实现 PocketVina 统一接口与 QuickVina2 三种子初筛

**Files:**
- Create: `src/airti_tf/screening/quickvina.py`
- Create: `tests/fixtures/quickvina_output.txt`
- Create: `tests/unit/test_quickvina.py`
- Create: `tests/integration/test_quickvina_adapter.py`

- [ ] **Step 1: 写命令生成、输出解析和失败分类测试**

~~~python
def test_command_contains_locked_seed_and_box(job):
    cmd = build_quickvina_command(job, seed=29)
    assert cmd[0] == "qvina2"
    assert ["--seed", "29"] == cmd[cmd.index("--seed"):cmd.index("--seed") + 2]
    assert "--center_x" in cmd and "--size_z" in cmd


def test_three_seed_summary_uses_median():
    records = [
        docking(seed=11, affinity=-8.0),
        docking(seed=29, affinity=-7.0),
        docking(seed=47, affinity=-9.0),
    ]
    result = summarize_seeds(records)
    assert result.affinity_median == -8.0
    assert result.seed_success_count == 3
~~~

Run: `pytest tests/unit/test_quickvina.py tests/integration/test_quickvina_adapter.py -q`
Expected: FAIL。

- [ ] **Step 2: 实现标准适配器**

每个 `ligand_state × pocket` 固定运行种子 11、29、47：

~~~bash
qvina2 \
  --receptor receptor.pdbqt \
  --ligand ligand.pdbqt \
  --center_x X --center_y Y --center_z Z \
  --size_x SX --size_y SY --size_z SZ \
  --exhaustiveness 8 --num_modes 9 --seed 11 \
  --out poses.seed11.pdbqt --log docking.seed11.log
~~~

`pocketvina-screen` 接口只暴露标准 `DockingJob -> DockingRecord`，后端命令和输出解析封装在适配器内。三种子至少成功 2 个才形成分数；1 个或 0 个成功标为技术失败并进入一次资源加倍重试。超时、PDBQT 错误、网格错误和非零退出码使用不同错误码。

- [ ] **Step 3: 添加 Vina 对照后端测试**

同一 10 例 smoke 集上运行 AutoDock Vina 作为对照，只比较候选稳定性和技术成功率，不混入 QuickVina2 生产分数。适配器契约测试确保两个后端产生同一 schema。

- [ ] **Step 4: 测试与提交**

Run:

~~~bash
.venv/bin/pytest tests/unit/test_quickvina.py tests/integration/test_quickvina_adapter.py -q
~~~

Expected: 单元测试全部通过；无 qvina2 二进制时集成测试明确 skip，而不是伪造通过。

~~~bash
git add src/airti_tf/screening/quickvina.py tests/fixtures/quickvina_output.txt tests/unit/test_quickvina.py tests/integration/test_quickvina_adapter.py
git commit -m "feat: add seeded QuickVina2 screening adapter"
~~~

### Task 10: 建立 100 个背景探针、跨蛋白校准与 Top 300 路由

**Files:**
- Create: `data/reference/background_probes_v1.smi`
- Create: `scripts/build_background_panel.py`
- Create: `src/airti_tf/screening/calibration.py`
- Create: `tests/unit/test_calibration.py`
- Create: `tests/integration/test_screen_routing.py`

- [ ] **Step 1: 写经验分位数与多口袋聚合失败测试**

~~~python
def test_empirical_percentile_is_directionally_correct():
    background = [-10.0, -8.0, -7.0, -6.0, -5.0]
    assert empirical_percentile(query=-9.0, background=background) == pytest.approx(5 / 6)


def test_top300_is_union_across_ligand_states_and_pockets():
    routed = route_screen_candidates(records, top_n=300)
    assert len(routed) <= 300
    assert routed["target_id"].is_unique
    assert set(routed.columns) >= {"target_id", "best_pocket_id", "best_state_id", "calibrated_score"}
~~~

Run: `pytest tests/unit/test_calibration.py tests/integration/test_screen_routing.py -q`
Expected: FAIL。

- [ ] **Step 2: 构建固定背景探针**

`scripts/build_background_panel.py` 从许可明确的公开分子集合中按 ECFP4 MaxMin 选择 100 个代表分子，并固定：

- 10 个分子量分层，每层 10 个；
- cLogP、HBD、HBA、可旋转键和形式电荷覆盖；
- 排除 PAINS、反应性基团、金属和未定义立体化学；
- 保存 canonical SMILES、来源 ID、描述符、选择脚本版本和输入集合哈希。

一旦用于生产，`background_probes_v1.smi` 只读；任何改变提升版本号并重建所有口袋背景分布。

- [ ] **Step 3: 实现每口袋经验校准**

每个合格口袋预计算 100 个背景探针的 QuickVina2 三种子中位数，保存经验分布。查询分数计算：

~~~python
percentile = (count(background_affinity >= query_affinity) + 1) / (len(background_affinity) + 1)
calibrated_score = percentile
~~~

其中 Vina 亲和力越负越好，因此 `query=-9` 应优于多数 `-6`。同时输出 z-score 仅作诊断，最终路由使用经验分位数。背景不足 95 个有效值的口袋不进入生产筛选。

- [ ] **Step 4: 实现靶点层聚合与多样性约束**

同一蛋白跨口袋、微观状态取最高校准分位数；并用第二好口袋、三种子离散度和姿势聚类一致性作次级排序。先取全局 Top 240，再按蛋白家族保留最多 15 个/家族并从未覆盖家族补足至最多 300；输出选择原因，禁止为了凑满 300 给 unsupported 蛋白赋分。

- [ ] **Step 5: 测试与提交**

Run:

~~~bash
.venv/bin/pytest tests/unit/test_calibration.py tests/integration/test_screen_routing.py -q
python scripts/build_background_panel.py --check-only data/reference/background_probes_v1.smi
~~~

Expected: 测试通过；探针正好 100 个、唯一、校验和稳定。

~~~bash
git add data/reference/background_probes_v1.smi scripts/build_background_panel.py src/airti_tf/screening/calibration.py tests/unit/test_calibration.py tests/integration/test_screen_routing.py
git commit -m "feat: calibrate cross-target docking and route top candidates"
~~~

### Task 11: 接入 Boltz-2 复合物与亲和力精评

**Files:**
- Create: `src/airti_tf/refinement/boltz2.py`
- Create: `tests/fixtures/boltz2_output/`
- Create: `tests/unit/test_boltz2.py`
- Create: `tests/integration/test_boltz2_adapter.py`

- [ ] **Step 1: 写 YAML、MSA、输出解析与路由失败测试**

~~~python
def test_boltz_yaml_contains_affinity_and_pocket_constraint(job):
    payload = build_boltz_yaml(job)
    assert payload["properties"][0]["affinity"]["binder"] == "B"
    assert payload["constraints"][0]["pocket"]["binder"] == "B"
    assert payload["sequences"][0]["protein"]["msa"].endswith(".a3m")


def test_missing_cached_msa_blocks_production(job_without_msa):
    with pytest.raises(MissingMSAError):
        build_boltz_yaml(job_without_msa, profile="production")
~~~

Run: `pytest tests/unit/test_boltz2.py tests/integration/test_boltz2_adapter.py -q`
Expected: FAIL。

- [ ] **Step 2: 生成官方 schema 兼容输入**

每个一级 Top 300 靶点以蛋白链 A、配体链 B 建立 YAML，包含：

~~~yaml
version: 1
sequences:
  - protein:
      id: A
      sequence: MPEPTIDE
      msa: /cache/msa/P00533.a3m
  - ligand:
      id: B
      smiles: "CCO"
constraints:
  - pocket:
      binder: B
      contacts: [[A, 718], [A, 719]]
properties:
  - affinity:
      binder: B
~~~

蛋白 MSA 按 UniProt accession、序列哈希和 MSA 数据库版本缓存。生产运行禁止 `msa: empty`；缓存缺失先进入独立 MSA 构建任务。

- [ ] **Step 3: 执行三种子精评并解析质量指标**

每个候选固定运行三个随机种子，单次命令：

~~~bash
boltz predict input.yaml \
  --out_dir output \
  --diffusion_samples 3 \
  --diffusion_samples_affinity 3 \
  --use_potentials
~~~

解析复合物置信度、链间置信度、口袋约束满足度、碰撞、亲和力预测和不确定性。至少 2/3 种子成功且无严重碰撞才形成 Boltz 共识；Top 300 路由到 Top 30 时同时考虑 docking 校准分位数、Boltz 置信度、亲和力与种子一致性。Top 30 再依据结构证据完整性选择 Top 10 进入 MD。

- [ ] **Step 4: 添加缓存和 OOM 重试**

失败分类包括 `missing_msa`、`invalid_yaml`、`ligand_too_large`、`cuda_oom`、`nan_output` 和 `constraint_violation`。`cuda_oom` 仅重试一次，降低并发而不改变模型参数；再次失败保留为缺失证据，不赋零分。

- [ ] **Step 5: 测试与提交**

Run:

~~~bash
.venv/bin/pytest tests/unit/test_boltz2.py tests/integration/test_boltz2_adapter.py -q
~~~

Expected: 夹具测试通过；无 Boltz/GPU 时真实集成测试明确 skip。

~~~bash
git add src/airti_tf/refinement/boltz2.py tests/fixtures/boltz2_output tests/unit/test_boltz2.py tests/integration/test_boltz2_adapter.py
git commit -m "feat: add Boltz-2 refinement with cached MSA support"
~~~

### Task 12: 实现 GROMACS 100 ns MD、续跑和轨迹分析

**Files:**
- Create: `data/forcefields/amberff19sb.ff/README.provenance.md`
- Create: `src/airti_tf/simulation/gromacs.py`
- Create: `src/airti_tf/simulation/analysis.py`
- Create: `tests/fixtures/gromacs/`
- Create: `tests/unit/test_gromacs.py`
- Create: `tests/unit/test_md_analysis.py`
- Create: `tests/integration/test_gromacs_resume.py`

- [ ] **Step 1: 写拓扑、时长、续跑和指标失败测试**

~~~python
def test_md_protocol_is_exactly_100_ns():
    mdp = render_production_mdp()
    assert mdp["dt"] == pytest.approx(0.002)
    assert mdp["nsteps"] == 50_000_000


def test_resume_uses_checkpoint_when_present(tmp_path):
    checkpoint = tmp_path / "md.cpt"
    checkpoint.write_bytes(b"checkpoint")
    cmd = build_mdrun_command(tmp_path)
    assert cmd[cmd.index("-cpi") + 1] == str(checkpoint)
    assert "-append" in cmd


def test_unfinished_trajectory_has_no_final_md_score():
    result = analyze_trajectory(frames=20_000, expected_frames=50_001)
    assert result.status == "failed"
    assert result.md_score is None
~~~

Run: `pytest tests/unit/test_gromacs.py tests/unit/test_md_analysis.py tests/integration/test_gromacs_resume.py -q`
Expected: FAIL。

- [ ] **Step 2: 固化力场来源与体系构建**

蛋白使用 Amber ff19SB，配体使用 GAFF2，水模型 TIP3P。`README.provenance.md` 必须记录 ff19SB GROMACS port 的上游来源、许可证、版本、下载 SHA-256、转换过程和与参考能量的验证结果。配体 AM1-BCC 电荷和 GAFF2 参数生成失败时标为 `ligand_parameterization_failed`，不改用未经确认的替代力场。

体系步骤固定为：复合物清理 → 配体参数 → dodecahedron 周期盒且溶质到边界 ≥ 1.0 nm → 加水 → 0.15 M NaCl 并中和 → steepest descent 最小化 → 100 ps NVT → 1 ns NPT → 100 ns production。

- [ ] **Step 3: 生成精确 MDP 与 checkpoint 命令**

production 采用 2 fs 步长、50,000,000 步，每 10 ps 保存坐标与能量，温度 300 K、压力 1 bar。执行：

~~~bash
gmx grompp -f md.mdp -c npt.gro -p topol.top -o md.tpr
gmx mdrun -deffnm md -cpi md.cpt -append
~~~

Top 10 各运行 1 个 100 ns 轨迹；按前两级临时共识排序的 Top 3 追加两个独立速度种子，总计每个 Top 3 有 3 条轨迹。不得以单条延长至 300 ns 替代独立重复。

- [ ] **Step 4: 实现轨迹完整性与证据指标**

完整性先检查时间范围、帧数、能量连续性和 PBC 处理，再计算：

- 蛋白 Cα RMSD、配体对齐口袋后的 RMSD；
- 口袋残基 RMSF；
- 蛋白–配体氢键占有率、疏水接触占有率、关键盐桥；
- 配体与口袋中心距离、接触原子数；
- MM/GBSA 仅作为可选解释性指标，不与实验自由能等同。

轨迹完成 < 95 ns 或后 80% 轨迹中配体持续离开口袋时，标记 `unstable`；仍保留轨迹与诊断，不删除候选。

- [ ] **Step 5: 测试续跑与提交**

Run:

~~~bash
.venv/bin/pytest tests/unit/test_gromacs.py tests/unit/test_md_analysis.py tests/integration/test_gromacs_resume.py -q
~~~

Expected: 夹具和命令测试通过；无 GROMACS 时真实二进制测试明确 skip。

~~~bash
git add data/forcefields/amberff19sb.ff/README.provenance.md src/airti_tf/simulation tests/fixtures/gromacs tests/unit/test_gromacs.py tests/unit/test_md_analysis.py tests/integration/test_gromacs_resume.py
git commit -m "feat: add resumable 100ns molecular dynamics stage"
~~~

### Task 13: 实现分阶段共识排序并正确处理缺失证据

**Files:**
- Create: `src/airti_tf/ranking/consensus.py`
- Create: `tests/unit/test_consensus.py`
- Create: `tests/integration/test_stage_routing.py`

- [ ] **Step 1: 写缺失值、风险惩罚和稳定性失败测试**

~~~python
def test_missing_md_is_not_converted_to_zero():
    ranked = rank_targets([target(md_score=None, md_status="failed")], stage="final")
    assert ranked[0].md_score is None
    assert "missing_md" in ranked[0].uncertainty_flags


def test_unsupported_target_is_excluded_but_counted():
    result = rank_targets(
        [target(status="ready"), target(status="unsupported", score=None)],
        stage="screen",
    )
    assert len(result.ranked) == 1
    assert result.coverage.unsupported == 1


def test_frozen_weights_reject_runtime_change(frozen_config):
    with pytest.raises(FrozenWeightError):
        frozen_config.final_weights.vina = 0.9
~~~

Run: `pytest tests/unit/test_consensus.py tests/integration/test_stage_routing.py -q`
Expected: FAIL。

- [ ] **Step 2: 实现阶段特异排序**

所有输入先在同一查询分子批次内转为 `[0, 1]` 百分位。回溯开发阶段的初始权重固定为：

~~~yaml
screen:
  vina: 0.65
  docking_consistency: 0.20
  structure_quality: 0.15
boltz:
  vina: 0.35
  boltz: 0.40
  docking_consistency: 0.10
  structure_quality: 0.15
final:
  vina: 0.25
  boltz: 0.30
  md: 0.30
  structure_quality: 0.15
risk_penalty_max: 0.15
~~~

风险项包括结构低置信度、Boltz 种子分歧、严重碰撞、配体大于 56 原子和技术缺失。某层缺失时不得用零填充：候选保留，`evidence_tier` 降级，并单列缺失原因。只有拥有完整 MD 的候选参与“完整三级证据”排序；其他候选进入“部分计算证据”附表。

- [ ] **Step 3: 固化确定性 tie-break**

同分时依次使用：证据层级、技术成功种子数、结构质量、目标 UniProt accession。路由输出必须包含 `rank_reason`、每项标准化分、原始分、风险项和权重版本。

- [ ] **Step 4: 测试与提交**

Run: `.venv/bin/pytest tests/unit/test_consensus.py tests/integration/test_stage_routing.py -q`
Expected: 所有测试通过，固定夹具的排名哈希不变。

~~~bash
git add src/airti_tf/ranking/consensus.py tests/unit/test_consensus.py tests/integration/test_stage_routing.py
git commit -m "feat: rank targets with stage-aware evidence handling"
~~~

### Task 14: 生成标准报告、证据卡和禁止性措辞检查

**Files:**
- Create: `src/airti_tf/reporting/render.py`
- Create: `templates/target_fishing_report.md.j2`
- Create: `tests/fixtures/report_context.json`
- Create: `tests/unit/test_report.py`
- Create: `tests/integration/test_report_traceability.py`

- [ ] **Step 1: 写报告完整性和措辞失败测试**

~~~python
def test_every_report_metric_has_artifact_provenance(report_context):
    report = render_report(report_context)
    assert find_untraced_metrics(report) == []


@pytest.mark.parametrize("phrase", ["确认靶点", "已证实直接结合", "真实结合概率"])
def test_prohibited_claims_block_release(report_context, phrase):
    report_context["conclusion"] = f"本研究{phrase} EGFR"
    with pytest.raises(ProhibitedClaimError):
        render_report(report_context, release=True)
~~~

Run: `pytest tests/unit/test_report.py tests/integration/test_report_traceability.py -q`
Expected: FAIL。

- [ ] **Step 2: 实现标准交付目录**

~~~text
final_report/
├── report.md
├── report_manifest.json
├── tables/
│   ├── proteome_coverage.csv
│   ├── top100_targets.csv
│   ├── boltz_top30.csv
│   ├── md_top10.csv
│   └── unsupported_targets.csv
├── evidence_cards/top01_to_top05/
├── structures/
├── pymol/
├── trajectories/
└── audit/
    ├── run_manifest.json
    ├── tool_versions.json
    ├── parameter_diff.json
    └── failure_catalog.csv
~~~

报告至少包含：项目输入、覆盖率、Top 100、Top 30 Boltz、Top 10 MD、Top 5 证据卡、结构/口袋质量、失败与 unsupported 清单、局限性、参数版本和下一步湿实验建议。

- [ ] **Step 3: 添加结论边界与追溯链接**

允许表述为“候选靶点”“结构支持的优先候选”“建议进一步验证”；正式发布模式禁止将纯计算结果写成“确认靶点”“直接结合已证实”或“真实结合概率”。固定免责声明不参与误报扫描，但其模板哈希必须锁定。每个表格数值包含 `task_id` 或 `artifact_id`，报告 manifest 验证链接存在且 SHA-256 一致。

- [ ] **Step 4: 测试与提交**

Run:

~~~bash
.venv/bin/pytest tests/unit/test_report.py tests/integration/test_report_traceability.py -q
~~~

Expected: 测试通过；任一缺失工件、失配哈希或禁止性结论均阻断 release。

~~~bash
git add src/airti_tf/reporting/render.py templates/target_fishing_report.md.j2 tests/fixtures/report_context.json tests/unit/test_report.py tests/integration/test_report_traceability.py
git commit -m "feat: render traceable target fishing reports"
~~~

### Task 15: 用 Nextflow DSL2 串联模块并验证断点续跑

**Files:**
- Create: `workflow/main.nf`
- Create: `workflow/nextflow.config`
- Create: `workflow/modules/target_library.nf`
- Create: `workflow/modules/ligand_prep.nf`
- Create: `workflow/modules/screen.nf`
- Create: `workflow/modules/refine.nf`
- Create: `workflow/modules/md.nf`
- Create: `workflow/modules/report.nf`
- Create: `tests/integration/test_nextflow_mock.py`

- [ ] **Step 1: 写模拟流程失败测试**

~~~python
def test_mock_workflow_reaches_report(tmp_path, run_nextflow):
    result = run_nextflow(
        "workflow/main.nf",
        profile="test",
        params={"queries": "tests/fixtures/ligands.smi", "outdir": tmp_path},
    )
    assert result.returncode == 0
    assert (tmp_path / "final_report/report.md").exists()
    assert (tmp_path / "job_status.sqlite").exists()
~~~

Run: `pytest tests/integration/test_nextflow_mock.py -q`
Expected: FAIL；Nextflow 不存在的本机应明确 skip，CI 容器内应先失败于流程文件缺失。

- [ ] **Step 2: 建立模块化 DAG**

`workflow/main.nf` 的数据流固定为：

~~~text
TARGET_LIBRARY ───────────────┐
LIGAND_PREP ─> SCREEN ─> CALIBRATE ─> BOLTZ ─> RANK_TOP30
                                           └────> SELECT_TOP10 ─> MD ─> FINAL_RANK ─> REPORT
~~~

每个 process 使用元组携带 `run_id`、`ligand_state_id`、`target_id`、`pocket_id` 与输入哈希。process 输出文件名可读，但缓存身份由值通道和工件哈希确定。

- [ ] **Step 3: 定义 profile 和资源标签**

`nextflow.config` 至少包含：

~~~groovy
profiles {
  test {
    params.mock_tools = true
    process.executor = 'local'
  }
  local {
    process.executor = 'local'
    docker.enabled = true
  }
  production {
    process.executor = 'local'
    docker.enabled = true
    process.errorStrategy = { task.exitStatus in [137, 143] ? 'retry' : 'terminate' }
    process.maxRetries = 1
  }
  hpc {
    process.executor = 'slurm'
    singularity.enabled = true
  }
}
~~~

使用 `withLabel` 分配 `cpu_small`、`cpu_screen`、`gpu_boltz`、`gpu_md`；MD 设置最长 walltime 和独占 GPU。镜像以 immutable digest 引用。

- [ ] **Step 4: 验证缓存续跑**

Run:

~~~bash
nextflow run workflow/main.nf -profile test --queries tests/fixtures/ligands.smi --outdir work/smoke
nextflow run workflow/main.nf -profile test -resume --queries tests/fixtures/ligands.smi --outdir work/smoke
~~~

Expected: 第一次完成；第二次所有输入未变的 process 显示 cached。修改一个 ligand 后，只重跑该 ligand 的下游任务。

- [ ] **Step 5: 提交**

~~~bash
git add workflow tests/integration/test_nextflow_mock.py
git commit -m "feat: orchestrate pipeline with resumable Nextflow DSL2"
~~~

### Task 16: 建立三层 benchmark、消融和发布判定

**Files:**
- Create: `configs/benchmark.yaml`
- Create: `src/airti_tf/benchmark.py`
- Create: `tests/fixtures/benchmark_predictions.csv`
- Create: `tests/unit/test_benchmark.py`
- Create: `tests/integration/test_release_gate.py`

- [ ] **Step 1: 写指标和数据泄漏失败测试**

~~~python
def test_success_at_k_uses_any_known_human_target():
    truth = {"drug1": {"P00533", "P04637"}}
    ranked = {"drug1": ["P11111", "P00533", "P22222"]}
    assert success_at_k(truth, ranked, k=2) == 1.0
    assert reciprocal_rank(truth, ranked) == pytest.approx(0.5)


def test_blind_set_cannot_be_used_to_fit_weights():
    with pytest.raises(DataLeakageError):
        fit_weights(dataset_role="blind")
~~~

Run: `pytest tests/unit/test_benchmark.py tests/integration/test_release_gate.py -q`
Expected: FAIL。

- [ ] **Step 2: 固化数据集角色和真值规则**

`configs/benchmark.yaml` 将分子明确分为：

- `smoke`：10 个，覆盖工具链和失败恢复；
- `retrospective_train` / `retrospective_validation`：合计约 100 个，按 Bemis–Murcko 骨架和靶点家族成组划分；
- `blind`：20–30 个，在所有权重、阈值和报告模板冻结后一次性解盲。

真值仅纳入人源、直接结合且有可追溯实验依据的蛋白；通路关联、表达变化和纯预测结果不得作为正例。多靶点药物按任一合格直接靶点命中计 Success@k，同时报告 target-level Recall@k。

- [ ] **Step 3: 实现指标与置信区间**

输出 Success@10/50/100、Recall@k、MRR、技术成功率、可计算蛋白比例、Top 20 Jaccard 稳定性、家族/结构来源/化学空间分层结果。使用按分子 bootstrap 10,000 次的 95% 置信区间；所有随机数种子固定并写入 benchmark manifest。

- [ ] **Step 4: 执行必做对照和消融**

固定比较：

1. QuickVina2 only；
2. QuickVina2 + Boltz-2；
3. 完整三级流程；
4. 去除结构质量权重；
5. 去除背景校准；
6. 去除家族多样化；
7. 随机排序与仅口袋可成药性排序；
8. AutoDock Vina 初筛对照。

权重仅在 retrospective_train 拟合，在 retrospective_validation 选择；随后写入 `frozen_weights.yaml` 并计算 SHA-256。盲测执行时拒绝未冻结权重。

- [ ] **Step 5: 实现发布门**

`release_decision.json` 只有同时满足以下条件才为 `pass`：

- 盲测 Success@100 ≥ 0.30；
- 技术任务成功率 ≥ 0.95；
- Boltz-2 加入后 Success@k 不低于 QuickVina2 基线；
- 至少三个主要靶点家族有成功案例；
- 相同输入重复运行的 Top 20 Jaccard 中位数 ≥ 0.70；
- 所有失败均有错误码和可重跑任务 ID；
- 报告完整性检查全部通过。

任一条件失败均输出 `fail` 与逐项证据，禁止用人工修改 JSON 放行；重新发布必须生成新的 benchmark 版本。

- [ ] **Step 6: 测试与提交**

Run:

~~~bash
.venv/bin/pytest tests/unit/test_benchmark.py tests/integration/test_release_gate.py -q
~~~

Expected: 测试通过；边界值 0.2999 失败，0.3000 通过。

~~~bash
git add configs/benchmark.yaml src/airti_tf/benchmark.py tests/fixtures/benchmark_predictions.csv tests/unit/test_benchmark.py tests/integration/test_release_gate.py
git commit -m "feat: benchmark target retrieval and enforce release gates"
~~~

### Task 17: 构建固定版本容器并做硬件 smoke 验证

**Files:**
- Create: `containers/base.Dockerfile`
- Create: `containers/screening.Dockerfile`
- Create: `containers/boltz2.Dockerfile`
- Create: `containers/gromacs.Dockerfile`
- Create: `tests/integration/test_container_contracts.py`
- Modify: `workflow/nextflow.config`

- [ ] **Step 1: 写镜像契约失败测试**

~~~python
@pytest.mark.parametrize(
    ("image", "command"),
    [
        ("screening", ["qvina2", "--help"]),
        ("boltz2", ["boltz", "--help"]),
        ("gromacs", ["gmx", "--version"]),
    ],
)
def test_image_exposes_required_command(container_runner, image, command):
    result = container_runner.run(image, command)
    assert result.returncode == 0
~~~

Run: `pytest tests/integration/test_container_contracts.py -q`
Expected: FAIL，因为镜像尚未构建。

- [ ] **Step 2: 构建分层镜像**

要求：

- 基础镜像固定 SHA-256 digest；
- 筛选镜像包含 fpocket、Meeko、QuickVina2 与 AutoDock Vina；
- Boltz-2 镜像固定模型和依赖版本，模型权重单独校验哈希；
- GROMACS 镜像启用 CUDA、PLUMED 不作为首版依赖；
- 镜像内 `airti-tf` 与宿主代码版本匹配；
- 生成 SPDX SBOM、许可证清单和镜像签名。

Run:

~~~bash
docker build -f containers/screening.Dockerfile -t airti-tf-screening:0.1.0 .
docker build -f containers/boltz2.Dockerfile -t airti-tf-boltz2:0.1.0 .
docker build -f containers/gromacs.Dockerfile -t airti-tf-gromacs:0.1.0 .
~~~

Expected: 三个镜像构建成功，并输出不可变 digest。

- [ ] **Step 3: 执行 CPU/GPU 工具 smoke**

CPU：

~~~bash
docker run --rm airti-tf-screening:0.1.0 qvina2 --help
~~~

GPU 节点：

~~~bash
docker run --rm --gpus all airti-tf-boltz2:0.1.0 boltz predict tests/fixtures/boltz/input.yaml --out_dir /tmp/boltz-smoke
docker run --rm --gpus all airti-tf-gromacs:0.1.0 gmx mdrun -s /opt/smoke/md.tpr -nsteps 1000 -deffnm /tmp/md-smoke
~~~

Expected: 命令成功；日志记录 GPU 型号、驱动、CUDA、峰值显存和运行时间。当前主机 GPU 不可通信，所以该步骤只能在满足 M0 的 GPU 节点执行。

- [ ] **Step 4: 固化 Docker digest，并保留 Apptainer 扩展验证**

Nextflow production profile 只引用已签名 OCI 镜像的不可变 digest。未来启用 HPC profile 时，再从相同 OCI digest 构建 SIF、保存 SIF SHA-256，并验证容器命令契约一致。

- [ ] **Step 5: 测试与提交**

Run: `.venv/bin/pytest tests/integration/test_container_contracts.py -q`
Expected: 所有容器命令契约通过。

~~~bash
git add containers workflow/nextflow.config tests/integration/test_container_contracts.py
git commit -m "build: package reproducible screening and GPU containers"
~~~

### Task 18: 完成 10 例端到端验收、运维文档和 v0.1.0 候选版

**Files:**
- Create: `scripts/run_smoke.sh`
- Create: `tests/integration/test_end_to_end_smoke.py`
- Modify: `README.md`
- Create: `docs/operations/runbook.md`
- Create: `docs/operations/failure-catalog.md`
- Create: `docs/operations/data-retention.md`
- Create: `docs/validation/v0.1.0-smoke-report.md`

- [ ] **Step 1: 写端到端验收测试**

验收断言：

~~~python
def test_smoke_delivery_is_complete(smoke_delivery):
    assert smoke_delivery.query_count == 10
    assert smoke_delivery.technical_success_rate >= 0.95
    assert smoke_delivery.report_manifest_valid
    assert smoke_delivery.all_metrics_traceable
    assert smoke_delivery.unsupported_targets_have_no_numeric_score
    assert smoke_delivery.resume_recomputed_task_count == 0
~~~

Run: `pytest tests/integration/test_end_to_end_smoke.py -q`
Expected: FAIL，直到完整流程和 10 例 smoke 产物存在。

- [ ] **Step 2: 编写可重复 smoke 命令**

`scripts/run_smoke.sh` 必须启用严格 shell 选项、先运行 production preflight、固定配置和镜像 digest，并调用：

~~~bash
nextflow run workflow/main.nf \
  -profile production \
  -params-file configs/smoke-frozen.yaml \
  --queries data/benchmark/smoke_v1.sdf \
  --outdir results/smoke_v0.1.0 \
  -with-report results/smoke_v0.1.0/nextflow-report.html \
  -with-trace results/smoke_v0.1.0/nextflow-trace.tsv
~~~

脚本在当前主机预期于 preflight 阶段安全退出；满足 M0 的生产节点才可继续。

- [ ] **Step 3: 进行故障注入和恢复验收**

在 smoke 运行中分别注入一次 QuickVina2 非零退出、Boltz OOM、MD 中断和报告工件哈希失配。断言前三类按限定策略恢复且有状态历史，哈希失配阻断报告发布。随后用 `-resume` 运行，已完成任务不得重算。

- [ ] **Step 4: 完成运维与边界文档**

`README.md` 和 runbook 必须覆盖：输入格式、1–5 分子限制、配置冻结、preflight、启动、恢复、错误码、人工复核、工件定位、版本升级和报告边界。首版明确不包含湿实验、客户 Web UI、多物种与一次 >5 个查询分子的生产支持。

- [ ] **Step 5: 生成 smoke 验证报告**

`docs/validation/v0.1.0-smoke-report.md` 记录 10 个分子、靶点库版本、可计算覆盖率、任务数、成功率、各阶段耗时、GPU/CPU 使用、失败注入结果、Top 20 稳定性和全部工件哈希。未运行的项目不得写为通过。

- [ ] **Step 6: 全量质量检查**

Run:

~~~bash
.venv/bin/ruff check src tests
.venv/bin/mypy src/airti_tf
.venv/bin/pytest -q --cov=airti_tf --cov-report=term-missing
nextflow config workflow/main.nf -profile production
airti-tf preflight --profile production
~~~

Expected: 代码检查和测试通过；生产节点 preflight 返回 0；核心 Python 模块分支覆盖率 ≥ 85%。

- [ ] **Step 7: 创建候选版本**

只有 Task 16 的 `release_decision.json` 为 `pass` 后才执行：

~~~bash
git add README.md scripts tests docs/operations docs/validation
git commit -m "release: validate reverse target fishing v0.1.0"
git tag -a v0.1.0-rc1 -m "Human proteome reverse target fishing v0.1.0 release candidate"
~~~

Expected: 标签指向包含锁定配置、容器 digest、参考库版本、benchmark 判定和 smoke 报告的提交。

## 3. 实施顺序与并行边界

严格依赖链为：

~~~text
Task 1 → 2 → 3 → 4
             ├→ 5 → 6 → 7 → 8 → 9 → 10
             └──────────────────────→ 11 → 12 → 13 → 14
Task 10 + 11 + 12 + 14 → 15 → 16 → 17 → 18
~~~

可安全并行的工作仅限接口已由 Task 2 固定之后：

- Task 5–7 的靶点库构建与 Task 8 的配体准备；
- Task 11 的 Boltz 适配器与 Task 12 的 GROMACS 适配器；
- Task 14 的报告模板与 Task 15 的 Nextflow mock 编排；
- 容器构建与纯 Python benchmark 指标测试。

任何并行分支不得自行改变 contracts、锁定配置、候选规模或发布门；需要改变时先合并契约变更并重新运行所有下游测试。

## 4. 关键验收清单

- [ ] 仅人源 taxonomy 9606 canonical 蛋白进入主索引；
- [ ] 全部蛋白都有 `ready`、`unsupported` 或 `failed` 状态，缺失不得伪装为零分；
- [ ] 查询批次严格限制 1–5 个小分子；
- [ ] QuickVina2 固定三种子，背景探针固定 100 个/口袋；
- [ ] 一级 Top 300 全部进入 Boltz-2，精评后保留 Top 30；
- [ ] Top 10 各 100 ns，Top 3 各追加两条独立重复；
- [ ] 最终交付 Top 5 证据卡、Top 100 表、覆盖率和不可判定清单；
- [ ] 10/约100/20–30 三层 benchmark 无数据泄漏；
- [ ] Success@100 ≥ 30% 与技术成功率 ≥ 95% 只作为内部发布门；
- [ ] 报告不把计算优先级解释为实验确认或真实结合概率；
- [ ] 所有数值、结构、轨迹、参数和结论均可回溯；
- [ ] 未通过生产 preflight 时不能提交全蛋白组任务。

## 5. 官方接口依据

- Boltz-2 输入、亲和力属性、口袋约束、MSA 与预测选项：[Boltz prediction documentation](https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md)
- fpocket 命令行与官方容器：[fpocket repository](https://github.com/Discngine/fpocket)
- Meeko 受体 PDBQT 和 Vina 网格准备：[Meeko receptor preparation](https://rwxmeeko.readthedocs.io/en/latest/cli_rec_prep.html)
- Nextflow process、容器、重试与资源指令：[Nextflow process reference](https://github.com/nextflow-io/nextflow/blob/master/docs/reference/process.md)
- GROMACS 运行、checkpoint 与 GPU 行为：[GROMACS user guide](https://manual.gromacs.org/documentation/2025-current/user-guide/index.html)

上述接口在实现时必须再次锁定精确工具版本和镜像 digest；官方文档用于接口契约，项目内测试用于防止版本漂移。
