# AIRTI 反向钓靶现有计算平台审计

**日期**：2026-08-12  
**结论**：现有蛋白质/多肽设计平台的 GPU、Docker、CPU、内存和大容量数据盘可直接复用。反向钓靶特有的软件栈已于 2026-08-13 封装为 AIRTI 统一镜像，并通过三条真实工具链 smoke；全人蛋白组参考库与检索基准仍需完成后方可进入生产任务。

## 1. 已验证基础设施

| 项目 | 实测结果 | AIRTI 用途 |
|---|---|---|
| GPU | NVIDIA GeForce RTX 4090；容器可见 48507 MiB；Torch CUDA 12.8 | Boltz-2 与 GROMACS GPU 任务 |
| GPU 容器 | Docker NVIDIA runtime 已注册 | 复用既有容器运行方式 |
| JAX GPU | 既有 AF3 smoke 返回 `devices [CudaDevice(id=0)]` | 证明容器到 GPU 的计算路径可用 |
| CPU | Intel Xeon Platinum 8347C，72 logical CPUs | 建库、口袋发现、QuickVina2 并行筛选 |
| 内存 | 503 GiB，总可用约 481 GiB | 全蛋白组清单、任务编排和批量筛选 |
| 大数据盘 | `/data` 可用约 12 TiB | 受体库、任务工件、轨迹与归档 |
| SSD | `/mnt/ssd4t` 可用约 2.5 TiB | 模型、数据库、缓存和高 I/O 临时文件 |
| 容器编排基础 | Docker Compose v5.1.4 | 延续既有平台的服务与 smoke-test 习惯 |

验证命令包括：

```bash
nvidia-smi
docker info
/data/protein-design/scripts/smoke-test.sh compose
/data/protein-design/scripts/smoke-test.sh host-gpu
/data/protein-design/examples/af3/run-check-or-full.sh
```

AF3 检查只执行环境和模型挂载 smoke，没有启动完整 AF3 推理。

## 2. 可复用的既有蛋白设计工作台

平台根目录为 `/data/protein-design`，版本标记为 `v0.5.0-dev`。已存在：

- `pd-af3-gpu:v3.0.2`；
- `pd-af2multimer-gpu:fixed`；
- `pd-bindcraft-gpu:0.19.1`；
- `pd-foundry-gpu:latest`；
- `pd-pepmimic-gpu:latest`；
- `pd-rfpeptide-gpu:fixed`；
- `pd-rosetta-cpu-parallel:latest`；
- 多方法 benchmark 与 PyRosetta 镜像；
- AF3 模型、627 GiB 公共数据库和 JAX 缓存；
- Docker Compose profile、smoke-test、输入/输出挂载和镜像外模型管理模式。

AIRTI 应复用这些工程约定，而不是修改或复制既有模型资产：

1. 轻量代码进入 AIRTI Git 仓库；
2. 模型、数据库和大结果留在仓库外；
3. AIRTI 科学工具使用一个固定版本统一镜像，模型与大数据外置；
4. 输入、输出、模型和缓存使用显式只读/读写挂载；
5. 统一镜像提供命令合同与三引擎真实 smoke，再开放完整计算开关；
6. 运行日志、参数、镜像 digest 和结果哈希进入审计清单。

## 3. AIRTI 统一镜像

宿主和既有设计镜像不直接安装反向钓靶工具。AIRTI 使用 `airti-tf:0.1.0-gpu` 统一提供：

- Boltz 2.2.1 与 Torch 2.7.1+cu128；
- GROMACS 2025.4 CUDA 构建、AmberTools 24.8 和 ParmEd 4.3.1；
- QuickVina2 2.1.0 与 AutoDock Vina 1.2.7；
- fpocket 4.2.3 与 Meeko 0.7.1；
- AIRTI Python 包 0.1.0。

选择单镜像是针对当前单节点、单 GPU 和串行 GPU 阶段做出的部署决策。其优势是 Nextflow 各阶段共享同一可审计环境，减少镜像间版本漂移；代价是本地镜像约 16.8 GB。Boltz 模型与 CCD 缓存约 7.6 GB，保留在 `/mnt/ssd4t/airti-target-fishing/boltz`，不随镜像复制。

最终本地内容标识、模型哈希与 SPDX SBOM 分别记录在 `containers/images.lock.yaml`、`containers/models.lock.yaml` 和 `docs/sbom/`。本地内容标识不是注册表 RepoDigest；推送 OCI 注册表后仍需记录远端不可变 digest。

## 4. 修订后的部署判断

首版不依赖尚未存在的 SLURM/Apptainer 集群，采用单机 Docker/NVIDIA 生产 profile：

- CPU 阶段最多使用 64 threads，预留系统资源；
- GPU 阶段默认并发为 1，Boltz-2 与 GROMACS 不同时占用 GPU；
- 大工件根目录设为 `/data/airti-target-fishing`；
- 模型与高 I/O 缓存根目录设为 `/mnt/ssd4t/airti-target-fishing`；
- 源码、测试和轻量清单保留在 `/home/a/Data/AIRTI`；
- 未来增加 GPU 节点时再启用 SLURM/Apptainer profile，不影响数据契约。

当前环境已通过 GPU、Docker、CPU、内存、容量、统一镜像命令合同和三引擎真实 smoke。启动全人蛋白组生产任务前仍需通过参考库完整性、背景校准库、10 例端到端检索和发布阈值门；现阶段结果不得表述为已完成全蛋白组验证。
