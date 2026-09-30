# 项目说明：SOFA 实时针–组织交互仿真

目标：用 SOFA v25.12 做实时的针–组织交互仿真，先做刚性针，再做多层组织和柔性针；最终服务于术前规划和闭环机器人控制。导师的原始要求：`~/sofa/papers/original/Manish对于项目的原始总体要求.docx`。

## 知识库（遇到问题先查这里）
| 内容 | 位置 |
|---|---|
| 文献笔记索引（20 篇） | `docs/literature/README.md` |
| 方法手册：统一符号的推导，§12 是理论和 SOFA 实现的逐项对照 | `docs/literature/handbook.md` |
| 文献综合分析、矛盾点、待和导师确认的问题 | `docs/literature/00_synthesis.md` |
| 插件逐行分析和已知问题（CollisionAlgorithm / ConstraintGeometry） | `docs/code/01_needle_plugins_v25.12.md` |
| SOFA 本体机制（时间步、λ 的单位、求解器、FEM、SofaPython3） | `docs/code/02_sofa_core_v25.12.md` |
| 其他插件（Cosserat、SoftRobots.Inverse、MOR、BeamAdapter） | `docs/code/03_other_plugins.md` |
| 验证实验汇总，以及可复现的脚本和原始数据 | `docs/code/04_experiments_summary.md`，`docs/code/experiments/` |
| 资源清单、学习计划 | `docs/resources/` |

## 查证规则
1. 回答问题时先读知识库中的相关笔记，**不要凭记忆**。
2. 笔记不够、有疑问或者要引用具体数值时，**回到原文或源码核对**：
   - 论文 PDF：`~/sofa/papers/original/`、`~/sofa/papers/extra/`（每篇笔记写了文件名和页码对应关系）；
   - 源码：`~/sofa/resources/`（各插件、SOFA v25.12 的部分源码），`~/sofa/reference_sources/`（我们编译用的 v25.12 插件）。代码笔记里的结论都标了 `文件:行号`。
3. 区分依据类型：〔原文〕〔代码〕〔运行〕〔推导〕〔推论〕；发现笔记有误时，改正并写明改正记录。

## 本机环境
- SOFA 二进制包：`~/sofa/SOFA_v25.12.00_Linux`。用 Python 驱动时：
  ```bash
  export SOFA_ROOT=~/sofa/SOFA_v25.12.00_Linux
  export PYTHONPATH=$SOFA_ROOT/plugins/SofaPython3/lib/python3/site-packages:$SOFA_ROOT/plugins/STLIB/lib/python3/site-packages
  export LD_LIBRARY_PATH=$SOFA_ROOT/lib
  ```
- 编译好的插件：`~/sofa/reference_install`（v25.12 版 CollisionAlgorithm / ConstraintGeometry）、`~/sofa/resources/install_master`（上游 master 版）、`~/sofa/resources/install_cosserat`。
- 编译带 Python 绑定的插件：用本地的 pybind11 2.12（`~/sofa/resources/pybind11-install`），加 `-DPYBIND11_BUILD_ABI=\"_cxxabi1016\"`（SOFA 二进制包是 GCC 11 编译的），必要时 `-DCMAKE_DISABLE_FIND_PACKAGE_SoftRobots=TRUE`。

## 注意事项
- `~/sofa/reference_*`、`~/sofa/diagnostic_*`、`~/sofa/resources/` 下的克隆：**只读，不要修改**。实验脚本要设 `sys.dont_write_bytecode = True`，避免往这些目录里写缓存。
- `NeedleInsertion_*.json`、`diagnostics/`：只保留在本地，**不提交、不删除**（已写进 `.git/info/exclude`）。
- 提交 git 前要征得同意；提交时只包含相关文件。
- 当前阶段：资源学习已完成（2026-09-30），**项目的实施要等用户批准后再开始**。
