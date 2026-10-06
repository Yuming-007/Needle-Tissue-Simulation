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
- 本项目的修补插件：源码 `plugins/NeedleSimFixes/`（FTC 刚度矩阵 bug 的修复），编译到 `~/sofa/needle_build/NeedleSimFixes/`（`cmake <src> -DCMAKE_PREFIX_PATH=$SOFA_ROOT/lib/cmake -DCMAKE_BUILD_TYPE=Release && make`）；`scenes/common.py` 自动加载。
- 实时配置（0d 冻结；**第 1–4 步的软实时基线，不是闭环控制的硬实时保证**）：组织用 `AsyncSparseLDLSolver` + `LinearSolverConstraintCorrection`，dt = 0.01 s，节点预算约 870（最重负载下）；验证实验改用 `SparseLDLSolver` 离线跑（确定性）；R3 和 R1 的差别中位数很小，但接触中偶发单步可达约 2.3%。图形界面：`RUNSOFA_OPTS="-g glfw" tools/run_gui.sh <场景>`（真实速度）；ImGui 界面在核显上每帧多约 12 ms，只用于查看。
- gmsh 4.15.2（官方 Linux64 SDK，用户 2026-10-03 同意安装；本机没有 pip）：`~/sofa/tools/gmsh-4.15.2-Linux64-sdk`，Python 接口在其 `lib/`（`gmsh.py` + `libgmsh.so`），使用时加进 `PYTHONPATH`。删除该目录即卸载。
- 编译带 Python 绑定的插件：用本地的 pybind11 2.12（`~/sofa/resources/pybind11-install`），加 `-DPYBIND11_BUILD_ABI=\"_cxxabi1016\"`（SOFA 二进制包是 GCC 11 编译的），必要时 `-DCMAKE_DISABLE_FIND_PACKAGE_SoftRobots=TRUE`。

## 注意事项
- `~/sofa/reference_*`、`~/sofa/diagnostic_*`、`~/sofa/resources/` 下的克隆：**只读，不要修改**。实验脚本要设 `sys.dont_write_bytecode = True`，避免往这些目录里写缓存。
- `NeedleInsertion_*.json`、`diagnostics/`：只保留在本地，**不提交、不删除**（已写进 `.git/info/exclude`）。
- 提交 git 前要征得同意；提交时只包含相关文件。
- 协作方式（用户 2026-09-30 确定）：Claude 是主要执行者和技术判断者；ChatGPT 作为第二技术视角，直接读 GitHub 仓库做关键节点检查。GitHub 仓库是双方共同依据的项目事实。每个阶段或一组有意义的结果完成后，先给用户一份简短总结（关键结论、实验结果、修改内容），用户确认后再 commit + push。文档要写到外部读者不看对话也能看懂（结果、脚本、数据路径）。
- 当前阶段：第 0 步、第 1 步已完成；第 1.5 步进行中。针尖表示（V0–V4）已完成并合并到 main：平底圆盘 + 分布点（插件多点接触）作为机制开发替身，冻结的是接口（针尖合力 + 合力矩 → 刺穿判据），a_eff ≠ R_phys（`docs/plan/step1_5_tip_survey.md` §5，结果 `docs/plan/step1_5_results.md`）。**非均匀工作网格**（方案 `docs/plan/step1_5_mesh_design.md`）：M0 生成（`scenes/mesh/make_graded_mesh.py`，网格在 `meshes/step1_5/`）、M1 导入、M2 计时、M3 筛选已完成（工作分支 `work/step1_5_mesh`，本地提交，未推送）。**当前停在**：推荐 A+（h_local 5 mm、加密区 R = D = 20 mm、远场 12 mm、971 节点，无界面约 75 帧/秒、d_0.1 = 10.75 mm），**等用户在 GUI 实测流畅度**（`RUNSOFA_OPTS="-g glfw" tools/run_gui.sh scenes/step1_5/m3_gui_scene.py A+`，用户把终端 `[流畅度]` 行贴回）；之后冻结 A+，在 A+ 上完成剩余 M3（① 侧向漂移、⑤ 包络已有筛选数据、⑥ 流畅度、⑦ R3 vs R1；② 已完成），总结 → 复核分支 → ChatGPT 复核 → 合并，然后进入第 2 步（先用单点调通刺穿机制，再切换圆盘；人工阈值仅用于验证）。粗规则网格（swapping=True）只作机制测试台。每个机制经用户确认后才进入下一个。
- 用户决定（2026-09-30）：速度 = 动画不卡顿、接近真实穿刺速度（每步计算时间 ≤ dt）。**2026-10-03 修改**：只要求 SOFA 动画播放流畅不卡，**不需要按真实速度播放**（可以是慢动作）；网格按精度选择，dt = 0.01、R3 不变；3D 组织；先刚性针；组织用共旋线弹性；参数先用文献值；"破裂"定义和柔性针模型暂缓。
- 用户决定（2026-09-30）：组织力场 = `FastTetrahedralCorotationalForceFieldFixed`（`plugins/NeedleSimFixes`），`method="polar"`；第 1–4 步机制验证用 ν = 0.45，这是**验证用参数，不是最终组织参数**，最终参数和接近不可压缩问题留到物理验证 / 参数辨识阶段。
- 原则：在能实现目标机制的前提下尽量复用已有资源，但不为了靠近某个示例而偏离自己的路线。
