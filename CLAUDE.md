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
- 实时配置（0d 冻结；**第 1–4 步的软实时基线，不是闭环控制的硬实时保证**）：组织用 `AsyncSparseLDLSolver` + `LinearSolverConstraintCorrection`，dt = 0.01 s（0d 的"节点预算约 870"是旧的严格实时判据下的结论，2026-10-03 起只要求动画流畅，见下）；验证实验改用 `SparseLDLSolver` 离线跑（确定性）；R3 和 R1 的差别中位数很小，但接触中偶发单步可达约 2.3%。图形界面：`RUNSOFA_OPTS="-g glfw" tools/run_gui.sh <场景> [场景参数...]`；ImGui 界面在核显上每帧多约 12 ms，只用于查看。
- gmsh 4.15.2（官方 Linux64 SDK，用户 2026-10-03 同意安装；本机没有 pip）：`~/sofa/tools/gmsh-4.15.2-Linux64-sdk`，Python 接口在其 `lib/`（`gmsh.py` + `libgmsh.so`），使用时加进 `PYTHONPATH`。删除该目录即卸载。
- 编译带 Python 绑定的插件：用本地的 pybind11 2.12（`~/sofa/resources/pybind11-install`），加 `-DPYBIND11_BUILD_ABI=\"_cxxabi1016\"`（SOFA 二进制包是 GCC 11 编译的），必要时 `-DCMAKE_DISABLE_FIND_PACKAGE_SoftRobots=TRUE`。

## 注意事项
- `~/sofa/reference_*`、`~/sofa/diagnostic_*`、`~/sofa/resources/` 下的克隆：**只读，不要修改**。实验脚本要设 `sys.dont_write_bytecode = True`，避免往这些目录里写缓存。
- `NeedleInsertion_*.json`、`diagnostics/`：只保留在本地，**不提交、不删除**（已写进 `.git/info/exclude`）。
- 提交 git 前要征得同意；提交时只包含相关文件。
- **工作方式（用户 2026-10-09 调整）**：目标是"机制正确、表现合理"的模型（课程项目尺度），不追求逐步的数值收敛验证。每个机制：实现 → 演示场景（用户 GUI 目测）→ 一条力–深度 / 力–时间曲线和文献（如 Okamura）定性对比 → 基本合理性检查（不崩溃、方向和量级、作用力 = 反作用力）→ 简短总结 → 提交。网格收敛、精确解对照、参数扫描、多轮敏感性测试只在用户要求时做；ChatGPT 复核只在大的节点或用户要求时做；优先用插件 / SOFA 的默认做法；已知局限写进报告。
- 协作方式（用户 2026-09-30 确定）：Claude 是主要执行者和技术判断者；ChatGPT 作为第二技术视角，直接读 GitHub 仓库做关键节点检查。GitHub 仓库是双方共同依据的项目事实。每个机制完成后，先给用户一份简短总结（关键结论、实验结果、修改内容），用户确认后再 commit + push。文档要写到外部读者不看对话也能看懂（结果、脚本、数据路径）。
- 当前阶段：第 0、1、1.5 步已完成。工作网格 **W50**（gmsh 非均匀网格，`meshes/step1_5/W50.npz`：h_local 5 mm、加密区沿竖直针道到 50 mm 深、半径 20 mm、远场 12 mm、1242 节点；`contact.add_tissue_with_surface(mesh=common.load_mesh("W50"))`），针**竖直对准顶面中心 (40, 40) mm**，最大插入深度 50 mm。针尖**默认用插件原生的单点**，平底圆盘（`contact.tip_patch_points` / `add_tip_patch`）保留为可选。已知局限：非结构网格不对称，tenting 时侧向力约为轴向的 3%（5 mm），第 8 步柔性针之前必须检查。下一步：**第 2 步刺穿机制**（插件原生刺穿状态机；验证用阈值让刺穿在约 5–6 mm 压深内发生）。结果和决策记录：`docs/plan/step1_5_results.md`；方案 `docs/plan/00_start_plan.md`。
- 用户决定（2026-09-30）：速度 = 动画不卡顿、接近真实穿刺速度（每步计算时间 ≤ dt）。**2026-10-03 修改**：只要求 SOFA 动画播放流畅不卡，**不需要按真实速度播放**（可以是慢动作）；网格按精度选择，dt = 0.01、R3 不变；3D 组织；先刚性针；组织用共旋线弹性；参数先用文献值；"破裂"定义和柔性针模型暂缓。
- 用户决定（2026-09-30）：组织力场 = `FastTetrahedralCorotationalForceFieldFixed`（`plugins/NeedleSimFixes`），`method="polar"`；第 1–4 步机制验证用 ν = 0.45，这是**验证用参数，不是最终组织参数**，最终参数和接近不可压缩问题留到物理验证 / 参数辨识阶段。
- 原则：在能实现目标机制的前提下尽量复用已有资源，但不为了靠近某个示例而偏离自己的路线。
