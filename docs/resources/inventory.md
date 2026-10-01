# 资源清单（阶段 1：资源普查）

- 日期：2026-09-30
- 本地位置：只读克隆放在 `~/sofa/resources/`（不在项目仓库里）；之前编译好的插件在 `~/sofa/reference_*`，诊断版本在 `~/sofa/diagnostic_*`（都没有改动）。
- 标记：**〔已确认〕** 表示看过代码、仓库或运行确认；**〔待核实〕** 表示只看了描述或文件名，要在阶段 2 读代码确认。
- 相关度：⭐⭐⭐ 直接用于起点；⭐⭐ 某个阶段会用到；⭐ 参考。

---

## 0. 已有的工作（另一个会话里完成的）〔已确认〕
- `~/sofa/reference_sources/{CollisionAlgorithm,ConstraintGeometry}`：InfinyTech3D fork 的 **`v25.12.00`** 版，已经用本机 SOFA v25.12 **在源码树外编译成功**，安装在 `~/sofa/reference_install/`。
- `~/sofa/diagnostic_sources/CollisionAlgorithm`：同一版本，只在 `InsertionAlgorithm::puncturePhase()` 里加了 12 行 `[PUNCTURE_DIAG]` 日志（打印 λ 原值、λ/dt、阈值和判定结果）。
- `needle_project/diagnostics/NeedleInsertion_puncture_diag.py`：基于官方 `NeedleInsertion.py` 的诊断场景。`NeedleInsertion_*.json` 是 SOFA 计时器输出的每步耗时。
- **当时的结论没有记录** → 阶段 2 重新核实。

## 1. 本机 SOFA v25.12（二进制版，commit 4d6c2b8e）

### 1.1 使用环境〔已确认〕
- 用 Python 驱动 SOFA 需要设置：
  ```bash
  export SOFA_ROOT=~/sofa/SOFA_v25.12.00_Linux
  export PYTHONPATH=$SOFA_ROOT/plugins/SofaPython3/lib/python3/site-packages
  export LD_LIBRARY_PATH=$SOFA_ROOT/lib:$LD_LIBRARY_PATH
  ```
  之后系统的 Python 3.12 就能 `import Sofa`（之前 `import Sofa` 失败，就是没设这些变量）。
- 已经导出全部 549 个注册组件：`~/sofa/resources/sofa_v25.12_components.tsv`。

### 1.2 和项目相关的组件〔已确认存在；具体行为待核实〕
| 类别 | 组件 | 用途 | 相关度 |
|---|---|---|---|
| 动画循环 | `FreeMotionAnimationLoop`、`ConstraintAnimationLoop`、`MultiStepAnimationLoop` | 约束法三步流程（handbook §1.4） | ⭐⭐⭐ |
| 约束求解器 | `BlockGaussSeidelConstraintSolver`、`UnbuiltGaussSeidelConstraintSolver`、`NNCGConstraintSolver`、`ImprovedJacobiConstraintSolver`、`LCPConstraintSolver`（v25.12 没有 `ProjectedGaussSeidelConstraintSolver`，已核对组件列表） | 各种 GS 变体；NNCG 是非光滑的非线性共轭梯度 | ⭐⭐⭐ |
| 柔度计算 | `LinearSolverConstraintCorrection`、`PrecomputedConstraintCorrection`、`UncoupledConstraintCorrection`、`GenericConstraintCorrection` | handbook §6.4 的几种 W 计算方式 | ⭐⭐⭐ |
| 线性求解器 | `BTDLinearSolver`（块三对角，针）、`PrecomputedLinearSolver`、`PrecomputedWarpPreconditioner`、Eigen 系列 | | ⭐⭐ |
| 时间积分 | `EulerImplicitSolver`（可设 firstOrder）、`StaticSolver` | 动力学格式和准静态格式 | ⭐⭐⭐ |
| 组织 FEM | `TetrahedronFEMForceField`、`FastTetrahedralCorotationalForceField`、`TetrahedralCorotationalFEMForceField`、`TetrahedronHyperelasticityFEMForceField`、`StandardTetrahedralFEMForceField`（超弹）、`MJEDTetrahedralForceField`、`ParallelTetrahedronFEMForceField`（多线程） | 线弹性共旋 / 超弹 | ⭐⭐⭐ |
| 粘弹 | `TetrahedronViscoelasticityFEMForceField`、`TetrahedronViscoHyperelasticityFEMForceField`（SofaViscoElastic） | 可以描述 Okamura 实验里停针后的松弛 | ⭐ |
| 针 | `BeamFEMForceField`（Timoshenko，共旋）、`BeamLinearMapping` | 刚性针和柔性针 | ⭐⭐⭐ |
| Lagrange 约束 | `BilateralLagrangianConstraint`、`UnilateralLagrangianConstraint`、`SlidingLagrangianConstraint`、`StopperLagrangianConstraint`、`AugmentedLagrangianConstraint`、`FixedLagrangianConstraint` | SOFA 自带的约束类型；`SlidingLagrangianConstraint` 可能能用来实现针身约束 | ⭐⭐ |
| 其他 | `NeedleTracker`（InfinyToolkit，判断针尖在哪个网格里）、`CarvingManager`（SofaCarving）、`TearingEngine`（Tearing） | 多层组织：判断针尖所在的层；拓扑切割 | ⭐ |

### 1.3 已安装、之前没注意到的插件〔已确认存在〕
- **SoftRobots.Inverse**：`QPInverseProblemSolver` + `PositionEffector` / `PositionEquality` / `SlidingActuator` / `ForcePointActuator` 等。这是**基于 QP 的逆向 FE 控制框架**（Coevoet 2019），可以把"针尖到达靶点"写成 effector，把"针根运动"写成 actuator → **可能直接服务于闭环控制目标**。和 Adagolodjo 的数值 Jacobian 路线是两种不同的做法。⭐⭐
- **ModelOrderReduction**：超降阶力场（包括四面体超弹）、`ModelOrderReductionMapping`、`MORUnilateralInteractionConstraint`。⭐⭐（后期提速）
- **BeamAdapter**：Kirchhoff 杆 + `AdaptiveBeamSlidingConstraint`（梁在约束点上滑动，本来是给导管、导丝用的）。⭐（柔性针的备选）

## 2. 针插入插件（核心候选）

### 2.1 InfinyTech3D/CollisionAlgorithm + ConstraintGeometry ⭐⭐⭐〔已确认〕
- GPL-3.0；**只有这个 fork 包含针插入算法**（`InsertionAlgorithm`、`ConstraintInsertion`、`InsertionResolution`）。ICube 原版里没有（grep 确认）。
- 有 `v25.12.00` tag，已经编译过；master 在 v25.12 之后改了几处，**可能很关键**：

  | 提交 | 内容 | 影响 |
  |---|---|---|
  | CA `17e4fd2`（2026-05） | 拔出检测改为检查针尖速度方向；修正 `prunePointsUsingEdges` 返回值写反的问题 | v25.12 版插入过程中也可能误删约束点；返回值写反但调用处不用，**不影响结果**（阶段 2 已确认，见 `docs/code/01`） |
  | CA `e532b38`（2026-05） | 适配场景里的 `contactDistance` | |
  | CG `b575d9e`（2026-05） | `ConstraintUnilateral` 加了 `contactDistance` | |
  | CG `4c62532`（2026-07） | `InsertionResolution`：`invertMatrix` 在 \|det\| < 机器精度时**会放弃求逆，而且不报错**，只留下空矩阵；新增 `scaleComplianceMatrix` 选项，先缩放再求逆 | ⚠️ **SI 单位的软组织，柔度块的行列式本来就很小，v25.12 版的插入约束可能会悄悄失效**（阶段 2 要验证） |
- ⚠️ **代码确认的潜在问题**：`ConstraintGeometry/src/ConstraintGeometry/constraint/UnilateralResolution.h:61-90` 的 `UnilateralFrictionResolution::resolution()` 仍然是 ICube 在 2025-12-18 修复前的写法（修复者是 Claire Martin，提交 `3e83cc4`，标题"FIX bug in Gauss Seidel constraint resolution"）：力被截断之后，传给下一个分量的修正量没有更新；切向摩擦是在耦合传递之后才截断的。它影响**刺穿前针尖和表面的摩擦接触**。实际影响有多大，阶段 2 分析。
- 文档：`doc/doc.md`、`doc/diagrams/`（阶段 2 读）。

### 2.2 ICube 原版（forge.icube.unistra.fr/sofa）⭐⭐〔已确认〕
- `sofa/CollisionAlgorithm`、`sofa/ConstraintGeometry`：公开，维护者是 H. Courtecuisse，2026-09 还在更新。和 InfinyTech3D fork 从 2025 年起分开发展（CollisionAlgorithm：ICube 独有 35 个提交，fork 独有 58 个）。
- ICube 版的方向：代码重构（"new_design"）、**Hermite 边几何（曲线针）**、修复 GS bug、迁移到 Assist 框架。**不包含针插入算法**。
- 用处：**可以借鉴修过的 `UnilateralResolution` / `BilateralResolution`，以及 Hermite 曲线针的几何**。
- `sofa/msofaplugin`：矩阵组装的底层设施（IncomingSparseMatrix、UnbuiltMatrix、PCG 等），2021 年起公开。⭐
- `sofa/sofaros2plugin`：SOFA 和 ROS2 的接口（需要 ROS2 环境）。⭐⭐（以后接机器人）

### 2.3 SofaDefrost/Cosserat ⭐⭐〔已确认〕
- 有 `release-v25.12` tag，需要自己编译。
- 针插入示例（`examples/python3/NeedleInsertion.py` + `cosserat/needle/needleController.py`，Adagolodjo 2021 写的）：用 **Python 控制器**实现穿刺，读 `GenericConstraintSolver` 的约束力，超过阈值 3 就关掉碰撞，开始添加约束点；**插入方向写死为 x 轴**；单位是 cm；组织 E = 100。**属于原型级别**，但它是 Manish 提到的 "Cosserat" 方案的现成起点。

## 3. 文献对应代码的公开情况〔已确认〕
| 工作 | 代码 | 备注 |
|---|---|---|
| Martin 2023/2025（IsoDOF、混合求解器、CTSVD） | ❌ 没有找到 | Assist 框架已公开（ICube GitLab 的 `assist/*`，27 个模块，2026-09 很活跃），但 `assist_sofa` 里只有 IsoDOF 的性能统计数据文件，**没有针插入组件** |
| Adagolodjo / Baksic / Ha（逆向 FE 控制） | ❌ 没有公开的仓库 | 部分思想可以用 SoftRobots.Inverse 实现；Cosserat 插件里有 Adagolodjo 的针示例 |
| Vanneste 2024（分区降阶） | ❌ 公开的 ModelOrderReduction 仓库里**没有**分区相关的代码（grep 确认） | 只有标准的 MOR 和超降阶 |
| Wang 2024/2025（JHU 2D 仿真器） | ❌ arXiv 页面没有代码链接 | 同一个实验室，可以问 Manish |
| Baksic 2022 博士论文 | ✅ 可以下载 PDF | HAL tel-03881361（**还没读**；可能有 SOFA 实现细节） |
| Assist | ✅ 公开 | 文档在 assist.cnrs.fr；需要 SOFA master；不是我们的起点 |

## 4. SofaPython3（v25.12 分支）〔已确认〕
- **没有内置的状态快照或序列化功能**。要实现 CE 或 MPC 所需的"保存和恢复状态"，只能自己做：①读写所有 MechanicalObject 的 position、velocity 等 Data；② CollisionAlgorithm 内部的约束点列表（`m_couplingPts`）**是 C++ 私有成员，Python 访问不到** → 这是规划接口的一个**已知难点**（阶段 2 细看）。
- 有用的官方示例：`access_compliance_matrix.py`、`access_constraint_matrix.py`、`access_contact_forces.py`、**`access_energy.py`（应变能，对应 Wang 2025 的代价函数项）**、`advanced_timer.py`、`RPYC/`（远程服务，可能用于多进程并行仿真）、`jax/`。

## 5. 数据
| 数据 | 用途 | 获取方式 |
|---|---|---|
| JHU 多层 plastisol 仿体（Wang 2024 的配方、参数和 CT 重建的针形状；Wang 2025 的双层仿体参数） | "两种组织的立方体"的参数和物理验证 | **要问 Manish** |
| 3D-IRCADb-01/02（IRCAD） | 针对具体患者的肝脏模型（Vanneste、Ha 都用过） | 公开，需要注册 |
| 文献里的力曲线（Okamura、DiMaio 等） | 定性验证 | 已经整理在文献笔记里 |

## 6. 社区资源〔待核实，还没读〕
- SOFA 论坛："Needle insertion simulation problem"、"[SOLVED] Needle Insertion into Soft Tissue"。
- SOFA 插件页面："Cosserat beam: cable, needle"。
- SOFA GitHub Discussions（CollisionAlgorithm README 推荐的求助渠道）。

---

## 7. 阶段 1 的结论
1. **起点候选没有变**：InfinyTech3D 的 CollisionAlgorithm + ConstraintGeometry 仍然是唯一现成的、和 SOFA v25.12 兼容的针插入实现。**但是要注意 v25.12 tag 之后的修复**（柔度块求逆失败、拔出检测），以及**还没修的 GS 摩擦 bug**。阶段 2 要决定：用 v25.12 tag、用 master，还是在 v25.12 的基础上挑几个修复补丁。
2. **新发现的有用资源**：
   - SoftRobots.Inverse（QP 逆向控制，已安装）；
   - ICube 原版里修过的求解代码和 Hermite 曲线针几何；
   - `sofaros2plugin`（以后接机器人）；
   - SofaPython3 的矩阵、能量访问示例；
   - Baksic 2022 博士论文。
3. **确认没有公开的**：Martin 的 IsoDOF、混合求解器、CTSVD，Vanneste 的分区降阶，Wang 的仿真器，逆向控制系列的代码。→ 如果需要这些，只能按论文自己实现（handbook 里已经有推导）。
4. **已知的难点**：规划和 MPC 需要保存和恢复状态，但插件的约束点状态在 C++ 内部。

## 8. 阶段 2 的阅读清单（按优先级）
1. InfinyTech3D CollisionAlgorithm：`InsertionAlgorithm.cpp`（约束点怎么生成、穿刺判定、拔出）、`NeedleOperations`、`doc/doc.md`；
2. InfinyTech3D ConstraintGeometry：`ConstraintInsertion`、`InsertionResolution`（摩擦是哪种形式）、`UnilateralResolution` / `BilateralResolution`（bug 的影响）；
3. SOFA 约束流水线：`FreeMotionAnimationLoop`、`BlockGaussSeidelConstraintSolver`、`LinearSolverConstraintCorrection`、`PrecomputedConstraintCorrection`（对应 handbook）；
4. 结合 `diagnostics/` 的诊断场景做运行验证（穿刺判定、柔度块求逆会不会失败）；
5. 可选：Baksic 2022 博士论文里讲 SOFA 实现的章节；SoftRobots.Inverse 的针相关用法。

---
**阶段 2 更新（2026-09-30）**：详细的代码结论和运行验证见 `docs/code/01_needle_plugins_v25.12.md`。最重要的两条：①`frictionCoeff` 实际是 GS 的欠松弛因子，摩擦是求解器没收敛的数值副产品；②（已纠正）穿刺阈值换算本身没错；官方场景把针设成一阶、组织是二阶，同一个 λ 在针一侧是力、在组织一侧是冲量，所以针根的力只有组织受力的 dt 倍（作用力 ≠ 反作用力）。

---

## 补充搜索（2026-09-30，进入第 1 步之前）：加密网格、已验证的场景、公开的验证数据
用户要求：查一下有没有别人已经验证过、可以直接用的资源，特别是加密好的组织网格和验证数据。下面所有条目都核对过来源页面；没有下载任何文件。

### 1. 公开的验证数据（最有价值）
- **TU Delft：针插入两个离体人肝脏的力测量数据集**（de Jong, Dankelman, van den Dobbelsteen, *Data in Brief* 11:308–310, 2017，doi:10.1016/j.dib.2017.01.018；数据在 Mendeley Data，doi:10.17632/94s7xd9mzt.2，**CC BY 许可**）〔原文：PMC5498459 页面〕
  - 18 G 套管针（trocar），三角形针尖；**恒速 5 mm/s** 插入和拔出；
  - 插入深度约 70 mm（整个肝脏厚度）；共 39 次插入（两个肝脏各约 20 次，间距 ≥ 10 mm）；
  - **针座的轴向力**，1 kHz 采样，原始未滤波；每次插入的数据包括力 [N]、位置 [mm]、时间 [s]；MATLAB 格式（.mat/.m/.fig），附带基础分析代码；
  - 肝脏包埋在 10 m% 明胶里：针先穿过明胶、再进入肝脏，天然是两层组织。
  - 对我们的用途：速度、测量量（针座轴向力–位置）和我们的设计完全对应，可以作为第 2–5 步（刺穿、摩擦、切割）定性形状和量级比较的数据，也可以用于第 7 步（多层组织）。注意：肝脏的几何和材料参数需要另外估计，数据本身不含这些。
- 已在资源清单中的：JHU 多层 plastisol 仿体（Wang 2024 / 2025 的配方和参数）、Okamura 2004 / DiMaio 2003 的文献曲线。

### 2. 已验证的 SOFA 插针场景
- 没有找到公开的、和实验做过定量对比的 SOFA 插针场景。能找到的是：
  - Cosserat 插件的插针示例（学习阶段已编译运行过，`docs/code/03`）；
  - InfinyTech3D CollisionAlgorithm 的官方场景（学习阶段已逐行分析，`docs/code/01`）；
  - SOFA 论坛上一个插针帖子（罚函数接触 + `LocalMinDistance`，只调了接触参数，没有验证）。
- 结论：针–组织交互机制的验证仍然要我们自己做，这和此前的判断一致。

### 3. 加密网格的工具和做法
- **gmsh**（开源）的尺寸场（Distance + Threshold）：距离给定几何（点、线，例如针道）≤ DistMin 时单元尺寸为 SizeMin，≥ DistMax 时为 SizeMax，中间线性插值，过渡带自动生成。有 Python API（官方教程 t10）。本机没有安装（`pip install --user gmsh`）。SOFA 的 `MeshGmshLoader` 可以直接读 .msh。→ **第 4.5 步的首选工具**。
- SOFA 的 Tetrahedral Mesh Tool（QTetraMesher，GPL v3）：Delaunay（CGAL）和 isosurface stuffing 两种四面体化方法；**文档没有提到局部尺寸控制**，不适合按针道加密。
- 文献中的做法：
  - Chentanez 2009：针经过时局部重新划分网格，让节点落在曲线针道上（前列腺网格 13,375 个四面体、2,763 个顶点，25 Hz）；
  - Goksel 2005 / 2006：在粗网格上用节点重定位 / 节点增加实现针–组织耦合，用 Woodbury 公式加速，触觉频率 > 1 kHz。
  - 这两种都是动态改网格，工作量远大于静态加密，不在两个月的范围内。我们用的是重心坐标约束（不需要节点落在针道上），配合静态加密网格。
- 非实时的详细有限元研究（Abaqus 等）的通常做法：接触区附近应变能密度梯度大，所以近处细、远处粗，并做网格收敛测试（例如最小单元 ≈ 0.5 倍针尖尺寸）。这些模型的单元比我们实时预算允许的细得多，只能作为方法参考，不能照搬数值。

### 4. 对计划的影响
- 第 4.5 步（非均匀网格）：用 gmsh 的 Distance + Threshold 尺寸场生成，管半径、细单元尺寸、过渡距离由收敛测试和节点预算决定；需要先征得用户同意安装 gmsh。
- 验证数据：TU Delft 数据集可以在第 2 步之后用来做定性和量级比较；下载前征得用户同意。

### 补充（2026-10-01）：SOFA 官方的 probe–tissue interaction 教程，以及现成网格
- **教程**：Eleonora Tagliabue（Verona 大学 Altair 实验室），https://www.sofa-framework.org/applications/plugins/tutorial-probe-tissue-interaction/ ，代码 https://gitlab.com/altairLab/probe-tissue-simulation （已浅克隆到 `~/sofa/resources/altairLab_probe-tissue-simulation`，分支 sofapython3，commit dc8116e，2022-07；只读）。论文：Tagliabue et al., "Biomechanical modelling of probe to tissue interaction during ultrasound scanning", IJCARS 2020。
  - README：同一个"超声探头压乳房组织"的场景，原本有罚函数、约束（拉格朗日乘子）、直接给定表面位移三种实现；SofaPython3 分支只实现了约束法。测试过的版本是 SOFA 21.12 / 22.06（`GenericConstraintSolver` 在 v25.12 里已经不能直接创建，要换成 `BlockGaussSeidelConstraintSolver` 等）。**README 没有写许可证**：只读学习，不直接拷贝代码。
  - 关键设计〔代码〕：
    - 探头是没有求解器、没有 ConstraintCorrection 的 Rigid3 物体，由控制器设定位姿，和我们 0b 的纯运动学针相同；
    - 组织碰撞用单独的低分辨率表面（`breast_500.stl`）经 `BarycentricMapping` 挂在体网格上；
    - 接触用 SOFA 自带的碰撞流水线（`MinProximityIntersection` / `LocalMinDistance`，三角形、线、点碰撞模型），`alarm_distance` = 2 mm，`contact_distance` = 1 mm，dt = 0.02 s；
    - 材料 NeoHookean，E = 2850 Pa，ν = 0.49；`inData/ground_truth/` 里有两个肿瘤在 4 级压深下的实测位置，是一套实验对照数据。
  - 网格〔运行：读 .msh 统计〕：`breast_13k.msh` 有 2469 个节点、9295 个四面体，单元平均边长中位数 8.6 mm（5%–95% 分位 5–11 mm），基本均匀，不是加密网格。
  - **对我们的价值**：
    - ① 第 1 步 A2（SOFA 自带接触）交叉验证的现成参考配置：接触距离怎么取、碰撞表面怎么挂；
    - ② "约束 vs 罚函数 vs 直接给定位移"的对比思路，和我们"接触 vs 给定位移"的验证一致；
    - ③ 面接触（探头端面）的做法，可以参考用于以后的多点针尖。
  - **局限**：没有刺穿和插入；面积大的探头接触和针尖点接触不同；SOFA 版本较旧。
- **现成网格**：
  - 均匀的组织 / 器官网格是有的：SOFA 自带（`share/sofa/mesh` 下的 liver、cube 等）、上面这个乳房网格、各种器官或仿体的几何数据集。
  - **没有找到可以直接下载的沿针道加密的网格**。文献里这类网格都是针对具体针道生成或动态调整的（Chentanez 2009 局部重新划分网格；Bui 等用误差估计实时自适应加密，arXiv 1704.07636，和文献笔记 `10_Bui2018` 同一系列）。
  - 原因〔推论〕：加密网格取决于针道位置、组织几何和节点预算，本来就是"针对问题生成"的东西。可复用的是**工具（gmsh）和方法（距离 / 误差驱动的尺寸场）**，不是某个具体网格。
