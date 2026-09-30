# 代码笔记 03：其他插件（Cosserat、SoftRobots.Inverse、ModelOrderReduction、BeamAdapter、拓扑切割）

- 日期：2026-09-30。标记同笔记 01、02。
- 本地：克隆在 `~/sofa/resources/`；Cosserat 编译在 `~/sofa/resources/{build,install}_cosserat`。

---

## C1. Cosserat（SofaDefrost，release-v25.12）〔代码 + 编译 + 运行〕
### 编译（可以在 SOFA v25.12 二进制包上编译，但需要三个变通）
1. **pybind11 版本**：SOFA 二进制包要求 pybind11 ≥ 2.12（系统里是 2.11.1）→ 把 v2.12.0 的源码下载到 `~/sofa/resources/pybind11-2.12.0`，本地安装到 `pybind11-install`，然后 `-Dpybind11_DIR=…/share/cmake/pybind11`（系统没有 python3-venv，没用 pip）。
2. **SoftRobots 的 CMake 配置不完整**（引用了包里没有的文件）→ `-DCMAKE_DISABLE_FIND_PACKAGE_SoftRobots=TRUE`（SoftRobots 对 Cosserat 是可选依赖）。
3. **pybind11 的 ABI 标识不一致**：SOFA 二进制包是 GCC 11 编译的（`__pybind11_internals_v5_gcc_libstdcpp_cxxabi1016__`），本机是 GCC 13（`cxxabi1018`）→ Python 绑定导入时报"unknown base type"→ 用 `-DCMAKE_CXX_FLAGS=-DPYBIND11_BUILD_ABI=\"_cxxabi1016\"` 重新编译。**以后我们自己写的任何 Python 绑定都要这样处理。**
- 另外 `-DCOSSERAT_BUILD_TESTS=OFF`。运行时 `PYTHONPATH` 要加上 `$SOFA_ROOT/plugins/STLIB/lib/python3/site-packages`（`splib3`、`stlib3`）。

### 针插入示例的状态
- **上游示例本身是坏的**（release-v25.12 和 master 都一样）：`cosserat/cosseratObject.py` 从 `cosserat.utils` 导入 `addEdgeCollision`，而这个函数实际在 `useful/utils.py`；示例还用了 `ProjectedGaussSeidelConstraintSolver`，**SOFA v25.12 里没有这个组件**。在 scratchpad 的副本里改了这两处之后才能运行（原仓库没动）。
- 示例是**键盘驱动**的；不开界面运行时改成每步把针根沿 x 推进 0.05 cm。
- 模型〔代码〕：
  - 针：Cosserat 杆（E = 1.2e9，单位制 cm-kg-s，长 15 cm，16 段）；组织：规则网格 FEM（E = 100，ν = 0.48）。
  - **穿刺**：Python 控制器读取 GS 的约束力，第一个分量大于 3 时关掉针和表面的碰撞，把进针点设为第一个约束点。
  - **插入约束**：`CosseratNeedleSlidingConstraint` + `DifferenceMultiMapping`（组织里的约束点到针上投影点的差），**只约束两个横向方向**（`useDirections` = (0, 1, 1)），轴向**完全没有约束 → 没有摩擦**；`BilateralConstraintResolution`。
  - 约束点：针上的滑动点沿 x 方向超过 `constraintDistance` 就加一个；往回退就删一个。**插入方向写死为 x 轴**（`computePositiveAlongXDistanceBetweenPoints`）。
- 运行结果：推进约 2 cm 后刺穿；每 1.5–2 cm 加一个约束点（推进 20 cm 时共 14 个）；约束力先随深度增大，16 cm 以后**大幅振荡**（7.1 → 3.4 → 0.9 → 5.2 → 8.7 → 4.6 → 1.0）→ 原型级别，深插入时不稳定。
- 结论：如果 Manish 坚持用 Cosserat 杆，**要自己把 CollisionAlgorithm 那一套交互（或我们重写的版本）接到 Cosserat 杆上**；这个示例只能作为"怎么把约束点映射到 Cosserat 杆"的参考。

## C2. SoftRobots.Inverse（release-v25.12，已装在二进制包里）〔代码〕
- ⚠️ **许可证：AGPL-v3，附带 Inria 专利**（FR3002047 等，"Method for controlling of a deformable robot"），AGPL 第 11 条授予使用这项专利的许可 → 以后如果项目代码依赖它并对外发布，要遵守 AGPL。
- 模型：`QPInverseProblemSolver`（qpOASES 或 proxQP）在约束空间里建一个 QP：**执行器**（actuator，未知的力或位移，可以设上下限）、**效应器**（effector，希望到达的目标）、**等式约束**、**传感器**，以及"接触"。
- ⚠️ **分类规则**（`QPMechanicalSetConstraint.cpp:70-115`）：不是 SoftRobots 自己类型的约束，**一律当作"接触"**（带 `responseFriction` 的单边约束，按 LCP / QPCC 处理）→ ConstraintGeometry 的针插入约束（本应是双边的）**会被当成单边接触**。
- 可以借鉴的组件：`PositionEffector`（例如针尖到靶点）、`SlidingActuator` / `JointActuator` / `ForcePointActuator`（例如针根运动）、`BeamRestPositionActuator`。
- 结论：它可以作为"逆向 FE 控制"的现成框架（Coevoet 2019 的路线，和 Adagolodjo 的数值 Jacobian 路线不同），但**要和针插入约束配合使用，需要把针的约束改写成 SoftRobots 的等式约束类型，或者修改它的分类逻辑**。

## C3. ModelOrderReduction（release-v25.12，已装）〔代码〕
- 组件：`ModelOrderReductionMapping`、超降阶力场（`HyperReducedTetrahedronFEMForceField`、`HyperReducedTetrahedralCorotationalFEMForceField`、`HyperReducedTetrahedronHyperelasticityFEMForceField` 等）、`MORUnilateralInteractionConstraint`、`MORContactMapping`。
- Python 工具（二进制包里有：`plugins/ModelOrderReduction/lib/python3/site-packages/{mor,morlib}`）：定义训练动作 → 运行场景采集快照 → POD → ECSW 超降阶 → 自动生成降阶场景。
- **Vanneste 2024 的分区降阶在公开仓库里没有**。
- 对针插入：约束的 H 通过 `ModelOrderReductionMapping` 映射到降阶坐标，原理上自定义约束也能用；但被穿刺区域的局部大变形很难用少量模态描述（这正是分区降阶要解决的问题）〔推导〕。

## C4. BeamAdapter（release-v25.12，已装）〔代码〕
- 主要面向导管和导丝（长度可变、从进入点送入）：`InterventionalRadiologyController`、`WireRestShape`、`AdaptiveBeamForceFieldAndMass`（Kirchhoff 杆）、`BeamInterpolation`（沿梁的形函数插值）。
- 可以借鉴的构件：`BeamProjectionDifferenceMultiMapping`（点到梁的投影差，作用类似 Cosserat 的 DifferenceMultiMapping）、`AdaptiveBeamSlidingConstraint`（刚体附着在梁上，只约束位置，可以滑动）。**`BeamInterpolation` 让针身光滑**，可以避免约束方向在节点处跳变（和 ICube 的 HermiteEdge 是同一个思路）。
- 针的长度是固定的，所以主要是参考，不是替代方案。

## C5. 拓扑切割和 InfinyToolkit〔代码〕
- **Tearing**：`TearingEngine` 基于 `TriangularFEMForceField` 的应力，**只支持二维三角形网格** → 不适用于三维组织里的针。
- **SofaCarving**：`CarvingManager` 删除工具附近的单元（"切除组织"），可以作用在四面体上面的三角形表面上 → 和针刺穿（不删除组织）是不同的现象。
- **InfinyToolkit / NeedleTracker**：判断针尖在哪一个"slice"网格里（包围盒粗筛 + 射线奇偶判断）→ 可以用于多层组织里"针尖当前在哪一层"的判断（方案 C 的界面检测）。
- → 和文献一致（Perrusi 2021 等）：**实时针插入里的"破裂"都不做拓扑修改**。

## 汇总：各插件对项目的作用
| 需求 | 最合适的资源 | 状态 |
|---|---|---|
| 针插入交互框架 | InfinyTech3D CollisionAlgorithm + ConstraintGeometry | 框架可用，力学部分要重写（笔记 01） |
| 刚性针 / 柔性针 | `BeamFEMForceField`（实际是 E-B 梁）+ `BTDLinearSolver`；Cosserat（要接交互）；BeamAdapter（参考） | 可用 |
| 多层组织 | 逐单元设材料（SOFA 原生）；超弹要拆成多个力场；NeedleTracker 判断层 | 可用 |
| 逆向控制 | SoftRobots.Inverse（AGPL + 专利；要改约束分类）；或者按 Adagolodjo 的方法自己做 | 要评估 |
| 提速 | 预计算柔度、MOR（标准版）；分区降阶要自己做 | 可用 / 要自己做 |
| 规划采样 | 进程分叉快照（笔记 02 B7） | 已验证可行 |
