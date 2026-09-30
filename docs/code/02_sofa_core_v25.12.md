# 代码笔记 02：SOFA v25.12 本体中和针插入相关的机制

- 日期：2026-09-30
- 源码：`~/sofa/resources/sofa_v25.12`（tag v25.12.00，commit 4d6c2b8，和安装包一致；部分目录）；SofaPython3：`~/sofa/resources/sofa-framework_SofaPython3`（`release-v25.12`）。
- 标记：〔代码〕〔运行〕〔推导〕〔待查〕同笔记 01。实验脚本在 `experiments/`。

---

## B1. 一个时间步（`FreeMotionAnimationLoop::step`，`FreeMotionAnimationLoop.cpp:154-322`）〔代码〕
1. `AnimateBeginEvent` → `CollisionLoop` 调用 `InsertionAlgorithm::doDetection()`：用**步首位置**；穿刺判定读取**上一步存下的 λ**；
2. `BehaviorUpdatePosition`、`UpdateInternalData`、`MechanicalBeginIntegration`；
3. **λ ×= 1/dt**，用于映射的几何刚度（`lambdaMultInvDt`）；
4. **自由运动**：`SolveVisitor`（ODE 求解器）→ freePos = pos + freeVel·dt；接着碰撞检测（`CollisionVisitor` 里调用约束对象的 `processGeometricalData()` → 约束容器按检测结果重建）；
5. **约束求解**：`buildSystem`（H、δ_free 用自由位置、W、ConstraintResolution）→ GS → 运动修正 → 存 λ；
6. `AnimateEndEvent` → `UpdateMapping`。
- 可选并行：`parallelCollisionDetectionAndFreeMotion`、`parallelODESolving`。

## B2. ODE 求解器和 λ 的单位〔代码 + 运行〕⚠️
- `EulerImplicitSolver`（二阶）：系统矩阵 = M(1 + h·rm) − h·B − h(h + rs)·K_sofa（SOFA 的 K 是 ∂f/∂x，为负定），解出**速度增量**，右端项 = h(f + …)（`EulerImplicitSolver.cpp:136-171`）。
- `firstOrder=True`：系统矩阵 = M − h·K_sofa，右端项 = f（**不乘 h**）。
- `getPositionIntegrationFactor()` = **dt**，`getVelocityIntegrationFactor()` = 1（`EulerImplicitSolver.h:132-140`）；约束柔度 W = J A⁻¹ Jᵀ × `correctionFactor`，POS_AND_VEL 时 `correctionFactor` = 位置积分因子（`BaseConstraintCorrection.cpp:54-74`）。
- ⇒ **二阶物体：λ 是冲量（N·s）**，力 = λ/dt；**一阶物体：同一个 λ 被当作力**。混用时作用力 ≠ 反作用力，相差 dt 倍（运行验证见笔记 01 问题 2）。SOFA 自带的 `SlidingLagrangianConstraint` 文档里也写着 "force (impulse)"。
- `StaticSolver`（配合 `NewtonRaphsonSolver`，`maxNbIterationsNewton` 默认 1，带线搜索）：位置积分因子也是 dt，速度因子 1（`StaticSolver.h:60-125`）。静力求解的是 KΔx = f，按 W = dt·J K⁻¹ Jᵀ 推导，λ 相当于 f/dt〔推导，未运行验证〕→ **又是一种不同的约定，不能和 EulerImplicit 的物体放在同一个约束问题里**。`BDFOdeSolver`、`NewmarkImplicitSolver` 也有〔没细读〕。
- **规则**〔推导 + 运行〕：同一个约束问题里，所有物体的 ODE 求解器要用同一种约定（都用二阶 EulerImplicit 最稳妥），否则同一个 λ 在不同物体上代表不同的物理量。

- **`EulerImplicitSolver` 之后 `force` 里是什么**〔代码，2026-09-30 补充〕：每步开始时在 (xₙ, vₙ) 上计算的内力 f（`EulerImplicitSolver.cpp:92,129`），不被投影（被投影的是它的拷贝 b，:145,160），不含 Rayleigh 项（:150-152 单独加到 b 上）。所以第 k+1 步之后读 `force`，得到的是第 k 步结束时位置上的弹性内力（受约束节点上也有值）。
- **数值阻尼**〔推导 + 运行〕：对无阻尼振子，特征值 z = 1/(1 − iωh)，σ = ln(1 + ω²h²)/(2h)。实测和理论相差 1% 以内（`docs/plan/step0_results.md` 0a-3）。
- **Rayleigh 刚度阻尼会抬高准静态力**，比例约为 r_s·v/δ（实测吻合）；r_sω/2 > 1 时变成过阻尼蠕变，时间常数约为 r_s。
- **`PartialLinearMovementProjectiveConstraint::projectVelocity` 会设置整个速度向量**（`.inl:212-214`，不看 `movedDirections`），所以未被驱动的方向速度也会被清零。

## B3. ConstraintCorrection（W 的计算方式）〔代码〕
| 组件 | 做法 | 条件和代价 |
|---|---|---|
| `LinearSolverConstraintCorrection` | 每步用关联的线性求解器计算 J A⁻¹ Jᵀ（`addJMInvJt`）| 通用；每个约束要一次回代。`wire_optimization`：约束沿线状拓扑从尖端到根部重排（适合针）|
| `PrecomputedConstraintCorrection` | 初始化时对每个自由度施加单位力做一次 EulerImplicit 求解，得到**稠密的 N×N A⁻¹**；`rotations` 做旋转修正；可以存成文件（文件名含维数和 dt）| **必须用 EulerImplicit**；**dt 固定**；线性（或共旋小变形）；内存 O(N²) |
| `UncoupledConstraintCorrection` | 每个自由度一个对角柔度 | 很快但粗糙，没有耦合 |
| `GenericConstraintCorrection` | 任意线性求解器 + `complianceFactor` | 通用 |

## B4. 约束求解器和 SOFA 自带的约束〔代码〕
- 求解器（v25.12 注册的只有这 5 个，已核对组件列表）：`BlockGaussSeidel`（按块调用 `ConstraintResolution`，官方针场景用的就是它）、`UnbuiltGaussSeidel`（不组装完整的 W）、`NNCG`（非光滑非线性共轭梯度）、`ImprovedJacobi`（投影 Jacobi，可以并行）、`LCP`。**v25.12 里没有 `ProjectedGaussSeidelConstraintSolver`**（那是之后版本的名字；`GenericConstraintSolver` 在 v25.12 是基类，不能直接创建）。
- 共同参数：`tolerance`（默认 1e-3）、`maxIterations`（默认 1000）、`sor`、`regularizationTerm`、**`scaleTolerance`（默认 true：容差 × 约束数）**、`allVerified`、`computeConstraintForces`。
- **每步的约束力从 0 开始，没有热启动**（`FullVector::resize` 会清零）。
- 自带的约束：`Bilateral`（附着）、`Unilateral`（带 mu 的接触）、`Sliding`（一个点在两个**固定**节点之间的线段上滑动，不适合针的插入）、`Stopper`、`Fixed`、`Uniform`、`AugmentedLagrangian`。

## B5. 组织 FEM〔代码 + 运行〕
- **可以逐单元设材料**：`TetrahedronFEMForceField`、`TetrahedralCorotationalFEMForceField`、`FastTetrahedralCorotationalForceField` 都继承自 `BaseLinearElasticityFEMForceField`，`youngModulus` 和 `poissonRatio` 都可以是向量（`getVecRealInElement`：向量长度大于单元编号就取对应的值，否则取第一个值）→ **多层方案 B 原生支持**。`TetrahedronFEM` 和 `TetrahedralCorotational` 还有 `localStiffnessFactor`。
- 超弹：`TetrahedronHyperelasticityFEMForceField`，材料有 ArrudaBoyce / Costa / MooneyRivlin / NeoHookean / **Ogden** / StVenantKirchhoff / VerondaWestman / StableNeoHookean；**参数是全局的一组**（`ParameterSet`），多层要拆成多个力场。
  - **Ogden 的参数换算**〔推导〕：SOFA 的 W_iso = (μ₁/α²)(J^(−α/3) tr C^(α/2) − 3)，W_vol = k₀(ln J)²/2，参数为 [μ₁, α₁, k₀]；Wang 2024/2025 用的是 W = (2μ/α²)(Σλ^α − 3) ⇒ **μ₁(SOFA) = 2 μ(Wang)**。
- 粘弹：`TetrahedronViscoelasticityFEMForceField`、`TetrahedronViscoHyperelasticityFEMForceField`（SofaViscoElastic 插件）。
- ⚠️ **弹性应变能读不到**：`FastTetrahedralCorotational` 和 `TetrahedralCorotational` 的 `getPotentialEnergy()` 没有实现（返回 0 并警告）；`TetrahedronFEM` 只在 `method="small"` 时计算，`large` / `polar` 都返回 0。`Node.computeEnergy()` 返回的势能因此只包含重力和弹簧〔运行确认〕。→ 如果需要应变能（比如 Wang 2025 的代价函数），要自己算。
- ⚠️ **v25.12 的 `FastTetrahedralCorotationalForceField::buildStiffnessMatrix` 有 bug**〔代码 + 运行，2026-09-30 补充〕：两个非对角块都写成了 −M（应为 −M 和 −Mᵀ，`FastTetrahedralCorotationalForceField.inl:551-552`），组装出的刚度矩阵在 ν ≠ 0.25 时是错的。凡是用组装矩阵的线性求解器（`EigenSimplicialLDLT`、`SparseLDLSolver` 等，包括插件官方针场景）都受影响；CG（走 `addDForce`）不受影响。静力测试中直接求解器的 Newton 发散，CG 收敛。上游 PR #6154（2026-06-26，v26.12）已修。**本项目的修复**：`plugins/NeedleSimFixes` 中的 `FastTetrahedralCorotationalForceFieldFixed`（继承原类，只覆盖 `buildStiffnessMatrix`）。`addKToMatrix`（旧接口）是正确的。详见 `docs/plan/step0_results.md` 发现 1。
- **旋转提取方法的差别**〔代码 + 运行〕：只有 FTC 的 `method="polar"` 对真正的变形梯度做极分解，均匀拉伸时结果精确；FTC `qr`（默认）/ `polar2` 和 TFEM `large` / `polar` 都会把拉伸误算成 O(ε) 的转动（单轴测试中横向应变误差约为 −2ε ～ −3ε）。
- 可以输出 von Mises 应力（`TetrahedronFEMForceField` 的 `computeVonMisesStress`）。

## B6. 针〔代码〕
- `BeamFEMForceField`：共旋（`large`）梁单元，按 Przemieniecki 的形式组装刚度，带剪切参数 φ；**但初始化时有效剪切面积 `_Asy` = `_Asz` = 0 ⇒ φ = 0 ⇒ 实际上是 Euler-Bernoulli 梁**（`BeamFEMForceField.inl:794-810`，注释里说 Timoshenko 用 10/9）。支持**空心截面**（`radiusInner`）；截面参数：A = π(r² − rᵢ²)，I = π(r⁴ − rᵢ⁴)/4，J = 2I，G = E/(2(1 + ν))。可以给子集（`listSegment`）。
- `BTDLinearSolver`：块三对角（6×6），Thomas 算法；`subpartSolve` 可以只算一部分（给柔度用）。

## B7. SofaPython3〔代码 + 运行〕
- 环境变量见资源清单 §1.1；不开界面运行：`Sofa.Simulation.initRoot(root)`，然后循环 `Sofa.Simulation.animate(root, dt)`。官方针场景 100 步约 0.6 s。
- 控制器事件：`onAnimateBeginEvent`、`onAnimateEndEvent`、`onBuildConstraintSystemEndEvent`（在这里可以用 `solver.W()` 读 W，需要 `from Sofa import SofaConstraintSolver`）、`onKeypressedEvent` 等。
- 读约束力：`computeConstraintForces=True` 之后读 `solver.constraintForces`（GS 的 f 向量，单位约定见 B2）。
- 能量：`node.computeEnergy()` → (K, U)，但见 B5 的限制。
- ⭐ **状态保存和恢复：用进程分叉快照（`os.fork()`）**〔运行，`experiments/fork_snapshot.py`〕：插入到 0.8 s 后分出 3 个子进程，各跑 20 步；设置相同的两个子进程和父进程结果**逐位相同**，改了参数的那个不同；**插件内部的约束点状态也一起复制**（子进程接着插入）；3 个分支并行共 0.18 s。→ 可以作为 CE / MPC 前向采样的基础（Linux 上可用；要注意各进程的内存）。
- 其他官方示例：矩阵访问（`access_*_matrix.py`）、接触力、计时器（`SofaRuntime.Timer`）、`RPYC/`（远程服务）、`jax/`。

## B8. 性能〔运行〕
- 官方基础场景（6×6×6 网格 → 四面体，针 20 个单元，插入阶段约 39 个约束），用你之前那次运行的计时文件（第 100–149 步）统计：**一步 8.8 ms**。
  - **Get Compliance（组装 W）：5.1 ms（58%）**
  - FreeMotion：2.5 ms（29%），其中矩阵分解 1.3 ms
  - GS 求解：0.4 ms（4%）
  - 碰撞 / 宽相位：约 0.3–0.4 ms
- → 和文献（Ha 2024：W 占 74–82%）一致：**提速的重点在 W**（预计算柔度、IsoDOF 思路、降阶），GS 不是瓶颈。
- **直接求解器的分解时间**〔运行，2026-09-30，第 0 步；静力问题，FTCfix，规则网格切四面体，每次 Newton 迭代都重新分解〕：

  | 自由度 | 组装矩阵 | `SparseLDLSolver`（默认 AMD 重排序）分解 |
  |---|---|---|
  | 2 187（n = 9） | 3 ms | 16 ms |
  | 6 591（n = 13） | 10 ms | 230 ms |
  | 14 739（n = 17） | 27 ms | 1.67–2.1 s |
  | 20 577（n = 19） | — | 6.4 s |
  | 27 783（n = 21） | — | 17 s |

  - 分解时间约按自由度的 **3.3 次方**增长（改正记录：最初只用前三个点估计为 2.4 次方，增加 n = 19、21 两个点后改正）。外推到 10.8 万自由度约 25 分钟一次。瓶颈是分解，不是组装。
  - `EigenSimplicialLDLT` 和 `SparseLDLSolver` 差不多（都是 simplicial，没有 supernodal）。
  - 可选的重排序组件：`AMDOrderingMethod`（默认，最快）、`NaturalOrderingMethod`（慢约 1.7 倍）、**`COLAMDOrderingMethod`（n = 17 时要 90 s，不要用）**。二进制包里**没有 `MetisOrderingMethod`**。
  - → 0d 要处理：实时网格的自由度上限，或者换不需要每步重新分解的方案。
- 多线程组件：`ParallelTetrahedronFEMForceField`、`BeamLinearMapping_mt`、动画循环的并行选项；`MultiThreading` 插件〔没测〕。
