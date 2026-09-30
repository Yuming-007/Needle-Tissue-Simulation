# 代码笔记 01：CollisionAlgorithm + ConstraintGeometry（InfinyTech3D，v25.12.00）与 SOFA 约束流水线

- 日期：2026-09-30
- 读的版本：`~/sofa/reference_sources/{CollisionAlgorithm,ConstraintGeometry}`（tag `v25.12.00`，也就是我们编译安装的版本）；SOFA 本体 v25.12.00（commit 4d6c2b8，和安装包一致；源码在 `~/sofa/resources/sofa_v25.12`，只拉了部分目录）。
- 标记：**〔代码〕** 表示读代码确认；**〔推导〕** 表示从代码逻辑推出；**〔运行〕** 表示用 `experiments/run_exp.py` 实测验证；**〔待查〕** 表示还没确认。
- 时间步里各环节的顺序（`FreeMotionAnimationLoop.cpp:154-322`）：①`AnimateBeginEvent`：`CollisionLoop` 在这里调用 `InsertionAlgorithm::doDetection()`（用步首位置，读上一步的 λ）；②更新行为模型，λ 乘以 1/dt 用于几何刚度；③自由运动，freePos = pos + freeVel·dt，接着碰撞检测，在 `CollisionVisitor` 里调用 `processGeometricalData()` 重建约束容器；④约束求解：W、δ_free（用自由位置）、GS、修正运动、存 λ；⑤`AnimateEndEvent`。〔代码〕
- 实验：官方场景 `scenes/NeedleInsertion.py`（组织 E = 1 MPa，6×6×6 规则网格转四面体；针 E = 1e12，20 个梁单元；dt = 0.01 s；针根用 ReadState 回放轨迹，前 1 s 以 90 mm/s 竖直插入，最终深度约 4 cm）。不开界面运行，100 步约 0.6 s。原始数据在 `experiments/results_2026-09-30.tar.gz`。

---

## 1. 算法流程（`InsertionAlgorithm.cpp`）〔代码〕
| 阶段 | 条件 | 做什么 | 行号 |
|---|---|---|---|
| 穿刺检测 | `m_couplingPts` 为空 | 针尖 ↔ 组织外表面最近点，输出一对接触点（交给 `ConstraintUnilateral`）；如果针尖上的约束力超过阈值，就把这个**表面点**存为第一个约束点 | 150–205 |
| 针身碰撞 | 还没刺穿 | 针身各边 ↔ 表面最近点 | 207–254 |
| 插入 | 已经刺穿 | ①删掉针尖"前方"的约束点；②针尖离最后一个约束点超过 `tipDistThreshold` 时，沿**针身方向**每隔 `tipDistThreshold` 取候选点，找最近的**组织四面体**点，只有候选点真的在四面体内才接受 | 256–349 |
| 重投影 | 每一步 | 每个约束点（组织上的物质点，重心坐标固定）↔ 针身上的最近点，交给 `ConstraintInsertion` | 351–378 |

要点：
- **约束点是组织里的物质点**，针那一侧的点每步重新投影，相当于针在固定的针道里滑动（Duriez 2009 的模型）。〔代码〕
- 约束点间距固定为 `tipDistThreshold`（不是网格相交法），沿针身方向生成（L316）。〔代码〕
- **刺穿之后针尖没有任何约束**：碰撞输出被清空（L110），针身碰撞也跳过（L124）→ **没有切割力**。轴向阻力只来自针身约束点。〔代码〕
- **只检测一个外表面（`surfGeom`）和一个体网格（`volGeom`）的穿刺**，不检测层间界面。一个体网格里分区域设材料是可以的（约束点可以落在任意四面体上），但层间界面不会产生穿刺事件（handbook §5.3 方案 B 的现状）。〔代码〕
- 拔出：只靠"删除针尖前方的约束点"。约束点全部删完后，回到穿刺检测状态。〔代码〕

## 2. 发现的问题（按严重程度排序）

### 问题 1 ⚠️⚠️ `frictionCoeff` 不是摩擦系数；摩擦是 GS 没收敛的数值副产品〔代码 + 推导 + 运行〕
- 代码：`InsertionResolution.h` 的 `resolution()`：
  - `corr[0] = frictionCoeff × (−invW·d)[0]`：先算出三行全部满足所需的修正，取轴向分量，乘以 `frictionCoeff`；
  - `corr[1]`、`corr[2]`：在已知 corr[0] 的条件下精确求解横向两行。
- 推导：收敛时所有修正为 0 ⇒ d₁ = d₂ = 0，并且（只要 frictionCoeff > 0）d₀ = 0 ⇒ **完全粘住**。frictionCoeff 只是轴向更新的**欠松弛因子**。
- 加上 SOFA GS 的三个特点：①**每步的力从 0 开始**（`FullVector::resize` 会清零，`FullVector.inl:100-111`；`GenericConstraintSolver::buildSystem` 调用 `clear`）；②停止判据是 Σ‖W_block·Δf‖ < tol（`BlockGaussSeidelConstraintSolver.cpp`，`gaussSeidel_increment`）；③默认 `scaleTolerance=true`，容差 = tol × 约束个数（`GenericConstraintSolver.cpp` 构造函数和 `doSolve`）⇒ 每步迭代 k 次后，轴向力 ≈ (1−(1−μ)^k) × 粘住所需的力 ≈ kμ × 粘住所需的力。**有效摩擦 = η_eff × 粘住力，η_eff 取决于迭代次数**（本质上就是 handbook §3 的 η 型摩擦，但 η 也不固定）。
- 运行结果（1.0 s 时，插入约 4 cm，针根轴向阻力，已扣除针的重量）：

  （以下都是**针根一侧**的力，受问题 2 影响，组织一侧 = 表中数值 × 100）

  | frictionCoeff | tol 1e-5，maxIt 5000（场景默认） | tol 1e-8（收敛） | tol 1e-5，maxIt 2 |
  |---|---|---|---|
  | 0 | 0.00 N | 0.01 N | 0.00 N |
  | 0.002（场景默认） | 0.08 N（每步只迭代约 4 次） | **51 N，完全粘住** | 0.08 N |
  | 0.02 | 35 N | **51 N** | 0.9 N |
  | 0.2 | 48 N | **51 N** | 21 N |
  | 1.0 | 51 N | 51 N | 51 N |

  "完全粘住"的表现：刺穿后约束数始终是 3（只有第一个约束点，不再新增），力每步线性增加约 1.2 N（整块组织被针推着走）。
- **后果**：①这个参数**无法用实验数据标定**（没有物理量纲，还依赖求解器设置）；②官方场景的默认值 0.002 实际上**几乎没有摩擦**（插入 4 cm 只有约 0.08 N，而 Okamura 的数据是约 0.5 N 以上）；③改 tolerance 或 maxIt 会让力学行为发生质变。
- 对照：handbook §2 的写法 (c)（η 型）。**要换成有物理意义的摩擦律**（handbook §4 的建议：μ × (径向预压力 + 横向约束力)，或者固定阈值 N/m × l_k）。

### 问题 1b ⚠️ 问题 1 在反复插拔中的失效方式：粘住后的恶性循环，组织最后挂在针上〔运行 + 推导〕
- 场景：官方 `NeedleInsertionCycles.py`（针 E = 1e11，组织 E = 8e5、8×8×8 网格，`punctureForceThreshold` = 200，`tipDistThreshold` = 3 mm，**场景自己的 `frictionCoeff` = 0**），两次插入和拔出。脚本 `experiments/run_exp_scene.py`。
- **frictionCoeff = 0（场景原始设置）**：两次插拔都正常（约束点最多 24 个，拔出时逐个删除，组织位移在 ±3 mm 以内）；只在两次重新刺穿的时刻，GS 达到 maxIt（5001 次）还没收敛。
- **frictionCoeff = 0.002**：第一次插拔正常；**第二次插入时针和组织粘住**：针尖插深 35 mm，约束点始终只有 3 个，组织被拖下去 32 mm（针根一侧的力约 21 N）；拔出时这 3 个约束点被组织带着走，始终在针尖后方，**删不掉**；针尖回到组织外（−92 mm，表面在 −100 mm）之后，**组织还被往上拽了 16 mm**，挂在针上。
- 机理：粘住 → 不再新增约束点 → 约束数少 → 容差（tol × 约束数）变严 → GS 迭代更多 → η_eff = 1−(1−μ)^k 趋近 1 → 更粘，形成恶性循环；删除约束点只看"是否在针尖前方"，约束点跟着组织走时永远满足不了这个条件。
- **结论：官方的反复插拔示例能正常工作，靠的是把摩擦关掉。**

### 另外两个官方场景（A5）〔代码 + 运行〕
- `NeedleInsertionCycles.py`：见问题 1b。**frictionCoeff = 0**。
- `NeedleInsertionHaptics.py`：针长 20 cm、40 个单元、E = 2e13、一阶；组织 E = 4e4（软）、8×8×8 网格；没有重力；`punctureForceThreshold` = 100，`tipDistThreshold` = 10 mm；**frictionCoeff = 0**；表面接触 `ConstraintUnilateral` 设了 **mu = 0.001** → 会触发问题 5 的 GS bug；带 `LCPForceFeedback`，可以接 Geomagic（默认关闭，用 ReadState 回放）。运行 501 步（`experiments/run_generic.py`）：最多 123 个约束（41 个约束点），组织最大位移 68 mm（大变形），4 步 GS 达到迭代上限；插入和拔出都正常。
- **三个官方场景里，只有基础场景设了非零的 frictionCoeff（0.002），而它实际上几乎没有摩擦**（问题 1）→ 可以认为官方示例**都是在"无摩擦"的状态下工作的**。

### 问题 2 ⚠️⚠️（2026-09-30 已纠正）穿刺阈值和 λ 的单位：真正的问题是场景里一阶针和二阶组织混用，导致作用力 ≠ 反作用力〔代码 + 运行〕
- **纠正**：这一节原来的结论是"λ 是力，插件多除了一次 dt"，**这是错的**。
- 代码：
  - 组织用二阶 `EulerImplicitSolver`：A = M(1+h·rm) + h(h+rs)K，解出速度增量，右端项 = h·(f + …)（`EulerImplicitSolver.cpp:136-171`）；
  - 针用 **`firstOrder=True`**：A = M + hK，右端项 = f（不乘 h）；
  - W 的积分因子 = `getPositionIntegrationFactor()` = **dt**（`BaseConstraintCorrection.cpp:54-74`，`EulerImplicitSolver.h:140`）。
  - ⇒ 对二阶物体，W = dt·J A⁻¹ Jᵀ，**λ 被当作冲量**（N·s）；对一阶物体，同一个 W 下 **λ 被当作力**（N）。
- 运行（`experiments/check_newton3.py`，刺穿前的准静态阶段）：

  | 设置 | 针根弹簧力 | GS 的 λ | 组织底部反力 − 自重 | 压入深度 |
  |---|---|---|---|---|
  | 官方场景（针一阶，dt = 0.01） | 16.4 N | 16.47 | ≈ 1647 N（= λ/dt） | 12.1 mm |
  | 针一阶，dt = 0.005 | 8.2 N | 8.24 | ≈ 1643 N（= λ/dt） | 12.1 mm |
  | 针二阶，dt = 0.01 | 232.5 N | 2.326 | 233 N（= λ/dt） | 2.4 mm |

- **结论**：
  1. 按 SOFA 的约定（二阶积分器 + `LinearSolverConstraintCorrection`），**GS 的 λ 是冲量，力 = λ/dt**；`getLambda()` 存的也是冲量 → **插件 `InsertionAlgorithm.cpp:191` 除以 dt 是对的**，从组织一侧看阈值就是力（1600 N），刺穿时的压入深度不随 dt 变化。
  2. **官方场景的一阶针把 λ 当作力** → 针受到的力 = 组织受到的力 × dt，**违反作用力与反作用力定律**（相差 1/dt 倍）。之前测到的"刺穿力 ≈ 阈值 × dt"其实是**针根一侧**的读数。
  3. **后果**：①针根的力（机器人或力传感器测到的力）被缩小了 dt 倍，不能直接和实验比较；②针的变形也是在被缩小的力下算出来的（一阶针 E = 1e12，看不出区别；柔性针就会出错）；③本节以及问题 1 的表格里，力都是**针根一侧**的数值。
  4. **修正方法**：针和组织用一致的积分器（都用二阶，或者都用准静态、一阶，并相应处理 W 的因子）。二阶针要注意稳定性（针很刚、质量很小）。〔待验证：哪种组合最稳定〕
- 另外两个和文献不一致的地方仍然成立〔代码〕：①判定用的是针尖约束力各分量**范数之和**（包括切向摩擦），不是法向力；②用的是**上一步**的 λ（在 `AnimateBeginEvent` 时读取，见 §3），阈值没有进入约束律 → 刺穿的那一步力可能超过阈值。

### 问题 3 ⚠️ 3×3 柔度块求逆失败时的未定义行为（上游 master 已部分修复）〔代码 + 运行〕
- 代码：`InsertionResolution.h` 的 `init()`：`SOFA_UNUSED(invertMatrix(invW, temp))`。SOFA 的 `invertMatrix` 在 **|det| ≤ 2.2e-16（绝对值）**时直接返回 false（`Mat.h:1078-1085`，`equalsZero` 用的是机器精度）；`Mat` 的元素类型是 `VecNoInit`（`Mat.h:66,90`；`Vec.h:708-713`），**不做初始化** ⇒ 失败时 invW 是未初始化的内存 ⇒ 行为未定义（可能轴向摩擦消失，也可能是随机值）。
- 运行：官方场景里插入约束块的 |det| 约为 2e-9（W 对角元约 6e-4 到 2e-3 m/N）→ **不会触发**。块的对角元小于约 6e-6 m/N 时会触发（组织更硬、单元更小、时间步更小）。〔推导〕
- 刺穿前的单边约束块 |det| 约为 1e-22，但单边约束的求解代码不求逆，不受影响。〔运行 + 代码〕
- 上游 master 的 `scaleComplianceMatrix` 选项先缩放再求逆，但**失败时仍然没有处理未初始化的问题**（`InsertionResolution.h` master 版）〔待查细节〕。

### 问题 4 拔出检测：删除约束点时没检查针的运动方向（上游 master 已修复）〔代码〕
- `NeedleOperations.cpp`（v25.12）：只要约束点在针尖边方向的前方就删掉，**插入过程中也会删**（比如针弯曲、组织变形导致约束点跑到前方）。返回值写反了，但调用处（L271）不用返回值，**不影响结果**。master 改为只在针尖速度背离针尖边方向时才删。

### 问题 5 带摩擦的单边约束 GS 局部求解有 bug（ICube 原版已修，fork 没修）〔代码〕
- `UnilateralResolution.h:61-90`：截断后传播的修正量不对；切向截断发生在耦合传递之后。和 ICube 提交 `3e83cc4`（Claire Martin，2025-12-18）对比确认。
- 在官方场景里**不会产生影响**：`ConstraintUnilateral` 没有设 mu，mu = 0 时切向力都被截断为 0；maxForce0 默认是无穷大，法向不会被截断。**只要设了 mu > 0 就会有影响。**

### 细节（不是 bug）
- 轴向方向 = 针身投影点所在边的方向（`FirstDirection` + `EdgeNormalHandler`），**沿针身分段为常数** → 投影点跨过针节点时方向会跳变（Martin 第 5 章指出会引起力的尖峰）。〔代码〕
- 横向两个方向用 Gram-Schmidt 从 z 轴生成（`ConstraintNormal.h:57-75`）；针沿 z 方向插入时，它们在平面内可能随机旋转，但两个横向约束都是双边的，张成的平面不变，**不影响结果**。〔代码〕
- `getSize()` 返回节点数，所以 `begin(getSize()-2)` 取到的是针尖所在的边（L268 正确）。〔代码〕

## 2b. 上游 master 和 ICube 原版（A6）〔代码 + 编译 + 运行〕
- **master（InfinyTech3D，2026-08）可以在 SOFA v25.12 上编译**：`~/sofa/resources/build_master` → `install_master`（Ninja，Release；CollisionAlgorithm 3 个警告，ConstraintGeometry 0 个）。
- master 相对 v25.12 的源码改动很小：拔出方向检查、`contactDistance`、柔度块缩放求逆、去掉几个 proximity 类的导出宏；其余是文档。
- **运行对比**（`run_exp4.py`，master 的 examples 场景 + master 插件）：基础场景的摩擦实验（μ = 0.002 / 0.02，tol = 1e-5 / 1e-8）**结果和 v25.12 完全相同**；带摩擦（0.002）的反复插拔**同样失败**（组织被拖下 33 mm、往上拽 16 mm，最后剩 3 个约束点）→ master 的拔出修复解决不了问题 1b。
- ICube 原版：通用工具包，不含针插入；值得借鉴的是 ① `UnilateralResolution` / `BilateralResolution` 的 GS 修复（问题 5）；② **HermiteEdgeGeometry**（2025-07，半成品）：用梁节点的 Rigid 朝向构造沿边的三次 Hermite 曲线，中心线和切向都连续 → 可以消除"约束方向在针节点处跳变"的问题（Martin 第 5 章）；`getVelocity` 还没实现。

## 2c. 上游 issue 和 PR 里的相关信息（D2，2026-09-30）〔网页〕
- **#60（开着）"NEEDLE-3 针插入算法的状态和功能"**：开发者自己勾选完成的包括"刺穿阈值""反复插拔""穿过多个边界或层，或者碰到不可穿透的层""多层插入场景（多个物体作为层）"；**没完成的**：力反馈超过阈值时关闭碰撞、算法说明、类关系、针模型说明、代码注释。
- **#68（分层场景的 PR）：没有合并**（改动文件数为 0，分支已经指向 master 的历史提交）→ **公开仓库里没有多层示例**。
- **#56（开着）**：刺穿后不再做三角形碰撞检测 → **针不能从组织另一侧穿出**。
- **#57（已关闭）**：刺穿第二层或碰到不可穿透物体——但 v25.12 代码里一个算法实例仍然只有一个 `surfGeom` 和一个 `volGeom` → 多个物体要用多个 `InsertionAlgorithm` 实例〔推断〕。
- **#109**：针弯曲时碰撞边不跟着弯 → 开发者回复是场景配置问题（碰撞子节点里用了 `BeamLinearMapping`）；官方示例用的是 `IdentityMapping`（Rigid → Vec3）。
- **CG#8**："区分有摩擦和无摩擦的插入求解"（已关闭，没有正文）。
- **#106（开着）**：考虑把两个仓库合并。
- SOFA 论坛的两个针插入帖子（2020、2021）：官方人员都推荐 Cosserat 插件，没有更多技术细节。Cosserat issue #32（视网膜手术针插入）：作者 Adagolodjo 确认可以用，但当时不支持 SOFA 19。

## 3. SOFA 约束流水线（和 handbook 的对应）
| handbook 里的步骤 | SOFA 组件和代码 | 结论 |
|---|---|---|
| §1.4 第 3 步，组装 W | `LinearSolverConstraintCorrection::addComplianceInConstraintSpace`：W += J A⁻¹ Jᵀ × `correctionFactor`（`LinearSolverConstraintCorrection.inl:196-207`） | 积分因子 = dt（EulerImplicit）→ 对二阶物体，W = dt·JA⁻¹Jᵀ 的量纲是 m/(N·s)，**λ 是冲量**（N·s），力 = λ/dt〔代码 + 运行，见问题 2〕；一阶物体会把同一个 λ 当作力 |
| §6.1 GS | `BlockGaussSeidelConstraintSolver::doSolve`：按块调用 `ConstraintResolution::resolution()`；停止判据 Σ‖W_block Δf‖ < tol（默认 × 约束个数）；**没有热启动** | 〔代码〕 |
| GS 参数 | `tolerance`（默认 1e-3）、`maxIterations`（默认 1000）、`sor`、`regularizationTerm`、`scaleTolerance`（默认 true）、`allVerified`、`computeConstraintForces` | 〔代码〕 |
| 读取 W 和 λ | Python：`from Sofa import SofaConstraintSolver`，在 `onBuildConstraintSystemEndEvent` 里调用 `solver.W()`；`computeConstraintForces=True` 之后读 `constraintForces` | 〔运行〕 |

## 4. 对项目的影响
1. **插件的算法框架可以用**（穿刺 → 插入的状态机、固定间距的物质约束点、重投影、单个体网格内可以跨区域），和文献里的 Duriez / Martin 模型一致。
2. **插件的力学参数化不能直接用**：摩擦（问题 1）和物理量对不上，而且依赖求解器设置；官方场景的一阶针会让针根的力缩小 dt 倍（问题 2）。这正好对应 Manish 的"组织破裂"和我们对"能用实验数据标定"的要求。**至少要改两处**：
   - 穿刺判定：保留 `/ dt`（符合 SOFA 约定），改成用法向力判断；场景里针和组织用一致的积分器；
   - 摩擦律：换成有物理量纲的形式（handbook §4），在 `InsertionResolution::resolution()` 里实现。
3. **还缺的功能**：切割力（刺穿后针尖的阻力）、层间界面穿刺、斜面针尖，都要在这个框架上加（handbook §5、§9）。
4. **版本选择建议**：以 v25.12 tag 为基础，挑选上游的修复（拔出方向检查、求逆缩放）以及 ICube 的 GS 修复，再加上我们自己的修改。**不建议直接用 master**：master 是为 SOFA v26.06 准备的，和我们的 v25.12 二进制可能不兼容〔待查〕。

## 5. 还没读的部分
- `TBaseConstraint` 里 H 的组装和 `getConstraintViolation` 的符号约定（问题 1 的结论不依赖这些细节）。
- `FindClosestProximity`、AABB 宽相位的性能特性（网格变大以后才重要）。
- `PrecomputedConstraintCorrection`、`EulerImplicitSolver` 的积分因子（性能阶段再看）。
