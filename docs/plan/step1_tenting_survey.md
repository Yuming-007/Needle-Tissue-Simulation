# 第 1 步之前：tenting 的已有实现方式调研与比较（2026-09-30）

- 起因：用户要求在动手写代码之前，先调研别人怎样实现 tenting，比较哪些和本项目更符合，并考虑和后面几步的衔接（reuse 原则）。
- 来源：本项目 20 篇文献笔记（`docs/literature/`，都已对照 PDF 核实）、本机的代码（InfinyTech3D 插件、Cosserat 插件、SOFA 本体），以及补充网页搜索（只看了摘要和页面描述，没有逐页读全文的标为〔摘要〕）。

## 1. 已有的实现方式

### A. 约束法：针尖一个点 ↔ 表面最近点，单边约束（Signorini）〔原文 + 代码〕
- 做法：针尖节点 Q 和它在组织表面上的最近点 P（用三角形的重心坐标插值到 3 个节点）之间建一个单边约束，法向取三角形法向；通常再加 2 个切向 Coulomb 摩擦约束；λ_n 超过刺穿阈值 f_p 时切换为"已刺穿"，不改网格拓扑。
- 文献：Duriez 2009（`01`）、Adagolodjo 2019（`14`，Signorini + Coulomb μ_s，阈值 p_f）、Perrusi 2021（`15`）、Baksic 2022（`20`）、Martin 2025（`04`，"针尖节点和它的表面投影，方向取三角形法向"）、Bui 2019（`11`，接触点用三角形形函数插值，Signorini + Coulomb，阈值 λ_p0）。**所有约束法的实时插针工作都用这种方式**（`00_synthesis.md`）。
- 本机代码，有两种现成实现：
  - **A1 InfinyTech3D CollisionAlgorithm + ConstraintGeometry**：`InsertionAlgorithm` 找针尖到外表面的最近点，交给 `ConstraintUnilateral`（mu 可选），超过阈值就把表面点存为第一个插入约束点，进入插入状态（`docs/code/01` §1）。**刺穿前 → 刺穿 → 插入 → 拔出的状态机是现成的**。已知问题：摩擦不是物理的（问题 1）、带摩擦的单边求解有 bug（问题 5，mu = 0 时不影响）、`ContactDirection` 在两点重合时方向可能不稳定（`step1_tenting.md` §1.4）。
  - **A2 SOFA 自带的碰撞流水线**（Cosserat 插针示例的做法）：`PointCollisionModel`（针尖）对 `TriangleCollisionModel`（组织表面），`LocalMinDistance` 做邻近检测（`contactDistance` 相当于接触厚度），`FrictionContactConstraint` 生成带摩擦的拉格朗日接触；刺穿由 Python 控制器判断 `constraintForces[0] > 阈值` 后关闭接触、手动加插入约束点（`Cosserat/examples/python3/NeedleInsertion.py`、`cosserat/needle/needleController.py`）〔代码〕。注意：这个示例直接拿 λ（冲量）和阈值比较，和我们"力 = λ/dt"的约定不一致。
- 优点：
  - 精确满足"不穿透、只推不拉"，没有人为的接触刚度；
  - 多点接触时耦合正确（Duriez 2006）；
  - 和我们已经验证过的约束框架、纯运动学针、读力方法（0b）完全一致；
  - 后面的刺穿、插入约束、摩擦都在同一个框架里。
- 缺点：针尖是一个点，刚度由网格决定（`step1_tenting.md` §1.3）。

### B. 罚函数接触（penalty）〔原文 + 摘要〕
- 做法：穿透深度 × 人为的接触刚度 = 接触力。例如 SOFA 论坛上那个插针例子（`LocalMinDistance` + 罚函数响应，只调了接触参数）；Bui 2018 用罚权重施加约束（`10`）。
- 缺点：力取决于人为刚度，刚度大时方程变刚、不稳定（Duriez 2006、Martin 2025 都指出了这一点）；会有微小穿透。和我们"不引入人为参数"的原则不符。

### C. 让网格节点落在针上（贴合网格、节点吸附、重新划分网格）〔原文 + 摘要〕
- 做法：针尖或针身直接和网格节点重合。DiMaio 2003（2D，`08`）；Goksel 2005 / 2006（节点重定位、增加节点，用 Woodbury 公式加速，触觉频率 > 1 kHz）；Chentanez 2009（局部重新划分网格、节点吸附，`09`）。
- 优点：针尖处的变形在节点上精确表示，不需要重心插值。Chentanez 认为针周围的变形梯度不连续，非贴合方法会把力"抹开"（`09`）。
- 缺点：要在运行中修改网格（拓扑、刚度矩阵、分解都要更新），实现量很大；和我们的异步分解（0d）冲突；插件没有这个功能。

### D. 离线的详细有限元：针尖几何 + 断裂模型〔摘要〕
- 做法：Abaqus 等软件，针尖几何（直径、斜面）完整建模，组织用内聚力单元（cohesive zone）或断裂力学（Mahvash & Dupont 2010：J 积分，能量释放率超过断裂韧性时裂纹扩展）描述刺穿。
- 结论（摘要层面）：针的直径对反力影响最大；刺穿力、刺穿时的位移、刺穿后的力跳变都随直径增大而增大。Mahvash & Dupont 还指出插入速度越快，刺穿力越小。
- 优点：物理上最完整（针尖面积、刺穿机理）。
- 缺点：网格要细到针尖尺寸以下，远不能实时。可以作为物理参考，不能直接复用。

### E. 针尖用多个点表示一个小面积（文献里没有看到，是我们可能的扩展）〔推论〕
- 做法：在针尖端面（例如半径 0.6 mm 的圆）上放一圈点，用 `RigidMapping` 挂在针上，每个点各建一个单边约束，GS 联合求解。相当于一个小的平底压头。
- 只有当刺入点附近的表面单元和针尖差不多大（≲ 1 mm）时才有意义。单元远大于针尖时，这一圈点都落在同一个三角形上，雅可比几乎线性相关，效果和一个点没有区别，还会让 W 接近奇异（`docs/code/04`：约束间距 ≪ 单元尺寸时 GS 不收敛）。

## 2. 比较（针对本项目）

| 方式 | 物理正确性（刺穿前） | 实时 | 能复用的现成实现 | 和我们已验证的框架是否一致 | 和第 2–8 步的衔接 | 网格依赖 |
|---|---|---|---|---|---|---|
| **A1 插件（单边约束）** | 精确的接触律 | 是 | **有，含完整的刺穿 / 插入状态机** | 一致 | **最好**：第 2–6 步直接沿用它的状态机，物理部分（阈值判据、摩擦、切割）替换成我们的 | 有（点接触） |
| A2 SOFA 自带碰撞流水线 | 精确的接触律 | 是 | 有（只有接触；刺穿和插入要自己写，如 Cosserat 示例的 Python 控制器） | 一致 | 一般：刺穿以后的部分都要自己写 | 有 |
| B 罚函数 | 依赖人为刚度 | 是 | 有 | 不一致 | 差 | 有 |
| C 贴合网格 / 重新划分 | 针尖处最精确 | 难（要更新分解） | 没有 | 冲突（异步分解） | 工作量大 | 小 |
| D 离线详细有限元 | 最完整 | 否 | 没有 | — | 只作物理参考 | 小 |
| E 多点针尖 | 网格够细时更接近真实 | 是 | 可在 A1 基础上加 | 一致 | 好 | 网格够细时减小 |

## 3. 建议
1. **第 1 步采用 A1**（InfinyTech3D 插件的 `InsertionAlgorithm` + `ConstraintUnilateral`，刺穿阈值设为无穷大）。理由：
   - 物理上和所有约束法文献一致；
   - 和第 0 步验证过的框架一致；
   - **它的状态机正是第 2–6 步要用的**，第 1 步就用它，后面不用换实现。
2. **用 A2 作为交叉验证**：在 P1（顶点接触）里，同一个工况同时用 SOFA 自带的接触流水线跑一遍。两种独立实现应该给出相同的力–深度曲线，这比只和"给定位移"比较多一层独立检查。
3. **针尖面积（E）暂不加**：在 h = 10 mm（甚至加密后的 2–3 mm）上，针尖仍比单元小很多，多点针尖和一个点没有区别。到加密网格的步骤再评估，并在那时决定"针尖刚度的网格依赖"怎么处理（加密 + 按实验数据标定，或者用能量 / 断裂判据定义刺穿，参考 Mahvash & Dupont 2010），这直接影响第 2 步。
4. **注意的已知问题**：`ContactDirection` 的方向稳定性（P6），法向朝向和正负号（P3），λ 和力的换算（Cosserat 示例的单位问题提醒我们：阈值必须按力 = λ/dt 定义）。

## 4. 来源
- 文献笔记：`docs/literature/01, 02, 04, 08, 09, 10, 11, 14, 15, 20, 00_synthesis.md`。
- 代码：`~/sofa/reference_sources/ConstraintGeometry/.../ConstraintUnilateral.h, ContactDirection.h, SecondDirection.h, UnilateralResolution.h`；`~/sofa/resources/SofaDefrost_Cosserat/examples/python3/NeedleInsertion.py`、`cosserat/needle/needleController.py`。
- 网页（摘要层面）：Chentanez 2009（https://people.eecs.berkeley.edu/~jrs/papers//needlesim.pdf）；Goksel 2005（https://people.ece.ubc.ca/orcung/pubs/Goksel05_miccai.pdf）；Mahvash & Dupont 2010（https://dash.harvard.edu/bitstreams/7312037e-7e37-6bd4-e053-0100007fdf3b/download）；内聚力有限元相关（https://www.researchgate.net/publication/221733179 、https://link.springer.com/content/pdf/10.1007/s10237-020-01310-x.pdf）；Misra 2008（http://reedlab.eng.usf.edu/publications/misra2008needle.pdf）。
