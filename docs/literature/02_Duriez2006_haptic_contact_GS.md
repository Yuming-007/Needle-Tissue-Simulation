# Duriez et al. 2006 — Realistic Haptic Rendering of Interacting Deformable Objects in Virtual Environments

- **出处**：IEEE TVCG 12(1), 2006, pp. 36–47。作者：C. Duriez, F. Dubois, A. Kheddar, C. Andriot（CEA/LIST, LMGC Montpellier, AIST）
- **本地文件**：`~/sofa/papers/original/Realistic_haptic_rendering_of_interacting_deformable_objects_in_virtual_environments.pdf`
- **阅读状态**：全文精读（含算法 1、图 1–18）；**2026-09-29 对照原文逐条核对，没有发现错误**
- **在本项目中的地位**：**数值基础论文**。给出了 SOFA 约束求解的基本框架：Signorini + Coulomb 摩擦律、Delassus 算子、按接触分块的非线性 Gauss–Seidel 求解器。Duriez 2009 的所有针约束都建立在它的上面。

---

## 1. 要解决的问题
- 多个可变形体之间**多点接触 + 干摩擦**时的实时求解，要求能稳定地输出 haptic 力（500 Hz – 1 kHz）。
- 当时常用的方法各有缺陷：
  - penalty 法：力取决于人为设定的刚度，刚度大时方程变刚、不稳定；多个接触时是分段线性的组合问题（每个接触两种状态，n 个接触就有 2ⁿ 种组合）。
  - imposed motion（双边约束 / 主从）：会"粘住"，也难处理摩擦。
  - 摩擦锥用 k 边棱锥近似成 LCP：问题规模变大，精度也差。

## 2. 核心贡献
1. 在实时可变形体仿真里首次使用 **Signorini 接触律 + Coulomb 摩擦律**，不做近似。
2. 在接触空间里把力学行为线性化 → **Delassus 算子 W**。
3. 用**按接触分块的非线性 Gauss–Seidel** 求解（每块 3×3：一个法向 + 两个切向），**不需要**用棱锥近似摩擦锥。
4. 用**全局共旋**方法把质量比和刚度比从时间步长里解耦出来，haptic 时可以用任意时间步长保持稳定。

## 3. 模型与公式

### 3.1 形变模型与线性化（§2.1–2.2）
- 离散化后：$M\ddot u + D\dot u + K(u) = f$（FEM 或质点弹簧都是这个形式）。
- 线性化 + 隐式积分后：$\tilde K u = \tilde f$（式 2、3）。$\tilde K$ 依赖于材料、积分格式和线性化方式。
- **柔度** $C = \tilde K^{-1}$（也叫 compliance 或 capacitance）。线性模型可以预先算好，并压缩到只含接触节点。

### 3.2 Signorini 律（§3）
- 物体 D1 上的点 P 和 D2 上的点 Q，法向 n，间隙 $\delta_n = \vec{QP}\cdot n$。
- 互补条件：$0\le\delta_n \perp f_n \ge 0$（式 6、7）。
- P、Q 都在三角形上，用线性插值：$U_P=\sum\Psi_\alpha(P)U_1(\alpha)$；力也用同一组形函数分配到节点（式 8–12）。
- 碰撞检测只需要提供 **P、Q、它们的重心坐标、法向**（§3.2）。
- 一次三角形–三角形碰撞可能产生多个接触点（Fig. 3）。

### 3.3 LCP 形式（§3.4）
- **自由运动**：忽略接触力得到的位置，对应间隙记为 $\delta^{free}$。
- $\delta_n = n^T(U_P-U_Q)+\delta_n^{free}$，写成矩阵：$\delta_n=[H_1]U_1-[H_2]U_2+\delta^{free}$（式 14）；力 $F_1=H_1^Tf_n$，$F_2=-H_2^Tf_n$（式 15）。
- $P = C_iF + P^{free}$（式 16），于是
  $$\delta_n = \Big(\sum_i H_iC_iH_i^T\Big)f_n + \delta_n^{free},\qquad W := \sum_i H_iC_iH_i^T\ \text{（Delassus 算子）}$$
  没有摩擦时就是一个 LCP（式 18）。
- **接触组（contact group）**：通过接触连在一起的物体一起求解，不相连的组可以分开解。

### 3.4 Coulomb 摩擦（§4）
- $\delta_t = 0 \Rightarrow \|f_t\| < \mu\|f_n\|$（粘着，stick）；$\delta_t \ne 0 \Rightarrow f_t = -\mu\|f_n\|\,\delta_t/\|\delta_t\|$（滑动，slip）（式 19）。
- δ_t 是一个时间步内切向速度的积分，也就是切向位移。
- Fig. 4 是**纯 Coulomb**：粘着时力在 ±μf_n 之间，滑动时力恒为 ±μf_n，**没有斜坡**。
  → 所以 Duriez 2009 Fig. 4 里滑动段的斜坡是另外加的（大概是和速度相关的项），不来自这篇文章。
- 接触空间的完整方程（式 20）：3×3 块，$[\delta_n;\delta_t] = [W_{nn}\,W_{nt};W_{tn}\,W_{tt}][f_n;f_t]+\delta^{free}$。
- 3D 里真正的难点是：**滑动方向事先不知道**，而且法向和切向是耦合的（式 21）。
- **k 边棱锥近似**（式 22–24）：把摩擦锥换成 k 边多面锥 → 得到规模为 m×(k+2) 的 LCP，可以用 Lemke 解。k = 8 时规模是无摩擦时的 10 倍，要 k ≥ 16 误差才小于 5%（Fig. 7）。

### 3.5 Gauss–Seidel 类算法（§5，算法 1）
- 对接触 α：$\delta_\alpha - [W_{\alpha\alpha}]f_\alpha = \sum_{\beta<\alpha}[W_{\alpha\beta}]f_\beta+\sum_{\beta>\alpha}[W_{\alpha\beta}]f_\beta+\delta^{free}_\alpha$（式 25），其他接触的力"冻结"，$[W_{\alpha\beta}]$ 是 3×3 块。
- 局部问题是非线性的（Signorini + Coulomb）。可以把 $W_{\alpha\alpha}$ 换成对角近似，再用"图求交"来解（式 25'）。
- **算法 1 的要点**：
  1. 预处理：$\Lambda_i = (\lambda_{min}+\lambda_{max})/2$，取 $[W_{ii}]_{tt}$ 特征值的平均，作为切向的标量柔度。
  2. 每次迭代，对每个接触 i：算出 $\delta^{free}+\sum_j W_{ij}f_j$（用最新的 f）。
  3. **法向**：$f_n \leftarrow f_n - \delta_n/W_{ii}(1,1)$。
  4. 如果 $f_n>\epsilon_1$：**切向** $f_t \leftarrow f_t - \delta_t/\Lambda_i$；若 $\|f_t\|>\mu f_n$，就把 f_t 投影回锥面：$f_t \leftarrow \mu f_n\,f_t/\|f_t\|$。
  5. 否则 $f=0$（脱离接触）。
  6. 收敛判据：$\sum\|f^{k}-f^{k-1}\|/\|f^k\| < \epsilon_2$。
- **复杂度**：每次迭代 O(m·m')，m' 是非零接触数。棱锥 LCP 最坏是 O(k·m²)。
- **性能**（Fig. 8，Matlab）：接触多了以后，GS 远快于 8 边棱锥 LCP（50 个接触：约 0.5 s 对约 9.5 s）。
- **收敛性**：GS 的收敛速度取决于 W 的**对角占优**程度。3D 线弹性 FEM 通常满足这个条件。tolerance ε₂ 是速度和精度之间的调节旋钮。

### 3.6 接触之间的力学耦合（§5.3）
- W 的非对角块体现接触之间通过物体传递的耦合（Fig. 9、11：桌子同一条腿上的接触耦合强，不同腿之间弱；加了横撑之后耦合变强）。
- **penalty 法和"单接触"方法都忽略了这种耦合**，所以处理不了真正的多点接触。
- 针插入时沿针身有几十个约束点，它们通过组织和针强烈耦合 → 正是 W 必须建准确的原因（Martin 论文后面的重点）。

### 3.7 W 的计算（§5.4）
- 线弹性（小位移）：$\tilde K$ 压缩到表面节点，预先分解。运行时只需要把坐标系从物体转到接触点，并做插值，复杂度 O(m²)。
- 如果 W 算起来太慢：可以只算**块对角**部分，非对角的"冻结贡献"用 $\tilde K^{-1}H^TF$ 在物体空间里算完再映射回来。这就是 SOFA 里 `LinearSolverConstraintCorrection` 和 `UncoupledConstraintCorrection` 等不同实现的思路来源。

### 3.8 全局共旋与时间步长（§6）
- 隐式 Euler（式 27）：$(M/\Delta t^2 + D/\Delta t + K)u = \ldots$。刚度越大，就必须用越小的时间步才能保持质量项占优，这和 haptic 要求的固定小时间步相冲突。
- **全局共旋**：物体运动 = 刚体运动 + 局部线性变形（Fig. 12）。刚体部分用 $J_c A^{-1} J_c^T$ 映射到接触空间（式 29），整体为
  $$\delta = \Big[\sum_i H_iC_iH_i^T + J_c\big(\tfrac{A}{\Delta t^2}\big)^{-1}J_c^T\Big]f+\delta^{free}\quad(30)$$
- 刚度趋于无穷时，模型自然退化成刚体运动，因为质量和刚度在 W 里是解耦的，所以 haptic 可以稳定输出。
- "准刚体"：刚体摩擦接触的解不唯一；加一点 FEM 柔度之后解变得唯一、平滑（Fig. 13）。
  → **对本项目的启示**："刚性针"最好不要用纯刚体，而是用刚度很大的可变形梁（CollisionAlgorithm 的例子里 E = 1e12，就是这个思路）。这样约束问题更良态。这是我的推论，要在实验里验证。

### 3.9 Haptic 耦合与实验（§7）
- 6D 虚拟耦合（Adams & Hannaford）：设备一端是阻抗（impedance）控制，仿真一端是导纳（admittance）控制，中间用 6 自由度的弹簧阻尼连接。
- 卡扣（snap-in）任务：推入 → 不稳定平衡 → 卡入，三个阶段。夹子材料 E = 700 MPa，时间步 3 ms。
- 性能：30 个无摩擦接触约 3 ms；20 个有摩擦接触约 4 ms（μ = 0.2 或 0.8 差别不大）。haptic 刷新率最高 250 Hz。

## 4. 局限
- 只处理了表面接触（Signorini），没有体内约束，也没有穿刺。
- 要快，就得依赖线性或共旋模型加预计算。
- Gauss–Seidel 在 W 对角占优差、约束很多、双边约束和摩擦混在一起时收敛慢（Martin 2023 改进的正是这一点）。
- 性能数据来自 2006 年的 Matlab 实现，只能看趋势。

## 5. 和本项目的关系
- **针刺穿之前**的针尖–表面接触，就是 Signorini + Coulomb，外加穿刺阈值 f_p。本文的算法 1 就是 CollisionAlgorithm 里 `ConstraintUnilateral` 的理论来源（待看代码确认）。
- 理解 SOFA 里 `FreeMotionAnimationLoop` → `ConstraintCorrection`（算 W）→ `ConstraintSolver`（GS）这三件事分别在做什么，这篇是必读。
- **调参时要直接用到**：GS 的 `tolerance`（ε₂）、`maxIt`、μ，以及 W 的对角占优程度（针和组织刚度相差悬殊时，W_n 和 W_t 的量级差别很大 → 可能影响收敛）。
- 刚性针用"刚度很大的可变形体"比纯刚体更稳定，原因见 §3.8。

## 6. 回答 Duriez 2009 笔记里的疑问
- Q2（摩擦图滑动段的斜率）：本文是纯 Coulomb，没有斜坡 → 斜坡是 Duriez 2009 另加的，最可能是速度相关（粘性）项。还要看 Martin 论文里有没有进一步说明。

## 7. 待解决的问题
- ~~SOFA 的 `BlockGaussSeidelConstraintSolver` 和 `ProjectedGaussSeidelConstraintSolver` 与算法 1 有什么区别？~~（已回答：v25.12 只有 BlockGS 等 5 个求解器；BlockGS 按块调用 `ConstraintResolution::resolution()`，结构和算法 1 相同，具体的局部求解由约束类型决定，见 `docs/code/01`、`02`）
- 针（梁）和组织（四面体）刚度相差几个数量级，Λ_i 用特征值平均来近似，在这种情况下是否还合理？

---
## 核对记录（2026-09-29）
- 核对了式 (1)–(30)、算法 1 的每一行、图 7/8/17/18 的数值、§7.2 的参数（E = 700 MPa，ν = 0.35，15 g，3 ms）、性能（30 个无摩擦接触约 3 ms，20 个有摩擦接触约 4 ms，最高 250 Hz）：**都和原文一致**。
- 补充一个细节：原文 Signorini 律的原始形式是用应力写的，$0\le\delta_n(P)\perp\sigma_n^{(1)}(P)\ge0$（式 6）；用线性形函数插值之后，才等价于用节点力写的形式（式 7）。
- 补充：式 (27) 的隐式 Euler 写的是 $(M/\Delta t^2+D/\Delta t+K)\,u_t=\ldots$，**未知量是位移 u**；而 Martin 论文和 SOFA 的写法以**速度增量**为未知量，$(M+hB+h^2K)\,\Delta v=\ldots$。两者相差一个 h² 因子，W 的量纲也就不同（前者是柔度 m/N，后者要再乘 h² 才是位移/力）。**以后比较 W 的数值时要注意是哪一种约定。**
