# Duriez et al. 2009 — Interactive Simulation of Flexible Needle Insertions Based on Constraint Models

- **出处**：MICCAI 2009（LNCS），8 页。作者：C. Duriez, C. Guébert, M. Marchal, S. Cotin, L. Grisoni（INRIA Lille / Rennes，SOFA 核心团队）
- **本地文件**：`~/Downloads/Interactive Simulation of Flexible Needle Insertions Based on Constraint Models.pdf`
- **阅读状态**：全文精读（含图 1–8）；**2026-09-29 对照原文逐条核对**，改正 2 处（见文末核对记录）
- **在本项目中的地位**：**核心方法论文**。SOFA 中针–组织交互的"标准模型"，CollisionAlgorithm / ConstraintGeometry 插件与 Martin 2023/2025 都是它的延续。

---

## 1. 要解决的问题

- 细长柔性器械（针、电极、活检针）插入软组织的交互式仿真，用于训练和**术前规划**。
- 以往 FEM 方法（DiMaio、Alterovitz）需要沿针道**重划网格**，难以实时。
- 需同时表达的非光滑现象：穿刺（puncture）、切割（cutting）、静/动摩擦、针尖导向（bevel steering）、非均质组织。

## 2. 核心思想

**用互补约束（complementarity constraints）+ Lagrange 乘子描述所有针–组织交互，且约束点可放在组织体内任意位置 → 不需要重划网格。**
方法与针、组织的具体力学模型无关（只要能提供映射 J 和柔度 W）。

## 3. 模型与公式

### 3.1 约束定位（§2.1，Fig. 1）
- 一条约束由两点定义：组织体内点 P、针曲线上点 Q，以及约束方向 **n**；δ 为 P、Q 沿 n 的距离。
- P 在四面体内用重心坐标插值：$\mathbf u_t = \mathbf J_t \Delta\mathbf q_t$；Q 用针模型插值：$\mathbf u_n = \mathbf J_n \Delta\mathbf q_n$。
- 约束违背量变化（式 1）：
  $$\Delta\delta = \mathbf n^T(\mathbf u_n-\mathbf u_t) = \mathbf H_n\Delta\mathbf q_n + \mathbf H_t\Delta\mathbf q_t$$
- 由虚功原理得到约束力映射（式 2）：$\mathbf f_t = \mathbf H_t^T\lambda,\ \mathbf f_n=\mathbf H_n^T\lambda$。λ 为组织作用于针的力。
- ⚠️ 原文式 (1) 写 $-\mathbf J_t$ 后得 $+\mathbf H_t$，说明 $\mathbf H_t$ 已吸收负号（$\mathbf H_t = -\mathbf n^T\mathbf J_t$）。

### 3.2 穿刺约束（§2.2，Fig. 2）——三状态互补律
Q=针尖，P=组织表面接触点，n=表面法向，f_p=穿刺阈值：

| 状态 | δ_p | λ_p | 物理含义 |
|---|---|---|---|
| step 1 接近 | δ_p ≥ 0 | λ_p = 0 | 未接触 |
| step 2 接触 | δ_p = 0 | 0 ≤ λ_p ≤ f_p | 表面被推动变形，未破 |
| step 3 穿刺 | δ_p ≤ 0 | λ_p = f_p | 针尖进入，力等于阈值 |

- 在 (δ, λ) 图上这是一条**单调**的阶梯图（δ 减小时 λ 不减），所以与 Gauss–Seidel 局部直线的交点唯一（见 §4）。
- **每层组织可以设不同的 f_p**，同一个约束在穿过多层时可以多次触发。
- f_p 取极大值 → 永远刺不穿 → 针沿表面滑动（模拟骨头）。
- 耦合：step 2 时在两个切向加**摩擦约束**；step 3 后切换成**针尖路径约束**（§3.4）。
- **对本项目的意义**：Manish 说的"组织破裂（rupture）"在这个框架里就是 step 2 → step 3 的状态切换，不涉及拓扑改变。

### 3.3 切割约束（§2.3）
- f_c：穿过一种组织结构需要的力，每层可以不同。
- 形式和穿刺约束类似，只是 δ_c 衡量的是**针尖相对"已有切割路径末端"的位移**。
- **沿原路径重新插入时没有切割力**（路径已经存在）→ 解释了 Fig. 6 中 (5) 重新插入时力比较小。
- 〔我的推论，原文没有写，Fig. 6 里也看不到〕如果 f_c < f_p，从穿刺切换到切割时针尖力应该会突然下降，对应 Okamura 2004 观察到的"穿刺峰值后骤降"。原文既没有给出 f_p 和 f_c 的数值，也没有讨论这个切换。
- ⚠️ 原文只有一段文字，**没给出 δ_c 的精确定义、方向，也没给出和摩擦、路径约束的组合方式**。实现细节要看代码（CollisionAlgorithm 是否实现了切割力，待核实）。

### 3.4 针尖路径与导向（§2.4，Fig. 3）
- 针尖有一个方向，用来约束针尖的**横向运动**：针尖只能沿"切割方向"前进（δ_t = 0，双边等式约束）。
- **斜面针尖导向（bevel steering）**：在针尖坐标系里把切割方向设成倾斜的，于是推针和拉针时约束方向不同（Fig. 3a/b）→ 针走曲线。
- **针身跟随针尖轨迹**：沿针身，在垂直于针轴的方向上令针和组织之间的相对位移为零（Fig. 3c），λ_t 就是维持这个等式所需的力。
- ⚠️ 原文说的是 "along the tangential directions of the needle curve"，但从 Fig. 3(c) 看 δ_t 是**垂直于针**的。理解为：针轴的两个法向方向用双边约束，轴向交给摩擦（和 Martin 论文 §3.3.4.3 的描述一致）。

### 3.5 摩擦（§2.5，Fig. 4）
- 两个状态：粘着（stick，δ_f = 0）和滑动（slip，δ_f ≠ 0）。
- 粘着的阈值是 μ·p，其中 μ 是摩擦系数，p 是组织对针的压力。**p 由组织刚度估算，是常数**（作者自己承认这是局限，计划以后用 FEM 应力来算）。
- Fig. 4 的滑动段是倾斜的（标了角度 α）：滑动时阻力随滑动量线性变化。因为每步的 δ_f 近似等于速度乘时间步，所以这相当于**"静摩擦 + 和速度相关的动摩擦"**。原文没解释 α，这是我从图上推出来的。
- 把阻力 r 沿已插入针身积分：每个约束点负责一段长为 l 的针身，
  $$\lambda_f = l\,\pi d\, r\quad(d\text{ 为针直径})$$
  → 摩擦力和插入深度成正比（Fig. 6 step 2 的线性上升）。

### 3.6 动力学与求解（§3）
- 针和组织各自满足（式 3、4）：$\mathbf M\dot{\mathbf v} = \mathbf p - \mathbb F(\mathbf q,\mathbf v) + \mathbf H^T\lambda$
- 时间积分用隐式后向 Euler，每步分三步走（也就是 SOFA 的 `FreeMotionAnimationLoop`）：
  1. **自由运动（free motion）**：先忽略约束力，求出每个模型的运动；
  2. **求约束力 λ**；
  3. **修正运动**：$\Delta\mathbf q = h\,\mathbf A^{-1}\mathbf H^T\lambda$。
- 约束空间方程（式 5）：
  $$\boldsymbol\delta = \Big[\underbrace{\mathbf H_n\mathbf A_n^{-1}\mathbf H_n^T}_{\mathbf W_n} + \underbrace{\mathbf H_t\mathbf A_t^{-1}\mathbf H_t^T}_{\mathbf W_t}\Big]\lambda + \boldsymbol\delta^{free},\quad \mathbf A = \tfrac{\mathbf M}{h^2}+\tfrac{d\mathbb F}{h\,d\mathbf v}+\tfrac{d\mathbb F}{d\mathbf q}$$
  W 就是 Delassus 算子（柔度矩阵）。
- **Gauss–Seidel 类求解**（式 6）：对每个约束 α，固定其他 λ_β，得到一条直线 $\delta_\alpha = W_{\alpha\alpha}\lambda_\alpha + \delta_\alpha^-$，和该约束的特征图求交点得到新的 λ_α（Fig. 5）。因为 $W_{\alpha\alpha}>0$，交点唯一。当 δ 的误差低于阈值时停止迭代。
- 这种写法的好处：**每种约束律只需要提供"特征图 + 求交"**，所以穿刺、切割、摩擦都能放进同一个求解器。

### 3.7 约束采样（§3 末）
- 约束点太多（超过相关自由度数）→ 过约束，收敛差；约束点太稀 → 精度差。
- 实际做法：约束点沿针曲线规则采样，和针的离散化一致；插入和拔出时**动态增删约束点**。

## 4. 实现细节与性能（§4.1）
- **针**：50 个梁单元，支持大位移（参数 E、ν、截面积 A，可以模拟空心针）。用带状三对角求解器算 W_n，只要几毫秒。
- **组织**：线弹性 Hooke 定律 + 共旋（大位移、小变形），ν = 0.49，Rayleigh 阻尼。
- **W_t 的快速计算**：用 compliance warping（Saupin 2008），也就是**预先算好柔度矩阵，运行时乘上旋转**。对应 SOFA 的 `PrecomputedConstraintCorrection`。前提是线性材料 + 小变形。
- **性能**：导向实验平均 28 FPS（2009 年的硬件）。

## 5. 实验（§4.2）
1. **反复插入和拔出**（Fig. 6，虚线是施加给针的运动）：(1) 刺穿表面时力**陡升**（图上是一段近乎竖直的上升，**没有出现"峰值后骤降"**）→ (2) 力随插入深度线性增加（摩擦）→ (3) 停止插入后力松弛下降 → (4) 部分拔出，力变为负值 → (5) 沿原路径再插入，力比第一次小（没有切割力）→ 反复几次 → (6) 完全拔出，负向摩擦力逐渐减小到零。作者说"和 Dehghan 2008 一致"，但**只是定性比较，没有数值**。
2. **避障导向**（Fig. 7）：先插入一半，旋转 180°，再插完。对应 Webster 2006，只是这里组织是可变形的。
3. **非均质组织**（Fig. 8）：两种刚度的组织。软针刺不穿硬区，沿表面滑开；加大针的刚度就能刺穿；也可以先插几根针让软区变硬，柔性针就能刺穿。**这就是 Manish 要的"两种组织的立方体"的原型实验**。

## 6. 局限（原文 + 我的分析）
- **没有定量验证**，也没给具体参数值（f_p、f_c、μ、p 都没有）。"用实验数据参数化"只是一句话。
- 组织只用线弹性共旋模型，限于小变形。预计算柔度要求刚度矩阵固定，**和非线性超弹材料、拓扑变化都不兼容**。
- 摩擦压力 p 是常数，和组织的实际应力无关。
- 切割约束的描述太简略，没法照着复现。
- 约束点只沿针身规则采样，**不一定和组织网格对齐**。Martin 2023 就是在这一点上改进的（用网格相交来定位约束点）。
- Gauss–Seidel 在约束很多、双边约束和摩擦混在一起时收敛慢。Martin 2023 用混合求解器解决这个问题。
- 没讨论多层之间的**内部界面**怎么检测：step 2 的"表面"在多层组织里包括层与层之间的交界面。

## 7. 和本项目的关系

| 项目需求 | 这篇论文给出的方案 | 备注 |
|---|---|---|
| 组织破裂 | 穿刺三状态律（f_p） | 力学层面的建模，不改拓扑 |
| 多层组织 | 每层不同的 f_p、f_c（μ、p 也应该分层） | 需要检测层间界面，原文没讲 |
| 刚性针 → 柔性针 | 约束和针模型无关，换成梁就行 | 刚性针可以用 E 很大的梁，或者单个 Rigid3 |
| 速度 | 预计算柔度 + 三对角求解 + Gauss–Seidel | 预计算柔度限制了材料模型 |
| 两种组织的立方体 | Fig. 8 就是原型 | 可以作为第一个复现目标 |
| 术前规划 / 闭环控制 | 结论里提到会用作规划工具 | Wang 2025 是这个方向的现代实现 |

**在 SOFA 里对应的组件（待核实）**：
- `FreeMotionAnimationLoop`
- `BlockGaussSeidelConstraintSolver`（v25.12 没有 ProjectedGaussSeidel）
- `LinearSolverConstraintCorrection` 或 `PrecomputedConstraintCorrection`
- `BeamFEMForceField`
- `TetrahedronFEMForceField` / `FastTetrahedralCorotationalForceField`
- CollisionAlgorithm 插件：`ConstraintUnilateral`（对应穿刺 step 2 + 摩擦）、`ConstraintInsertion`（对应横向双边约束 + 轴向摩擦）
- 待确认：切割约束、针尖路径约束、斜面导向有没有实现

## 8. 待解决的问题
1. 切割约束的 δ_c 到底怎么定义？CollisionAlgorithm 里有没有等价实现？
2. 摩擦图里滑动段的斜率 α 是什么物理量？（看 Duriez 2006 和 Martin 的论文能不能解释）
3. 多层组织里，穿过第一层之后，第二层的表面接触是怎么检测的？
4. 从穿刺（f_p）到切割（f_c）的切换，在数值上会不会引起振荡？

---
## 核对记录（2026-09-29，对照原文 PDF 逐页核对）
- 式 (1)–(6)、Fig. 2/4/5 的约束律、§4.1 的模型参数（50 个梁单元、ν = 0.49、Rayleigh 阻尼、compliance warping）、28 FPS：**和原文一致**。
- **改正 1**：Fig. 6 的描述。原来写的"出现穿刺峰值"不准确。原图在刺穿时是陡升，接着线性增长，没有峰值后骤降；力下降出现在 (3) 停止插入后的松弛阶段。
- **改正 2**：§3.3 里"f_p → f_c 切换时力会突然下降"是我的推论，原文没有写，现在已经标明。
- 补充的交叉对照：式 (5) 里的 $\mathbf A=\mathbf M/h^2+d\mathbb F/(h\,d\mathbf v)+d\mathbb F/d\mathbf q=(\mathbf M+h\mathbf B+h^2\mathbf K)/h^2$，所以本文的 $\mathbf W=\mathbf H\mathbf A^{-1}\mathbf H^T$ 等于 $h^2\mathbf H(\mathbf M+h\mathbf B+h^2\mathbf K)^{-1}\mathbf H^T$，**和 Martin 论文式 3.24 的 $\mathbf W=h^2\sum\mathbf H\mathbf A^{-1}\mathbf H^T$ 是同一个量**，只是两边的 A 定义不同（差一个 h² 因子）。
