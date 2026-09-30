# Martin 2025 博士论文 — Modeling and Resolution of Needle-Tissue Interactions for a Fast and Stable Haptic Rendering: Application to Hepatic Percutaneous Procedures

- **出处**：Université de Strasbourg（ICube / INRIA MIMESIS），2025 年 12 月 4 日答辩。作者 Claire Martin，导师 Hadrien Courtecuisse。答辩委员：M. Marchal、F. Bello（报告人）、F. Zara、A. Lelevé、Y. Adagolodjo。HAL tel-05470170
- **本地文件**：`~/sofa/papers/original/Modeling and Resolution of Needle-Tissue Interactions for a Fast and Stable Haptic Rendering.pdf`（196 页；PDF 页码 = 印刷页码 + 17）
- **阅读状态**：第 1–8 章正文逐页精读（含全部公式、算法、图表）。第 9 章法语摘要和参考文献只浏览。**2026-09-29 对照原文核对了第 3–6 章；2026-09-30 补核第 2 章 §2.3.3–2.4.2（摩擦、穿刺、刚度力模型和求解器综述）**；第 7 章"组件私有"的结论是直接引用的原文，补充 4 处（见文末）
- **在本项目中的地位**：Duriez 2009 这条路线**最新、最完整**的发展（2025 年）。对"约束怎么放、怎么解、为什么不稳定"讲得最透彻。**但核心代码没有开源**（§7.3.3）。

---

## 0. 总览：四项贡献

| 章 | 贡献 | 场景 | 发表 |
|---|---|---|---|
| 4 | 异构模型组合 + 用网格相交定位约束点 + 混合求解器（直接解双边块、迭代解摩擦） | 肝脏 + 周围器官，低频仿真 | EG 2023，MICCAI 2024（Vanneste） |
| 5 | 异步 haptic 框架：针和局部组织在 1 kHz 下变形，约束方向在高频下更新 | 单个体，无摩擦 | IROS 2024 |
| 6 | 非线性粘滑摩擦的求解：先估计状态 → 线性化 → 截断 SVD 正则化 → 直接求解 | 单个体，有摩擦 | — |
| 7 | 模拟器集成：SOFA + Assist + ROS2 + Unity VR + 透视成像 | 肝脏 RFA | — |

## 1. 问题与挑战（第 1 章）
- 应用场景：肝癌射频消融（RFA）的穿刺训练模拟器，带 haptics。
- 三类挑战：
  - ① 计算效率：视觉 60 Hz，haptic 最好 1 kHz，最低 600 Hz；
  - ② 稳定性：约束问题病态；
  - ③ **针插入是非完整（non-holonomic）过程**：每一步的约束解都会写进针道，误差会沿轨迹不断累积。

## 2. 文献综述要点（第 2 章）

### 2.1 形变模型
- 质点弹簧：快，但不准。
- FEM：金标准，可以配合模型降阶（MOR，Goury & Duriez 2018）、区域分解等提速。
- PBD / XPBD（Macklin；Guo 2025 用于肝脏）：快且稳定，但物理参数难以对应。
- 可微分仿真。

### 2.2 交互建模
- penalty 法：效果依赖刚度参数，容易不稳定。
- 能量障碍法：IPC（Li 2020）；也有用于杆件的版本（Choi 2021）。
- **约束法 + Lagrange 乘子**（Signorini 接触、Coulomb 摩擦、等式约束）。
- Erleben 2020 提出**各向异性摩擦锥**，这很适合针：沿轴向和横向可以用不同的摩擦（这是我的联想，论文没用）。

### 2.3 针插入力模型（对参数化非常重要）
- **三类力**（Okamura 2004 分类）：摩擦、切割、刚度（穿刺前的组织形变）。
- **摩擦模型**：
  - Coulomb（Kikuuwe 2006）；
  - Karnopp（零速附近 ±Δv 的带内为静摩擦）；
  - **修正 Karnopp**（Okamura 2004：带粘性斜率 b_n、b_p）；
  - LuGre（有 Stribeck 效应，6 个参数，很难辨识；Khalaji 2013、Alamilla 2022）；
  - 傅里叶级数形式；
  - Kobayashi 2009：速度仿射项 + 低速时的对数项。
  - **Martin 的结论：针的摩擦强烈依赖速度，所以纯 Coulomb 不合适。**
- **穿刺**：
  - Okamura：针尖位置的二次多项式（穿刺前的刚度），切割力对某种组织是常数；
  - Misra 2008：针尖力对**断裂韧性**比对弹性更敏感；
  - Hing 2006：速度越低，穿刺力越大、切割力越小。
- **横向力**：
  - 弹簧模型：Roesthuis 2011，Liu 2020（非线性弹簧）；
  - 约束模型：Chentanez 2009、Duriez 2009、Wang 2012。
- **其他模型**：
  - Asadian 2010：用深度代替速度的 LuGre 式模型；
  - Pepley 2018：分段指数模型（尸体颈部数据）；
  - Dehghan 2007：三参数模型（摩擦力密度 + 切割力密度）。

### 2.4 求解器
- 直接法（LU、LDLᵀ）和迭代法（Jacobi、GS，约束的排序会影响收敛，Andrews 2017）。
- LCP 主元法（Murty、块主元法）、PGS。
- **直接和迭代混合**（Silcowitz 2010；Lacoursière & Linde 2011）→ 这是第 4 章混合求解器的来源。
- QP 求解器（OSQP、PIQP）、非光滑 Newton（Alart 1997；Bertails-Descoubes 2011；Macklin 2019）。
- **稳定化**：Tikhonov 正则化、截断 SVD（TSVD；截断阈值的选法见 Gavish & Donoho 2014）→ 这是第 6 章的来源。

### 2.5 Haptics 和已有模拟器
- 最低 600 Hz。Phantom Omni 最常用。
- 异步框架：Courtecuisse 2015（约束问题在低频下更新）、Peterlik 2011。
- 已有模拟器：Sutherland 2013（FEM 2000 节点，视觉 20 Hz / haptic 1 kHz）、Chan 2010、Alamilla 2022。
- 临床医生的建议：**针身必须沿针尖轨迹走**；要能在插入过程中改变轨迹；力要随器官刚度变化。

## 3. 数值基础（第 3 章）——整篇论文和 SOFA 的共同框架

### 3.1 模型
- **器官**：四面体，每节点 3 个自由度，**线弹性 + 共旋**（Felippa 2000）：$\tilde{\mathbf K}_e = \mathbf R_e\mathbf K_e\mathbf R_e^T$，每个单元一个 3×3 旋转。因为计算成本，**不用超弹性**。
- **针**：Timoshenko 梁（考虑剪切），两节点单元，每节点 6 个自由度；$\mathbf K_e = \Lambda\mathbf K_0\Lambda^T$。**针尖是对称的，没有斜面效应**（留作以后的工作）。

### 3.2 时间积分（式 3.4–3.10）
$$\underbrace{(\mathbf M + h\mathbf B + h^2\mathbf K)}_{\mathbf A}\underbrace{h\mathbf a^{t+h}}_{\mathbf x} \approx \underbrace{h(\mathbf g-\mathbf f_t) - h^2\mathbf K\mathbf v^t}_{\mathbf b} + h\mathbf c$$
隐式 Euler，只做一步 Newton 线性化，x 是速度增量。碰撞检测用邻近检测，不用连续碰撞检测。

### 3.3 约束系统（式 3.11–3.16、3.23–3.26）
- 约束点用重心坐标插值：$\mathbf q=\sum\phi_j\mathbf p_{k,j}$；$\mathbf c=\mathbf H^T\lambda$，$\delta=h\mathbf H\mathbf x$。
- **H 在一个时间步内保持不变**，在步首由碰撞检测给出。
- KKT 系统：$\mathbf A_i\mathbf x_i=\mathbf b_i+h\mathbf H_i^T\lambda$，$\delta^{t+h}=\delta^t+h\sum\mathbf H_i\mathbf x_i$。
- 每步 5 个步骤：
  1. 自由运动 $\mathbf x^{free}=\mathbf A^{-1}\mathbf b$；
  2. 组装约束系统：$\delta^{t+h} = \delta^{free} + \mathbf W\lambda$，其中 $\mathbf W = h^2\sum \mathbf H_i\mathbf A_i^{-1}\mathbf H_i^T$，$\delta^{free}=\delta^t+h\sum\mathbf H_i\mathbf x_i^{free}$；
  3. 分块 GS 求解（算法 1）；
  4. 运动修正 $\mathbf x^{corr}=h\mathbf A^{-1}\mathbf H^T\lambda$；
  5. 时间积分。
- **GS 的收敛判据用 $\|\lambda^{k+1}-\lambda^k\|$**，不用残差，因为摩擦滑动时 δ ≠ 0，残差没有意义。

### 3.4 四种约束律（§3.3.3）
1. **双边约束**：δ = 0。
2. **单边约束（Signorini）**：$0\le\delta_n\perp\lambda_n\ge0$。
3. **Coulomb 摩擦**（只用于表面接触）：粘着时 $\|\lambda_t\|<\mu\|\lambda_n\|$，滑动时 $\lambda_t=-\mu\|\lambda_n\|\delta_t/\|\delta_t\|$。静、动摩擦用同一个 μ，两个切向。
4. **针摩擦**（每个约束点一个，沿针轴方向）：
   $$\lambda_f = \eta\,\lambda_b,\quad \eta\in[0,1)\qquad(3.20)$$
   λ_b 是"假如沿针轴也是双边约束"时的力。η 是比例：η < 1 就允许滑动，而且 λ_b 和 δ^free（相对速度 × h）成正比 → 摩擦力**随速度增大**。这个版本**没有静摩擦**。
   - ⚠️ **我的分析**：标量情况下 $\lambda_b\approx-\delta_b^{free}/W_{bb}$，所以 $\lambda_f\approx\eta\,v_{rel}h/W_{bb}$。**等效的粘性系数取决于 h 和局部柔度 W**（也就是网格和刚度）→ η **不是物理参数**，换了时间步长或网格就要重新标定。这对"用实验数据定参数"影响很大。

### 3.5 三种交互模型（§3.3.4，图 3.9–3.10）
1. **器官之间**：LDI 体积接触（Allard 2010）+ Coulomb 棱锥摩擦。
2. **穿刺前**：一对约束点，一个是针尖节点，一个是它在组织表面的投影；方向取三角形法向。约束律是 Signorini + 2 个 Coulomb 摩擦。**穿刺阈值 λ_p**：0 < λ_n < λ_p 时表面变形；λ_n ≥ λ_p 时判定刺穿，切换到内部交互。
3. **穿刺后**：沿针尖轨迹，每隔固定距离 d 生成一个约束点：$\mathbf q_{k+1}=\mathbf q_k+d\,(\mathbf p_{tip}-\mathbf q_k)/\|\cdot\|$（式 3.22）。每个点：**2 个垂直于针的双边约束 + 1 个沿针轴的针摩擦约束**。
   - ⚠️ **穿刺后不再有切割力（f_c）**，这一点和 Duriez 2009 不同。第 6 章的力曲线里，针身摩擦占主导。

### 3.6 性能手段（§3.5）
- 瓶颈：自由运动和构造 W 都需要 A⁻¹。
- **针**：A_n 是块三对角（BTD）→ 用 Thomas 算法（SOFA 里有 `BTDLinearSolver`）。
- **MOR**（§3.5.2）：
  - snapshot POD：$\mathbf p\approx\mathbf p_0+\Phi\alpha$，截断误差 $\nu^2=\sum_{k+1}\sigma_i^2/\sum\sigma_i^2$；
  - 加上超降阶（hyperreduction，Ryckelynck 2005）：只在一部分单元（RID）上组装；
  - 然后 $\mathbf W_r=\mathbf H_r\tilde{\mathbf A}_r^{-1}\mathbf H_r^T$ 很便宜。
  - **只用于不被穿刺的器官**，因为穿刺路径没法事先知道。
  - SOFA 的 `ModelOrderReduction` 插件已经装好了。
- **异步预条件器**（Courtecuisse 2015）：A⁻¹ 在另一个线程里做 LDLᵀ 分解，要花几步；这期间仿真用上一次的分解。
- **IsoDOF**（Zeng 2022）：$\mathbf H=\hat{\mathbf H}\bar{\mathbf I}$，$\mathbf W_v=\hat{\mathbf H}(\bar{\mathbf I}\mathbf A^{-1}\bar{\mathbf I}^T)\hat{\mathbf H}^T$。中间那一项在相邻时间步之间变化很小，可以增量更新。

### 3.7 Haptics（§3.6）
- God-object 方法：设备和针根之间一根弹簧，$\mathbf g=-\mathbf K_{hg}(\mathbf p_g-\mathbf p_h)$。
- Courtecuisse 2015 的异步方案：W 在低频下更新，δ^free 在 haptic 频率下更新。

## 4. 第 4 章：仿真环境与交互处理

### 4.1 按用途选不同模型（表 4.3，图 4.2）

| 类别 | 用在哪里 | A⁻¹ 的处理方式 | 在 SOFA 里大致对应 |
|---|---|---|---|
| 针 n | 针 | BTD + Thomas | `BTDLinearSolver` |
| 精细体 v | **被穿刺的**：肝脏（1197 节点）、皮肤（142） | 异步预条件器 + IsoDOF | `LinearSolverConstraintCorrection`（没有异步） |
| 降阶体 r | 大肠、胃+食管 | MOR + 超降阶 | ModelOrderReduction 插件 |
| 粗糙体 c | 膈肌、小肠（各 60 节点） | 预先算好 A⁻¹ | `PrecomputedConstraintCorrection` |

- 以上都统一到 $\mathbf W=\sum\mathbf H_i\mathbf A_i^{-1}\mathbf H_i^T$ 里一起求解。
- **性能**（Ryzen 9 5900X，RTX 3070 Ti）：
  - 80–143 个约束；自由运动 7.5–9.6 ms，组装 W 2.5–6 ms，整步 19–28 ms；
  - **时间步 25 ms → 刚好实时**；
  - 注意**网格很粗**（肝脏才约 1200 个节点）。

### 4.2 用网格相交定位约束点（§4.3.1，图 4.3）
- 固定间距 d 的两难：d 大 → 组织和针之间出现不真实的相对滑动；d 小 → 过约束。
- 新做法：**针的边和组织三角面的每一个交点**都放一个约束点，不需要参数，密度自动匹配网格。
- 靠近顶点时约束点会挨得很近 → 依然需要一个稳健的求解器。

### 4.3 混合求解器（§4.3.2，算法 2）
- 约束排序：所有双边约束放在前面，摩擦约束放在后面 → $\mathbf W=\begin{pmatrix}\mathbf W_b&\mathbf W_c\\\mathbf W_c^T&\mathbf W_f\end{pmatrix}$（式 4.8）。
- 每次迭代：
  - 双边块直接求：$\lambda_b^{k+1}=\mathbf W_b^+(-\delta_b^{free}/h^2-\mathbf W_c\lambda_f^k)$（式 4.9），用**伪逆**；
  - 摩擦约束逐个 GS 更新，再乘以 η。
- $\mathbf W_b^+$ 每次求解只算一次。
- 可以扩展到表面接触：在后面追加（1 个法向 + 2 个切向）块，按 GS 顺序处理（图 4.5）。

### 4.4 结果（§4.4.2）
- 网格相交 + IsoDOF：1072 节点的组织，16 节点的针，109–171 个约束；每步新增的 IsoDOF 不到 1 个；组装 W_v 约 1.0–1.25 ms。
- 混合求解器 vs GS：数据和 Martin 2023 的表 2 相同。**约束间距一旦小于单元尺寸（约 1 cm），GS 就急剧恶化**；混合求解器始终保持在 130–140 次迭代左右。
- **Vanneste et al. 2024（MICCAI）**：**沿预先确定的针道**对肝脏做局部降阶——针道附近用细网格，其余部分降阶。Martin 没采用，因为训练时用户的轨迹不固定。
  - ⭐ **对本项目很有价值**：术前规划时，候选轨迹是已知或有限的，这个方法可以直接用来提速。要找来读（TODO）。

## 5. 第 5 章：高频下的形变和约束更新
- **准静态阻尼格式**（式 5.1–5.4）：令 a ≈ 0，得 $(\tfrac1h\mathbf D+\mathbf K)\Delta\mathbf p\approx\mathbf g-\mathbf f_t+\mathbf c$，$\mathbf p^{t+h}=\mathbf p^t+\Delta\mathbf p$。这样对 h 的依赖更弱，便于两个不同频率的回路共享数据。
  - ⭐ **对本项目的意义**：机器人慢速插针本来就是准静态过程，这种不带惯性的格式很合适。SOFA 里大致对应 `EulerImplicitSolver firstOrder=True`（CollisionAlgorithm 例子里的针正是这么设的）或者 `StaticSolver`。
- **两个回路**（图 5.2）：
  - 仿真回路（50–100 Hz）：组织的自由运动、碰撞检测、H、异步 A⁻¹、W̄_v、整个组织的运动修正（式 5.10）；
  - haptic 回路（约 1 kHz）：只处理针（BTD），组装 $\mathbf W^{HR}=\hat{\mathbf H}_v\bar{\mathbf W}_v\hat{\mathbf H}_v^T+\mathbf H_n\mathbf A_n^{-1}\mathbf H_n^T$（式 5.12），GS 求解，针的运动修正，**组织只在 IsoDOF 上做局部修正**（式 5.11）。这样既是局部形变，又考虑了整个组织的柔度。
- **双边约束的方向取为垂直于"针轨迹"**（约束点连成的折线），而不是垂直于针的边（图 5.3、5.4）。否则当穿刺点落在弯曲针的节点上时，约束方向会在相邻两条边的方向之间来回跳，力曲线出现振荡和尖峰。
- **结果**：
  - haptic 一步平均 0.39 ms（最大 0.97 ms），其中求解 0.22 ms（12–16 个约束）；
  - 在 haptic 频率下更新 H，可以**消除针节点经过穿刺点时的力尖峰**（图 5.7、5.8）；
  - 组织从 1000 到 5832 个节点，IsoDOF 只从 36 增到 63，haptic 回路的耗时基本不变。

## 6. 第 6 章：非线性摩擦的高频求解

### 6.1 粘滑摩擦模型（式 6.1）
$$\lambda_f=\begin{cases}\lambda_{f,b} & \|\lambda_{f,b}\|<\|\lambda_{max\,stick}\|\ \text{（粘着）}\\ \eta\,\lambda_{f,b} & \text{否则（滑动）}\end{cases}$$
- 作者自己说 **λ_maxstick 在约束空间里"没有明确的物理量纲"**，取值 150 或 250 是凭"力的感觉合理"定的。**参数没有经过物理辨识。**
- 滑动时的力是 η × 粘着等效力，**没有上界**（不像 Coulomb 那样有上限），仍然依赖 δ^free 和 W。

### 6.2 为什么 GS 不行
- GS 提前终止时，**最后处理的约束解得准，最前面的误差大** → 新约束点加进来后，前面约束点的解会前后不一致 → 针道被扰动、出现漂移（非完整过程会把误差累积下来）。
- 摩擦状态频繁切换也会拖慢收敛。

### 6.3 方法（图 6.1，算法 4、5）
1. **估计摩擦状态**：先做 10 次 GS 迭代（GS1），判断每个摩擦约束是粘着还是滑动，然后固定下来。
2. **线性化**：$\mathbf W_l=\mathbf W-\mathbf G\mathbf W\mathbf G^T+\tfrac1\eta\mathbf G\mathbf W\mathbf G^T$（式 6.2），G 是滑动约束的指示对角阵。也就是把滑动约束对应的行和列都乘以 1/η（增加柔度），然后令 $\delta_f^{t+h}=0$ → 得到线性方程 $\mathbf W_l\lambda=-\delta^{free}$。
   - **我验证过**：标量情况下 $\lambda=-\delta/(W/\eta)=\eta\lambda_b$，和式 6.1 一致。
3. **正则化**：用截断 SVD（Spectra 库），$\lambda\approx-\mathbf W^+_{l,\bar k}\delta^{free}$（式 6.7）。
   - $k_t$：在时间预算内最多能算多少个奇异值（按系统规模 3N+1 查表）；
   - $\sigma_s$：低于这个阈值的奇异值都去掉，用来滤掉高频抖动；
   - 用上一步的 $k_{prev}+m$ 做预估（CTSVD）。
4. **直接求解**。

### 6.4 参数与结果
- **σ_s = 1e-4**：1e-2 会剧烈抖动；1e-5 和 1e-6 在硬组织上偶尔出现很大的偏差；1e-1 很平滑，但轨迹不真实。奇异值的范围大约是 1e-7 到 3。
- 实验用的组织杨氏模量 **E_v = 5.5 kPa（接近肝脏）和 20 kPa**。
- **性能**（表 6.3）：34 到 61 个约束时，求解 0.15–0.41 ms（GS1 100 次迭代要 0.26–0.66 ms）；55 个约束以内能到 1 kHz，61 个约束时还在 600 Hz 以上。20 个约束点对应 16 cm 长的针。
- 病态测试（约束点偏离针轴 1、2、5 mm）：GS1 的误差不均匀，偏向前面的约束；CTSVD 的平均误差接近 0，STD 和 MAX 也更小。
- 插入测试：在 20 kPa 的硬组织上，GS 的平均双边约束误差达到 1.5 mm；CTSVD 接近 0。
- **力曲线**（图 6.10）：以 7 mm/s 直线插入再拔出，针根的力随深度大致线性增长（摩擦），呈锯齿状（粘滑）；λ_maxstick 取 strong 时最大约 4 N，取 low 时约 2.5 N；拔出时对称变为负值。**没有穿刺峰值或切割段的描述。**

## 7. 第 7 章：模拟器
- 机器 1（Ubuntu）：SOFA 仿真（30–50 Hz）+ haptic（约 1 kHz），用第 4 章的混合求解器；haptic 用的是 §3.6 的传统异步框架（**第 5、6 章的方法没集成进来**）。
- 机器 2（Windows）：Unity VR，由 InfinyTech3D 开发；透视图像用 Beer–Lambert 光线投射（Vidal 2009）计算。
- 两台机器之间用 ROS2 通信；整个流程由 **Assist**（assist.cnrs.fr）驱动。
- **"这篇论文里开发的主要组件目前是私有的；Assist 正在逐步开源。"**

## 8. 局限（作者自述 + 我的分析）
**作者自述**：
- 第 5、6 章的框架只在单个组织体、没有表面接触的场景下验证过；
- 准静态假设在有呼吸运动时不成立；
- 摩擦状态的估计还比较粗；
- λ_maxstick 是常数，以后应该改成随法向压力变化。

**我的分析**：
1. **参数没有物理意义**：η、λ_maxstick 都依赖 h 和 W，没法直接用 Okamura 这类实验数据来标定。**这是我们做"物理可信的规划仿真"时必须解决的问题。**
2. **只有线弹性共旋材料**，没有超弹性；Manish 要求考虑非线性弹性。
3. **没有斜面针尖**，而规划和转向都需要。
4. **没有切割力，也没有多层组织**：实验都只有一种组织（肝脏加皮肤是两个独立的物体，而不是一个体里分层）。
5. 性能数据都来自粗网格。

## 9. 和本项目的关系

| 项目需求 | 论文给出的内容 | 能不能直接用 |
|---|---|---|
| 刚性针 → 柔性针 | Timoshenko 梁 + BTD 求解 | SOFA 有 `BeamFEMForceField` + `BTDLinearSolver` ✅ |
| 组织破裂 | 穿刺阈值 λ_p，单边约束切换成内部约束 | 已经在 CollisionAlgorithm 里（待核实）✅ |
| 多层组织 | 皮肤和肝脏是两个独立物体，各有自己的 λ_p | 要么做成多个物体，要么扩展成一个体内分区 ⚠️ |
| 速度 | IsoDOF、异步 → 私有 ❌；MOR 和预计算 → SOFA 有 ✅；混合求解器 → 算法简单，可以自己实现 ✅ | |
| 稳定性 | CTSVD、约束方向沿轨迹 | 算法清楚，可以复现（需要 Eigen / Spectra）⚠️ |
| 规划 / 控制 | 准静态格式 ✅；沿预定路径局部降阶（Vanneste 2024）⭐ | |
| 参数标定 | η、λ_maxstick 是经验值 | 需要我们自己建立和物理量的对应关系 ❗ |

**选约束点间距和求解器时要记住的教训**：
- 约束点间距 ≥ 组织单元尺寸，或者用网格相交法 + 稳健的求解器；
- 双边约束的方向应该垂直于"约束点轨迹"，而不是垂直于针的边；
- GS 提前终止的误差会沿针道累积。

## 10. 待解决的问题
1. CollisionAlgorithm 用的是固定间距还是网格相交？双边约束方向是按针的边还是按轨迹定的？摩擦是哪种模型（η 型、粘滑型，还是 Coulomb）？
2. Vanneste et al. 2024 MICCAI 的具体方法，以及有没有开源代码。
3. 怎么把 η、λ_maxstick 换算成物理量（单位长度的摩擦力 N/m、和速度相关的系数），让它们能用 Okamura 这类数据来标定？一个可能的办法是在约束律里直接用物理的摩擦力密度乘以约束点负责的针长 l，就像 Duriez 2009 的 $\lambda_f=l\pi d\,r$ 那样。

---
## 核对记录（2026-09-29，第 3–6 章逐页对照）
- 核对了式 3.1–3.35、4.1–4.9、5.1–5.14、6.1–6.8，算法 1–5，表 4.1–4.5、5.1–5.2、6.1–6.4，图 5.7、6.10 的数值：**原笔记没有发现错误**。
- **补充 1（约定差异，很重要）**：第 3、4 章是**动力学格式**：$\delta^{t+h}=\delta^t+h\sum\mathbf H\mathbf x$，$\mathbf W=h^2\sum\mathbf H\mathbf A^{-1}\mathbf H^T$，其中 $\mathbf A=\mathbf M+h\mathbf B+h^2\mathbf K$，x 是速度增量。第 5、6 章换成**准静态格式**：x = Δp（位置增量），$\mathbf A=\mathbf D/h+\mathbf K$，$\delta^{t+h}=\delta^t+\sum\mathbf H\mathbf x$，**$\mathbf W=\sum\mathbf H\mathbf A^{-1}\mathbf H^T$，没有 h²**（式 5.5–5.9）。所以同一篇论文里 W 就有两种定义。
- **补充 2**：第 6 章系统规模是 **3N+1**：每个约束点 3 个约束（2 个双边约束 + 1 个摩擦约束），**+1 是针尖处的一个摩擦约束**（原文 p.110）。
- **补充 3**：图 5.7（无摩擦、只有双边约束，343 个节点的组织）里，针根 Y 方向的反馈力大约在 **0.28–0.53 N** 之间波动。图 6.10（有摩擦，插入行程约 140 mm，7 mm/s）：strong 时插入末端约 4 N，拔出最低约 −3 N；low 时约 2.5 N 和 −2 N。**拔出时的力并不完全对称**（比插入时的峰值小）。
- **补充 4**：第 6 章 σ_s 的研究（表 6.2），针有 9 个节点，N0 是针根、N4 是中点、N8 是针尖；指标是节点在"垂直于相邻两节点连线的平面"内的最大偏移（mm）。σ_s = 1e-4 时，最大值出现在 S3（20 kPa）的 N8，为 0.676 mm。
- **补核第 2 章（2026-09-30）**：§2.3.4 的摩擦模型（Coulomb / Karnopp / 修正 Karnopp / LuGre / Fourier / Kobayashi）、Martin 的结论"Coulomb 不合适"（依据是文献模型里都有速度依赖）、穿刺和切割（Okamura、Misra 2008、Hing 2006）、刚度力、其他模型（Asadian、Pepley、Dehghan）、求解器综述（Andrews 2017、Silcowitz 2010、Lacoursière & Linde 2011、Erleben 2020 的各向异性摩擦锥）：**和笔记一致**。
