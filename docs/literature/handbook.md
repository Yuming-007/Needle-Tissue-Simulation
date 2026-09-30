# 方法手册：约束法针–组织交互仿真（统一符号与推导）

> 这份手册用一套统一的符号，从头推导 20 篇论文里的核心方法，并标出各论文写法和这里的对应关系。
> 标记约定：**〔原文〕** 表示结论直接来自论文；**〔推导〕** 表示我自己推导的，已检查过量纲和极限情况；**〔推论 / 提议〕** 表示我的推断或设计建议，还没有被文献或实验验证。
> 笔记编号见 `README.md`。

---

## 0. 符号

| 符号 | 含义 | 单位 |
|---|---|---|
| $i\in\{n,t\}$ | 物体下标：n 为针，t 为组织（多个组织体时记 $t_1,t_2,\dots$） | |
| $\mathbf q_i,\ \mathbf v_i$ | 物体 i 的节点位置和速度 | m, m/s |
| $\mathbf M_i,\ \mathbf B_i,\ \mathbf K_i$ | 质量、阻尼、切线刚度 | kg, N·s/m, N/m |
| $\mathbf D_i$ | 准静态格式里的（对角）阻尼 | N·s/m |
| $h$ | 时间步长 | s |
| $m$ | 约束的个数 | |
| $\mathbf H_i\in\mathbb R^{m\times N_i}$ | 约束 Jacobian：$\Delta\delta=\sum_i\mathbf H_i\Delta\mathbf q_i$ | 无量纲 |
| $\delta\in\mathbb R^m$ | 约束违背量（间隙或相对位移） | m |
| $\lambda\in\mathbb R^m$ | **约束力**（约定是力，不是冲量） | N |
| $\mathbf W$ | 约束空间里的柔度（Delassus 算子） | m/N |
| $\mathbf K^{\rm eff}_i$ | 等效刚度（见 §1） | N/m |
| $l_k$ | 约束点 k 负责的针身长度 | m |
| $d$ | 针直径 | m |

---

## 1. 时间离散和约束空间：统一三种约定〔推导〕

### 1.1 动力学格式（隐式 Euler，一步线性化）
未知量是 $\mathbf v^+$：
$$\underbrace{(\mathbf M+h\mathbf B+h^2\mathbf K)}_{\mathbf A}\mathbf v^+=\mathbf M\mathbf v+h(\mathbf g-\mathbf f(\mathbf q))-h^2\mathbf K\mathbf v\ \ (\text{按线性化方式可有差别})\ +\ h\mathbf H^T\lambda$$
所以 $\mathbf v^+=\mathbf v^{free}+h\mathbf A^{-1}\mathbf H^T\lambda$，$\mathbf q^+=\mathbf q+h\mathbf v^+$，
$$\delta^+=\delta+\mathbf H(\mathbf q^+-\mathbf q)=\underbrace{\delta+h\mathbf H\mathbf v^{free}}_{\delta^{free}}+\underbrace{h^2\mathbf H\mathbf A^{-1}\mathbf H^T}_{\mathbf W}\lambda .$$
**关键观察**：$h^2\mathbf A^{-1}=(\mathbf M/h^2+\mathbf B/h+\mathbf K)^{-1}$，所以
$$\boxed{\mathbf W=\sum_i\mathbf H_i(\mathbf K_i^{\rm eff})^{-1}\mathbf H_i^T,\qquad \mathbf K^{\rm eff}=\tfrac{\mathbf M}{h^2}+\tfrac{\mathbf B}{h}+\mathbf K\quad(\text{动力学})}$$

### 1.2 准静态格式（Martin 第 5、6 章；Adagolodjo 2019）
令 $\mathbf a\approx0$，得 $(\mathbf D/h+\mathbf K)\Delta\mathbf q=\mathbf g-\mathbf f+\mathbf H^T\lambda$，
$$\mathbf W=\sum\mathbf H(\mathbf K^{\rm eff})^{-1}\mathbf H^T,\qquad\mathbf K^{\rm eff}=\tfrac{\mathbf D}{h}+\mathbf K\quad(\text{准静态})$$

### 1.3 结论
- **所有论文里的 W 都是"约束空间里的柔度"**（量纲 m/N），区别只在 $\mathbf K^{\rm eff}$ 里包含了哪些项。论文写法对照：

  | 论文 | 写法 | 和本手册的关系 |
  |---|---|---|
  | Duriez 2009 式 (5) | $\mathbf H(\mathbf M/h^2+\partial\mathbb F/(h\partial\mathbf v)+\partial\mathbb F/\partial\mathbf q)^{-1}\mathbf H^T$ | 相同（**直接就是 K^eff**） |
  | Duriez 2006 式 (27)–(30) | 未知量是位移 u，$(M/\Delta t^2+D/\Delta t+K)$ | 相同 |
  | Martin 第 3、4 章式 3.24 | $h^2\sum\mathbf H\mathbf A^{-1}\mathbf H^T$，其中 $\mathbf A=\mathbf M+h\mathbf B+h^2\mathbf K$ | 相同（h² 吸收进去） |
  | Martin 2023 式 (3) | $\sum\mathbf H\mathbf A^{-1}\mathbf H^T$，**A 把 h 吸收了** | 相同 |
  | Martin 第 5、6 章，Adagolodjo 2019 | 准静态 $\mathbf A=\mathbf D/h+\mathbf K$（Adagolodjo 用 $\mathbf M+\mathbf K$，M 是只在受约束自由度上的正则化项） | 相同 |
- **物理含义**：$W_{\alpha\alpha}$ 是"在约束 α 上施加 1 N 的力，这一步里约束方向上产生多少位移"。动力学格式里质量项让 $\mathbf K^{\rm eff}$ 变大（$M/h^2$），**时间步越小，系统越"硬"**。准静态格式里，D/h 起同样的正则化作用。
- **刚体极限**：$\mathbf K\to\infty$ 时，动力学格式里 $\mathbf K^{\rm eff}\to\infty$、$\mathbf W\to0$；但 Duriez 2006 的全局共旋把刚体运动单独拿出来，所以 W 退化为刚体的 $\mathbf J_c(\mathbf M/h^2)^{-1}\mathbf J_c^T$，不会趋于零。**这是"刚性针用很硬的梁、而不是纯刚体"时问题依然良态的原因**〔推论〕。

### 1.4 一步的完整流程（SOFA 的 `FreeMotionAnimationLoop` 做的就是这些）〔原文 Duriez 2006/2009，Martin §3.4〕
1. 自由运动：$\lambda=0$，每个物体各自求 $\mathbf v^{free}$（或 $\Delta\mathbf q^{free}$）→ 需要解 $\mathbf A_i$；
2. 碰撞检测：确定约束点、方向 → 得到 $\mathbf H_i$ 和 $\delta$；**H 在这一步里固定不变**；
3. 组装 $\mathbf W=\sum\mathbf H_i(\mathbf K^{\rm eff}_i)^{-1}\mathbf H_i^T$ → 要对 $\mathbf H_i^T$ 的每一列求解一次，**一共约 m 次回代**（主要开销）；
4. 在约束空间里解非线性互补问题，得到 λ；
5. 修正运动：$\mathbf q_i\mathrel{+}=(\mathbf K^{\rm eff}_i)^{-1}\mathbf H_i^T\lambda$（针对各自的约定乘相应的 h 因子）。

---

## 2. 约束律：用"特征图 + 局部求解"的统一写法〔推导，依据 Duriez 2006/2009〕

GS 的每次局部更新：其他约束的 λ 固定不动，
$$\delta_\alpha=W_{\alpha\alpha}\lambda_\alpha+\tilde\delta_\alpha,\qquad\tilde\delta_\alpha=\delta^{free}_\alpha+\textstyle\sum_{\beta\ne\alpha}W_{\alpha\beta}\lambda_\beta .$$
记 **"粘着试探力"** $\lambda^{\star}_\alpha=-\tilde\delta_\alpha/W_{\alpha\alpha}$（也就是让 $\delta_\alpha=0$ 所需的力）。各种约束律的局部解如下：

| 约束 | 特征图 | 局部解 | 来源 |
|---|---|---|---|
| 双边约束 | δ = 0 | $\lambda=\lambda^\star$ | 全部 |
| 单边约束（Signorini） | $0\le\delta\perp\lambda\ge0$ | $\lambda=\max(0,\lambda^\star)$ | Duriez 2006 |
| **穿刺（三状态）** | δ≥0, λ=0 / δ=0, 0≤λ≤f_p / δ≤0, λ=f_p | $\lambda=\mathrm{clip}(\lambda^\star,0,f_p)$；**λ 达到 f_p → 状态切换为"已刺穿"** | Duriez 2009 |
| 切割（针尖前进） | 类似穿刺，阈值 f_c，δ 相对于已切路径的末端 | $\lambda=\mathrm{clip}(\lambda^\star,0,f_c)$ | Duriez 2009（描述很简略） |
| Coulomb 块（1 个法向 + 2 个切向） | 锥 | 法向：$\max(0,\cdot)$；切向：先按粘着求，超出 $\mu\lambda_n$ 就投影回锥面 | Duriez 2006 算法 1 |

**针身摩擦的四种写法**（同一个约束行，δ 用这一步内的相对位移，也就是 $\delta^t=0$）：

| 写法 | 局部解 | 参数和量纲 | 来源 |
|---|---|---|---|
| (a) 固定阈值（粘滑） | $\lambda=\mathrm{clip}(\lambda^\star,-F_s,F_s)$；滑动时 $\lambda=\mp F_d$ | $F_s=D\,l_k$，$F_d=(C+b|v|)\,l_k$（N） | DiMaio 2003，Chentanez 2009，Duriez 2009（$l\pi d\cdot r$） |
| (b) Coulomb（法向力来自横向约束） | $|\lambda_f|\le\mu\|\lambda_\perp\|$ | μ，无量纲 | Bui 2018/2019 |
| (c) η 型 | $\lambda=\eta\lambda^\star$ | η ∈ [0, 1)，**无量纲，而且没有物理意义**（见 §3） | Martin 第 3、4 章；Adagolodjo 的 μ_n |
| (d) 粘滑 + η | $|\lambda^\star|<\lambda_{ms}$ 时 $\lambda=\lambda^\star$，否则 $\lambda=\eta\lambda^\star$ | η、λ_ms（作者自己说"没有明确的物理量纲"） | Martin 第 6 章 |

---

## 3. η 型摩擦的等效物理参数〔推导〕

对单个摩擦约束（忽略耦合）：$\lambda^\star=-\tilde\delta/W_{ff}$。δ 取这一步内的相对位移，$\tilde\delta\approx-h\,v_{rel}$（$v_{rel}$ 是自由运动下针相对组织的轴向速度）。于是
$$\lambda_f=\eta\,\frac{h}{W_{ff}}\,v_{rel}\ \Rightarrow\ c_{\rm eff}=\frac{\eta\,h}{W_{ff}}\ [\mathrm{N\cdot s/m}],\qquad c_{\rm eff}/l_k=\frac{\eta\,h}{W_{ff}\,l_k}\ [\mathrm{N\cdot s/m^2}].$$
- 量纲检查：$h$ [s] / $W$ [m/N] → N·s/m ✓。
- **这是纯粘性摩擦**：没有静摩擦、没有 Coulomb 项，力和速度成正比 → 和 Martin 的"摩擦强烈依赖速度"的说法**自洽**（这是模型本身的性质，不一定是物理结论）。
- **要让它和 Okamura 的粘性系数 b 一致**：$\eta=b\,l_k\,W_{ff}/h$。$W_{ff}$ 取决于网格、材料和 $\mathbf K^{\rm eff}$（也就是也取决于 h）→ **换了网格、时间步或者材料，η 都得重新标定**。粗略估计 $W_{ff}\sim1/(E\,a)$（a 是单元的特征尺寸）〔推论〕→ $\eta\sim b\,l_k/(h\,E\,a)$。
- 数值例子〔推导〕：b = 212 N·s/m²（Okamura），l_k = 5 mm，h = 0.02 s，E = 5.5 kPa，a = 5 mm → η ≈ 212 × 0.005 / (0.02 × 5500 × 0.005) ≈ 1.9 > 1 → **超出了 η 的允许范围**。说明用软组织加这个时间步时，η 型摩擦**达不到**实验里的粘性系数；反过来说，只有当 W 很大（组织很软或单元很大）时 η 才能小于 1。这是一个量级估计，要用仿真验证。

---

## 4. 1D 针模型的一个根本局限，以及统一的摩擦律〔推导 + 提议〕

**问题**：在约束法里，针是一条 1D 曲线，组织对它施加的横向力来自横向双边约束的 $\lambda_\perp$。**针笔直插入、组织没有被侧向推开时，$\lambda_\perp\approx0$**。所以 Bui 的 Coulomb 形式（b）在直线插入时**给出的摩擦接近于零**。
- 但真实的摩擦主要来自**组织包裹有直径的针所产生的径向预压力**（针把组织撑开了一个直径为 d 的孔），在 1D 模型里**完全没有**。
- **Duriez 2009 的"压力 p 由组织刚度估计、取常数"，正好就是这个径向预压力**：$\lambda_f\le\mu\,p\,\pi d\,l_k$。
- **Bui 2018 的强烈网格依赖**可以用这一点解释：Bui 的摩擦完全来自 $\lambda_\perp$，而粗网格下组织更"硬"、针道更容易被扭曲 → $\lambda_\perp$ 更大 → 摩擦更大。这和他们观察到的"**网格越细力越小**"，以及"**μ 越小，网格的影响越小**"（图 11）都一致〔推论，和原文结果吻合，但原文没有给出这个解释〕。

**提议的统一摩擦律**（把两种来源都包括进去）〔提议〕：
$$|\lambda_{f,k}|\le\mu\,\big(\underbrace{p_0\,\pi d\,l_k}_{\text{径向预压力，按层取}}+\|\lambda_{\perp,k}\|\big)\ \ (+\ \text{可选的粘性项}\ b\,l_k|v|)$$
- 直线插入时退化为 Duriez 2009 或 Okamura 的"单位长度摩擦"：$\mu p_0\pi d\leftrightarrow C$（N/m）→ **可以直接用 Okamura 的 C、D 来标定 $\mu p_0\pi d$**；
- 针弯曲或受横向推动时，横向力会自然增加摩擦（Bui 的机理）；
- 参数都有物理量纲，而且和网格基本无关（$\lambda_\perp$ 部分仍然依赖网格）。
- **要确认**：CollisionAlgorithm 的 `ConstraintInsertion` 现在用的是哪种形式，能不能改成这种形式。

---

## 5. 穿刺、切割和多层组织：怎样再现实验里的力曲线〔推导 + 推论〕

### 5.1 Okamura 力曲线的约束法解释
| Okamura 的阶段 | 约束法里的机制 |
|---|---|
| 刺穿前二次上升 $a_1z+a_2z^2$ | 针尖和表面之间的单边约束 + 组织 FEM（线性材料下曲线形状主要来自几何非线性：共旋大转动，以及接触区域变化）→ **线弹性共旋未必能给出正确的二次形状**，可能需要超弹材料或调参 |
| 峰值（刺穿） | 单边约束力 λ_n 达到穿刺阈值 $f_p$ |
| 骤降 | 刺穿之后，针尖的阻力换成更低的切割阈值 $f_c$，**同时**表面回弹，储存的弹性能释放出来 |
| 线性增长 | 针身摩擦 ∝ 插入长度（§4） |
| 拔出变负 | 摩擦反向，没有切割 |

**定量推导**〔推导，数据来自 Okamura 2004〕：刺穿时（约 16.7 mm 深）力 ≈ 2.30 N；骤降 0.66 N → 刺穿后的总力 ≈ 1.64 N。这时针身只进去了很短一段，摩擦 ≈ C·l ≈ 10.6 N/m × 几 mm ≈ 0.03–0.05 N；切割 ≈ 0.94 N。
→ 1.64 − 0.94 − 0.05 ≈ **0.65 N 是剩余的组织弹性力**（表面没有完全回弹，z₃ > z₁）。
→ **仿真里要再现这条曲线，除了 $f_p>f_c$ 之外，组织在刺穿后还要保持一部分"被推进去"的变形**（约束法里会自然出现：已刺穿的约束点把组织钉在针上）。
→ 所以**不能简单地设 $f_p$ = 2.3 N、$f_c$ = 0.94 N 就指望力曲线和实验一致**；力曲线是 f_p、f_c、组织刚度、边界条件共同决定的。

### 5.2 Duriez 2009 和 Bui 2019 的差别
- Bui 2019：表面穿刺（λ_p0）→ 针尖切割（**针尖约束**，$\lambda_n\ge\mu\lambda_t+\lambda_{c0}$ 才前进）→ 力呈锯齿状（每次前进之后松弛）。图 25 中有明显的"刺穿后骤降"。
- Duriez 2009：Fig. 6 里看不到骤降，也没给 f_p、f_c 的数值。
- **锯齿**〔推论〕：针尖是"粘着 → 力积累 → 超过阈值 → 前进一段 → 松弛"的离散过程。锯齿的周期取决于约束点的间距（每新增一个约束点，就有一次状态切换）。对规划来说它是**数值噪声**；锯齿幅度可以通过减小约束间距或时间步来降低，但要受过约束的限制（§6.3）。

### 5.3 多层组织的三种实现方式〔推论 / 提议〕
| 方式 | 做法 | 优点 | 缺点 | 文献依据 |
|---|---|---|---|---|
| A. 多个独立组织体 | 每层一个网格，各带一个穿刺算法；层与层之间用约束或共享节点连起来 | 每层的穿刺检测都是现成的"表面检测" | 层间界面需要额外的耦合；网格生成更麻烦 | Martin（皮肤 + 肝脏是两个物体）、Chentanez（贴合网格的多个区域） |
| B. 一个网格、按单元设材料参数 | 每个四面体按所在区域取 E，界面不必和网格对齐 | 网格简单，改层厚不用重新划网格 | **层间的"表面"不存在** → 穿刺第二层时需要另外判断（按针尖所在单元的区域标签切换 $f_p$ 或 $f_c$） | Bui 2019（CutFEM 的简化版）、Chentanez（f_cut 按区域查表，界面附近取大值） |
| C. B + 界面上的虚拟穿刺 | 针尖从区域 1 的单元进入区域 2 的单元时，**临时加一个穿刺约束**（阈值取区域 2 的 f_p），刺穿后恢复为区域 2 的 f_c | 可以再现 van Gerwen 说的"每穿过一层断裂韧性更低的组织，出现一次峰值和骤降" | 需要自己实现界面检测和状态管理 | Duriez 2009（同一个穿刺约束可以多次触发） |

**建议**〔提议〕：第一版用 **B**（最简单，看插件的支持情况），同时让摩擦参数（§4 的 $p_0$、μ）和切割阈值按区域取值；如果验证时需要再现层间的峰值，再升级到 C。

---

## 6. 求解器〔推导，依据 Duriez 2006、Martin〕

### 6.1 GS 的收敛
- GS 对对称正定的 W 收敛，速度取决于**对角占优**的程度。
- **过约束**：约束数超过相关的自由度数 → W 奇异（秩 ≤ 相关自由度数）→ λ 不唯一，GS 在零空间里漂移，收敛非常慢。这就是"约束间距小于单元尺寸就出问题"的原因（Martin 2023 表 2、Adagolodjo 2019）。
- **提前终止的偏差**（Martin 第 6 章）：GS 按顺序更新，最后更新的约束残差最小 → 误差集中在前面的约束上；而针道上的约束是按插入顺序排列的 → 误差集中在针道**起点**附近，并且每一步都会被写进针道（非完整过程）。
  - 缓解〔推论〕：打乱或交替扫描顺序（对称 GS），用上一步的 λ 热启动。

### 6.2 混合求解器（Martin 2023）= 两块的块 GS
把 λ 分成双边部分 b 和摩擦部分 f，
$$\begin{pmatrix}\mathbf W_b&\mathbf W_c\\\mathbf W_c^T&\mathbf W_f\end{pmatrix}\begin{pmatrix}\lambda_b\\\lambda_f\end{pmatrix}+\delta^{free}=\begin{pmatrix}0\\\delta_f\end{pmatrix}.$$
每次迭代：$\lambda_b\leftarrow-\mathbf W_b^+(\delta_b^{free}+\mathbf W_c\lambda_f)$（精确地解双边块），然后对摩擦逐个做 GS。
- 这就是**块 GS**，其中一块精确求解 → 双边块内部的耦合一步就解决了，所以迭代次数和约束间距基本无关（Martin 2023 表 2：约 130–142 次）。
- 用伪逆是为了处理双边块的**秩亏**（过约束）：伪逆给出最小范数解。
- 实现成本低：$\mathbf W_b^+$ 每一步只算一次（规模 = 双边约束数）。

### 6.3 CTSVD（Martin 第 6 章）
- 固定摩擦状态后，滑动约束的特征图变成直线 $\lambda=\eta\lambda^\star$ ⇔ $\delta=\tilde\delta+(W/\eta)\lambda$……〔推导〕：局部有 $\lambda=-\eta\tilde\delta/W$，即 $(W/\eta)\lambda+\tilde\delta=0$，也就是**把这一行和这一列的柔度除以 η**（对应式 6.2 的 $\mathbf W_l$）→ 整个问题变成线性的 $\mathbf W_l\lambda=-\delta^{free}$。
- 再用截断 SVD 求解：去掉小奇异值 ↔ 去掉病态方向（相邻约束几乎共线时产生的方向），相当于正则化。
- 局限：摩擦状态估计错了就得不到真正的解；只适用于 η 型摩擦（粘滑 + η）。

### 6.4 W 的计算方式（性能的关键）〔原文；SOFA 组件已核实，见 `docs/code/02` B3〕
| 方式 | 适用条件 | SOFA 组件（待核实） |
|---|---|---|
| 直接分解 A、对 Hᵀ 的每一列回代 | 通用；网格小 | `LinearSolverConstraintCorrection` + 直接求解器 |
| 预先算好 A⁻¹（compliance warping） | 线性材料、网格固定、小变形 | `PrecomputedConstraintCorrection` |
| 只用对角近似 | 粗糙，但很快 | `UncoupledConstraintCorrection` |
| 块三对角（针） | 1D 梁 | `BTDLinearSolver` |
| 异步预条件器 + IsoDOF | 大网格 | 私有 |
| 降阶模型 | 非穿刺区域；或已知针道 | ModelOrderReduction 插件 |

**开销量级**〔原文〕：Ha 2024：组装 W（MM）占整步的 74–82%；Adagolodjo 2019：自由运动 26%、W 23%、7 次约束求解 42%。

---

## 7. 在约束空间里求 Jacobian（逆向 FE 控制）〔推导，依据 Adagolodjo 2019、Baksic 2020、Ha 2024〕

设控制输入 X（针根位姿）只通过针根的双边约束进入问题：$\delta^{free}(X)=\delta^{free}_0+\mathbf G\,\Delta X$。在同一步里 W、H 固定：
1. 如果所有约束都是线性的（全是双边约束，没有状态切换）：$\lambda=-\mathbf W^{-1}\delta^{free}$ → $\Delta\mathbf q=(\mathbf K^{\rm eff})^{-1}\mathbf H^T\lambda$ **对 X 是线性的** → Jacobian 可以精确求出：
   $$\frac{\partial\mathbf q}{\partial X}=-(\mathbf K^{\rm eff})^{-1}\mathbf H^T\mathbf W^{-1}\mathbf G .$$
2. 有摩擦或单边约束时 λ(X) 是分段线性的 → 用**有限差分**：每个扰动只需要重新解一次约束问题（小系统），不必重新组装 W → Adagolodjo 的"快 4.2 倍"就来自这里。中心差分（Baksic）可以扩大分段线性区域的覆盖范围。
3. **目标函数作为"虚拟约束行"**（Baksic 2020）：若 $e=\mathbf H_o\mathbf q$（位置的线性函数），则
   $$e=\underbrace{\mathbf H_o\mathbf q^{free}}_{e^{free}}+\underbrace{\mathbf H_o(\mathbf K^{\rm eff})^{-1}\mathbf H^T}_{\mathbf W_{oc}}\lambda .$$
   → 只需要 $\mathbf W_{oc}$ 的**行**（目标对约束力的敏感度），**不需要列**（目标行本身不施加力）→ 这正是原文"这些行的列设为 0、行保留"的含义。
4. **孤立目标约束**（Ha 2024）：不再对 e 线性化，而是算出目标所涉及自由度的**真实位置** $\mathbf x_{iso}=\mathbf x^{free}_{iso}+\mathbf W_O\lambda$，再计算任意非线性的 e（比如角度）。$\mathbf W_O$ 就是 $\mathbf W_{oc}$ 推广到"所有目标自由度"的版本。
- **对规划和 MPC 的意义**〔推论〕：CE 采样需要**多步**前向仿真，而上述技巧只在**一步之内**有效（W 固定）。用于多步时要么每步重算 W，要么在短时域内冻结 W（近似）。

---

## 8. 针模型〔原文 + 推论〕
| 模型 | 自由度 | 特点 | 论文 | SOFA |
|---|---|---|---|---|
| Timoshenko 梁（共旋） | 每节点 6 | 考虑剪切；针很细长时和 Euler-Bernoulli 几乎一样 | Martin、Adagolodjo、Baksic、Ha | `BeamFEMForceField` |
| Euler-Bernoulli 梁 | 每节点 6（3D）/ 3（2D） | 忽略剪切 | Bui 2019、Wang 2024/2025 | — |
| 离散弹性杆 + CORDE | 位置 + 扭转 | 不可伸长 | Chentanez 2009 | — |
| Cosserat 杆 | 应变参数化 | 能描述弯、扭、拉、剪 | （Manish 提到过） | Cosserat 插件（需编译） |
| Kirchhoff 杆 | | 导丝、导管 | | BeamAdapter（已装） |

- **刚性针**：用 E 很大的梁（CollisionAlgorithm 例子用 E = 1e12；Perrusi 的刚性针 E = 100–200 GPa，半径 8 mm），比纯刚体更良态（§1.3）。
- **离散化**：15 cm 针用约 50 个单元（Duriez 2009、Wang 2025）；12 cm 针用 28 个单元时收敛（Adagolodjo 2019）；13–16 个单元也常用（Baksic、Ha、Perrusi）。
- 钢针 E = 200 GPa（大多数论文）；Wang 2024 调参用 80 GPa；Baksic 2022 用挂重法标定活检针得到 77.94 GPa。〔推导〕这是因为**空心针按实心截面建模**：E_eff = E(1 − (rᵢ/rₒ)⁴)，反推 rᵢ/rₒ ≈ 0.88。→ SOFA 的 `BeamFEMForceField` 有 `radiusInner`，可以直接按真实截面建模；否则用约 80 GPa 的等效模量。
- **SOFA v25.12 的 `BeamFEMForceField` 实际上是 Euler-Bernoulli 梁**（有效剪切面积设为 0，见 `docs/code/02` B6）。

## 9. 斜面针尖：三种写法其实是同一件事〔推导〕
| 论文 | 写法 |
|---|---|
| Chentanez 2009 | 新的针尖节点在物质坐标系里横向偏移 $h=d\tan\psi$（d 是这一步的前进距离） |
| Wang 2024/2025 | 新的约束点带一个偏移 b（"预压缩弹簧"） |
| Duriez 2009 | 针尖路径约束的方向相对针轴倾斜 |
- 共同点：**每前进单位长度，就在针尖处给组织一个横向的"预应变"**，偏移率约为 tan ψ（或 b/d）。之后组织的弹性力把针推弯，**曲率自然取决于组织和针的相对刚度**。
- 在约束法里的实现〔提议〕：新建内部约束点时，把组织一侧的点沿斜面方向偏移 $d_k\tan\psi$（或让针尖路径约束的方向倾斜角度 ψ）→ 和 CollisionAlgorithm 的约束点生成逻辑结合（要看代码）。

## 10. 一个网格里的多种材料〔依据 Bui 2019，简化〕
- 被切单元的刚度：$\mathbf K_e=\sum_{r}\mathbf B^T\mathbf E_r\mathbf B\,V_{e,r}$（线性四面体的 B 是常数）→ **按体积分数加权** $E_e=\sum_r E_r V_{e,r}/V_e$ 在线性单元里是**精确等价**的〔推导：B 为常数时，$\int\mathbf B^T\mathbf E\mathbf B=\mathbf B^T(\sum\mathbf E_rV_r)\mathbf B$，只要 ν 相同，就等价于用加权后的 E〕。
  → 对"两种组织的立方体"，**每个单元按体积分数设 E 就等于 CutFEM 的精度**（前提是 ν 相同、线性四面体）。
- 更简单的近似：按质心所在区域取 E（阶梯近似），误差集中在界面附近的一层单元里。

## 11. 规划和控制接口需要什么（从 Wang 2025、Adagolodjo、Ha 总结）〔推论〕
- 状态 = 节点的位置和速度 + **约束点集合**（包括每个点的重心坐标和状态），两者都要能保存和恢复（Baksic 2020 的教训：约束位置不一致会让控制器发散）。
- 输入：针根位姿（插入、横移、旋转）；输出：针尖位姿、针的形状、针根的力、各层应变能、到关键结构的 SDF 距离。
- 敏感度：用 §7 的方法（单步），或者用有限差分（多步）。
- 速度：单步耗时 × 每次决策的前向仿真次数（CE：N 个样本 × 时域步数）必须满足控制频率。

---

## 12. 理论和 SOFA v25.12 实现的逐项对照（2026-09-30，依据 `docs/code/01–03`）

| handbook 里的概念 | SOFA / 插件里的实际实现 | 注意事项 |
|---|---|---|
| §1 W = Σ H (K^eff)⁻¹ Hᵀ（柔度，m/N），λ 是力 | W = J A⁻¹ Jᵀ × dt（二阶 EulerImplicit，A = M + h(h+rs)K + …）→ **λ 是冲量（N·s）**，力 = λ/dt〔代码 + 运行 + Baksic 2022 p.39〕 | 只是约定不同（λ_SOFA = h·λ_handbook，W_SOFA = W_handbook / h）；**一阶积分器会把同一个 λ 当作力** → 针和组织的积分器必须一致 |
| §1.4 一步五个环节 | `FreeMotionAnimationLoop::step`：`AnimateBeginEvent`（CollisionLoop 检测）→ 自由运动 → `processGeometricalData`（重建约束）→ 组装 W、GS → 修正 | 检测用步首位置，穿刺判定读上一步的 λ |
| §2 GS 局部求解 | `BlockGaussSeidelConstraintSolver` 按块调用 `ConstraintResolution::resolution()`；**每步从 λ = 0 开始**；停止判据 Σ‖W_block Δλ‖ < tol × 约束数（默认） | 默认容差很松，会掩盖过约束和数值摩擦 |
| §2 穿刺三状态律 | **事后检测**：上一步针尖约束冲量各分量范数之和 / dt > 阈值 → 切换状态；阈值没有进入约束律 | 用的是合力（含切向），不是法向力；滞后一步 |
| §2 切割约束 f_c | **没有实现** | 要自己加 |
| §2 针身摩擦 (a)/(b)/(c) | `InsertionResolution`：frictionCoeff 是**轴向欠松弛因子** → 收敛时完全粘住，实际是 η_eff ≈ k·μ 的 η 型 | 要换成 (a) 或 §4 的统一摩擦律 |
| §3 η 型摩擦依赖 h 和 W | 问题 1 的实验：同样的 frictionCoeff，阻力随 tol、maxIt 相差约 40 倍 | 比 §3 预测的还不稳定（η 还依赖迭代次数） |
| §4 1D 针没有径向预压力 | 插件的横向约束只把组织点拉回针轴，没有预压力项 | 和 §4 的分析一致 |
| §5.3 多层方案 B | SOFA 线弹性四面体力场的 `youngModulus` 可以逐单元设置；插件只检测一个外表面 | 界面穿刺（方案 C）要自己加 |
| §6.1 约束间距 ≥ 单元尺寸 | 实验：50 mm 单元、3 mm 间距时 GS 每步都达不到收敛（5001 次迭代）；25 mm 时约 15 次 | 验证了 |
| §6.2 混合求解器 | v25.12 的 5 个求解器里都没有 | 要在 `ConstraintResolution` 层面或新求解器里实现 |
| §6.4 W 的计算方式 | `LinearSolver` / `Precomputed`（稠密 N×N，EulerImplicit，dt 固定）/ `Uncoupled` / `Generic` | 官方场景里组装 W 占每步 58% |
| §7 约束空间 Jacobian | 需要保存和恢复状态 → **进程分叉快照可行**（插件内部状态一起复制） | 单步内还可以用 `solver.W()` 读 W |
| §8 Timoshenko 梁 | `BeamFEMForceField` 实际上是 **Euler-Bernoulli**（有效剪切面积 = 0）；支持空心截面 | 细长针可以忽略差别 |
| §9 斜面针尖 | 插件没有 | 在新建约束点时加横向偏移（要改 `InsertionAlgorithm`） |
| §10 逐单元材料 | 原生支持（线弹性）；超弹只能全局一组参数；Ogden 的 μ₁(SOFA) = 2μ(Wang) | |
| §11 应变能输出 | 共旋四面体力场的 `getPotentialEnergy()` **返回 0**（只有 `small` 方法实现了） | 要自己算 |
