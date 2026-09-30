# Vanneste, Martin, Goury, Courtecuisse, Pernod, Cotin, Duriez 2024 — Towards Realistic Needle Insertion Training Simulator Using Partitioned Model Order Reduction

- **出处**：MICCAI 2024（LNCS，doi 10.1007/978-3-031-72089-5_62）。HAL hal-04717755，CC BY 4.0。INRIA Lille（DEFROST）+ INRIA Strasbourg + **InfinyTech3D**（Erik Pernod，也就是 CollisionAlgorithm 插件的维护方）
- **本地文件**：`~/Downloads/extra_papers/A4_Vanneste2024_MICCAI_partitioned_MOR.pdf`
- **阅读状态**：全文精读；**2026-09-30 用原文文本核对了关键数值和结论**（见文末）
- **在本项目中的地位**：**在 SOFA + 开源 ModelOrderReduction 插件里实现的**"针道附近用全阶精细网格、其他地方降阶"。是本项目**提速的主要候选方案**之一，特别适合规划（候选轨迹已知）。

---

## 1. 问题
- 针插入需要细网格才能准确描述针道附近强烈的局部非线性形变，但细网格太慢。
- 标准 MOR 在针插入的区域不适用：**降阶模型描述不了局部的强非线性**。

## 2. 方法

### 2.1 模型（§2.1）
- 针：梁单元，每节点 6 个自由度；器官：线性四面体 + **共旋线弹性**。
- $\mathbf M\dot{\mathbf v}=\mathbf P-\mathbb F(\mathbf q,\mathbf v)+\mathbf H^T\lambda$；隐式积分，线性化，忽略粘性：$\mathbf A_n d\mathbf v=\mathbf b_n+h\mathbf H^T\lambda$（式 1–3）。

### 2.2 标准 MOR（§2.2）
- snapshot POD：$\mathbf q(t)\approx\mathbf q(0)+\Phi\alpha(t)$（式 4）；截断误差 $\nu^2=\sum_{p+1}\sigma_i^2/\sum\sigma_i^2$（式 5）。
- 降阶方程：$\Phi^T\mathbf A_n\Phi\,d\alpha=\Phi^T\mathbf b_n+h(\mathbf H\Phi)^T\lambda$（式 6）。
- **超降阶（hyperreduction）用 ECSW**（energy-conserving sampling and weighting）：只在降阶积分域（RID）上组装。

### 2.3 分区降阶（§2.3）⭐
- **离线**把肝脏的自由度分成两组：F（**针可能插入的区域，全阶精细网格，不降阶**）和 R（其余部分，降阶）。
- 分块方程（式 7）：$\begin{bmatrix}\mathbf A_{RR}&\mathbf A_{RF}\\\mathbf A_{FR}&\mathbf A_{FF}\end{bmatrix}\begin{bmatrix}d\mathbf v_R\\d\mathbf v_F\end{bmatrix}=\begin{bmatrix}\mathbf b_R\\\mathbf b_F\end{bmatrix}$。
- 只对 R 降阶（式 8）：$\begin{bmatrix}\phi^T\mathbf A_{RR}\phi&\phi^T\mathbf A_{RF}\\\mathbf A_{FR}\phi&\mathbf A_{FF}\end{bmatrix}\begin{bmatrix}d\alpha\\d\mathbf v_F\end{bmatrix}=\begin{bmatrix}\phi^T\mathbf b_R\\\mathbf b_F\end{bmatrix}$ → 混合了降阶坐标和全阶坐标。
- R 部分还可以做超降阶，但**RID 必须包含两个区域交界上的所有单元**，以保证耦合正确（式 9）。
- 思路来自计算力学里的损伤问题（Kerfriden 2012/2013，局部–整体降阶）。

## 3. 实验（§3）
- **实现在 SOFA 和开源的 ModelOrderReduction 插件里**（github.com/SofaDefrost/ModelOrderReduction），用 SOFA 的 CUDA 实现。
- 解剖结构：肝脏、胃、大肠（精细网格，用弹簧固定在肝脏上）、膈肌（模拟呼吸）、小肠（粗六面体网格）；数据来自 IRCAD 3D-IRCADb-01 公开数据集。
- 网格：粗网格约 1k / 4k / 5k 个四面体，细网格约 9k / 19k / 40k 个四面体（大肠 / 胃 / 肝脏）。
- **生成降阶模型**：离线做一系列插针（针对某个肿瘤选定进针区域、针的方向和深度），针完全插入后**做一个小圆周运动**，尽量激发附近单元的变形；拔针复位；慢速进行（240 步），同时膈肌模拟呼吸。
- **性能**（表 1，i7-7820HQ + Quadro M1200，旧硬件）：只用 CPU 0.0046 fps；只用 CUDA 3.7 fps（240 ms）；**CUDA 碰撞 + MOR 16.1 fps（51 ms），比全 CUDA 快 4.35 倍**。
- **精度**（图 2，以细网格全阶仿真为基准）：只有呼吸时，精细降阶模型比粗全阶模型误差低 50%（肝脏低 71%）；**呼吸 + 插针时，精细降阶模型的肝脏误差比粗全阶模型低 65%，平均误差约 2–3 mm**。

## 4. 局限
- 插入区域 F 要**离线指定**，还要针对它训练降阶模型 → 训练时插针位置不固定就不好用（Martin 论文因此没采用）；
- F 区域不降阶，F 越大越慢；
- 只有线弹性共旋；
- 测试硬件是旧的笔记本 GPU。

## 5. 和本项目的关系
1. ⭐ **对"术前规划"的场景非常合适**：规划时候选的进针区域和方向是已知的或有限的 → 可以离线对这个区域生成分区降阶模型 → 在线（CE 采样、MPC）时快得多。
2. **全部在开源的 SOFA 生态里**：ModelOrderReduction 插件**已经装在我们的 SOFA 里**（`plugins/ModelOrderReduction`）。**要确认已装的版本里有没有"分区"功能**（这篇论文的方法是不是已经合并到插件的主分支）。
3. 和 InfinyTech3D（CollisionAlgorithm 的维护方）有合作 → **CollisionAlgorithm 和 MOR 大概率可以一起用**。
4. 对我们的"立方体"场景，降阶的收益可能不大（网格本来就小）；**等做到"针对具体患者的解剖模型"阶段，它才是关键的提速手段**。所以它是**后期的优化选项**，不是起点。
5. 训练数据的生成方法值得借鉴：插针 + 在最深处做小圆周运动，用来激发附近区域的变形模态。

---
## 核对记录（2026-09-30）
- 用原文文本逐项核对了：ECSW、ModelOrderReduction 插件的 GitHub 地址、IRCAD 数据、240 步、i7-7820HQ + Quadro M1200、表 1（0.0046 / 3.7 / 16.1 fps；216975 / 240 / 51 ms）、快 4.35 倍、误差下降（只有呼吸时整体 50%、肝脏 71%；插针时肝脏 65%）：**和原文一致**。
- **补充**：原文还给了一个综合结论：**整体误差下降 68%，同时提速 4.35 倍**。
