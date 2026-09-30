# 文献笔记索引

PDF 原文位于 `~/sofa/papers/original/`，补充论文在 `~/sofa/papers/extra/`（都没有放进仓库）。每篇笔记都是逐页读原始 PDF 页面（含公式、图表）后写成的，并在 2026-09-29/30 对照原文逐篇核对过（每篇末尾有核对记录）。

| # | 笔记 | 论文 | 类别 |
|---|---|---|---|
| 00 | [00_synthesis.md](00_synthesis.md) | **综合分析、对比、待确认的问题**（§14 是核对后的更正） | — |
| H | [handbook.md](handbook.md) | **方法手册：统一符号，从头推导核心方法** | — |
| Q | [open_questions.md](open_questions.md) | 各篇笔记里待解决问题的处理情况 | — |
| T | [selftest.md](selftest.md) | 应用自测（18 题） | — |
| 01 | [01_Duriez2009_constraint_needle.md](01_Duriez2009_constraint_needle.md) | Duriez et al., MICCAI 2009 | 核心：约束法针插入 |
| 02 | [02_Duriez2006_haptic_contact_GS.md](02_Duriez2006_haptic_contact_GS.md) | Duriez et al., TVCG 2006 | 核心：Signorini + Coulomb + GS |
| 03 | [03_Martin2023_hybrid_solver_isodofs.md](03_Martin2023_hybrid_solver_isodofs.md) | Martin, Zeng, Courtecuisse, EG 2023 | 核心：网格相交 + 混合求解器 |
| 04 | [04_Martin2025_thesis.md](04_Martin2025_thesis.md) | Martin, 博士论文 2025 | 核心：完整体系 |
| 05 | [05_Okamura2004_force_modeling.md](05_Okamura2004_force_modeling.md) | Okamura et al., TBME 2004 | 参数：三部分力模型 |
| 06 | [06_Simone2002_insertion_forces.md](06_Simone2002_insertion_forces.md) | Simone & Okamura, ICRA 2002 | 参数（前身） |
| 07 | [07_vanGerwen2012_force_survey.md](07_vanGerwen2012_force_survey.md) | van Gerwen et al., MEP 2012 | 参数：综述 |
| 08 | [08_DiMaio2003_needle_modeling.md](08_DiMaio2003_needle_modeling.md) | DiMaio & Salcudean, TRA 2003 | 参数 + 早期仿真 |
| 09 | [09_Chentanez2009_remeshing_rod.md](09_Chentanez2009_remeshing_rod.md) | Chentanez et al., SIGGRAPH 2009 | 另一条路线：重划网格 |
| 10 | [10_Bui2018_error_control.md](10_Bui2018_error_control.md) | Bui et al., TBME 2018 | 另一条路线：自适应加密（SOFA） |
| 11 | [11_Bui2019_CutFEM.md](11_Bui2019_CutFEM.md) | Bui et al., CMAME 2019 | 另一条路线：CutFEM 多材料（SOFA） |
| 12 | [12_Wang2025_CE_MPC_flexible_needle.md](12_Wang2025_CE_MPC_flexible_needle.md) | Wang et al., RA-L 2025 | 应用：规划和控制（JHU） |
| 13 | [13_Ishida2024_situational_force_control.md](13_Ishida2024_situational_force_control.md) | Ishida, …, Sahu, IJCARS 2024 | 应用：数字孪生（Manish） |
| — | **补充阅读（2026-09-29）** | | |
| 14 | [14_Adagolodjo2019_inverseFE_robotic_needle.md](14_Adagolodjo2019_inverseFE_robotic_needle.md) | Adagolodjo et al., T-RO 2019 | ⭐ SOFA 逆向 FE + 机器人闭环（真实机器人） |
| 15 | [15_Perrusi2021_needle_laceration.md](15_Perrusi2021_needle_laceration.md) | Perrusi, Baksic, Courtecuisse, EG 2021 | 横向撕裂（rupture 的第三种含义） |
| 16 | [16_Wang2024_bevel_multilayer.md](16_Wang2024_bevel_multilayer.md) | Wang et al., ICRA 2024 | JHU 多层仿体 + 2D 弹簧模型，批评 SOFA 缺乏验证 |
| 17 | [17_Vanneste2024_partitioned_MOR.md](17_Vanneste2024_partitioned_MOR.md) | Vanneste et al., MICCAI 2024 | 分区降阶（SOFA MOR 插件） |
| 18 | [18_Baksic2020_moving_tissue_inverseFE.md](18_Baksic2020_moving_tissue_inverseFE.md) | Baksic et al., ICRA 2020 | 约束型目标函数，呼吸运动 |
| 19 | [19_Ha2024_2026_isoconstraint_learning.md](19_Ha2024_2026_isoconstraint_learning.md) | Ha, Bert, Courtecuisse, IROS 2024 + ICRA 2026 | 孤立目标约束 + 神经网络代替逆向仿真 |

| 20 | [20_Baksic2022_thesis.md](20_Baksic2022_thesis.md) | Baksic, 博士论文 2022 | SOFA 实现细节 + 针和仿体的参数标定方法 |

未精读：`B4_Adagolodjo2016_IROS_inverseFE.pdf`（T-RO 2019 的仿真版前身，内容被 T-RO 2019 覆盖）。
