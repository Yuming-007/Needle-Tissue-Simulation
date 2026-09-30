# 资源学习清单（达到"完全掌握、灵活运用"的标准）

状态：☐ 未开始 ◐ 进行中 ☑ 完成（完成后注明笔记位置）

## A. 针插入插件（InfinyTech3D CollisionAlgorithm + ConstraintGeometry）
- ☑ A1 InsertionAlgorithm、InsertionResolution、穿刺和摩擦的机制（`docs/code/01`，已运行验证）
- ☑ A2 约束组装：TBaseConstraint、InternalConstraint、H 的构建、violation 的符号约定、storeLambda
- ☑ A3 几何层：Geometry / Element / Proximity、FindClosestProximity、过滤器、AABB 宽相位、ContainsPoint、法向处理器（Phong、Gouraud）
- ☑ A4 CollisionLoop / CollisionPipeline 的调度时序（在一个时间步里什么时候执行）
- ☑ A5 其他官方场景：NeedleInsertionCycles（反复插拔）、NeedleInsertionHaptics；运行并分析
- ☑ A6 上游 master 和 ICube 原版的改动细节（new_design、Hermite 曲线针、GS 修复），和 v25.12 的兼容性

## B. SOFA 本体（v25.12）
- ☑ B1 FreeMotionAnimationLoop 的完整时间步（自由运动、碰撞、约束、修正的顺序和事件）
- ☑ B2 ODE 求解器：EulerImplicitSolver（积分因子、firstOrder、Rayleigh 阻尼）→ 发现一阶针 + 二阶组织导致作用力≠反作用力（`docs/code/01` 问题 2）；StaticSolver 已看（`docs/code/02`）
- ☑ B3 ConstraintCorrection：LinearSolver / Precomputed / Uncoupled / Generic 的实现与适用条件
- ☑ B4 其他约束求解器：ProjectedGS、UnbuiltGS、NNCG、ImprovedJacobi；SOFA 自带的 Lagrange 约束（Bilateral、Sliding、Unilateral）
- ☑ B5 组织 FEM：TetrahedronFEMForceField（能否逐单元设 E）、FastTetrahedralCorotational、超弹性 FEM、粘弹性
- ☑ B6 针：BeamFEMForceField（Timoshenko 实现细节）、BTDLinearSolver
- ☑ B7 SofaPython3：控制器、事件、数据读写、矩阵和能量访问、状态保存与恢复的可行方案、多进程
- ☑ B8 性能工具：AdvancedTimer、多线程组件

## C. 其他插件
- ☑ C1 Cosserat：理论、组件、针插入示例（尝试编译并运行）
- ☑ C2 SoftRobots.Inverse：QP 逆向问题的建模方式，能否用于针插入控制
- ☑ C3 ModelOrderReduction：工作流程、和约束的兼容性
- ☑ C4 BeamAdapter：和针相关的部分（滑动约束）
- ☑ C5 拓扑切割（SofaCarving、Tearing）和 InfinyToolkit（NeedleTracker）：和"组织破裂"的关系

## D. 文献和社区里的实现细节
- ☑ D1 Baksic 2022 博士论文（SOFA 里的针插入和逆向控制实现）
- ☑ D2 SOFA 论坛和 GitHub Discussions 里的针插入讨论
- ☑ D3 SOFA 官方文档：约束、碰撞、Lagrange 乘子相关页面

## E. 综合
- ☑ E1 用实验把代码笔记里的结论串起来（例如网格收敛、约束间距的影响、性能量级）
- ☑ E2 更新 handbook：SOFA 实现和理论的逐项对照
- ☑ E3 应用自测（代码层面），自评是否达到标准
