# 第 1.5 步：非均匀工作网格的复用调研（2026-10-02）

- 用户要求：在正式写网格方案、安装 gmsh 之前，做一次范围受控的复用调研。原则是复用优先、逐步增加复杂度、用证据验证。
- 不展开的内容：自适应重划分、切割 / 重划分、病人个体化网格、高阶单元、无网格 / XFEM、网格优化文献、柔性针沿路径变化的网格。
- 本文**只做调研和方案判断**：没有安装 gmsh，没有写网格生成代码。所有参数（a_eff、h_local、加密半径、过渡宽度、节点预算）都**不在这里冻结**。

## 1. 查了哪些资源
| 资源 | 内容 | 依据 |
|---|---|---|
| gmsh 官方手册（4.15，Mesh size fields 一节） | 各尺寸场的定义和参数；哪些 3D 算法支持尺寸场；质量指标；输出格式和物理组 | 〔原文〕https://gmsh.info/doc/texinfo/gmsh.html |
| gmsh 官方教程 t10（Python） | Distance + Threshold + MathEval + Box + Min 组合成背景尺寸场的完整示例 | 〔原文〕github.com/live-clones/gmsh `tutorials/python/t10.py` |
| SOFA v25.12 `MeshGmshLoader` 源码 | 支持的格式、节点编号、单元类型、物理组 | 〔代码〕github.com/sofa-framework/sofa（v25.12 标签）`Sofa/Component/IO/Mesh/.../MeshGmshLoader.cpp` |
| SOFA 官方示例 `share/sofa/examples/Component/IO/Mesh/MeshGmshLoader.scn` 和 `share/sofa/mesh/msh4_cube.msh` | 读取 msh 2 / msh 4 四面体网格接 `TetrahedronSetTopologyContainer` 和 FEM 的写法 | 〔代码〕本地二进制包 |
| SOFA `TetrahedronSetGeometryAlgorithms.h` | 现成的单元检查函数 | 〔代码〕本地头文件 |
| Sandia CUBIT 手册（四面体质量指标表） | 常用指标和"可接受范围" | 〔原文〕sandia.gov CUBIT 15.7 tetrahedral metrics |
| 项目知识库 `docs/resources/inventory.md` §3，文献笔记 03 / 04（Martin）、10（Bui）、14（Adagolodjo），probe–tissue 教程 | 已有针 / 探头仿真的网格做法 | 〔原文〕〔代码〕已有笔记 |
| 检索 Chentanez 2009（Berkeley） | 沿曲线针道局部重划分 | 〔摘要〕people.eecs.berkeley.edu/~jrs/papers/needlesim.pdf |

## 2. gmsh 的局部加密机制
- **Distance**：到给定点（`PointsList`）或曲线（`CurvesList`，按 `Sampling` 采样）的距离。〔原文〕
- **Threshold**：输入一个场（通常是 Distance）。距离 ≤ DistMin 时尺寸为 SizeMin，≥ DistMax 时为 SizeMax，中间**线性插值**；`Sigmoid` 改成平滑过渡，`StopAtDistMax` 在 DistMax 以外不起作用。〔原文〕
- **Box / Cylinder / Ball**：区域内 VIn、区域外 VOut。Box 和 Ball 有 `Thickness`（过渡层厚度）；Cylinder 是突变。〔原文〕
- **MathEval**：任意表达式，可以直接写距离公式，不依赖几何实体。〔原文〕
- **Min**：取多个场的最小值，用来组合。〔原文〕
- **限制**：一般尺寸场只被 3D 的 **Delaunay** 和 **HXT** 算法支持，Frontal（Netgen）不支持。〔原文〕
  - t10 还关掉了 `MeshSizeExtendFromBoundary`、`MeshSizeFromPoints`、`MeshSizeFromCurvature`，否则边界上的尺寸会干扰尺寸场。〔原文〕
- **判断**〔推论〕：
  - **最适合的组合是 Distance（到一段线段）+ Threshold**，必要时用 Min 再叠一个 Box（例如表面附近的加密层）。这正是 t10 的写法，只需要改参数。
  - 加密区 = 到"针道线段"距离 ≤ DistMin 的区域。**只在针尖附近加密**就是"线段很短"的特例，**沿针道的走廊**就是"线段更长"，用的是同一套代码。
  - 过渡宽度 = DistMax − DistMin；相邻单元尺寸的增长率约为 (SizeMax − SizeMin) / (DistMax − DistMin)。
  - 线段不一定要作为几何实体嵌入体内：可以用 MathEval 直接写"到线段的距离"，避免嵌入实体带来的网格约束。这一点在实际生成时验证。
  - 组织顶面是 2D 网格，同样服从尺寸场，所以**接触表面会一起被加密**（这是我们需要的）。

## 3. SOFA 对 gmsh 四面体网格的支持
- `MeshGmshLoader`〔代码〕：
  - 支持 msh 1 / 2 / 4，**只支持 ASCII**，不支持二进制；
  - 节点编号转成从 0 开始。msh 2 用映射表处理编号不连续（"in case of hole or switch"）；**msh 4 默认假定每个块内的编号是连续的**；
  - 读取点、线、三角形、四边形、四面体、六面体（含二阶），物理标签变成 `trianglesGroups` / `tetrahedraGroups` 等分组。
- 官方示例：`MeshGmshLoader → MechanicalObject + TetrahedronSetTopologyContainer（src=@loader）→ FEM`。和我们现在"拓扑容器 + FEM + `Tetra2TriangleTopologicalMapping` 提取表面 + 插件几何"的结构兼容，**只换掉网格来源**。
- 可以直接复用的现有部分：
  - `FastTetrahedralCorotationalForceFieldFixed`（我们的修复针对刚度矩阵组装，和网格是否规则无关）；
  - `Tetra2TriangleTopologicalMapping` + IdentityMapping 提取表面；插件的 `TriangleGeometry`、`PhongTriangleNormalHandler`、`ConstraintUnilateral`、重心接触；
  - 针、驱动、读力、多点针尖、R1 / R3 求解路线，都和网格来源无关。
- **已知的坑 / 必须检查的项目**：
  1. **输出格式**：用 **msh 2.2 ASCII**（`Mesh.MshFileVersion = 2.2`），避开 msh 4 的"编号连续"假设和二进制；或者直接用 gmsh 的 Python API 取节点和单元数组传给 SOFA（我们的场景本来就是 Python 写的）。〔代码 + 推论〕
  2. **只导出四面体**：定义物理体后，gmsh 默认只导出属于物理组的单元〔原文〕。只设物理体、不导出表面三角形，表面统一由 `Tetra2TriangleTopologicalMapping` 提取，避免和读取器读入的三角形混在一起。〔推论〕
  3. **四面体的方向**：表面法向朝内还是朝外取决于四面体节点的顺序。第 1 步在规则网格上确定 `flipNormals=False` 时法向朝内；换成 gmsh 网格后要**重新检查**：有符号体积是否全为正，表面法向是否朝内。SOFA 有 `checkNodeSequence`。〔代码 + 推论〕
  4. **对称性没有了**：非结构网格不对称。V3-0 的"四分之一模型 = 全模型"不再成立；局部加载下的**人为侧向漂移**（第 1 步发现规则网格默认切分有约 30%）和**落点依赖**（V2）都必须在新网格上重测。〔推论〕
  5. probe 教程用的是单独的低分辨率碰撞表面 + BarycentricMapping，我们是直接从体网格提取表面，两者不同。对针尖接触来说，碰撞表面要和加密的体网格一致，所以**沿用我们的做法**。〔推论〕
- **没有找到**：SOFA 官方或开源的、在非均匀四面体网格上做探头 / 针接触的场景。probe 教程的乳房网格是均匀的（边长中位数 8.6 mm，见 inventory）。

## 4. 已有针 / 探头仿真的网格做法（只看网格）
| 工作 | 网格做法 | 规模 | 和我们的关系 |
|---|---|---|---|
| Chentanez 2009 | 针经过时**局部重划分**，让节点落在曲线针道上 | 13,375 四面体 / 2,763 节点，25 Hz | 动态改网格，和我们"不重划分、重心约束"的路线不同，**不照搬** |
| Goksel 2005 / 2006 | 粗网格上节点重定位 / 增加节点 | — | 同上，不照搬 |
| Bui 2018 | 误差驱动的**自适应**六面体加密（模板细分 + 悬挂节点） | — | 自适应，超出范围；可借鉴的是"在误差大的地方加密"的思路 |
| Martin 2023 / 2025 | 静态网格，沿针方向较密（15 cm 50 个节点），横向较粗；肝脏约 1200 节点 | 约 1000 节点量级 | 实时约束法的网格规模参考；**约束间距不要小于单元尺寸**（否则 GS 不收敛） |
| Adagolodjo 2019 | 由 CT 分割的均匀四面体网格，2592 个四面体 | — | 没有局部加密 |
| 离线有限元（Abaqus 等） | 针尖附近细、远处粗，做收敛测试 | 远超实时预算 | 只借鉴方法，不借鉴数值 |
- **证据边界**（ChatGPT 复核 2026-10-03，Claude 同意）：我们的测试台只能证明"在刺穿前的点接触中，单点刚度强烈依赖 h，且大压深下的稳健性随加细变差"。**不能**进一步推出"别人因为用单点所以不加密"：Bui 2018 同样属于约束法（刺穿前沿用 Duriez 2009 的单点），却做了自适应加密。很多实时工作不做局部加密，可能和实时成本、研究目标、交互架构都有关，目前没有证据把主要原因归到单点针尖。（修正记录：Claude 曾在讨论中口头这样说过，没有写进文档。）
- **结论**〔推论〕：实时的约束法插针工作基本用**均匀或接近均匀的粗网格**（约 1000 节点量级），或者动态改网格。**"静态的、针尖 / 针道附近细、远处粗"在实时插针里没有找到可以直接借用的公开实现。** 可复用的是工具（gmsh 尺寸场）和方法（距离驱动的尺寸场），具体参数要自己定。

## 5. 网格质量指标
- **现成工具**：
  - gmsh：`Mesh.QualityType`（gamma / SICN / SIGE），`getElementQualities`、`Plugin(AnalyseMeshQuality)`，`Mesh.Optimize` / `Mesh.OptimizeNetgen` 优化。〔原文〕
  - SOFA：`TetrahedronSetGeometryAlgorithms::computeBadTetrahedron`（二面角 20–160°、最长边 / 最短边 ≤ 10、节点顺序），`checkNodeSequence`，`TetrahedronSetTopologyContainer::checkTopology`。〔代码〕这些是 C++ 函数，不一定有 Python 绑定，我们可以用 numpy 自己算（V4 已经在算有符号体积）。
- **对低阶四面体有限元、大变形和接触最关键的指标**〔原文 + 推论〕：
  1. **有符号体积 / Jacobian**：必须全为正（硬条件）。V4 的翻转判据就是它。
  2. **二面角**：小二面角使刚度矩阵病态，大二面角损害精度〔原文：Shewchuk 一系的结论〕。最小二面角是初筛的首选指标。
  3. **形状质量**（gamma、SICN、scaled Jacobian、radius ratio，本质相近）：衡量单元离正四面体有多远。
  4. **相邻单元的尺寸比 / 尺寸梯度**：过渡太陡会造成局部各向异性和人为刚度跳变，对接触区的力读数影响大。
- **经验阈值**（**只能作初筛，不是硬标准**）：
  - CUBIT：长宽比 1–3，scaled Jacobian 0.2–1，形状质量 0.2–1 为"可接受"〔原文〕；
  - SOFA 自带检查：二面角 20–160°、边长比 ≤ 10〔代码〕；
  - 二面角能做到约 12–160° 已属较好的网格〔原文：检索结果中的质量改进文献〕。
  - 判断〔推论〕：这些阈值来自一般有限元经验，**真正的标准是我们自己的验证**：在新网格上重跑小压深接触对照（V1 类）、落点依赖（V2 类）、侧向漂移，以及大压深翻转压深（V4 类）。

## 6. 初步倾向（待用户确认）
1. **gmsh 尺寸场**：Distance（到一段竖直线段，或者用 MathEval 写距离公式）+ Threshold（线性或 Sigmoid）+ 必要时 Min / Box；3D 算法用 Delaunay 或 HXT；生成后用 `Mesh.Optimize` 优化；输出 msh 2.2 ASCII，只导出四面体。
2. **先做"以针尖为中心的加密体"，代码按"线段"写**：第 1.5 / 2 步只需要覆盖 tenting 和刺穿附近；线段加长就是第 3 步的针道走廊，参数化写法不变。现在就做完整走廊会浪费节点。
3. **SOFA 侧新增代码很少**：一个"从 gmsh 网格建组织"的函数（读取器或数组 → 拓扑容器），替换 `common.add_tissue` 里的规则网格部分，加上方向和法向检查；表面提取、接触、FEM、求解器、针、读力全部复用。估计几十行。
4. **验证复用**：V1（LCP 参照）、V2（落点）、V4（翻转）的脚本可以直接换网格重跑；四分之一模型的捷径不能用了，LCP 工具要用全模型（节点数少，应当可以接受）。

## 7. 必须靠我们自己实验决定的问题
- **h_local / a_eff**：按 V3-3 曲线取舍，但曲线要在新网格上重新验证。
- **加密半径（DistMin）和深度**：至少覆盖针尖支撑和 tenting 变形最大的区域（V4 显示圆盘可用压深约 1.1 a_eff）。具体多大由"加大半径时力是否还在变化"决定。
- **过渡宽度（DistMax − DistMin）和远场尺寸**：由尺寸梯度、质量指标和对力的影响决定。
- **节点数和实时成本**：见下面的风险。

## 8. 可能改变路线的风险
- **实时节点预算和加密需求之间的矛盾**（最重要）〔推论，粗估〕：
  - 0d 的节点预算约 870（规则网格、最重负载下）。非结构四面体网格在相同单元尺寸下，节点数大约是规则网格的 1.5 倍。
  - 粗估：a_eff = 10 mm，加密区取半径约 2a、深约 2a 的圆柱（约 2.5×10⁴ mm³）。h_local = 5 mm 时加密区约 300 个节点；h_local = 2.5 mm 时约 2400 个，**远超预算**。远场 80 mm 立方体用 10–15 mm 的单元约 250–800 个节点。
  - 所以在现有预算下，h_local 很可能只能取到约 5 mm（h/a ≈ 0.5，V3-3 中偏硬约 35–40%），除非：缩小加密区、缩小组织块、增大 a_eff，或者在新网格上**重新测实时预算**（870 是规则网格、R3 的结果，非结构网格的稀疏结构不同，需要重测）。
  - 这需要用户在"精度（偏硬）和实时"之间做取舍，网格方案里会给出几组候选。
- **非结构网格的不对称**：可能重新带来侧向漂移和落点依赖，必须重测（第 1 步曾在规则网格上发现约 30% 的漂移，是改用 swapping=True 才解决的）。
- 没有发现必须引入自适应网格、高阶单元或重划分的证据：现有路线（静态网格 + 重心约束）仍然可行。
