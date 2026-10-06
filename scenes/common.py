"""项目公共工具：SOFA 环境、组织立方体构建。

用法：在脚本开头 `import common`（scenes/ 在 sys.path 中）。
约定（0c）：SI 单位（m, kg, s, N, Pa）。
"""
import sys
sys.dont_write_bytecode = True
import os

SOFA_ROOT = os.path.expanduser("~/sofa/SOFA_v25.12.00_Linux")
os.environ.setdefault("SOFA_ROOT", SOFA_ROOT)

import numpy as np
import Sofa
import Sofa.Core
import Sofa.Simulation
import SofaRuntime

SofaRuntime.importPlugin("Sofa.Component")

# 本项目的修补插件（plugins/NeedleSimFixes，编译到仓库外）
NEEDLESIM_PLUGIN_DIR = os.path.expanduser("~/sofa/needle_build/NeedleSimFixes")
SofaRuntime.PluginRepository.addFirstPath(NEEDLESIM_PLUGIN_DIR)
SofaRuntime.importPlugin("NeedleSimFixes")


def _vec3(v):
    return np.array(v if np.ndim(v) else [v, v, v], dtype=float)


def grid_nodes(L, n, origin=(0.0, 0.0, 0.0)):
    """RegularGridTopology 的节点顺序：x 最快，其次 y，最后 z。返回 (nx*ny*nz, 3) 坐标。

    L、n 可以是标量（立方体）或三元组（长方体：各方向尺寸 / 节点数）。
    """
    L3, n3, o = _vec3(L), _vec3(n).astype(int), _vec3(origin)
    sx, sy, sz = (np.linspace(o[i], o[i] + L3[i], n3[i]) for i in range(3))
    z, y, x = np.meshgrid(sz, sy, sx, indexing="ij")
    return np.stack([x.ravel(), y.ravel(), z.ravel()], axis=1)


def add_tissue(parent, L=0.08, n=9, E=5000.0, nu=0.3, rho=1000.0, method="polar",
               name="Tissue", origin=(0.0, 0.0, 0.0), ff_type="FTCfix", swapping=False, mesh=None):
    """8 cm 立方体（默认），n×n×n 节点的规则网格切成四面体，共旋线弹性。

    mesh：给定时改用这套四面体网格（dict，X：节点坐标 (N, 3)，tets：四面体 (M, 4)，有符号体积为正），
    例如第 1.5 步 gmsh 生成的非均匀网格（scenes/mesh/make_graded_mesh.py 输出的 npz）；此时忽略 L、n、origin、swapping。
    parent 节点下不放求解器，由调用者决定（静力 / 动力）。
    返回 (tissue 节点, MechanicalObject, 力场)。
    """
    if mesh is not None:
        t = parent.addChild(name)
        t.addObject("TetrahedronSetTopologyContainer", name="topo",
                    position=np.asarray(mesh["X"], float).tolist(), tetrahedra=np.asarray(mesh["tets"], int).tolist())
        t.addObject("TetrahedronSetTopologyModifier", name="modifier")
    else:
        # 规则网格节点必须先于组织节点创建（初始化按场景图顺序进行，Hexa2Tetra 需要六面体拓扑已初始化）
        o, L3, n3 = _vec3(origin), _vec3(L), _vec3(n).astype(int)
        topo = parent.addChild(name + "Grid")
        topo.addObject("RegularGridTopology", name="hexa", n=n3.tolist(),
                       min=o.tolist(), max=(o + L3).tolist())
        t = parent.addChild(name)
        t.addObject("TetrahedronSetTopologyContainer", name="topo",
                    position="@../%sGrid/hexa.position" % name)
        t.addObject("TetrahedronSetTopologyModifier", name="modifier")
        # swapping=True：相邻六面体的切分方向交替（见 step1 结果：swapping=False 时网格不对称，
        # 中心节点受竖直点载荷会侧向漂移约 24%）
        t.addObject("Hexa2TetraTopologicalMapping", input="@../%sGrid/hexa" % name,
                    output="@topo", swapping=swapping)
    mo = t.addObject("MechanicalObject", name="dofs", template="Vec3d")
    # 质量：静力测试不用，但动力测试需要；密度 × 体积
    t.addObject("MeshMatrixMass", name="mass", massDensity=rho, topology="@topo")
    cls = {"FTCfix": "FastTetrahedralCorotationalForceFieldFixed",
           "FTC": "FastTetrahedralCorotationalForceField", "TFEM": "TetrahedronFEMForceField"}[ff_type]
    ff = t.addObject(cls, name="fem", youngModulus=E, poissonRatio=nu, method=method, topology="@topo")
    return t, mo, ff


def load_mesh(name_or_path):
    """读取 scenes/mesh/make_graded_mesh.py 输出的 npz（给名字时到 meshes/step1_5/ 下找）。返回 dict（X、tets、参数）。"""
    p = name_or_path if name_or_path.endswith(".npz") else os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "../meshes/step1_5", name_or_path + ".npz")
    d = np.load(p)
    return {k: d[k] for k in d.files}


def add_static_solver(node, linear="cg", newton_iters=50, abs_tol=1e-20):
    """静力求解：StaticSolver + NewtonRaphsonSolver + 线性求解器。

    linear: "cg"（不组装矩阵）/ "ldl"（SparseLDLSolver）/ "eigen"（EigenSimplicialLDLT）。
    原版 FTC 的组装矩阵有 bug，只能用 "cg"；修复后的 FTCfix 三种都可以（见 docs/plan/step0_results.md）。
    abs_tol：残差范数（N）的绝对停止阈值；按问题的力的量级设（例如 1e-9 × 典型力）。
    阈值过严（低于舍入误差下限）时 Newton 会跑满 newton_iters 次。
    """
    node.addObject("NewtonRaphsonSolver", name="newton",
                   maxNbIterationsNewton=newton_iters,
                   absoluteResidualStoppingThreshold=abs_tol,
                   relativeInitialStoppingThreshold=1e-14,
                   relativeSuccessiveStoppingThreshold=1e-14,
                   relativeEstimateDifferenceThreshold=1e-14,
                   absoluteEstimateDifferenceThreshold=1e-16,
                   warnWhenLineSearchFails=False, warnWhenDiverge=False)
    node.addObject("StaticSolver", name="static", newtonSolver="@newton")
    if linear == "cg":
        node.addObject("CGLinearSolver", name="linsolver", iterations=5000,
                       tolerance=1e-24, threshold=1e-30)
    elif linear == "eigen":
        node.addObject("EigenSimplicialLDLT", name="linsolver",
                       template="CompressedRowSparseMatrixMat3x3d")
    else:
        node.addObject("SparseLDLSolver", name="linsolver",
                       template="CompressedRowSparseMatrixMat3x3d")


def add_dynamic_solver(node, rayleigh_mass=0.0, rayleigh_stiffness=0.0, linear="ldl"):
    """动力学：EulerImplicitSolver（后向 Euler，每步一次线性化）+ 直接线性求解器。

    注意〔代码，EulerImplicitSolver.cpp:92,129,145,160〕：每步开始时在当前 (x_n, v_n) 上算内力 f，
    存进 MechanicalObject 的 force（不被投影、不含 Rayleigh 项）。所以第 n+1 步 animate 之后读 force，
    得到的是 x_{n+1}（第 n 步结束时的位置）上的内力。
    """
    node.addObject("EulerImplicitSolver", name="odesolver",
                   rayleighMass=rayleigh_mass, rayleighStiffness=rayleigh_stiffness)
    if linear == "eigen":
        node.addObject("EigenSimplicialLDLT", name="linsolver",
                       template="CompressedRowSparseMatrixMat3x3d")
    elif linear == "async":   # 0d 的 R3：另一个线程做分解，求解用上一次完成的分解
        node.addObject("AsyncSparseLDLSolver", name="linsolver",
                       template="CompressedRowSparseMatrixMat3x3d")
    else:
        node.addObject("SparseLDLSolver", name="linsolver",
                       template="CompressedRowSparseMatrixMat3x3d")


def internal_force(root, node, constraints):
    """读当前构形下的节点内力（只用于静力测试的最后一步，会拆掉约束）。

    投影约束会把受约束自由度上的力清零，所以先移除 constraints，再让 Newton 迭代 0 次：
    只计算一次内力、不移动节点。返回 (N,3) 数组：f = 内力（外载为 0 时，反力 = -f）。
    """
    for c in constraints:
        node.removeObject(c)
    root.newton.maxNbIterationsNewton = 0
    mo = node.dofs
    x0 = mo.position.array().copy()
    Sofa.Simulation.animate(root, root.dt.value)
    if np.all(np.isfinite(x0)):  # 发散（NaN）的对照工况跳过检查
        assert np.abs(mo.position.array() - x0).max() == 0.0
    return mo.force.array().copy()


class ForceProbe:
    """独立的"力探针"场景：和组织同样的网格和材料，没有约束，Newton 迭代 0 次。

    force(x) 返回位置 x 上的节点内力（N,3），不改变任何状态。用来在动力学中读一致状态下的力，
    以及用有限差分算 K·v ≈ -(f(x + εv) - f(x))/ε（Rayleigh 刚度阻尼力 = -r_s K v）。
    """

    def __init__(self, **tissue_kw):
        self.root = make_root()
        add_static_solver(self.root, linear="cg", newton_iters=0)
        self.node, self.mo, _ = add_tissue(self.root, **tissue_kw)
        Sofa.Simulation.init(self.root)

    def force(self, x):
        self.mo.position.value = x
        Sofa.Simulation.animate(self.root, self.root.dt.value)
        return self.mo.force.array().copy()

    def K_times(self, x, v, eps=1e-4):
        return -(self.force(x + eps * v) - self.force(x)) / eps


def quat_from_axis_angle(axis, angle):
    """SOFA 的四元数顺序是 [x, y, z, w]。"""
    a = np.asarray(axis, float); a = a / np.linalg.norm(a)
    return np.r_[np.sin(angle / 2) * a, np.cos(angle / 2)]


def quat_mul(q1, q2):
    """四元数乘法 q1*q2（[x, y, z, w] 顺序），表示先转 q2 再转 q1。"""
    x1, y1, z1, w1 = q1; x2, y2, z2, w2 = q2
    return np.array([w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
                     w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
                     w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2,
                     w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2])


def quat_to_mat(q):
    x, y, z, w = q
    return np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                     [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                     [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])


def add_kinematic_needle(parent, length=0.10, n_body=21, name="Needle", pose0=None, local_points=None):
    """纯运动学的刚性针（0b 方案 C）：一个 Rigid3 节点，没有 ODE 求解器，也没有 ConstraintCorrection。

    局部坐标：针根在原点，针轴沿局部 -z，针尖在 (0, 0, -length)。
    子节点 Body：n_body 个针身点（含针根和针尖），Tip：针尖点；都用 RigidMapping 挂在针根刚体上。
    pose0 = (p[3], q[4])：初始位姿，应等于轨迹在 t = 0 的位姿（否则第一步会有一次姿态跳变，
    映射的雅可比在那一步是错的，见 0b-1）。子节点坐标按局部坐标给出（globalToLocalCoords=False）。
    位姿由 NeedleDriver 控制器每步设定。返回 (needle 节点, 针根 MO, 针身 MO, 针尖 MO, 局部坐标数组)。
    """
    p0, q0 = (np.zeros(3), np.array([0, 0, 0, 1.0])) if pose0 is None else (np.asarray(pose0[0], float), np.asarray(pose0[1], float))
    nd = parent.addChild(name)
    base = nd.addObject("MechanicalObject", name="base", template="Rigid3d", position=[np.r_[p0, q0].tolist()])
    if local_points is None:
        local = np.stack([np.zeros(n_body), np.zeros(n_body), -np.linspace(0, length, n_body)], 1)
    else:                                   # 自定义针身点（局部坐标），例如和组织里的约束点一一对应
        local = np.asarray(local_points, float)
    body_node = nd.addChild("Body")
    body = body_node.addObject("MechanicalObject", name="dofs", template="Vec3d", position=local.tolist())
    body_node.addObject("RigidMapping", input="@../base", output="@dofs", globalToLocalCoords=False)
    tip_node = nd.addChild("Tip")
    tip = tip_node.addObject("MechanicalObject", name="dofs", template="Vec3d", position=[[0, 0, -length]])
    tip_node.addObject("RigidMapping", input="@../base", output="@dofs", globalToLocalCoords=False)
    return nd, base, body, tip, local


def quat_inv(q):
    return np.array([-q[0], -q[1], -q[2], q[3]]) / np.dot(q, q)


def rotvec_from_quat(q):
    """单位四元数 → 转动向量（轴 × 角）。"""
    q = q / np.linalg.norm(q)
    if q[3] < 0:
        q = -q
    s = np.linalg.norm(q[:3])
    if s < 1e-15:
        return 2 * q[:3]
    return 2 * np.arctan2(s, q[3]) * q[:3] / s


def add_needle_visual(nd, length, radius=0.0008, segments=12, color=(0.25, 0.25, 0.3, 1.0)):
    """给纯运动学针加显示模型：一根圆柱（局部坐标从针根 (0,0,0) 到针尖 (0,0,-length)），
    用 RigidMapping 挂在针根刚体上。只用于显示，不参与力学和接触。
    放在单独的子节点里（一个节点只允许一个状态对象 / 映射，见 step1_results.md 发现 8）。"""
    ang = np.linspace(0, 2 * np.pi, segments, endpoint=False)
    ring = np.stack([radius * np.cos(ang), radius * np.sin(ang)], 1)
    pts = [[x, y, 0.0] for x, y in ring] + [[x, y, -length] for x, y in ring] + [[0, 0, 0.0], [0, 0, -length]]
    tris = []
    for i in range(segments):
        j = (i + 1) % segments
        tris += [[i, j, segments + j], [i, segments + j, segments + i]]          # 侧面
        tris += [[2 * segments, j, i], [2 * segments + 1, segments + i, segments + j]]  # 两个端面
    v = nd.addChild("Visual")
    v.addObject("OglModel", name="visual", position=pts, triangles=tris, color=list(color))
    v.addObject("RigidMapping", input="@../base", output="@visual", globalToLocalCoords=False)
    return v


class NeedleDriver(Sofa.Core.Controller):
    """每步开始时（AnimateBeginEvent）设定纯运动学针的位姿。

    trajectory(t) 返回 (位置 p[3], 四元数 q[4], ...)，世界坐标；速度由相邻两步的位姿差分得到（步内平均速度）。
    mode="end"：  pos = 指令(t_{n+1})，freeVel = 0       → freePos = 指令(t_{n+1})（精确）
    mode="start"：pos = 指令(t_n)，  freeVel = 步内平均速度 → freePos = pos + freeVel·dt（转动时是一阶近似）
    依据〔代码〕：FreeMotionAnimationLoop.cpp:417 对所有状态（包括映射的子节点）计算 freePos = pos + freeVel·dt；
          没有求解器的物体，freeVel 不会被自由运动更新（:406-411）。
    子节点必须同时写入〔运行，0b-1〕：子节点的 position 要到一步末尾（:291 propagateXAndV）才由映射更新，
          而 :417 用的是子节点当时的 position，所以只改针根时，子节点的 freePos 会落后一步。
          因此这里按刚体公式 x_i = p + R r_i 直接写入 children = [(子节点 MO, 局部坐标), ...]。
          映射仍然负责把力从子节点传回针根（applyJT）。
    时间用步数 k·dt 计算（避免累加误差）。
    """

    def __init__(self, base, trajectory, dt, children=(), mode="end", **kw):
        super().__init__(**kw)
        self.base, self.traj, self.dt, self.mode, self.k = base, trajectory, dt, mode, 0
        self.children = list(children)

    def _set(self, p, q, v, w):
        self.base.position.value = [np.r_[p, q].tolist()]
        self.base.velocity.value = [np.r_[v, w].tolist()]
        self.base.free_velocity.value = [np.r_[v, w].tolist()]
        R = quat_to_mat(q)
        for mo, local in self.children:
            x = p + local @ R.T
            vx = v + np.cross(w, local @ R.T)
            mo.position.value = x
            mo.velocity.value = vx
            mo.free_velocity.value = vx
            mo.free_position.value = x + vx * self.dt

    def onAnimateBeginEvent(self, event):
        t_n = self.k * self.dt
        p0, q0 = self.traj(t_n)[:2]
        p1, q1 = self.traj(t_n + self.dt)[:2]
        if self.mode == "end":
            self._set(p1, q1, np.zeros(3), np.zeros(3))
        else:
            v = (p1 - p0) / self.dt
            w = rotvec_from_quat(quat_mul(q1, quat_inv(q0))) / self.dt
            self._set(p0, q0, v, w)
        self.k += 1


def needle_wrench(child_mos, lam, dt, p_base):
    """纯运动学针受到的约束合力和力矩（作用在针上，世界坐标，力矩以针根原点 p_base 为参考点）。

    child_mos：挂在针根下的 Vec3 子节点 MO（针身点、针尖），约束建在它们上面；
    lam：约束求解器的 constraintForces（λ，二阶积分器下是冲量 N·s），力 = λ / dt（0c 约定）；
    每个子节点的 constraint 数据（scipy 稀疏矩阵，行 = 约束编号，列 = 3·点号 + 分量）就是 J，
    点上的力 f_i = (Jᵀ λ)_i / dt；F = Σ f_i，τ = Σ (x_i - p_base) × f_i，x_i 用 free_position（约束作用的构形）。
    不读针根刚体自己的 constraint 数据：SofaPython3 把 Rigid3 的约束矩阵转成 scipy 时出错（得到 0×0 并报错），见 0b-2。
    返回 (F[3], τ[3], 各子节点的点力列表)。
    """
    lam = np.asarray(lam, float)
    F, tau, per = np.zeros(3), np.zeros(3), []
    for mo in child_mos:
        J = mo.constraint.value
        x = mo.free_position.array()
        n, ncol = J.shape
        if n == 0 or ncol == 0 or J.nnz == 0:          # 这个子节点上没有约束
            per.append(np.zeros_like(x)); continue
        lam_pad = np.zeros(n); m = min(n, len(lam)); lam_pad[:m] = lam[:m]
        g = np.zeros(3 * len(x)); g[:ncol] = J.T @ lam_pad   # 列数可能只到最后一个被约束的点
        f = g.reshape(-1, 3) / dt
        F += f.sum(0)
        tau += np.cross(x - p_base, f).sum(0)
        per.append(f)
    return F, tau, per


def make_root(gravity=(0.0, 0.0, 0.0), dt=0.01, loop="default", root=None):
    """loop="default"：DefaultAnimationLoop（静力、无约束动力学）；
    loop="free"：FreeMotionAnimationLoop + BlockGaussSeidelConstraintSolver（拉格朗日约束，第 1 步要用的配置）。
    root：给定时（例如 runSofa 的 createScene 传进来的根节点）就在它上面设置，否则新建。"""
    if root is None:
        root = Sofa.Core.Node("root")
    root.gravity = list(gravity)
    root.dt = dt
    if loop == "free":
        root.addObject("FreeMotionAnimationLoop")
        root.addObject("BlockGaussSeidelConstraintSolver", name="csolver", tolerance=1e-12, maxIterations=1000,
                       computeConstraintForces=True)
    else:
        root.addObject("DefaultAnimationLoop")
    return root
