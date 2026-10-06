"""第 1.5 步 V3 的线性 LCP 工具：刚性针尖（多点）压平整组织表面的小压深静力接触。见 step1_5_tip_survey.md §7、
step1_5_results.md V3-0。

和 V1 参照（v1_multipoint_lcp.py）的区别只在于计算方式，参照问题相同：
  - 刚度：直接取 SOFA 组装的静力切线刚度 K（StaticSolver + SparseLDLSolver，未变形构形，一次 Newton 迭代；
    MatrixLinearSystem.get_system_matrix()）。V3-0 已验证 K⁻¹f 和 SOFA 一次迭代的位移一致到 1e-12。
    固定自由度（底面；四分之一模型的对称面法向）从系统中去掉。
  - 接触柔度 W = Bᵀ K⁻¹ B：B 的每一列对应一个接触点（法向 n = +z，权重 = 该点在未变形顶面三角形中的重心坐标），
    K 只分解一次（scipy splu），每个接触点一个右端项（V1 是每个节点分量两次静力求解）。
  - 间隙：表面平整、未变形，刚性平移 δ 时 q_i = z_tip,i − L（圆盘 = −δ）。
  - LCP：μ ≥ 0，g = q + Wμ ≥ 0，μ·g = 0；等价于 min ½μᵀWμ + qᵀμ（μ ≥ 0）。点距细于 h 时 W 可能半正定（冗余），
    加极小的正则 εI（ε = 1e-12·迹/N）后用 Cholesky + NNLS；冗余时只比较合力（各点 μ 不唯一）。
四分之一模型（x ≥ L/2、y ≥ L/2，对称面上固定法向位移）：只保留四分之一区域的点；对称面上的点载荷按一半（两个面
交点按四分之一）施加，记 μ_i = φ_i λ_i；全模型合力 F = 4 Σ μ_i（每个四分之一模型的点代表 4φ_i 个全模型的点：内部 4、对称面上 2、交点 1；所以 Σ_全 λ = Σ 4φ_i λ_i = 4Σμ_i）。
"""
import sys
sys.dont_write_bytecode = True
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
from scipy.optimize import nnls
import Sofa, Sofa.Simulation, Sofa.SofaLinearSystem
import common

L, E, NU = 0.08, 5000.0, 0.45


class Mesh:
    """一个规则网格（全模型或四分之一模型）：节点、四面体、顶面三角形、K 的分解。"""

    def __init__(self, h, quarter=False):
        N = int(round(L / h))
        self.h, self.quarter = h, quarter
        if quarter:
            size, nn, org = (L / 2, L / 2, L), (N // 2 + 1, N // 2 + 1, N + 1), (L / 2, L / 2, 0.0)
        else:
            size, nn, org = L, N + 1, (0.0, 0.0, 0.0)
        X = common.grid_nodes(size, nn, org)
        bot = np.where(np.isclose(X[:, 2], 0))[0]
        root = common.make_root(); common.add_static_solver(root, linear="ldl", abs_tol=1e-30, newton_iters=1)
        t, mo, _ = common.add_tissue(root, L=size, n=nn, E=E, nu=NU, origin=org, swapping=True)
        fixed = [3 * i + c for i in bot for c in range(3)]
        t.addObject("FixedProjectiveConstraint", indices=bot.tolist())
        if quarter:
            symx = np.setdiff1d(np.where(np.isclose(X[:, 0], L / 2))[0], bot)
            symy = np.setdiff1d(np.where(np.isclose(X[:, 1], L / 2))[0], bot)
            t.addObject("PartialFixedProjectiveConstraint", indices=symx.tolist(), fixedDirections=[1, 0, 0])
            t.addObject("PartialFixedProjectiveConstraint", indices=symy.tolist(), fixedDirections=[0, 1, 0])
            fixed += [3 * i for i in symx] + [3 * i + 1 for i in symy]
        top = np.where(np.isclose(X[:, 2], L))[0]
        t.addObject("ConstantForceField", indices=[int(top[0])], forces=[[0, 0, -1e-6]])   # 触发组装
        Sofa.Simulation.init(root)
        Sofa.Simulation.animate(root, 0.01)
        K = sp.csr_matrix(root.MatrixLinearSystem1.get_system_matrix())
        tets = np.array(t.topo.tetrahedra.value)
        Sofa.Simulation.unload(root)
        self.X, self.tets = X, tets
        self.free = np.setdiff1d(np.arange(3 * len(X)), np.unique(fixed))
        self.lu = sla.splu(K[self.free][:, self.free].tocsc(), permc_spec="COLAMD")
        self.ndof = len(self.free)
        # 顶面三角形：四面体中三个节点都在 z = L 的面
        faces = set()
        ontop = np.isclose(X[:, 2], L)
        for tt in tets:
            s = [v for v in tt if ontop[v]]
            if len(s) == 3:
                faces.add(tuple(sorted(s)))
        self.tris = np.array(sorted(faces))
        A = X[self.tris][:, :, :2]
        self._T = np.stack([A[:, 0] - A[:, 2], A[:, 1] - A[:, 2]], axis=2)     # (nt, 2, 2)
        self._c = A[:, 2]

    def locate(self, xy):
        """点 xy 所在的顶面三角形和重心权重（在边 / 顶点上时取任一包含它的三角形，权重相同）。"""
        l = np.linalg.solve(self._T, (xy - self._c)[..., None])[..., 0]
        w = np.column_stack([l, 1 - l.sum(1)])
        ok = np.where((w >= -1e-9).all(1))[0]
        if len(ok) == 0:
            raise ValueError(f"point {xy} not on top face")
        i = ok[0]
        return self.tris[i], np.clip(w[i], 0, None) / np.clip(w[i], 0, None).sum()

    def contact_compliance(self, pts_xy, phi):
        """W = Bᵀ K⁻¹ B（法向 +z）。pts_xy：接触点的 xy；phi：载荷系数（全模型 1；对称面 1/2；交点 1/4）。"""
        B = np.zeros((3 * len(self.X), len(pts_xy)))
        info = []
        for r, xy in enumerate(pts_xy):
            tri, w = self.locate(np.asarray(xy))
            for j, wj in zip(tri, w):
                B[3 * j + 2, r] += wj
            info.append((tri.tolist(), w.tolist()))
        Bf = sp.csc_matrix(B[self.free])
        W = np.zeros((len(pts_xy), len(pts_xy)))
        for c0 in range(0, len(pts_xy), 64):                 # 分块求解，避免 (自由度 × 点数) 的稠密矩阵过大
            c1 = min(c0 + 64, len(pts_xy))
            Y = self.lu.solve(Bf[:, c0:c1].toarray())
            W[:, c0:c1] = Bf.T @ Y
        return 0.5 * (W + W.T), info


class MeshArrays(Mesh):
    """任意四面体网格（全模型，例如 gmsh 非均匀网格），接口同 Mesh。M3 起使用（2026-10-03）。

    X：节点（米），tets：四面体（有符号体积为正）；底面 z = 0 的节点全固定；K 同样取 SOFA 组装的静力切线刚度。
    """

    def __init__(self, X, tets, h=None):
        X = np.asarray(X, float); tets = np.asarray(tets, int)
        self.h, self.quarter = h, False
        bot = np.where(np.isclose(X[:, 2], 0, atol=1e-9))[0]
        root = common.make_root(); common.add_static_solver(root, linear="ldl", abs_tol=1e-30, newton_iters=1)
        t, mo, _ = common.add_tissue(root, E=E, nu=NU, mesh=dict(X=X, tets=tets))
        t.addObject("FixedProjectiveConstraint", indices=bot.tolist())
        top = np.where(np.isclose(X[:, 2], L, atol=1e-9))[0]
        t.addObject("ConstantForceField", indices=[int(top[0])], forces=[[0, 0, -1e-6]])   # 触发组装
        Sofa.Simulation.init(root)
        Sofa.Simulation.animate(root, 0.01)
        K = sp.csr_matrix(root.MatrixLinearSystem1.get_system_matrix())
        Sofa.Simulation.unload(root)
        fixed = [3 * i + c for i in bot for c in range(3)]
        self.X, self.tets = X, tets
        self.free = np.setdiff1d(np.arange(3 * len(X)), fixed)
        self.lu = sla.splu(K[self.free][:, self.free].tocsc(), permc_spec="COLAMD")
        self.ndof = len(self.free)
        faces = set()
        ontop = np.isclose(X[:, 2], L, atol=1e-9)
        for tt in tets:
            s_ = [v for v in tt if ontop[v]]
            if len(s_) == 3:
                faces.add(tuple(sorted(s_)))
        self.tris = np.array(sorted(faces))
        A_ = X[self.tris][:, :, :2]
        self._T = np.stack([A_[:, 0] - A_[:, 2], A_[:, 1] - A_[:, 2]], axis=2)
        self._c = A_[:, 2]


def disk_stiffness_at(mesh, a, s, center, delta=1e-4):
    """刚性平底圆盘（中心 center 的 xy，半径 a，点距 s）压下 δ 的线性接触刚度（全模型）、合力矩、受力点数。"""
    pts = disk_points(a, s, center)
    W, info = mesh.contact_compliance(pts, np.ones(len(pts)))
    q = -delta * np.ones(len(pts))
    mu, g = solve_lcp(W, q)
    F = mu.sum()
    r = pts - np.asarray(center)
    M = np.array([(mu * r[:, 1]).sum(), -(mu * r[:, 0]).sum()])      # 力沿 +z 作用在针上：M = r × F
    act = mu > 1e-9 * mu.max()
    return dict(k=F / delta, M_xy=M.tolist(), cop=(np.array([-M[1], M[0]]) / F).tolist(), n_active=int(act.sum()),
                n_points=int(len(pts)), complementarity=float(np.abs(mu * g).max() / (mu.max() * delta)))


def solve_lcp(W, q):
    n = len(q)
    eps = 1e-12 * np.trace(W) / n
    Ws = W + eps * np.eye(n)
    R = np.linalg.cholesky(Ws).T
    mu, _ = nnls(R, -np.linalg.solve(R.T, q), maxiter=50 * n)
    g = q + W @ mu
    return mu, g


def disk_points(a, s, center):
    """平底圆盘的采样点（xy），同 contact.tip_patch_points(shape="disk") 的布点：圆心 + 同心圆环，点距约 s。"""
    pts = [[0.0, 0.0]]
    nr = max(1, int(round(a / s)))
    for k in range(1, nr + 1):
        r = a * k / nr
        m = max(4, int(round(2 * np.pi * r / s / 4)) * 4)
        for j in range(m):
            ph = 2 * np.pi * j / m
            pts.append([r * np.cos(ph), r * np.sin(ph)])
    return np.array(pts) + np.asarray(center)


def quarter_select(pts, center, tol=1e-12):
    """四分之一区域（x ≥ cx、y ≥ cy）内的点及其载荷系数 φ。"""
    d = pts - np.asarray(center)
    keep = (d[:, 0] >= -tol) & (d[:, 1] >= -tol)
    on_x = np.abs(d[:, 0]) <= tol
    on_y = np.abs(d[:, 1]) <= tol
    phi = np.where(on_x & on_y, 0.25, np.where(on_x | on_y, 0.5, 1.0))
    return pts[keep], phi[keep]


def disk_stiffness(mesh, a, s, delta=1e-4):
    """刚性平底圆盘（半径 a，点距 s）压下 δ 的线性接触刚度 k = F/δ（全模型等效）、受力点比例。"""
    c = np.array([L / 2, L / 2])
    pts = disk_points(a, s, c)
    if mesh.quarter:
        pts, phi = quarter_select(pts, c)
        mult = 4.0
    else:
        phi, mult = np.ones(len(pts)), 1.0
    W, info = mesh.contact_compliance(pts, phi)
    # 载荷 μ_i = φ_i λ_i 施加在四分之一模型上；间隙 g = q + W μ
    q = -delta * np.ones(len(pts))
    mu, g = solve_lcp(W, q)
    F = mult * mu.sum()
    act = mu > 1e-9 * mu.max()
    # 四分之一模型中每个点代表 4φ 个全模型的点（内部 4、对称面上 2、交点 1）
    nfull = (4 * phi).sum() if mesh.quarter else len(pts)
    n_act_full = (4 * phi[act]).sum() if mesh.quarter else act.sum()
    comp = float(np.abs(mu * g).max() / (mu.max() * delta))
    return dict(k=F / delta, n_points=int(round(nfull)), n_active=int(round(n_act_full)), n_rows=len(pts),
                min_gap_rel=float(g.min() / delta), complementarity=comp, W_cond=float(np.linalg.cond(W)))
