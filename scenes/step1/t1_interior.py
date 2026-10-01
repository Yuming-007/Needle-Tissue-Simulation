"""第 1 步 P2：三角形内部接触 vs 同一离散接触雅可比的参照。见 docs/plan/step1_tenting.md。

针尖对准顶面中心附近某个三角形内部的点 P（重心权重 w，三个节点 a、b、c）。
参照（静力、线性）：分别在 a、b、c 上施加 z 向单位力（x、y 自由），得到 3×3 的 z 向柔度矩阵 C，
  k_ref = 1 / (wᵀ C w)（接触力按 w 分到三个节点，约束的是加权位移 Σ w_i u_iz）。
接触仿真：和 t1_vertex 相同的测试台（A1，R1，dt = 0.01），针尖从 P 上方 1 mm 以 0.2 mm/s 压到 0.05 mm，
  用压深 0.01–0.05 mm 的线性拟合得到刚度 k；同时记录接触点权重的变化（无摩擦，接触点可以滑动）。
  第一版用 0.1–0.5 mm 拟合，接触 k 比线性参照低 1.0–1.55%：参照是线性（初始）柔度，而 F(δ) 已经开始变软
  （顶点静力曲线约每毫米降 5%）。缩小拟合区间检验这个解释。另加中心顶点作对照（w = (1,0,0)）。
运行：python3 t1_interior.py → results/step1/t1_interior.json
"""
import sys
sys.dont_write_bytecode = True
import os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
import Sofa, Sofa.Simulation
import common, contact
import t1_vertex as V

L, N, NU, E, h = V.L, V.N, V.NU, V.E, V.L / (V.N - 1)


def top_triangles():
    """顶面上的三角形（和场景里 Tetra2Triangle 提取的相同）及其节点编号。"""
    root = common.make_root(dt=0.01, loop="free")
    tis = contact.add_tissue_with_surface(root, n=N, L=L, E=E, nu=NU)
    Sofa.Simulation.init(root)
    tris = np.array(tis["surface"].tris.triangles.value)
    Sofa.Simulation.unload(root)
    top = np.all(np.isclose(V.X[tris][:, :, 2], L), axis=1)
    return tris[top]


def locate(tris, xy):
    for tri in tris:
        A = V.X[tri][:, :2]
        T = np.array([A[0] - A[2], A[1] - A[2]]).T
        l1, l2 = np.linalg.solve(T, xy - A[2])
        w = np.array([l1, l2, 1 - l1 - l2])
        if (w >= -1e-12).all():
            return tri, w
    raise ValueError("not on top face")


def compliance(nodes, f=1e-3):
    Cm = np.zeros((3, 3))
    for j, nj in enumerate(nodes):
        root = common.make_root(); common.add_static_solver(root, linear="ldl", abs_tol=1e-14)
        t, mo, _ = common.add_tissue(root, L=L, n=N, E=E, nu=NU, swapping=True)
        t.addObject("FixedProjectiveConstraint", indices=V.BOT.tolist())
        t.addObject("ConstantForceField", indices=[int(nj)], forces=[[0, 0, -f]])
        Sofa.Simulation.init(root); Sofa.Simulation.animate(root, 0.01)
        u = mo.position.array() - V.X
        Cm[:, j] = -u[nodes, 2] / f
        Sofa.Simulation.unload(root)
    return Cm


def contact_run(xy, depth=0.00005, gap=0.001, v=0.0002, dt=0.01):
    zf, T = V.piecewise(L + gap, [(gap / v, -v), (depth / v, -v)])
    root, hnd = V.build(dt, zf)
    # 把针移到 xy：重建针的轨迹（build 默认对准中心）
    q = V.Q
    drv = root.driver
    drv.traj = lambda t: (np.array([xy[0], xy[1], zf(t) + V.LEN]), q)
    p0, q0 = drv.traj(0.0)
    hnd["base"].position.value = [np.r_[p0, q0].tolist()]
    Sofa.Simulation.init(root)
    rec = {k: [] for k in ["t", "tip_z", "node", "lam", "F_needle", "tau", "F_tissue", "n", "contact_nodes", "F_bottom_prev"]}
    for k in range(int(round(T / dt))):
        Sofa.Simulation.animate(root, dt)
        V.record(root, hnd, dt, rec, None)
    Sofa.Simulation.unload(root)
    return rec


if __name__ == "__main__":
    out = []
    tris = top_triangles()
    for off in [(0.0, 0.0), (0.30, 0.15), (0.55, 0.25), (-0.40, 0.20)]:
        xy = np.array([L / 2 + off[0] * h, L / 2 + off[1] * h])
        tri, w = locate(tris, xy)
        Cm = compliance(tri)
        k_ref = 1.0 / (w @ Cm @ w)
        rec = contact_run(xy)
        tz = np.array(rec["tip_z"]); F = np.array(rec["F_needle"])[:, 2]
        dlt = L - tz
        m = (dlt >= 1e-5 - 1e-12) & (dlt <= 5e-5 + 1e-12)
        k = np.polyfit(dlt[m], F[m], 1)[0]
        # 接触点权重的变化
        ws = [sorted([abs(nw[1][2]) for nw in r]) for r in rec["contact_nodes"] if r]
        nodes_seen = sorted({nw[0] for r in rec["contact_nodes"] for nw in r})
        lat = float(np.abs(np.array(rec["F_needle"])[m][:, :2]).max() / np.abs(F[m]).max())
        d = dict(offset=off, xy=xy.tolist(), tri=tri.tolist(), w=w.tolist(), C=Cm.tolist(), k_ref=float(k_ref), k=float(k),
                 rel=float(k / k_ref - 1), nodes_seen=nodes_seen, lateral_over_axial=lat)
        out.append(d)
        print(f"偏移 ({off[0]:+.2f}h, {off[1]:+.2f}h)：三角形 {tri.tolist()}，权重 {np.round(w,4).tolist()}；"
              f"k_ref = {k_ref:.3f} N/m，接触 k = {k:.3f} N/m，相对差 {k/k_ref-1:+.2e}；接触过的节点 {nodes_seen}；"
              f"侧向力 / 轴向力 最大 {lat:.2e}", flush=True)
    json.dump(out, open(os.path.join(V.OUT, "t1_interior.json"), "w"), indent=1)
