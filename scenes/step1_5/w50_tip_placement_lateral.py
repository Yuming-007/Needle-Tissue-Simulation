"""单点针尖落点对侧向力和可压深度的影响（W50 / W50c，轻量检查，2026-10-09）。

三种落点（针竖直下压，5 mm/s，压到 8 mm，不刺穿；R1，dt = 0.01）：
  W50  中心：针对准顶面中心 (40, 40) mm，那里没有节点（落在三角形内，约 82% 的力压在 0.8 mm 外的节点上）；
  W50c 节点：网格在顶面中心嵌入了一个节点，针正好压在节点上；
  W50  重心：针对准"中心所在三角形"的重心（三个节点各 1/3）。
每个算例输出：压深 0.5–8 mm 时的轴向力、侧向力 / 轴向力、最差单元体积比（< 0 表示翻转）。
运行：python3 w50_tip_placement_lateral.py → results/step1_5/w50_tip_placement_lateral.log（终端输出）
"""
import sys
sys.dont_write_bytecode = True
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../step1"))
import numpy as np
import Sofa, Sofa.Simulation
import common, contact
import t1_vertex as V
from v4_large_depth import tet_volumes

contact.load_plugins()
L, LEN, dt = V.L, V.LEN, 0.01


def centroid_of_center_triangle(mesh):
    X, t = mesh["X"], mesh["tets"]
    c = np.array([L / 2, L / 2])
    top = np.isclose(X[:, 2], L)
    faces = set()
    for tt in t:
        s = [v for v in tt if top[v]]
        if len(s) == 3:
            faces.add(tuple(sorted(s)))
    for tri in faces:
        A = X[list(tri)][:, :2]
        T = np.array([A[0] - A[2], A[1] - A[2]]).T
        l = np.linalg.solve(T, c - A[2])
        if min(l[0], l[1], 1 - l.sum()) >= -1e-12:
            return A.mean(0)


def run(mesh, xy, depth=0.008, v=0.005, gap=0.0005):
    pts = np.array([[0.0, 0.0, -LEN]])
    zf, T = V.piecewise(L + gap, [(gap / v, -v), (depth / v, -v)])
    root = common.make_root(dt=dt, loop="free")
    traj = lambda tt: (np.array([xy[0], xy[1], zf(tt) + LEN]), V.Q)
    nd, base, body, tip, local = common.add_kinematic_needle(root, length=LEN, n_body=V.NB, pose0=traj(0.0))
    pmo, pg = contact.add_tip_patch(nd, pts)
    root.addObject(common.NeedleDriver(base, traj, dt, children=[(body, local), (tip, local[-1:]), (pmo, pts)], name="driver"))
    root.addObject("CollisionLoop")
    tis = contact.add_tissue_with_surface(root, L=L, E=V.E, nu=V.NU, linear="ldl", mesh=mesh)
    _, bg = contact.add_needle_geometry(nd, body, tip, V.NB)
    contact.add_tip_contact(root, pg, bg, tis, distance=0.02)
    Sofa.Simulation.init(root)
    tets = np.array(tis["node"].topo.tetrahedra.value)
    X0 = np.asarray(tis["X"])
    V0 = tet_volumes(X0, tets)
    out = []
    for k in range(int(round(T / dt))):
        Sofa.Simulation.animate(root, dt)
        lam = np.array(root.csolver.constraintForces.value)
        z = zf((k + 1) * dt)
        F, _, _ = common.needle_wrench([pmo], lam, dt, np.array([xy[0], xy[1], z]))
        vr = (tet_volumes(tis["dofs"].position.array(), tets) / V0).min()
        out.append(((L - z) * 1e3, F.copy(), vr))
        if not np.isfinite(F).all() or vr < -0.5:
            break
    Sofa.Simulation.unload(root)
    return out


if __name__ == "__main__":
    w50, w50c = common.load_mesh("W50"), common.load_mesh("W50c")
    cen = centroid_of_center_triangle(w50)
    cases = [("W50 中心（三角形内，靠近节点）", w50, np.array([L / 2, L / 2])),
             ("W50c 节点", w50c, np.array([L / 2, L / 2])),
             (f"W50 三角形重心 ({cen[0]*1e3:.2f}, {cen[1]*1e3:.2f}) mm", w50, cen)]
    for name, mesh, xy in cases:
        out = run(mesh, xy)
        print(f"== {name}", flush=True)
        for dd in [0.5, 1, 2, 3, 4, 5, 6, 7, 8]:
            d, F, vr = min(out, key=lambda o: abs(o[0] - dd))
            if abs(d - dd) > 0.06:
                continue
            print(f"   {d:.1f} mm：F_z {F[2]:.3f} N，侧向/轴向 {np.hypot(F[0], F[1])/abs(F[2]):.2e}，最差体积比 {vr:.2f}", flush=True)
        print(f"   （到达 {out[-1][0]:.1f} mm）", flush=True)
