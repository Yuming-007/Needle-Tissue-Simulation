"""第 1.5 步 M2：非均匀网格上的实时预算复测。见 docs/plan/step1_5_mesh_design.md §5。

求解路线 R3（AsyncSparseLDLSolver + LinearSolverConstraintCorrection），dt = 0.01，FTCfix polar，ν = 0.45。
两种负载：
  TIP：a_eff = 10 mm 圆盘针尖（点距 = 标称 h_local；规则网格用 h），从表面上方 0.5 mm 以 5 mm/s 压到 5 mm
       （第 2 步刺穿前的情形）。只统计接触阶段的步（压深 ≥ 0.5 mm）。
  W2：和 0d 相同的负载类型：组织内一条竖直线上 7 个嵌入材料点（BarycentricMapping），和纯运动学针身上的对应点用
       BilateralLagrangianConstraint 相连（21 行），针以 5 mm/s 沿轴向运动（相当于刺穿后完全粘住的针道）。
       线的位置 x = L/2 + 0.31·4 mm、y = L/2 + 0.17·4 mm，点距 10 mm（7 个点 = 0d n = 9 的 21 行）。
网格：gmsh 候选 A-、A_far16、A_g05、A、A+、B、E，以及规则网格 n = 9、11（和 0d 对照，校准机器状态）。
计时：Sofa.Simulation.animate 的墙钟时间（含 Python 控制器），中位数、p95、最大值。
运行：python3 m2_timing.py [网格名 ...] → results/step1_5/m2_timing.json
"""
import sys
sys.dont_write_bytecode = True
import os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../step1"))
import numpy as np
import Sofa, Sofa.Simulation
import common, contact
import t1_vertex as V

contact.load_plugins()
L, LEN, dt, A = V.L, V.LEN, 0.01, 0.010
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../results/step1_5")
MESHES = ["regular_n9", "regular_n11", "A_far16", "A-", "A_g05", "A", "A+", "B", "E"]


def tissue(root, name):
    if name.startswith("regular"):
        n = int(name.split("_n")[1])
        return contact.add_tissue_with_surface(root, n=n, L=L, E=V.E, nu=V.NU, linear="async"), (L / (n - 1))
    m = common.load_mesh(name)
    return contact.add_tissue_with_surface(root, L=L, E=V.E, nu=V.NU, linear="async", mesh=m), float(m["h_local"])


def run_tip(name, v=0.005, depth=0.005, gap=0.0005):
    zf, T = V.piecewise(L + gap, [(gap / v, -v), (depth / v, -v)])
    root = common.make_root(dt=dt, loop="free")
    traj = lambda t: (np.array([L / 2, L / 2, zf(t) + LEN]), V.Q)
    nd, base, body, tip, local = common.add_kinematic_needle(root, length=LEN, n_body=V.NB, pose0=traj(0.0))
    root.addObject("CollisionLoop")
    tis, h = tissue(root, name)
    pts = contact.tip_patch_points(LEN, A, h, "disk")
    pmo, pg = contact.add_tip_patch(nd, pts)
    root.addObject(common.NeedleDriver(base, traj, dt, children=[(body, local), (tip, local[-1:]), (pmo, pts)], name="driver"))
    _, bg = contact.add_needle_geometry(nd, body, tip, V.NB)
    contact.add_tip_contact(root, pg, bg, tis, distance=0.02)
    Sofa.Simulation.init(root)
    ts, rows = [], []
    for k in range(int(round(T / dt))):
        a = time.perf_counter(); Sofa.Simulation.animate(root, dt); te = time.perf_counter() - a
        if L - zf((k + 1) * dt) >= 0.0005:
            ts.append(te); rows.append(len(root.csolver.constraintForces.value))
    nodes = len(tis["dofs"].position.value)
    Sofa.Simulation.unload(root)
    return nodes, len(pts), np.array(ts), int(np.median(rows))


def run_w2(name, nsteps=200, skip=20):
    root = common.make_root(dt=dt, loop="free")
    tis, h = tissue(root, name)
    xl, yl = L / 2 + 0.31 * 0.004, L / 2 + 0.17 * 0.004
    pts = np.array([[xl, yl, L - (i + 0.5) * 0.010] for i in range(7)])
    base0 = np.array([xl, yl, L + 0.02])
    local = pts - base0
    traj = lambda tt: (base0 + np.array([0, 0, -0.005 * tt]), np.array([0, 0, 0, 1.0]))
    nd, base, body, tip, _ = common.add_kinematic_needle(root, pose0=traj(0.0), local_points=local)
    root.addObject(common.NeedleDriver(base, traj, dt, children=[(body, local), (tip, local[-1:])], name="driver"))
    emb = tis["node"].addChild("Embedded")
    emo = emb.addObject("MechanicalObject", name="dofs", template="Vec3d", position=pts.tolist())
    emb.addObject("BarycentricMapping", input="@../dofs", output="@dofs")
    root.addObject("BilateralLagrangianConstraint", object1=body.getLinkPath(), object2=emo.getLinkPath(),
                   first_point=list(range(7)), second_point=list(range(7)))
    Sofa.Simulation.init(root)
    ts = []
    for k in range(nsteps):
        a = time.perf_counter(); Sofa.Simulation.animate(root, dt); ts.append(time.perf_counter() - a)
    nodes = len(tis["dofs"].position.value)
    Sofa.Simulation.unload(root)
    return nodes, np.array(ts[skip:])


def stats(ts):
    return dict(med_ms=1e3 * float(np.median(ts)), p95_ms=1e3 * float(np.percentile(ts, 95)), max_ms=1e3 * float(ts.max()), n=int(len(ts)))


if __name__ == "__main__":
    out = []
    for name in (sys.argv[1:] or MESHES):
        nodes, npts, ts_tip, rows = run_tip(name)
        _, ts_w2 = run_w2(name)
        r = dict(name=name, nodes=nodes, tip_points=npts, tip_rows=rows, tip=stats(ts_tip), w2=stats(ts_w2))
        out.append(r)
        print(f"{name:12s} 节点 {nodes:5d}：TIP（{npts} 点，约 {rows} 行）中位数 {r['tip']['med_ms']:6.2f} ms，p95 {r['tip']['p95_ms']:6.2f}，"
              f"最大 {r['tip']['max_ms']:6.2f}；W2（21 行）中位数 {r['w2']['med_ms']:6.2f} ms，p95 {r['w2']['p95_ms']:6.2f}，最大 {r['w2']['max_ms']:6.2f}",
              flush=True)
        json.dump(out, open(os.path.join(OUT, "m2_timing" + ("_extra" if len(sys.argv) > 1 else "") + ".json"), "w"), indent=1)
