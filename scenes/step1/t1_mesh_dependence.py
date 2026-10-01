"""第 1 步 P9：点接触刚度的网格依赖（静态、固定压深）。见 docs/plan/step1_tenting.md。

P1 已证明"顶点接触 ≡ 给同一顶点施加位移（x、y 自由）"，所以这里直接给顶面中心节点施加 z 位移，做静力求解。
全模型：h = 10、5 mm（n = 9、17）；四分之一对称模型（x ≥ L/2、y ≥ L/2，对称面上固定法向位移，力 × 4）：
h = 10、5、2.5 mm。都用 swapping=True。压深：全模型 1、5、10 mm；四分之一模型 1、2、5 mm（细网格避免单元严重扭曲）。
记录：反力 F、割线刚度 F/δ、最差单元的体积比 V/V0（< 0 表示翻转）。
运行：python3 t1_mesh_dependence.py → results/step1/t1_mesh_dependence.json
"""
import sys
sys.dont_write_bytecode = True
import os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
import Sofa, Sofa.Simulation
import common

L, E, NU = 0.08, 5000.0, 0.45
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../results/step1")


def tet_volumes(x, tets):
    a, b, c, d = (x[tets[:, i]] for i in range(4))
    return np.einsum("ij,ij->i", np.cross(b - a, c - a), d - a) / 6.0


def run(h, depth, quarter):
    N = int(round(L / h))
    if quarter:
        size, nn, org = (L / 2, L / 2, L), (N // 2 + 1, N // 2 + 1, N + 1), (L / 2, L / 2, 0.0)
    else:
        size, nn, org = L, N + 1, (0.0, 0.0, 0.0)
    root = common.make_root(); common.add_static_solver(root, linear="ldl", abs_tol=1e-11, newton_iters=100)
    t, mo, _ = common.add_tissue(root, L=size, n=nn, E=E, nu=NU, origin=org, swapping=True)
    X = common.grid_nodes(size, nn, org)
    bot = np.where(np.isclose(X[:, 2], 0))[0]
    c = int(np.where(np.all(np.isclose(X, [L / 2, L / 2, L]), axis=1))[0][0])
    cons = [t.addObject("FixedProjectiveConstraint", indices=bot.tolist()),
            t.addObject("PartialFixedProjectiveConstraint", name="tip", indices=[c], fixedDirections=[0, 0, 1])]
    if quarter:
        symx = np.setdiff1d(np.where(np.isclose(X[:, 0], L / 2))[0], bot)
        symy = np.setdiff1d(np.where(np.isclose(X[:, 1], L / 2))[0], bot)
        cons.append(t.addObject("PartialFixedProjectiveConstraint", indices=symx.tolist(), fixedDirections=[1, 0, 0]))
        cons.append(t.addObject("PartialFixedProjectiveConstraint", indices=symy.tolist(), fixedDirections=[0, 1, 0]))
    Sofa.Simulation.init(root)
    tets = np.array(t.topo.tetrahedra.value)
    V0 = tet_volumes(X, tets)
    p = X.copy(); p[c, 2] -= depth
    mo.position.value = p
    a = time.perf_counter(); Sofa.Simulation.animate(root, 0.01); t_solve = time.perf_counter() - a
    status = str(root.newton.status.value)
    x = mo.position.array().copy()
    vr = (tet_volumes(x, tets) / V0).min()
    f = common.internal_force(root, t, cons)
    F = (4.0 if quarter else 1.0) * float(f[c, 2])
    Sofa.Simulation.unload(root)
    return dict(h_mm=h * 1e3, depth_mm=depth * 1e3, quarter=quarter, nodes=int(len(X)), F=F, k_secant=F / depth,
                min_volume_ratio=float(vr), status=status, t_solve=t_solve)


if __name__ == "__main__":
    out = []
    cases = [(h, d, False) for h in [0.01, 0.005] for d in [0.001, 0.005, 0.010]]
    cases += [(h, d, True) for h in [0.01, 0.005, 0.0025] for d in [0.001, 0.002, 0.005]]
    for h, d, q in cases:
        r = run(h, d, q); out.append(r)
        print(f"{'四分之一' if q else '全模型'} h={r['h_mm']:4.1f}mm 节点 {r['nodes']:6d} 压深 {r['depth_mm']:4.1f}mm：F = {r['F']:.4f} N，"
              f"割线刚度 {r['k_secant']:7.2f} N/m，最差单元体积比 {r['min_volume_ratio']:.3f}，{r['status']}，{r['t_solve']:.1f} s", flush=True)
        json.dump(out, open(os.path.join(OUT, "t1_mesh_dependence.json"), "w"), indent=1)
