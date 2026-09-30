"""0a-1 补充：组装矩阵（直接求解器）与不组装矩阵（CG）的一致性。

目的：验证 FastTetrahedralCorotationalForceFieldFixed（plugins/NeedleSimFixes）的 buildStiffnessMatrix。
工况：底面全固定；顶面绕 z 轴扭转 θ 并下压 c·L（大变形、大转动，所有单元 R ≠ I）。
读数：LDL 与 CG 最终位置的最大差 / L；Newton 迭代次数；残差。
对照：原版 FastTetrahedralCorotationalForceField（v25.12，组装矩阵有 bug）。
运行：python3 t0a1_solver_consistency.py → results/step0/t0a1_solver_consistency.json
"""
import sys
sys.dont_write_bytecode = True
import os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
import Sofa, Sofa.Simulation
import common

L, E = 0.08, 5000.0


def run(ff_type, method, nu, linear, twist_deg, comp, n=9):
    root = common.make_root()
    common.add_static_solver(root, linear=linear, newton_iters=100)
    t, mo, _ = common.add_tissue(root, L=L, n=n, E=E, nu=nu, method=method, ff_type=ff_type)
    X = common.grid_nodes(L, n)
    bot = np.where(np.isclose(X[:, 2], 0))[0]
    top = np.where(np.isclose(X[:, 2], L))[0]
    t.addObject("FixedProjectiveConstraint", indices=bot.tolist())
    t.addObject("FixedProjectiveConstraint", indices=top.tolist(), name="topfix")
    Sofa.Simulation.init(root)
    th = np.radians(twist_deg); c = np.array([L / 2, L / 2, 0])
    Rz = np.array([[np.cos(th), -np.sin(th), 0], [np.sin(th), np.cos(th), 0], [0, 0, 1]])
    p = X.copy()
    p[top] = (X[top] - c) @ Rz.T + c
    p[top, 2] = L * (1 - comp)
    mo.position.value = p
    Sofa.Simulation.animate(root, root.dt.value)
    # residualGraph 的字符串形式："residual r0 r1 ... rk"，r0 是初始残差 → 迭代次数 = k
    res = root.newton.residualGraph.getValueString().split()[1:]
    iters = len(res) - 1
    return mo.position.array().copy(), str(root.newton.status.value), iters


if __name__ == "__main__":
    out = []
    for ff_type in ["FTCfix", "FTC"]:
        for method in ["polar", "none"]:
            for nu in [0.3, 0.49]:
                for twist, comp in [(0, 0.01), (45, 0.10)]:
                    xl, sl, il = run(ff_type, method, nu, "ldl", twist, comp)
                    xc, sc, ic = run(ff_type, method, nu, "cg", twist, comp)
                    d = float(np.abs(xl - xc).max() / L) if np.all(np.isfinite(xl)) else float("inf")
                    r = dict(ff=ff_type, method=method, nu=nu, twist_deg=twist, comp=comp,
                             ldl_status=sl, ldl_iters=il, cg_status=sc, cg_iters=ic, max_diff_over_L=d)
                    out.append(r)
                    print(f"{ff_type:6s} {method:5s} nu={nu:.2f} twist={twist:2d} comp={comp:.2f}: "
                          f"LDL {sl:30s} it={il:3d} | CG {sc:30s} it={ic:3d} | max|x_ldl-x_cg|/L={d:.2e}", flush=True)
    res = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../results/step0/t0a1_solver_consistency.json")
    json.dump(out, open(res, "w"), indent=1)
