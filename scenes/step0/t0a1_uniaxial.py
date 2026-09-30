"""0a-1 单轴拉压测试（静力）。见 docs/plan/step0_foundation.md §2.1。

边界条件只消除刚体运动：底面 z 固定；角点 (0,0,0) 再固定 x,y；相邻角点 (h,0,0) 再固定 y。
顶面节点 z 向位移 δ = ε·L（先移动再用投影约束锁住 z），x、y 自由。
读数：总轴向力 F（顶面反力之和）、横向应变、整个位移场与解析解的最大偏差。
解析解（线弹性小变形）：u = (-νεx, -νεy, εz)，F = E·L²·ε。

运行：python3 t0a1_uniaxial.py  → results/step0/t0a1_uniaxial.json
"""
import sys
sys.dont_write_bytecode = True
import os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
import Sofa, Sofa.Simulation
import common

L, E = 0.08, 5000.0


def run(nu, n, method, eps, linear="cg", ff_type="FTC"):
    root = common.make_root()
    common.add_static_solver(root, linear=linear)
    t, mo, ff = common.add_tissue(root, L=L, n=n, E=E, nu=nu, method=method, ff_type=ff_type)
    X = common.grid_nodes(L, n)
    bot = np.where(np.isclose(X[:, 2], 0))[0]
    top = np.where(np.isclose(X[:, 2], L))[0]
    cons = [t.addObject("PartialFixedProjectiveConstraint", name=nm, indices=idx, fixedDirections=d)
            for nm, idx, d in [("bottomZ", bot.tolist(), [0, 0, 1]), ("topZ", top.tolist(), [0, 0, 1]),
                               ("cornerXY", [0], [1, 1, 0]), ("cornerY", [1], [0, 1, 0])]]
    Sofa.Simulation.init(root)
    assert np.abs(mo.position.array() - X).max() < 1e-12
    p = X.copy(); p[top, 2] += eps * L
    mo.position.value = p
    Sofa.Simulation.animate(root, root.dt.value)
    status = str(root.newton.status.value)
    x = mo.position.array().copy()
    u = x - X
    ua = np.stack([-nu * eps * X[:, 0], -nu * eps * X[:, 1], eps * X[:, 2]], 1)
    f = common.internal_force(root, t, cons)
    F = -f[top, 2].sum()                     # 施加在顶面上的总力（拉为正）
    Fref = E * L * L * eps
    lat = ((x[:, 0].max() - x[:, 0].min()) / L - 1) / eps   # 横向应变 / 轴向应变
    interior = np.setdiff1d(np.arange(len(X)), np.r_[top, bot])
    return dict(ff=ff_type, nu=nu, n=n, h_mm=1000 * L / (n - 1), method=method, eps=eps, linear=linear, status=status,
                F=F, F_ref=Fref, F_err=F / Fref - 1, lat_ratio=lat, lat_err=lat / (-nu) - 1,
                field_err=float(np.abs(u - ua).max() / (abs(eps) * L)),
                interior_resid=float(np.abs(f[interior]).max()), bottom_F=float(f[bot, 2].sum()))


def show(r, tag=""):
    print(f"{tag}{r['ff']:6s} {r['method']:6s} {r['linear']:3s} nu={r['nu']:.2f} h={r['h_mm']:4.1f}mm eps={r['eps']:+.0e} "
          f"F={r['F']:+.5f}N ({r['F_err']:+.2e}) lat={r['lat_ratio']:+.5f} ({r['lat_err']:+.2e}) "
          f"field={r['field_err']:.2e} {r['status']}", flush=True)


if __name__ == "__main__":
    out = []
    # 主力场：FastTetrahedralCorotationalForceFieldFixed + polar，直接求解器和 CG 都跑
    for lin in ["ldl", "cg"]:
        for nu in [0.3, 0.49]:
            for n in [5, 9, 17]:
                out.append(run(nu, n, "polar", 0.01, linear=lin, ff_type="FTCfix")); show(out[-1])
            for e in [1e-4, 1e-3, -0.01, 0.05]:
                out.append(run(nu, 9, "polar", e, linear=lin, ff_type="FTCfix")); show(out[-1])
    # 对照：旋转提取方法（qr = 原默认值，none = 线性），修好的力场 + LDL
    for m in ["qr", "none"]:
        for nu in [0.3, 0.49]:
            for e in [1e-4, 1e-3, 0.01, -0.01, 0.05]:
                out.append(run(nu, 9, m, e, linear="ldl", ff_type="FTCfix")); show(out[-1])
    # 对照：原版 FTC（v25.12）+ LDL，暴露 buildStiffnessMatrix 问题
    for nu in [0.25, 0.3, 0.45, 0.49]:
        for m in ["polar", "none"]:
            out.append(run(nu, 9, m, 0.01, linear="ldl", ff_type="FTC")); show(out[-1], "[orig] ")
    # 对照：TetrahedronFEMForceField
    for m in ["large", "polar"]:
        for nu in [0.3, 0.49]:
            for e in [0.01, 0.05]:
                out.append(run(nu, 9, m, e, linear="ldl", ff_type="TFEM")); show(out[-1])
    res = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../results/step0/t0a1_uniaxial.json")
    json.dump(out, open(res, "w"), indent=1)
