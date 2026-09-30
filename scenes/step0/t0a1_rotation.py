"""0a-1 整体转动测试。见 docs/plan/step0_foundation.md §2.2。

(a) 无载荷的立方体整体刚性转动 θ：内力应 ≈ 0（共旋）；method=none/small 会产生虚假内力。
(b) 单轴拉伸 1% 的平衡构形整体转动 θ：内力应随之转动，f' = R f（客观性）。
读数：(a) Σ|f_i| 与 0.32 N（1% 单轴拉伸的总力）之比；(b) ‖f' − R f‖ / ‖f‖。
运行：python3 t0a1_rotation.py → results/step0/t0a1_rotation.json
"""
import sys
sys.dont_write_bytecode = True
import os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
import Sofa, Sofa.Simulation
import common

L, E, EPS = 0.08, 5000.0, 0.01
F1 = E * L * L * EPS   # 0.32 N


def rotmat(axis, deg):
    a = np.asarray(axis, float); a /= np.linalg.norm(a)
    th = np.radians(deg); K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * K @ K


def force_at(x, nu, n, ff_type, method):
    """在给定节点位置 x 下（参考构形为规则网格）计算内力。"""
    root = common.make_root()
    common.add_static_solver(root)
    t, mo, _ = common.add_tissue(root, L=L, n=n, E=E, nu=nu, method=method, ff_type=ff_type)
    Sofa.Simulation.init(root)
    mo.position.value = x
    return common.internal_force(root, t, [])


def stretched_state(nu, n):
    """1% 单轴拉伸的均匀场（解析）：u = (-νεx, -νεy, εz)。"""
    X = common.grid_nodes(L, n)
    return X * np.array([1 - nu * EPS, 1 - nu * EPS, 1 + EPS])


if __name__ == "__main__":
    out = []
    n = 9
    X = common.grid_nodes(L, n)
    c = X.mean(0)
    rots = [((0, 0, 1), 30), ((1, 1, 1), 30), ((1, 0, 0), 90), ((1, 2, 3), 90)]
    models = [("FTCfix", "polar"), ("FTCfix", "qr"), ("FTCfix", "none"), ("FTC", "qr"), ("FTC", "polar"), ("FTC", "polar2"), ("FTC", "none"),
              ("TFEM", "large"), ("TFEM", "polar"), ("TFEM", "svd"), ("TFEM", "small")]
    for nu in [0.3, 0.49]:
        Xs = stretched_state(nu, n)
        for ff_type, m in models:
            f0 = force_at(Xs, nu, n, ff_type, m)      # 未转动的拉伸构形
            for axis, deg in rots:
                R = rotmat(axis, deg)
                fa = force_at((X - c) @ R.T + c, nu, n, ff_type, m)
                fb = force_at((Xs - c) @ R.T + c, nu, n, ff_type, m)
                ra = np.abs(fa).sum() / F1
                rb = np.linalg.norm(fb - f0 @ R.T) / np.linalg.norm(f0)
                out.append(dict(nu=nu, ff=ff_type, method=m, axis=axis, deg=deg, rigid_sum_over_F1=ra,
                                stretched_objectivity_err=rb))
                print(f"nu={nu} {ff_type:6s} {m:6s} axis={axis} {deg:3d}deg: "
                      f"(a) sum|f|/0.32N={ra:.2e}  (b) |f'-Rf|/|f|={rb:.2e}", flush=True)
    res = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../results/step0/t0a1_rotation.json")
    json.dump(out, open(res, "w"), indent=1)
