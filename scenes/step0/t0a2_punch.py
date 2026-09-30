"""0a-2 局部组织测试：无摩擦刚性平底圆柱压头，和 Sneddon 解比较。见 docs/plan/step0_foundation.md §3。

场景：边长 L 的立方体，底面全固定，侧面自由；顶面中心 r ≤ a 的节点给定 z 向位移 -δ，
x、y 自由（无摩擦；bonded=True 时三个方向都固定，作为对照）。静力求解（LDL，修复后的 FTC，polar）。
读数：压头节点总反力 F，k = F/δ；与 k_S = 2 a E / (1-ν²) 比较：
  - 用名义半径 a；
  - 用等效半径 a_eff = sqrt(A/π)，A = 被压节点的从属面积（每个节点 h²）。
运行：python3 t0a2_punch.py → results/step0/t0a2_punch.json
"""
import sys
sys.dont_write_bytecode = True
import os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
import Sofa, Sofa.Simulation
import common

E, A_PUNCH = 5000.0, 0.010


def run(L, h, nu, bonded=False, delta_over_a=0.01, linear="ldl", quarter=False):
    """quarter=True：只算 x ≥ L/2、y ≥ L/2 的四分之一，对称面上固定法向位移；力和面积 × 4。"""
    N = int(round(L / h))
    delta = delta_over_a * A_PUNCH
    if quarter:
        size, nn, org = (L / 2, L / 2, L), (N // 2 + 1, N // 2 + 1, N + 1), (L / 2, L / 2, 0.0)
    else:
        size, nn, org = L, N + 1, (0.0, 0.0, 0.0)
    root = common.make_root()
    common.add_static_solver(root, linear=linear, abs_tol=1e-9)   # F ≈ 0.01 N → 相对 1e-7
    t, mo, _ = common.add_tissue(root, L=size, n=nn, E=E, nu=nu, origin=org)
    X = common.grid_nodes(size, nn, org)
    bot = np.where(np.isclose(X[:, 2], 0))[0]
    top = np.isclose(X[:, 2], L)
    r = np.hypot(X[:, 0] - L / 2, X[:, 1] - L / 2)
    punch = np.where(top & (r <= A_PUNCH + 1e-9))[0]
    cons = [t.addObject("FixedProjectiveConstraint", name="bottom", indices=bot.tolist())]
    if bonded:
        cons.append(t.addObject("FixedProjectiveConstraint", name="punch", indices=punch.tolist()))
    else:
        cons.append(t.addObject("PartialFixedProjectiveConstraint", name="punch", indices=punch.tolist(),
                                fixedDirections=[0, 0, 1]))
    # 从属面积权重（四分之一模型中，对称面上的节点只有一半或四分之一属于本模型）
    w = np.ones(len(X))
    if quarter:
        symx = np.setdiff1d(np.where(np.isclose(X[:, 0], L / 2))[0], bot)
        symy = np.setdiff1d(np.where(np.isclose(X[:, 1], L / 2))[0], bot)
        cons.append(t.addObject("PartialFixedProjectiveConstraint", name="symX", indices=symx.tolist(),
                                fixedDirections=[1, 0, 0]))
        cons.append(t.addObject("PartialFixedProjectiveConstraint", name="symY", indices=symy.tolist(),
                                fixedDirections=[0, 1, 0]))
        w[np.isclose(X[:, 0], L / 2)] *= 0.5
        w[np.isclose(X[:, 1], L / 2)] *= 0.5
    mult = 4.0 if quarter else 1.0
    Sofa.Simulation.init(root)
    assert np.abs(mo.position.array() - X).max() < 1e-12
    p = X.copy(); p[punch, 2] -= delta
    mo.position.value = p
    t0 = time.perf_counter()
    Sofa.Simulation.animate(root, root.dt.value)
    t_solve = time.perf_counter() - t0
    status = str(root.newton.status.value)
    iters = len(root.newton.residualGraph.getValueString().split()) - 2
    x = mo.position.array().copy()
    f = common.internal_force(root, t, cons)
    F = mult * f[punch, 2].sum()             # 组织对压头的反力（向上为正），全模型等效
    F_bot = mult * f[bot, 2].sum()
    free = np.setdiff1d(np.arange(len(X)), np.r_[punch, bot])
    if quarter:   # 对称面上的节点受对称约束，只检查非约束方向；这里只检查完全自由的内部节点
        free = np.setdiff1d(free, np.r_[symx, symy])
    a_eff = np.sqrt(mult * w[punch].sum() * h * h / np.pi)
    kS = lambda a: 2 * a * E / (1 - nu ** 2)
    k = F / delta
    return dict(L=L, h_mm=h * 1000, n=[int(v) for v in np.atleast_1d(nn)], dofs=int(3 * len(X)), quarter=quarter, nu=nu,
                bonded=bonded, delta=delta, n_punch=int(len(punch)),
                a_eff_mm=a_eff * 1000, status=status, newton_iters=iters, F=F, F_bottom=F_bot, k=k,
                kS_nominal=kS(A_PUNCH), ratio_nominal=k / kS(A_PUNCH),
                kS_eff=kS(a_eff), ratio_eff=k / kS(a_eff),
                interior_resid=float(np.abs(f[free]).max()), t_solve=t_solve)


def show(d):
    print(f"L={d['L']*100:.0f}cm h={d['h_mm']:4.1f}mm {'Q' if d['quarter'] else 'F'} dofs={d['dofs']:6d} nu={d['nu']:.2f} "
          f"{'bonded' if d['bonded'] else 'frictionless':12s} "
          f"punch={d['n_punch']:3d} a_eff={d['a_eff_mm']:5.2f}mm k={d['k']:7.2f}N/m "
          f"k/kS(a)={d['ratio_nominal']:.3f} k/kS(a_eff)={d['ratio_eff']:.3f} "
          f"F+Fbot={d['F']+d['F_bottom']:+.1e} resid={d['interior_resid']:.1e} t={d['t_solve']:.1f}s it={d['newton_iters']} {d['status']}",
          flush=True)


if __name__ == "__main__":
    # 用法：python3 t0a2_punch.py coarse|quarter
    #   coarse： 全模型：8 cm h = 10/5 mm，16 cm h = 10 mm，粘结对照，线性性检查
    #   quarter：四分之一对称模型：先和全模型对照（h = 10/5 mm），再算
    #            细网格（8 cm h = 2.5 mm）和有限尺寸（h = 10 mm：8/16/24/32 cm；h = 5 mm：8/16 cm）
    # 为什么用四分之一模型：SparseLDL 分解时间约按自由度的 3.3 次方增长，全模型 h = 2.5 mm（10.8 万自由度）
    # 每次分解约 25 分钟（见 docs/plan/step0_results.md）。
    part = sys.argv[1] if len(sys.argv) > 1 else "coarse"
    res = os.path.join(os.path.dirname(os.path.abspath(__file__)), f"../../results/step0/t0a2_punch_{part}.json")
    out = []

    def rec(d):
        out.append(d); show(d)
        json.dump(out, open(res, "w"), indent=1)

    if part == "coarse":
        for nu in [0.3, 0.45, 0.49]:
            for L, h in [(0.08, 0.01), (0.08, 0.005), (0.16, 0.01)]:
                rec(run(L, h, nu))
        for nu in [0.3, 0.49]:
            rec(run(0.08, 0.005, nu, bonded=True))
        for da in [0.001, 0.05]:
            rec(dict(run(0.08, 0.005, 0.3, delta_over_a=da), delta_over_a=da))
    elif part == "quarter":
        for nu in [0.3, 0.49, 0.45]:
            for L, h in [(0.08, 0.01), (0.08, 0.005), (0.08, 0.0025),
                         (0.16, 0.01), (0.24, 0.01), (0.32, 0.01), (0.16, 0.005)]:
                rec(run(L, h, nu, quarter=True))
