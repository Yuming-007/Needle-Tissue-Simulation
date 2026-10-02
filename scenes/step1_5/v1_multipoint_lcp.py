"""第 1.5 步 V1：多点圆盘接触 vs 同一离散接触雅可比的 LCP 参照（第 1 步 P2 的多点推广）。见 step1_5_tip_survey.md §7。

动态接触：同 V0 / V5 的测试台（a_eff = 10 mm 圆盘，点距约等于 h，R1，dt = 0.01），针从表面上方 0.1 mm 以 0.2 mm/s
  压到 0.05 mm（准静态、线性范围）。在压深 0.02、0.05 mm 记录逐行数据（v5_inactive_points.rows_info）：
  针尖点、方向 n、投影节点和重心权重 w、λ、间隙。
参照（静力、线性、同一雅可比）：
  - 组织节点柔度：对所有投影节点 S 的每个分量施加 ±f 的静力（对称差分消去二阶项），得 3|S|×3|S| 的 C；
  - 接触柔度 W_rs = (n_r⊗w_r)ᵀ C (n_s⊗w_s)，n、w 用动态接触在该压深实际使用的值；
  - 未变形组织下的间隙 q_r = n_r·(x_tip,r − Σ w X0)，线性化间隙 g = q + Wλ（λ 为力，N）；
  - LCP：λ ≥ 0，g ≥ 0，λ·g = 0。W 对称正定 ⇒ 等价于 min ½λᵀWλ + qᵀλ（λ ≥ 0），用 Cholesky + NNLS 精确求解。
比较：受力点集合、合力 F_tip 和合力矩 M_tip、逐点 λ、不受力点的间隙。
反面对照："全部点接触"参照 λ = −W⁻¹q（不要求 λ ≥ 0），看会出现多少拉力、合力差多少。
判据：受力点集合一致；合力相对差 ≤ 1%；作用力 = 反作用力。
运行：python3 v1_multipoint_lcp.py → results/step1_5/v1_multipoint_lcp.json
"""
import sys
sys.dont_write_bytecode = True
import os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../step1"))
import numpy as np
from scipy.optimize import nnls
import Sofa, Sofa.Simulation
import common, contact
import t1_vertex as V
from v5_inactive_points import rows_info

contact.load_plugins()
L, LEN, dt, A_EFF = V.L, V.LEN, 0.01, 0.010
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../results/step1_5")
PROBE_MM = [0.02, 0.05]


def dynamic(n, gap=0.0001, v=0.0002, depth=0.00005):
    h = L / (n - 1)
    pts = contact.tip_patch_points(LEN, A_EFF, h, "disk")
    zf, T = V.piecewise(L + gap, [(gap / v, -v), (depth / v, -v)])
    root = common.make_root(dt=dt, loop="free")
    traj = lambda t: (np.array([L / 2, L / 2, zf(t) + LEN]), V.Q)
    nd, base, body, tip, local = common.add_kinematic_needle(root, length=LEN, n_body=V.NB, pose0=traj(0.0))
    pmo, pg = contact.add_tip_patch(nd, pts)
    root.addObject(common.NeedleDriver(base, traj, dt, children=[(body, local), (tip, local[-1:]), (pmo, pts)],
                                       name="driver"))
    root.addObject("CollisionLoop")
    tis = contact.add_tissue_with_surface(root, n=n, L=L, E=V.E, nu=V.NU, linear="ldl")
    _, bg = contact.add_needle_geometry(nd, body, tip, V.NB)
    contact.add_tip_contact(root, pg, bg, tis, distance=0.02)
    Sofa.Simulation.init(root)
    cs = root.csolver
    probes, todo = [], list(PROBE_MM)
    for k in range(int(round(T / dt))):
        Sofa.Simulation.animate(root, dt)
        z_tip = zf((k + 1) * dt)
        d_mm = (L - z_tip) * 1e3
        if todo and d_mm >= todo[0] - 1e-9:
            todo.pop(0)
            lam = np.array(cs.constraintForces.value)
            rows, _ = rows_info(lam, pmo, tis)
            p_ref = np.array([L / 2, L / 2, z_tip])
            F, M, _ = common.needle_wrench([pmo], lam, dt, p_ref)
            # 组织一侧受到的约束力（只读表面 MO，避免和体 MO 重复计数，见第 1 步）
            Js = tis["surf_dofs"].constraint.value
            lp = np.zeros(Js.shape[0]); m = min(len(lam), Js.shape[0]); lp[:m] = lam[:m]
            F_tis = (Js.T @ lp).reshape(-1, 3).sum(0) / dt
            probes.append(dict(depth_mm=d_mm, rows=rows, F=F.tolist(), M_tip=M.tolist(), F_tissue=F_tis.tolist(),
                               p_ref=p_ref.tolist(), gs_iter=int(cs.currentIterations.value),
                               gs_err=float(cs.currentError.value)))
    Sofa.Simulation.unload(root)
    return probes


def node_compliance(n, S, f=1e-3):
    """S 中各节点 3 个分量的静力柔度（m/N），对称差分 (u(+f) − u(−f)) / 2f。"""
    X = common.grid_nodes(L, n)
    bot = np.where(np.isclose(X[:, 2], 0))[0]
    root = common.make_root(); common.add_static_solver(root, linear="ldl", abs_tol=1e-14)
    t, mo, _ = common.add_tissue(root, L=L, n=n, E=V.E, nu=V.NU, swapping=True)
    t.addObject("FixedProjectiveConstraint", indices=bot.tolist())
    cff = t.addObject("ConstantForceField", indices=list(map(int, S)), forces=np.zeros((len(S), 3)).tolist())
    Sofa.Simulation.init(root)
    C = np.zeros((3 * len(S), 3 * len(S)))
    for j in range(len(S)):
        for c in range(3):
            u = []
            for sgn in (1.0, -1.0):
                fr = np.zeros((len(S), 3)); fr[j, c] = sgn * f
                cff.forces.value = fr.tolist()
                mo.position.value = X.tolist(); mo.velocity.value = np.zeros_like(X).tolist()
                Sofa.Simulation.animate(root, 0.01)
                u.append(mo.position.array()[S].copy() - X[S])
            C[:, 3 * j + c] = ((u[0] - u[1]) / (2 * f)).ravel()
    Sofa.Simulation.unload(root)
    return C


def lcp_reference(rows, S, C, X0):
    idx = {s: i for i, s in enumerate(S)}
    R = len(rows)
    B = np.zeros((3 * len(S), R))          # 第 r 列 = n_r ⊗ w_r（在 S 的 3|S| 空间）
    q = np.zeros(R)
    for r, x in enumerate(rows):
        nvec = np.array(x["n"])
        xp = sum(w * X0[j] for j, w in zip(x["nodes"], x["w"]))
        q[r] = nvec @ (np.array(x["xp"]) - xp)
        for j, w in zip(x["nodes"], x["w"]):
            B[3 * idx[j]:3 * idx[j] + 3, r] += w * nvec
    W = B.T @ C @ B
    Wsym = 0.5 * (W + W.T)
    Rc = np.linalg.cholesky(Wsym).T        # W = Rcᵀ Rc
    lam, _ = nnls(Rc, -np.linalg.solve(Rc.T, q))
    g = q + Wsym @ lam
    lam_all = -np.linalg.solve(Wsym, q)
    return dict(lam=lam, g=g, lam_all=lam_all, q=q, W_asym=float(np.abs(W - W.T).max() / np.abs(W).max()),
                W_cond=float(np.linalg.cond(Wsym)))


def run(n):
    X0 = common.grid_nodes(L, n)
    a = time.perf_counter()
    probes = dynamic(n)
    t_dyn = time.perf_counter() - a
    S = sorted({j for p in probes for x in p["rows"] for j in x["nodes"]})
    a = time.perf_counter()
    C = node_compliance(n, np.array(S))
    t_c = time.perf_counter() - a
    res = []
    for p in probes:
        rows = p["rows"]
        ref = lcp_reference(rows, S, C, X0)
        lam_dyn = np.array([x["lam"] for x in rows]) / dt
        g_dyn = np.array([x["gap"] for x in rows])
        nv = np.array([x["n"] for x in rows]); xp = np.array([x["xp"] for x in rows]); pr = np.array(p["p_ref"])
        # 作用在针上的力 = λ n（行 = n ⊗ e_point）
        F_ref = (ref["lam"][:, None] * nv).sum(0)
        M_ref = np.cross(xp - pr, ref["lam"][:, None] * nv).sum(0)
        F_all = (ref["lam_all"][:, None] * nv).sum(0)
        act_dyn = lam_dyn > 1e-12 * max(lam_dyn.max(), 1e-30)
        act_ref = ref["lam"] > 1e-12 * max(ref["lam"].max(), 1e-30)
        F_dyn = np.array(p["F"])
        d = dict(depth_mm=p["depth_mm"], n_rows=len(rows), points=[x["point"] for x in rows],
                 active_dyn=act_dyn.tolist(), active_ref=act_ref.tolist(), same_active_set=bool((act_dyn == act_ref).all()),
                 lam_dyn=lam_dyn.tolist(), lam_ref=ref["lam"].tolist(), lam_all=ref["lam_all"].tolist(),
                 gap_dyn=g_dyn.tolist(), gap_ref=ref["g"].tolist(),
                 F_dyn=F_dyn.tolist(), F_ref=F_ref.tolist(), F_all=F_all.tolist(),
                 M_dyn=p["M_tip"], M_ref=M_ref.tolist(), F_tissue=p["F_tissue"],
                 rel_F=float(F_dyn[2] / F_ref[2] - 1), rel_F_all=float(F_all[2] / F_ref[2] - 1),
                 lam_rel_max=float(np.abs(lam_dyn - ref["lam"])[act_ref].max() / ref["lam"].max()),
                 n_tension_all=int((ref["lam_all"] < 0).sum()),
                 action_reaction=float(np.linalg.norm(F_dyn + np.array(p["F_tissue"])) / np.linalg.norm(F_dyn)),
                 W_asym=ref["W_asym"], W_cond=ref["W_cond"], gs_iter=p["gs_iter"], gs_err=p["gs_err"])
        res.append(d)
    return dict(n=n, h_mm=L / (n - 1) * 1e3, S=S, t_dynamic=t_dyn, t_compliance=t_c, probes=res)


if __name__ == "__main__":
    out = []
    for n in [9, 17]:
        r = run(n)
        out.append(r)
        print(f"== h = {r['h_mm']:.1f} mm 圆盘：投影节点 {len(r['S'])} 个；动态 {r['t_dynamic']:.0f} s，柔度矩阵 {r['t_compliance']:.0f} s", flush=True)
        for d in r["probes"]:
            ina_dyn = [pt for pt, a in zip(d["points"], d["active_dyn"]) if not a]
            ina_ref = [pt for pt, a in zip(d["points"], d["active_ref"]) if not a]
            gi = [i for i, a in enumerate(d["active_ref"]) if not a]
            gd = ", ".join(f"{d['gap_dyn'][i]*1e3:.4e}/{d['gap_ref'][i]*1e3:.4e}" for i in gi)
            print(f"  压深 {d['depth_mm']:.3f} mm：不受力点 动态 {ina_dyn} / LCP {ina_ref}，集合一致 {d['same_active_set']}；"
                  f"F_z 动态 {d['F_dyn'][2]:.6e} N，LCP {d['F_ref'][2]:.6e} N，相对差 {d['rel_F']:+.2e}；"
                  f"逐点 λ 最大差 / λmax {d['lam_rel_max']:.2e}", flush=True)
            print(f"     不受力点间隙（mm，动态/LCP）：{gd}；|M| 动态 {np.linalg.norm(d['M_dyn']):.1e}，LCP {np.linalg.norm(d['M_ref']):.1e}；"
                  f"作用 = 反作用 {d['action_reaction']:.1e}；W 不对称 {d['W_asym']:.1e}，条件数 {d['W_cond']:.1e}；GS {d['gs_iter']}", flush=True)
            print(f"     反面对照（全部点接触）：拉力点 {d['n_tension_all']} 个，F_z {d['F_all'][2]:.6e} N，相对 LCP {d['rel_F_all']:+.2e}", flush=True)
        json.dump(out, open(os.path.join(OUT, "v1_multipoint_lcp.json"), "w"))
