"""第 3 步：插入导向（横向约束，不加摩擦）的验证运行。见 docs/plan/step3_results.md。

场景和第 2 步完全相同（s2_puncture.build：W50，FTCfix polar，E = 5 kPa，ν = 0.45，β = 0.1，阈值 0.22 N，耦合点间距 5 mm，
frictionCoeff = 0，刚性纯运动学针，竖直对准顶面中心，单点针尖），只换轨迹：
  deep   （3A）：5 mm/s 下插到 50 mm（W50 加密区末端），停 1 s，退回，停 0.5 s；
  lateral（3B）：5 mm/s 下插到 30 mm，停 0.5 s，整根针以 1 mm/s 沿 +x 平移 3 mm，停 2 s，以 1 mm/s 平移回原位，停 1 s，
                 退回，停 0.5 s。针是刚性的，刺入点也跟着平移。
每步记录：针上的约束合力（作用在针上）、组织一侧的合力（作用力 = 反作用力）、每个耦合点的深度 / 组织点 P 到针轴的横向距离 /
  该点作用在针上的力、几个观察节点的 x 位移（3B：15 mm 深处离针 5 / 10 / 20 / 35 mm，离针 20 mm 处深 40 / 70 mm）、GS 迭代次数和误差、最差单元体积比、每步计算时间。
底面反力（受力平衡检查）：每 REACT_EVERY 步用组织力场自己组装的矩阵估计内力。assembleKMatrix 返回 SOFA 约定的 df/dx = −K
  （当前构形的切线刚度，负定），所以 f = (df/dx)·(x − X0) + β·(df/dx)·v；底面固定节点上的内力之和取负 = 支座作用在组织上的力
  R_bottom。静止时应 = 针上的力 F（组织：针的推力 −F + 支座力 R = 0）。共旋单元用 K(x)(x − X0) 代替 R·K0·(Rᵀx − X0)，
  误差是转角的二阶小量（3B 中表现为随平移距离平方增长的假 R_z）。
运行：python3 s3_insertion.py <deep|lateral> [ldl|async] [平移幅度 mm] [耦合点间距 mm]   （默认 ldl = R1；async = R3，用于和 R1
  对照、测每步时间；幅度默认 3 mm；间距默认 5 mm，4 / 6 mm 是"单点力 vs 合力"的网格敏感性检查）
  → results/step3/s3_<mode>_<linear>[_amp<幅度>][_tip<间距>].npz + 终端摘要
"""
import sys
sys.dont_write_bytecode = True
import os, time
HERE = os.path.dirname(os.path.abspath(__file__))
for d in ["..", "../step1", "../step1_5", "../step2"]:
    sys.path.insert(0, os.path.join(HERE, d))
import numpy as np
import Sofa, Sofa.Simulation
import common
import t1_vertex as V
import s2_puncture as S
from v4_large_depth import tet_volumes

OUT = os.path.join(HERE, "../../results/step3")
SPEED, GAP = 0.005, 0.0005
LAT_DEPTH, LAT_AMP, LAT_SPEED = 0.030, 0.003, 0.001
PROBES = [(0.005, 0.015), (0.010, 0.015), (0.020, 0.015), (0.035, 0.015), (0.020, 0.040), (0.020, 0.070)]   # (离针 x 距离, 深度) m
NCP = 12                                                    # 记录时每步最多的耦合点数
REACT_EVERY = 5                                             # 每隔几步算一次底面反力（组装 K 较慢）


def trajectories(mode, amp=LAT_AMP):
    """返回 (zf, xf, 总时长)：针尖 z(t) 和针的横向平移 x(t)。amp：3B 的平移幅度（m）。"""
    if mode == "deep":
        d = 0.050
        zf, T = V.piecewise(S.L + GAP, [((GAP + d) / SPEED, -SPEED), (1.0, 0.0), ((GAP + d) / SPEED, SPEED), (0.5, 0.0)])
        return zf, (lambda t: 0.0), T
    d, a, v = LAT_DEPTH, amp, LAT_SPEED
    t_in = (GAP + d) / SPEED
    t_lat = [(t_in, 0.0), (0.5, 0.0), (a / v, v), (2.0, 0.0), (a / v, -v), (1.0, 0.0)]
    hold = sum(s[0] for s in t_lat)
    zf, T = V.piecewise(S.L + GAP, [(t_in, -SPEED), (hold - t_in, 0.0), (t_in, SPEED), (0.5, 0.0)])
    xf, _ = V.piecewise(0.0, t_lat)
    return zf, xf, T


def coupling_details(s, lam, dt):
    """每个耦合点：(深度 m，组织点 P 到针轴的横向距离 m，作用在针上的力 [3] N)。
    针身 MO 和组织 MO 的约束雅可比里，每个耦合点占连续 3 行（同一行号）；一行在某个 MO 上的各点系数 = 权重 × 方向，
    所以 P = Σ|系数块|·x / Σ|系数块|（组织 MO，重心坐标插值），该行的方向 = Σ 系数块（针身 MO，权重和为 1），
    该点作用在针上的力 = Σ_3 行 λ·方向 / dt。"""
    Jb = s["body"].constraint.value
    if Jb.shape[0] == 0 or Jb.nnz == 0:
        return []
    Jb, Jt = Jb.tocsr(), s["tis"]["dofs"].constraint.value.tocsr()
    xt = s["tis"]["dofs"].position.array()
    nb = 3 * len(s["body"].position.array())
    rows = sorted(set(Jb.nonzero()[0]))
    base = np.asarray(s["base"].position.value)[0][:3]          # 针竖直：针轴是过针根的竖直线
    out = []
    for i in range(0, len(rows) - 2, 3):
        r3 = rows[i:i + 3]
        f = np.zeros(3)
        for r in r3:
            blk = np.pad(Jb.getrow(r).toarray().ravel(), (0, nb - Jb.shape[1]))[:nb].reshape(-1, 3)
            f += lam[r] * blk.sum(0) / dt
        row = np.pad(Jt.getrow(r3[0]).toarray().ravel(), (0, 3 * len(xt) - Jt.shape[1]))
        w = np.linalg.norm(row.reshape(-1, 3), axis=1)
        P = w @ xt / w.sum()
        out.append((S.L - P[2], np.hypot(P[0] - base[0], P[1] - base[1]), f))
    return sorted(out, key=lambda o: o[0])


def run(mode, linear="ldl", dt=0.01, amp=LAT_AMP, tip_dist=S.TIP_DIST):
    zf, xf, T = trajectories(mode, amp)
    root = Sofa.Core.Node("root")
    s = S.build(root, dt=dt, linear=linear, zf=zf, xf=xf, tip_dist=tip_dist)
    Sofa.Simulation.init(root)
    tis = s["tis"]
    X0 = np.asarray(tis["X"])
    tets = np.array(tis["node"].topo.tetrahedra.value)
    V0 = tet_volumes(X0, tets)
    probes = [int(np.argmin(np.linalg.norm(X0 - [S.L / 2 + r, S.L / 2, S.L - dp], axis=1))) for r, dp in PROBES]
    rec = {k: [] for k in ["t", "depth", "dx", "F", "F_tissue", "n_cp", "cp_depth", "cp_gap", "cp_f",
                           "probe_dx", "iters", "err", "vol_min", "step_ms", "R_bottom"]}
    bot = np.asarray(tis["bottom"])
    for k in range(int(round(T / dt))):
        a = time.perf_counter()
        Sofa.Simulation.animate(root, dt)
        ms = 1e3 * (time.perf_counter() - a)
        t = (k + 1) * dt
        lam = np.array(root.csolver.constraintForces.value)
        p_base = np.asarray(s["base"].position.value)[0][:3]
        F, _, _ = common.needle_wrench([s["pmo"], s["body"]], lam, dt, p_base)
        Ft, _, _ = common.needle_wrench([tis["dofs"]], lam, dt, p_base)
        cps = coupling_details(s, lam, dt) if len(lam) else []
        pad = lambda v, w=1: np.pad(np.array(v, float).reshape(-1, w), ((0, NCP - len(v)), (0, 0)),
                                    constant_values=np.nan)[:NCP].squeeze()
        x = tis["dofs"].position.array()
        rec["t"].append(t); rec["depth"].append((S.L - zf(t)) * 1e3); rec["dx"].append(xf(t) * 1e3)
        rec["F"].append(F.copy()); rec["F_tissue"].append(Ft.copy()); rec["n_cp"].append(len(cps))
        rec["cp_depth"].append(pad([c[0] * 1e3 for c in cps]))
        rec["cp_gap"].append(pad([c[1] * 1e3 for c in cps]))
        rec["cp_f"].append(pad([c[2] for c in cps], 3) if cps else np.full((NCP, 3), np.nan))
        rec["probe_dx"].append((x[probes, 0] - X0[probes, 0]) * 1e3)
        rec["iters"].append(root.csolver.currentIterations.value); rec["err"].append(root.csolver.currentError.value)
        rec["vol_min"].append((tet_volumes(x, tets) / V0).min()); rec["step_ms"].append(ms)
        R = np.full(3, np.nan)
        if (k + 1) % REACT_EVERY == 0:
            dfdx = tis["fem"].assembleKMatrix()                 # df/dx = −K（当前构形，未投影固定约束）
            v = tis["dofs"].velocity.array()
            f = dfdx @ (x - X0).ravel() + S.BETA * (dfdx @ v.ravel())
            R = -f.reshape(-1, 3)[bot].sum(0)
        rec["R_bottom"].append(R)
        if not np.isfinite(F).all() or rec["vol_min"][-1] < -1:
            print(f"!! 第 {k+1} 步发散或单元严重翻转，停止", flush=True)
            break
    Sofa.Simulation.unload(root)
    rec = {k: np.array(v) for k, v in rec.items()}
    rec["probe_r"] = np.array([p[0] for p in PROBES]) * 1e3
    rec["probe_depth"] = np.array([p[1] for p in PROBES]) * 1e3
    rec["probe_X"] = X0[probes] * 1e3
    os.makedirs(OUT, exist_ok=True)
    tag = "" if amp == LAT_AMP or mode != "lateral" else f"_amp{amp*1e3:g}"
    tag += "" if np.isclose(tip_dist, S.TIP_DIST) else f"_tip{tip_dist*1e3:g}"
    np.savez(os.path.join(OUT, f"s3_{mode}_{linear}{tag}.npz"), **rec)
    return rec


def summary(mode, rec, linear):
    d, F, n = rec["depth"], rec["F"], rec["n_cp"]
    kp = int(np.argmax(n > 0))
    lat = np.hypot(F[:, 0], F[:, 1])
    post = slice(kp, None)
    print(f"== {mode}，{linear}：刺穿 第 {kp+1} 步，压深 {d[kp]:.2f} mm", flush=True)
    if mode == "deep":
        for dd in [10, 20, 30, 40, 50]:
            i = int(np.argmin(np.abs(d[:np.argmax(d) + 1] - dd)))
            cz = rec["cp_depth"][i][~np.isnan(rec["cp_depth"][i])]
            print(f"   下插 {d[i]:.1f} mm：耦合点 {n[i]} 个 {np.round(cz, 1).tolist()} mm；轴向力 {F[i, 2]:.1e} N，"
                  f"侧向 {lat[i]*1e3:.1f} mN；P 到针轴最大 {np.nanmax(rec['cp_gap'][i]):.1e} mm", flush=True)
        back = np.where(np.arange(len(n)) > np.argmax(d))[0]
        gone = back[n[back] == 0]
        print(f"   后退：耦合点在压深 {d[gone[0]]:.2f} mm 全部删除" if len(gone) else "   后退：耦合点没有删完", flush=True)
    else:
        dx = rec["dx"]
        i3 = int(np.argmax(dx))
        hold = np.where(np.isclose(dx, dx.max()))[0]
        i_end = hold[-1]                                        # 平移到位、保持 2 s 的最后一步
        back0 = np.where((np.arange(len(dx)) > i_end) & np.isclose(dx, 0) & (n > 0))[0]
        print(f"   平移到 {dx[i3]:.2f} mm 的第一步：针上的力 F = {np.round(F[i3]*1e3, 2).tolist()} mN；保持结束："
              f"{np.round(F[i_end]*1e3, 2).tolist()} mN（耦合点 {n[i_end]} 个）", flush=True)
        for frac in [1 / 3, 2 / 3]:
            j = int(np.argmin(np.abs(dx[:i3 + 1] - frac * dx.max())))
            print(f"   平移 {dx[j]:.2f} mm：F_x = {F[j, 0]*1e3:.2f} mN（F_x / dx = {F[j, 0]*1e3/dx[j]:.2f} mN/mm）", flush=True)
        print(f"   保持结束 F_x / dx = {F[i_end, 0]*1e3/dx[i_end]:.2f} mN/mm", flush=True)
        if len(back0):
            j = back0[0]
            print(f"   回到原位第一步：F = {np.round(F[j]*1e3, 2).tolist()} mN；回位后 1 s 内 |F| 最大 "
                  f"{np.abs(F[back0]).max()*1e3:.2f} mN", flush=True)
        cf, cd = rec["cp_f"][i_end], rec["cp_depth"][i_end]
        m = ~np.isnan(cd)
        print(f"   保持结束各耦合点（深度 mm：F_x mN）：" +
              "，".join(f"{a:.1f}: {b*1e3:.2f}" for a, b in zip(cd[m], cf[m, 0])), flush=True)
        print(f"   Σ 各点 F_x = {np.nansum(cf[m, 0])*1e3:.2f} mN（合力 F_x {F[i_end, 0]*1e3:.2f} mN）", flush=True)
        Rb = rec["R_bottom"]
        j = max(i for i in range(i_end + 1) if np.isfinite(Rb[i, 0]))
        print(f"   底面反力（支座作用在组织上）保持结束（第 {j+1} 步）：R = {np.round(Rb[j]*1e3, 2).tolist()} mN；"
              f"针上的力 F = {np.round(F[j]*1e3, 2).tolist()} mN", flush=True)
        print(f"   保持结束：P 到针轴最大 {np.nanmax(rec['cp_gap'][i_end]):.1e} mm；节点 x 位移（离针 mm, 深度 mm → mm）：" +
              "，".join(f"({r:.0f}, {dp:.0f}) → {u:.3f}" for r, dp, u in
                       zip(rec["probe_r"], rec["probe_depth"], rec["probe_dx"][i_end])), flush=True)
    print(f"   刺穿后：|轴向力| 最大 {np.abs(F[post, 2]).max():.1e} N；P 到针轴最大 {np.nanmax(rec['cp_gap'][post]):.1e} mm", flush=True)
    res = np.linalg.norm(F + rec["F_tissue"], axis=1)
    print(f"   作用力 + 反作用力：最大 {res.max():.1e} N；GS 迭代 中位数 {np.median(rec['iters']):.0f} / 最大 "
          f"{rec['iters'].max()}，最大误差 {rec['err'].max():.1e}；最差单元体积比 {rec['vol_min'].min():.2f}", flush=True)
    ms = rec["step_ms"][20:]
    print(f"   每步计算（含记录前的 animate）：中位数 {np.median(ms):.1f} / p95 {np.percentile(ms, 95):.1f} / 最大 {ms.max():.1f} ms",
          flush=True)


if __name__ == "__main__":
    mode = sys.argv[1]
    linear = sys.argv[2] if len(sys.argv) > 2 else "ldl"
    amp = float(sys.argv[3]) * 1e-3 if len(sys.argv) > 3 else LAT_AMP
    tip = float(sys.argv[4]) * 1e-3 if len(sys.argv) > 4 else S.TIP_DIST
    rec = run(mode, linear, amp=amp, tip_dist=tip)
    summary(mode, rec, linear)
