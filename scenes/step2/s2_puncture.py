"""第 2 步：插件原生刺穿机制的验证运行（W50，单点针尖，刚性纯运动学针）。见 docs/plan/step2_results.md。

机制（InsertionAlgorithm，docs/code/01）：刺穿前针尖–表面单边约束（tenting）；上一步针尖 Σ‖λ‖/dt 超过阈值 →
  针尖下的表面点成为第一个耦合点，撤掉针尖接触；之后针尖离最后一个耦合点 ≥ tipDistThreshold 时沿针身每 5 mm 加一个
  耦合点（只接受落在组织四面体内的候选点），每步把耦合点重投影到针身上，ConstraintInsertion 给 2 行横向约束 +
  1 行轴向（frictionCoeff = 0 → 轴向无力）；针后退时删掉针尖前方的耦合点，删光则回到刺穿前状态。
场景：W50（1242 节点），FTCfix polar，E = 5 kPa，ν = 0.45（验证用参数），底面固定，无重力；针竖直对准顶面中心。
  组织黏性：Rayleigh 刚度阻尼 β = 0.1 s（C = βK，α = 0），取 Bui 2019 的文献值，未标定（为什么加、影响多大：
  docs/plan/step2_results.md；不加阻尼时实时求解器在刺穿回弹时发散）。
  轨迹：表面上方 0.5 mm 出发，5 mm/s 下插到 30 mm 深，停 1 s，再以 5 mm/s 退回出发点，停 0.5 s。
  刺穿阈值 0.22 N（只用于验证，让刺穿发生在约 5–6 mm 压深；W50 上 5 mm 时约 0.20 N、6 mm 时约 0.23 N）。
  离线运行用 SparseLDLSolver（确定性）；GUI 场景（s2_gui_scene.py）用 AsyncSparseLDLSolver。
每步记录：针的约束合力（common.needle_wrench，作用在针上）、其中针尖 / 针身各自的部分、组织一侧的合力（作用力 = 反作用力检查）、
  耦合点数和位置（由针身约束的雅可比反推，= 耦合点在针轴上的投影）、顶面中心附近节点的 z（回弹）、GS 迭代次数、最差单元体积比。
运行：python3 s2_puncture.py [dt] [β]   （默认 dt = 0.01、β = 0.1；dt = 0.005 只用于比较刺穿时刻；β = 0 是不加阻尼的对照）
  → results/step2/s2_puncture_dt{dt}.npz（β ≠ 0.1 时为 s2_puncture_dt{dt}_beta{β}.npz）+ 终端摘要
"""
import sys
sys.dont_write_bytecode = True
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../step1"))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../step1_5"))
import numpy as np
import Sofa, Sofa.Simulation
import common, contact
import t1_vertex as V
from v4_large_depth import tet_volumes

contact.load_plugins()
L, LEN, NB, Q = V.L, V.LEN, V.NB, V.Q
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../results/step2")
THRESHOLD, TIP_DIST, SPEED, DEPTH, GAP = 0.22, 0.005, 0.005, 0.030, 0.0005
BETA = 0.1                                                   # 组织 Rayleigh 刚度阻尼 β（s），Bui 2019


def trajectory(v=SPEED, depth=DEPTH, gap=GAP, hold=1.0, hold_end=0.5):
    """针尖 z(t)：下插 → 停 → 退回 → 停。返回 (zf, 总时长)。"""
    return V.piecewise(L + gap, [((gap + depth) / v, -v), (hold, 0.0), ((gap + depth) / v, v), (hold_end, 0.0)])


def build(root, dt=0.01, linear="ldl", threshold=THRESHOLD, tip_dist=TIP_DIST, zf=None, rayleigh_stiffness=BETA, rayleigh_mass=0.0):
    """搭建第 2 步场景（离线脚本和 GUI 场景共用）。返回各部件的 dict。
    rayleigh_stiffness：组织的 Rayleigh 刚度阻尼 β（EulerImplicitSolver.rayleighStiffness），默认 BETA = 0.1；
    rayleigh_mass：Rayleigh 质量阻尼 α（1/s），默认 0（只在对照检查里用）。"""
    if zf is None:
        zf, _ = trajectory()
    common.make_root(dt=dt, loop="free", root=root)
    xy = np.array([L / 2, L / 2])
    traj = lambda t: (np.array([xy[0], xy[1], zf(t) + LEN]), Q)
    nd, base, body, tip, local = common.add_kinematic_needle(root, length=LEN, n_body=NB, pose0=traj(0.0))
    pts = np.array([[0.0, 0.0, -LEN]])                       # 单点针尖（插件原生，默认）
    pmo, pg = contact.add_tip_patch(nd, pts)
    root.addObject(common.NeedleDriver(base, traj, dt, children=[(body, local), (tip, local[-1:]), (pmo, pts)],
                                       name="driver"))
    root.addObject("CollisionLoop")
    tis = contact.add_tissue_with_surface(root, L=L, E=V.E, nu=V.NU, linear=linear, mesh=common.load_mesh("W50"),
                                        rayleigh_stiffness=rayleigh_stiffness, rayleigh_mass=rayleigh_mass)
    _, bg = contact.add_needle_geometry(nd, body, tip, NB)
    algo, cu, ci = contact.add_puncture_insertion(root, pg, bg, tis, threshold=threshold, tip_dist=tip_dist)
    X = tis["X"]
    top = np.where(np.isclose(X[:, 2], L))[0]
    probe = int(top[np.argmin(np.linalg.norm(X[top, :2] - xy, axis=1))])   # 顶面中心最近的节点（约 0.8 mm 外）
    return dict(root=root, base=base, body=body, pmo=pmo, tis=tis, algo=algo, zf=zf, xy=xy, probe=probe, dt=dt)


def coupling_points(body):
    """由针身 MO 上的约束雅可比反推每个耦合点在针身上的对应点（每 3 行一组，第 1 行的点权重 = 各点系数的模）。"""
    J = body.constraint.value
    x = body.free_position.array()
    if J.shape[0] == 0 or J.nnz == 0:
        return np.zeros((0, 3))
    J = J.tocsr()
    rows = sorted(set(J.nonzero()[0]))
    out = []
    for r in rows[::3]:
        row = J.getrow(r).toarray().ravel()
        row = np.pad(row, (0, 3 * len(x) - len(row)))
        w = np.linalg.norm(row.reshape(-1, 3), axis=1)
        out.append(w @ x / w.sum())
    return np.array(out)


def run(dt=0.01, beta=BETA, linear="ldl", tag=None, alpha=0.0):
    """tag：输出文件名的后缀（默认 dt{dt}）；对照算例（不同 β / 求解器）用别的 tag，不覆盖主结果。"""
    zf, T = trajectory()
    root = Sofa.Core.Node("root")
    s = build(root, dt=dt, zf=zf, linear=linear, rayleigh_stiffness=beta, rayleigh_mass=alpha)
    Sofa.Simulation.init(root)
    tis = s["tis"]
    tets = np.array(tis["node"].topo.tetrahedra.value)
    V0 = tet_volumes(np.asarray(tis["X"]), tets)
    rec = {k: [] for k in ["t", "depth", "F", "F_tip", "F_shaft", "F_tissue", "n_cp", "cp_z", "probe_z",
                           "iters", "err", "n_rows", "vol_min"]}
    for k in range(int(round(T / dt))):
        Sofa.Simulation.animate(root, dt)
        t = (k + 1) * dt
        z = zf(t)
        lam = np.array(root.csolver.constraintForces.value)
        p_base = np.array([s["xy"][0], s["xy"][1], z + LEN])
        F, _, per = common.needle_wrench([s["pmo"], s["body"]], lam, dt, p_base)
        # 只读组织主 MO：表面 MO 上的约束行经 IdentityMapping 已经累加到主 MO 的雅可比里（两个都读会重复计数）
        Ft, _, _ = common.needle_wrench([tis["dofs"]], lam, dt, p_base)
        cp = coupling_points(s["body"])
        x = tis["dofs"].position.array()
        rec["t"].append(t); rec["depth"].append((L - z) * 1e3)
        rec["F"].append(F.copy()); rec["F_tip"].append(per[0].sum(0)); rec["F_shaft"].append(per[1].sum(0))
        rec["F_tissue"].append(Ft.copy())
        rec["n_cp"].append(len(cp)); rec["cp_z"].append(np.pad((L - cp[:, 2]) * 1e3, (0, 12 - len(cp)),
                                                               constant_values=np.nan)[:12])
        rec["probe_z"].append((x[s["probe"], 2] - L) * 1e3)
        rec["iters"].append(root.csolver.currentIterations.value); rec["err"].append(root.csolver.currentError.value)
        rec["n_rows"].append(len(lam))
        rec["vol_min"].append((tet_volumes(x, tets) / V0).min())
        if not np.isfinite(F).all() or not abs(rec["probe_z"][-1]) < 100:
            print(f"!! 第 {k+1} 步发散（力不是有限值或顶面节点位移 > 100 mm），停止", flush=True)
            break
    Sofa.Simulation.unload(root)
    rec = {k: np.array(v) for k, v in rec.items()}
    os.makedirs(OUT, exist_ok=True)
    np.savez(os.path.join(OUT, f"s2_puncture_{tag or f'dt{dt:g}'}.npz"), threshold=THRESHOLD, tip_dist=TIP_DIST,
             beta=beta, alpha=alpha, linear=linear, **rec)
    return rec


def summary(rec, dt, beta=BETA):
    d, Fz, n = rec["depth"], rec["F"][:, 2], rec["n_cp"]      # 作用在针上的力，Fz > 0：组织把针往上推（阻力）
    kp = int(np.argmax(n > 0))                                 # 第一个有耦合点的步（刺穿后的第一步）
    print(f"== dt = {dt:g} s，β = {beta:g} s，阈值 {THRESHOLD} N，耦合点间距 {TIP_DIST*1e3:.0f} mm", flush=True)
    print(f"   刺穿：第 {kp+1} 步（t = {rec['t'][kp]:.2f} s），针尖压深 {d[kp]:.2f} mm；"
          f"刺穿前一步（判据读的那一步）针尖力 {rec['F_tip'][kp-1][2]:.4f} N，再前一步 {rec['F_tip'][kp-2][2]:.4f} N", flush=True)
    print(f"   刺穿后第 1 步轴向力 {Fz[kp]:.4f} N，第 10 步 {Fz[kp+9]:.4f} N；顶面中心节点 z：刺穿前 {rec['probe_z'][kp-1]:.2f} mm，"
          f"刺穿后 0.5 s {rec['probe_z'][kp+int(0.5/dt)]:.2f} mm，1 s {rec['probe_z'][kp+int(1.0/dt)]:.2f} mm", flush=True)
    for dd in [10, 15, 20, 25, 30]:
        i = int(np.argmin(np.abs(d[:np.argmax(d)+1] - dd)))
        cz = rec["cp_z"][i][~np.isnan(rec["cp_z"][i])]
        print(f"   下插 {d[i]:.1f} mm：耦合点 {n[i]} 个，深度 {np.round(cz, 2).tolist()} mm；轴向力 {Fz[i]:.2e} N，"
              f"侧向 {np.hypot(*rec['F'][i, :2]):.2e} N", flush=True)
    imax = int(np.argmax(d))
    post = slice(kp, None)
    lat = np.hypot(rec["F"][post, 0], rec["F"][post, 1])
    print(f"   刺穿后：|轴向力| 最大 {np.abs(Fz[post]).max():.2e} N，侧向力最大 {lat.max():.2e} N", flush=True)
    back = np.where(np.arange(len(n)) > imax)[0]
    gone = back[n[back] == 0]
    if len(gone):
        print(f"   后退：耦合点在压深 {d[gone[0]]:.2f} mm（t = {rec['t'][gone[0]]:.2f} s）全部删除，回到刺穿前状态；"
              f"之后的最大针尖力 {np.abs(Fz[gone[0]:]).max():.4f} N", flush=True)
    drop = [(rec['t'][i], d[i], n[i-1], n[i]) for i in back if n[i] < n[i-1]]
    print(f"   后退时耦合点减少的时刻（t, 压深, 前, 后）：{[(round(a,2), round(b,2), c, e) for a,b,c,e in drop]}", flush=True)
    res = np.linalg.norm(rec["F"] + rec["F_tissue"], axis=1)
    print(f"   作用力 + 反作用力：最大 |F_针 + F_组织| = {res.max():.2e} N（|F_针| 最大 {np.linalg.norm(rec['F'], axis=1).max():.2e} N）", flush=True)
    print(f"   GS：迭代次数 中位数 {np.median(rec['iters']):.0f}，最大 {rec['iters'].max()}，最大约束行数 {rec['n_rows'].max()}，"
          f"最大误差 {rec['err'].max():.1e}；最差单元体积比 {rec['vol_min'].min():.2f}", flush=True)


if __name__ == "__main__":
    dt = float(sys.argv[1]) if len(sys.argv) > 1 else 0.01
    beta = float(sys.argv[2]) if len(sys.argv) > 2 else BETA
    rec = run(dt, beta, tag=None if beta == BETA else f"dt{dt:g}_beta{beta:g}")
    summary(rec, dt, beta)
