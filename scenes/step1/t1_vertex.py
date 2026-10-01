"""第 1 步：顶点接触测试（P1、P3、P4、P5、P6、P7）。见 docs/plan/step1_tenting.md。

测试台：8 cm 立方体，n = 9，swapping=True（中心节点处镜像对称），FTCfix polar，ν = 0.45，E = 5 kPa，
底面固定，重力关闭；动力学 R1（SparseLDLSolver，确定性），dt = 0.01 s，不加 Rayleigh。
针：纯运动学，竖直（针轴 -z），针尖对准顶面中心节点。

接触运行（A1：插件 InsertionAlgorithm + ConstraintUnilateral + SecondDirection，不刺穿）：
  针尖从表面上方 gap 处出发，以 v 下压到压深 depth，保持 t_hold，再以 v_back 后退到原来的高度，再保持。
  gap 取 v·dt 的整数倍，使针尖恰好在某一步结束时到达表面。
精确参照（P1a）：同一组织，用 SOFA 自带的 BilateralLagrangianConstraint 把针尖直接连到中心节点，
  针尖从节点处出发，以同样的 v 下压同样的深度。和接触运行"到达表面之后"的阶段逐步比较 λ。
  双边约束是 3 维的；由于网格在中心对称，x、y 方向的 λ 应为 0。
物理参照（P1b）：同一网格的静力解 F_static(δ)（中心节点给定 z 位移，x、y 自由）。
读数（每步）：λ（力 = λ/dt）、针根力（common.needle_wrench）、约束作用在组织上的力（表面 MO 的 J）、
  约束方向 n（针尖 MO 的 J）、接触点所在节点和权重、中心节点位置、底面反力（组织 force，对应上一步构形）。
运行：python3 t1_vertex.py → results/step1/t1_vertex.json
"""
import sys
sys.dont_write_bytecode = True
import os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
import Sofa, Sofa.Simulation
import common, contact

contact.load_plugins()
L, N, NU, E, LEN, NB = 0.08, 9, 0.45, 5000.0, 0.10, 11
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../results/step1")
os.makedirs(OUT, exist_ok=True)
X = common.grid_nodes(L, N)
C = int(np.where(np.all(np.isclose(X, [L / 2, L / 2, L]), axis=1))[0][0])
BOT = np.where(np.isclose(X[:, 2], 0))[0]
Q = np.array([0, 0, 0, 1.0])


def piecewise(z0, segs):
    """针尖 z 的分段线性轨迹：segs = [(持续时间, 速度), ...]，速度向下为负。返回 f(t) -> z。"""
    ts = np.cumsum([0] + [s[0] for s in segs]); zs = [z0]
    for d, v in segs:
        zs.append(zs[-1] + v * d)

    def f(t):
        for i, (d, v) in enumerate(segs):
            if t <= ts[i + 1] + 1e-12:
                return zs[i] + v * (t - ts[i])
        return zs[-1]
    return f, ts[-1]


def build(dt, zf, linear="ldl", with_contact=True):
    root = common.make_root(dt=dt, loop="free")
    traj = lambda t: (np.array([L / 2, L / 2, zf(t) + LEN]), Q)
    nd, base, body, tip, local = common.add_kinematic_needle(root, length=LEN, n_body=NB, pose0=traj(0.0))
    root.addObject(common.NeedleDriver(base, traj, dt, children=[(body, local), (tip, local[-1:])], name="driver"))
    root.addObject("CollisionLoop")
    tis = contact.add_tissue_with_surface(root, n=N, L=L, E=E, nu=NU, linear=linear)
    tg, bg = contact.add_needle_geometry(nd, body, tip, NB)
    if with_contact:
        contact.add_tip_contact(root, tg, bg, tis)
    else:   # 精确参照：针尖直接用双边约束连到中心节点
        root.addObject("BilateralLagrangianConstraint", name="bilat", object1=tip.getLinkPath(),
                       object2=tis["dofs"].getLinkPath(), first_point=[0], second_point=[C])
    return root, dict(base=base, body=body, tip=tip, tis=tis)


def record(root, h, dt, rec, f_prev):
    lam = np.array(root.csolver.constraintForces.value)
    p_base = np.asarray(h["base"].free_position.value)[0][:3]
    F, tau, _ = common.needle_wrench([h["body"], h["tip"]], lam, dt, p_base)
    tis = h["tis"]
    # 约束作用在组织上的力：只取约束直接建立的那一层（接触 → 表面 MO；双边 → 体 MO）。
    # 表面 MO 经 IdentityMapping 挂在体 MO 上，约束行也会被映射到体 MO，两层都加会重复计数。
    Ft = np.zeros(3); nodes_w = []
    for mo in (tis["surf_dofs"], tis["dofs"]):
        J = mo.constraint.value
        if J.shape[0] == 0 or J.nnz == 0:
            continue
        lp = np.zeros(J.shape[0]); m = min(len(lam), J.shape[0]); lp[:m] = lam[:m]
        g = np.zeros(3 * len(mo.position.array())); g[:J.shape[1]] = J.T @ lp
        f = g.reshape(-1, 3) / dt
        Ft += f.sum(0)
        row = np.zeros(3 * len(mo.position.array())); row[:J.shape[1]] = J.getrow(0).toarray().ravel()
        row = row.reshape(-1, 3)
        nodes_w = [(int(i), row[i].tolist()) for i in np.where(np.abs(row).sum(1) > 0)[0]]
        break
    Jt = h["tip"].constraint.value
    n_dir = (Jt.getrow(0).toarray().ravel()[:3].tolist() if Jt.shape[0] and Jt.nnz else [0, 0, 0])
    x = tis["dofs"].position.array()
    fint = tis["dofs"].force.array()
    rec["t"].append(len(rec["t"]) * dt + dt)
    rec["tip_z"].append(float(h["tip"].position.array()[0][2]))
    rec["node"].append(x[C].tolist())
    rec["lam"].append(lam.tolist())
    rec["F_needle"].append(F.tolist()); rec["tau"].append(tau.tolist())
    rec["F_tissue"].append(Ft.tolist())
    rec["n"].append(n_dir); rec["contact_nodes"].append(nodes_w)
    rec["F_bottom_prev"].append(float(-fint[BOT, 2].sum()))   # 上一步结束时构形的底面反力（向上为正）


def run(dt, zf, t_end, with_contact=True, linear="ldl"):
    root, h = build(dt, zf, linear=linear, with_contact=with_contact)
    Sofa.Simulation.init(root)
    rec = {k: [] for k in ["t", "tip_z", "node", "lam", "F_needle", "tau", "F_tissue", "n", "contact_nodes", "F_bottom_prev"]}
    ts = []
    for k in range(int(round(t_end / dt))):
        a = time.perf_counter(); Sofa.Simulation.animate(root, dt); ts.append(time.perf_counter() - a)
        record(root, h, dt, rec, None)
    Sofa.Simulation.unload(root)
    rec["t_step"] = ts
    return rec


def static_curve(depths):
    """物理参照：中心节点给定 z 位移（x、y 自由）的静力解，返回压头反力（向上为正）。"""
    out = []
    for d in depths:
        root = common.make_root(); common.add_static_solver(root, linear="ldl", abs_tol=1e-10)
        t, mo, _ = common.add_tissue(root, L=L, n=N, E=E, nu=NU, swapping=True)
        cons = [t.addObject("FixedProjectiveConstraint", indices=BOT.tolist()),
                t.addObject("PartialFixedProjectiveConstraint", indices=[C], fixedDirections=[0, 0, 1])]
        Sofa.Simulation.init(root)
        p = X.copy(); p[C, 2] -= d; mo.position.value = p
        Sofa.Simulation.animate(root, 0.01)
        f = common.internal_force(root, t, cons)
        out.append(float(f[C, 2]))
        Sofa.Simulation.unload(root)
    return np.array(out)


if __name__ == "__main__":
    res = {}
    depth, gap = 0.010, 0.001
    dtab = np.linspace(0, depth, 11)
    Ftab = np.r_[0.0, static_curve(dtab[1:])]
    res["static_curve"] = dict(depth=dtab.tolist(), F=Ftab.tolist())
    print("静力曲线 F(δ)：" + "  ".join(f"{d*1e3:.0f}mm:{F:.4f}N" for d, F in zip(dtab, Ftab)), flush=True)

    # ---- P1a / P1b / P3 / P4 / P6 / P7：v = 1 mm/s，dt = 0.01
    v, dt = 0.001, 0.01
    t_gap, t_push = gap / v, depth / v
    zf, T = piecewise(L + gap, [(t_gap, -v), (t_push, -v), (1.0, 0.0), (t_push + t_gap, +v), (1.0, 0.0)])
    rc = run(dt, zf, T)
    zr, Tr = piecewise(L, [(t_push, -v), (1.0, 0.0)])
    rb = run(dt, zr, Tr, with_contact=False)
    res["contact_v1"] = rc; res["bilateral_v1"] = rb

    k0 = int(round(t_gap / dt))                       # 到达表面的那一步（之后的步和参照逐步对应）
    nb = len(rb["lam"])
    lam_c = np.array([l[0] if l else 0.0 for l in rc["lam"][k0:k0 + nb]])
    lam_b = np.array([l[2] for l in rb["lam"]])       # 双边约束的 z 分量：δ = x_节点 - x_针尖，节点受 +λ
    lam_bxy = np.abs(np.array([l[:2] for l in rb["lam"]])).max()
    # 接触：组织受到的力沿 -z（向下）；双边：节点受到 +λ_z，推压时 λ_z < 0。比较组织受力的大小
    Fc = np.array(rc["F_tissue"][k0:k0 + nb])[:, 2]; Fb = np.array(rb["F_tissue"])[:, 2]
    p1a = float(np.abs(Fc - Fb).max() / np.abs(Fb).max())
    print(f"P1a 接触 vs 双边（同一动力学，逐步）：组织受力最大相对差 {p1a:.2e}；双边 x/y 分量最大 |λ| = {lam_bxy:.1e}", flush=True)

    tipz = np.array(rc["tip_z"]); dlt = L - tipz
    Fn = np.array(rc["F_needle"])[:, 2]
    push = (np.arange(len(dlt)) >= k0 + int(0.2 / dt)) & (np.arange(len(dlt)) < k0 + int(t_push / dt))
    rel = Fn[push] / np.interp(dlt[push], dtab, Ftab) - 1
    print(f"P1b 准静态（1 mm/s，起步 0.2 s 之后）：针根力相对静力解 最大 {np.abs(rel).max():.2e}，平均 {rel.mean():+.2e}", flush=True)

    lam_all = np.array([l[0] if l else 0.0 for l in rc["lam"]])
    nodez = np.array(rc["node"])[:, 2]
    sep = tipz > L + 1e-9                               # 针尖在未变形表面之上
    print(f"P3 单边性：λ 最小值 {lam_all.min():.2e}；针尖在表面之上时 λ 最大 {lam_all[sep].max():.2e}；"
          f"结束时中心节点 z 偏离静止位置 {abs(nodez[-1]-L):.2e} m", flush=True)

    Ft = np.array(rc["F_tissue"]); Fnv = np.array(rc["F_needle"])
    p4a = float(np.abs(Fnv + Ft).max() / np.abs(Fnv).max())
    hold = (np.arange(len(dlt)) >= k0 + int(t_push / dt) + 50) & (np.arange(len(dlt)) < k0 + int((t_push + 1.0) / dt))
    Fbot = np.array(rc["F_bottom_prev"])
    p4b = float(np.abs(Fbot[1:][hold[:-1]] - Fn[:-1][hold[:-1]]).max() / Fn[hold].max())
    print(f"P4 作用力 = 反作用力：|F_针 + F_组织| / |F_针| 最大 {p4a:.2e}；保持阶段 底面反力 vs 针根力 相对差 {p4b:.2e}", flush=True)

    nd = np.array(rc["n"]); act = lam_all > 0
    ang = np.degrees(np.arccos(np.clip(np.abs(nd[act][:, 2]) / np.linalg.norm(nd[act], axis=1), -1, 1)))
    print(f"P6 约束方向和 z 轴的最大夹角 {ang.max():.2e} 度", flush=True)
    cn = [r for r, a in zip(rc["contact_nodes"], act) if a]
    only_c = all(len(r) == 1 and r[0][0] == C for r in cn)
    lat = np.abs(np.array(rc["node"])[:, :2] - L / 2).max()
    print(f"P7 接触时只作用在中心节点上：{only_c}；中心节点最大侧向位移 {lat:.2e} m", flush=True)
    res["summary_v1"] = dict(P1a=p1a, P1a_bilat_xy=float(lam_bxy), P1b_max=float(np.abs(rel).max()), P1b_mean=float(rel.mean()),
                             P3_lam_min=float(lam_all.min()), P3_lam_above=float(lam_all[sep].max()),
                             P3_end_offset=float(abs(nodez[-1] - L)), P4a=p4a, P4b=p4b, P6_deg=float(ang.max()),
                             P7_only_center=bool(only_c), P7_lateral=float(lat))

    # ---- P3（快速后退）：v_back = 1 m/s（两步退出 20 mm，到表面以上 10 mm），快于组织回弹（约 10 mm × 44 rad/s ≈ 0.44 m/s），
    # 组织应和针尖分离：针尖在未变形表面之上时 λ = 0，组织自由振动后回到静止位置。
    # （第一版用 50 mm/s：组织回弹比针尖后退快，一直追着针尖，λ > 0 是正确的，不能用来检验分离。）
    zf2, T2 = piecewise(L + gap, [(t_gap, -v), (t_push, -v), (0.5, 0.0), (0.02, +1.0), (2.0, 0.0)])
    rf = run(dt, zf2, T2)
    lam_f = np.array([l[0] if l else 0.0 for l in rf["lam"]])
    tz = np.array(rf["tip_z"]); nz = np.array(rf["node"])[:, 2]
    i_back = int(round((t_gap + t_push + 0.5) / dt))
    above = tz > L + 1e-9
    print(f"P3 快速后退（1 m/s）：针尖在未变形表面之上的步数 {int(above[i_back:].sum())}，这些步 λ 最大 {lam_f[i_back:][above[i_back:]].max():.2e}；"
          f"后退后中心节点最高到 {nz[i_back:].max()*1e3:.3f} mm（静止位置 {L*1e3:.1f} mm），结束时偏离 {abs(nz[-1]-L)*1e3:.4f} mm", flush=True)
    res["fast_retract"] = rf
    res["summary_v1"].update(P3_fast_lam_above=float(lam_f[i_back:][above[i_back:]].max()),
                             P3_fast_end_offset=float(abs(nz[-1] - L)))
    # 后退阶段接触点偏离节点（漏斗形凹坑里最近点投影落在侧壁上）
    t_arr = np.array(rc["t"]); retr = (t_arr > t_gap + t_push + 1.0) & (lam_all > 0)
    wmax = max([max(abs(nw[1][2]) for nw in rc["contact_nodes"][i] if nw[0] != C) for i in np.where(retr)[0]
                if not (len(rc["contact_nodes"][i]) == 1 and rc["contact_nodes"][i][0][0] == C)] or [0.0])
    pushm = (t_arr > t_gap) & (t_arr <= t_gap + t_push) & (lam_all > 0)
    only_push = all(len(rc["contact_nodes"][i]) == 1 and rc["contact_nodes"][i][0][0] == C for i in np.where(pushm)[0])
    angp = float(ang[pushm[act]].max()) if pushm.any() else 0.0
    print(f"P6/P7 按阶段：推压阶段只作用在中心节点 {only_push}，推压阶段方向最大夹角 {angp:.2e} 度；"
          f"后退阶段邻近节点最大权重 {wmax:.2e}", flush=True)
    res["summary_v1"].update(P7_push_only_center=bool(only_push), P6_push_deg=angp, retract_neighbor_weight=float(wmax))

    # ---- P5：dt 和速度
    res["P5"] = []
    for v5 in [0.001, 0.005]:
        for dt5 in [0.005, 0.01, 0.02]:
            tg5, tp5 = gap / v5, depth / v5
            zf5, T5 = piecewise(L + gap, [(tg5, -v5), (tp5, -v5)])
            r5 = run(dt5, zf5, T5)
            tz5 = np.array(r5["tip_z"]); F5 = np.array(r5["F_needle"])[:, 2]
            k5 = int(round(tg5 / dt5))
            m5 = (np.arange(len(tz5)) >= k5) & ((L - tz5) >= 0.002)        # 压深 ≥ 2 mm（起步瞬态之后）
            rel5 = F5[m5] / np.interp(L - tz5[m5], dtab, Ftab) - 1
            res["P5"].append(dict(v=v5, dt=dt5, max_rel=float(np.abs(rel5).max()), mean_rel=float(rel5.mean())))
            print(f"P5 v={v5*1e3:.0f}mm/s dt={dt5}: 压深 ≥ 2 mm 时针根力相对静力解 最大 {np.abs(rel5).max():.2e}，平均 {rel5.mean():+.2e}", flush=True)

    json.dump(res, open(os.path.join(OUT, "t1_vertex.json"), "w"))
