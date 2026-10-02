"""第 1.5 步 V5：多点针尖中"不受力的点"是真正脱离接触，还是约束冗余。见 docs/plan/step1_5_tip_survey.md §7。

测试台同 V0（v0_b2a_patch.py）：8 cm 立方体，h = 10 / 5 mm，swapping=True，R1，dt = 0.01，a_eff = 10 mm，
5 mm/s 压到 5 mm。在若干压深（小压深 0.05–0.5 mm 属线性范围，大压深 1–5 mm）逐行记录：
  - 该行属于哪个针尖点（针尖面片 MO 的 J 的非零列）、约束方向 n；
  - 投影到的组织表面节点和重心权重（表面 MO 的 J，行 = n ⊗ w）；
  - λ、求解后的间隙 g = n·(x_tip − Σ w x_node)（用步末的组织位置；针是运动学的，位置 = 自由位置）；
  - 每组相同节点（同一三角形 / 边 / 顶点）上的点数；表面 J 的秩（秩 < 行数 ⇒ 约束线性相关、W 奇异）。
判读：λ = 0 且 g 明显 > 0 ⇒ 真正脱离接触；λ = 0 且 g ≈ 0 ⇒ 处于接触但不受力（冗余 / 退化，力的分配不唯一）。
线性小变形理论下刚性平底压头整个接触面受压〔推导〕，所以小压深下的不受力点更可能是冗余。
运行：python3 v5_inactive_points.py → results/step1_5/v5_inactive_points.json
"""
import sys
sys.dont_write_bytecode = True
import os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../step1"))
import numpy as np
import Sofa, Sofa.Simulation
import common, contact
import t1_vertex as V

contact.load_plugins()
L, LEN, dt = V.L, V.LEN, 0.01
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../results/step1_5")
PROBE_MM = [0.05, 0.1, 0.25, 0.5, 1.0, 2.5, 5.0]


def rows_info(lam, pmo, tis):
    """逐行：针尖点号、方向、表面节点和权重、λ、间隙。"""
    Jp = pmo.constraint.value.tocsr()
    Js = tis["surf_dofs"].constraint.value.tocsr()
    xp = pmo.free_position.array()
    xs = tis["surf_dofs"].position.array()
    out = []
    for r in range(len(lam)):
        rp = Jp.getrow(r) if r < Jp.shape[0] else None
        if rp is None or rp.nnz == 0:
            continue
        cols, vals = rp.indices, rp.data
        ip = int(cols[0] // 3)
        n = np.zeros(3)
        for c, v in zip(cols, vals):
            n[c % 3] += v
        rs = Js.getrow(r)
        nodes = {}
        for c, v in zip(rs.indices, rs.data):
            nodes.setdefault(int(c // 3), np.zeros(3))[c % 3] += v
        # 表面一侧的行 = -n ⊗ w（作用力与反作用力），权重 w_j = -(row_j · n)/|n|²
        w = {j: float(-(vec @ n) / (n @ n)) for j, vec in nodes.items()}
        xproj = sum(wj * xs[j] for j, wj in w.items())
        g = float(n @ (xp[ip] - xproj) / np.linalg.norm(n))
        out.append(dict(row=r, point=ip, n=(n / np.linalg.norm(n)).tolist(), nodes=sorted(w), w=[w[j] for j in sorted(w)],
                        lam=float(lam[r]), gap=g, xp=xp[ip].tolist()))
    return out, Js


def surface_rank(Js, rows):
    if not rows:
        return 0, 0.0
    M = Js[[r["row"] for r in rows], :].toarray()
    s = np.linalg.svd(M, compute_uv=False)
    rank = int((s > 1e-10 * s[0]).sum())
    return rank, float(s[0] / s[-1]) if s[-1] > 0 else float("inf")


def run(n, shape, a_eff=0.010, depth=0.005, gap=0.0005, v=0.005):
    h = L / (n - 1)
    pts = contact.tip_patch_points(LEN, a_eff, h, shape)
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
    probes = []
    todo = list(PROBE_MM)
    for k in range(int(round(T / dt))):
        Sofa.Simulation.animate(root, dt)
        d_mm = (L - zf((k + 1) * dt)) * 1e3
        if todo and d_mm >= todo[0] - 1e-9:
            todo.pop(0)
            lam = np.array(cs.constraintForces.value)
            rows, Js = rows_info(lam, pmo, tis)
            rank, cond = surface_rank(Js, rows)
            groups = {}
            for r in rows:
                groups.setdefault(tuple(r["nodes"]), []).append(r["point"])
            F, M, _ = common.needle_wrench([pmo], lam, dt, np.array([L / 2, L / 2, zf((k + 1) * dt)]))
            probes.append(dict(depth_mm=d_mm, rows=rows, rank=rank, n_rows=len(rows), cond=cond,
                               groups=[dict(nodes=list(g), points=p) for g, p in groups.items()],
                               F=F.tolist(), M_tip=M.tolist(),
                               gs_iter=int(cs.currentIterations.value), gs_err=float(cs.currentError.value)))
    Sofa.Simulation.unload(root)
    return dict(n=n, h_mm=h * 1e3, shape=shape, a_eff_mm=a_eff * 1e3, points=pts.tolist(), probes=probes)


def summarize(r):
    print(f"== h={r['h_mm']:.1f} mm {r['shape']}，{len(r['points'])} 个点", flush=True)
    for p in r["probes"]:
        rows = p["rows"]
        lam = np.array([x["lam"] for x in rows]); g = np.array([x["gap"] for x in rows])
        act = lam > 1e-14
        lam_scale = lam.max() if len(lam) and lam.max() > 0 else 1.0
        ina = [(x["point"], x["gap"] * 1e3, len(x["nodes"])) for x in rows if x["lam"] <= 1e-14]
        maxgrp = max((len(gp["points"]) for gp in p["groups"]), default=0)
        print(f"  压深 {p['depth_mm']:4.2f} mm：行 {p['n_rows']:2d}，受力 {act.sum():2d}，表面 J 秩 {p['rank']:2d}，"
              f"同组最多 {maxgrp} 点；受力点间隙 |g|max {np.abs(g[act]).max()*1e3 if act.any() else 0:.1e} mm；"
              f"F_z {p['F'][2]:.4f} N，|M_tip| {np.linalg.norm(p['M_tip']):.1e} N·m，GS {p['gs_iter']}", flush=True)
        if ina:
            print("     不受力的点（点号, 间隙 mm, 投影节点数）：" + ", ".join(f"({i},{gg:+.3e},{nn})" for i, gg, nn in ina), flush=True)


if __name__ == "__main__":
    out = []
    for n, shape in [(9, "disk"), (17, "disk"), (17, "cap")]:
        r = run(n, shape)
        out.append(r)
        summarize(r)
        json.dump(out, open(os.path.join(OUT, "v5_inactive_points.json"), "w"))
