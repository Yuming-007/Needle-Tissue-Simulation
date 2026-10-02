"""第 1.5 步 V2：落点依赖。圆盘中心对准顶点 / 边中点 / 三角形内部不同位置时，刚度变化多少。见 step1_5_tip_survey.md §7。

第 1 步 P2 中单点针尖在 h = 10 mm 网格上：顶点 74.7 N/m，三角形内部 106–124 N/m（+43%–66%）。
这里在同一测试台（V0 / V1：a_eff = 10 mm 圆盘，点距约等于 h，R1，dt = 0.01）上把针轴平移到
xy = 中心 + (ox, oy)·h，偏移取第 1 步 P2 的 4 个位置，另加边中点 (0.5, 0) 和方格中心 (0.5, 0.5)；
每个位置同时跑单点针尖（同一脚本、同一轨迹）作对照。h = 10、5 mm。
加载：0.2 mm/s 准静态压到 0.05 mm（线性范围；V1 已证明此时刚接触的惯性过渡已衰减到约 1e-4），
k = F_z(0.05 mm) / 0.05 mm。另记录侧向力、压力中心偏移 = M_tip × ẑ / F_z（相对针轴）、受力点数。
判据：圆盘刚度在各落点之间的变化（max/min − 1）≤ 约 10%。
运行：python3 v2_placement.py → results/step1_5/v2_placement.json
"""
import sys
sys.dont_write_bytecode = True
import os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../step1"))
import numpy as np
import Sofa, Sofa.Simulation
import common, contact
import t1_vertex as V

contact.load_plugins()
L, LEN, dt, A_EFF = V.L, V.LEN, 0.01, 0.010
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../results/step1_5")
OFFSETS = [(0.0, 0.0), (0.30, 0.15), (0.55, 0.25), (-0.40, 0.20), (0.5, 0.0), (0.5, 0.5)]


def run(n, tip, off, gap=0.0001, v=0.0002, depth=0.00005):
    h = L / (n - 1)
    pts = contact.tip_patch_points(LEN, A_EFF, h, "disk") if tip == "disk" else np.array([[0.0, 0.0, -LEN]])
    xy = np.array([L / 2 + off[0] * h, L / 2 + off[1] * h])
    zf, T = V.piecewise(L + gap, [(gap / v, -v), (depth / v, -v)])
    root = common.make_root(dt=dt, loop="free")
    traj = lambda t: (np.array([xy[0], xy[1], zf(t) + LEN]), V.Q)
    nd, base, body, tip_mo, local = common.add_kinematic_needle(root, length=LEN, n_body=V.NB, pose0=traj(0.0))
    pmo, pg = contact.add_tip_patch(nd, pts)
    root.addObject(common.NeedleDriver(base, traj, dt, children=[(body, local), (tip_mo, local[-1:]), (pmo, pts)],
                                       name="driver"))
    root.addObject("CollisionLoop")
    tis = contact.add_tissue_with_surface(root, n=n, L=L, E=V.E, nu=V.NU, linear="ldl")
    _, bg = contact.add_needle_geometry(nd, body, tip_mo, V.NB)
    contact.add_tip_contact(root, pg, bg, tis, distance=0.02)
    Sofa.Simulation.init(root)
    a = time.perf_counter()
    for k in range(int(round(T / dt))):
        Sofa.Simulation.animate(root, dt)
    t_run = time.perf_counter() - a
    cs = root.csolver
    lam = np.array(cs.constraintForces.value)
    z_tip = zf(T)
    F, M, _ = common.needle_wrench([pmo], lam, dt, np.array([xy[0], xy[1], z_tip]))
    Sofa.Simulation.unload(root)
    d_real = L - z_tip
    cop = np.array([-M[1], M[0]]) / F[2]        # 压力中心相对针轴的偏移：M = r × F，F ≈ F_z ẑ ⇒ r_xy = (−M_y, M_x)/F_z
    return dict(n=n, h_mm=h * 1e3, tip=tip, offset=list(off), xy=xy.tolist(), depth_mm=d_real * 1e3,
                F=F.tolist(), M_tip=M.tolist(), k=float(F[2] / d_real), lateral_over_axial=float(np.abs(F[:2]).max() / abs(F[2])),
                cop_mm=(cop * 1e3).tolist(), n_active=int((lam > 1e-12 * max(lam.max(), 1e-30)).sum()), n_points=int(len(pts)),
                t_run=t_run)


if __name__ == "__main__":
    out = []
    for n in [9, 17]:
        for tip in ["point", "disk"]:
            ks = []
            for off in OFFSETS:
                r = run(n, tip, off)
                out.append(r); ks.append(r["k"])
                print(f"h={r['h_mm']:4.1f}mm {tip:5s} 偏移 ({off[0]:+.2f}h,{off[1]:+.2f}h)：k = {r['k']:8.3f} N/m，"
                      f"受力点 {r['n_active']}/{r['n_points']}，侧向/轴向 {r['lateral_over_axial']:.1e}，"
                      f"压力中心偏移 ({r['cop_mm'][0]:+.3f}, {r['cop_mm'][1]:+.3f}) mm，{r['t_run']:.0f} s", flush=True)
                json.dump(out, open(os.path.join(OUT, "v2_placement.json"), "w"), indent=1)
            ks = np.array(ks)
            print(f"  → h={L/(n-1)*1e3:.0f}mm {tip}：k 范围 {ks.min():.3f}–{ks.max():.3f} N/m，max/min − 1 = {ks.max()/ks.min()-1:+.1%}，"
                  f"相对顶点落点 {ks.min()/ks[0]-1:+.1%} ~ {ks.max()/ks[0]-1:+.1%}", flush=True)
