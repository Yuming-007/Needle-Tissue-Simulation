"""第 1.5 步 V4：大压深下的数值稳健性包络（numerical robustness envelope）。见 step1_5_tip_survey.md §7、
step1_5_results.md V4。

问题：有限支撑针尖压得很深时，组织单元会不会翻转（数值崩溃）、从多深开始。
  - 只判断数值稳健性，不判断力是否符合物理（共旋线弹性在大压缩下本身不符合物理）。
  - 输出"安全压深"包络，不用来迁就刺穿阈值：阈值由物理模型 / 数据决定，如果合理阈值超出包络，应回头改网格 /
    针尖 / 材料表示（ChatGPT 复核 2026-10-02）。20 mm 只是压力测试的上限，不是真实针的验证目标（a_eff ≠ R_phys）。
预测〔推导，粗估〕：平底圆盘边缘应力奇异，边缘单元剪切最大；边缘一个单元宽度内的位移差约 0.3–0.4 d，和 h 相当时
  翻转 ⇒ d_inv ≈ 3h，且网格越细越早翻转（d_inv/h 大致不变）；半球没有尖锐边缘，应明显更稳。
测试台：全模型（插件不能给对称面上的点施加一半的力，所以不用四分之一模型），R1，dt = 0.01，swapping=True，ν = 0.45；
  针尖最低点从表面上方 0.5 mm 出发，匀速压到 20 mm。
每步记录：压深、F、M_tip（针尖参考点）、受力点数、GS 迭代 / 残差、最差单元的有符号体积比 V/V0 及其位置
  （单元形心到针轴的距离、深度）、每步耗时、是否出现非有限值。
d_0.1 = V/V0 首次 < 0.1 的压深（严重扭曲预警），d_inv = 首次 ≤ 0 的压深（翻转）。翻转后再跑 10 步后停止。
运行：python3 v4_large_depth.py <组名>   组名：A（h = 10 全部 + h = 5 全部）、B（h = 4）、C（速度检查 1 mm/s）、M（M3 筛选：A+、B+）
     → results/step1_5/v4_large_depth_<组名>.json
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
L, LEN, dt = V.L, V.LEN, 0.01
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../results/step1_5")

# (n, 针尖, a_eff, 速度, 备注)
GROUPS = {
    "A": [(9, "point", 0.0, 0.005, ""), (9, "disk", 0.005, 0.005, "h/a = 2：超出适用范围的负对照"),
          (9, "disk", 0.010, 0.005, "h/a = 1：适用范围边界"), (9, "cap", 0.010, 0.005, "h/a = 1"),
          (17, "point", 0.0, 0.005, ""), (17, "disk", 0.005, 0.005, "h/a = 1：适用范围边界"),
          (17, "disk", 0.010, 0.005, ""), (17, "cap", 0.010, 0.005, "")],
    "B": [(21, "disk", 0.010, 0.005, ""), (21, "cap", 0.010, 0.005, "")],
    "C": [(17, "disk", 0.010, 0.001, "速度检查：对照 A 组 5 mm/s")],
    "M": [("A+", "disk", 0.010, 0.005, "M3 筛选：非均匀网格 A+"), ("B+", "disk", 0.010, 0.005, "M3 筛选：非均匀网格 B+")],
}


def tet_volumes(x, tets):
    a, b, c, d = (x[tets[:, i]] for i in range(4))
    return np.einsum("ij,ij->i", np.cross(b - a, c - a), d - a) / 6.0


def run(n, tip, a, v, depth=0.020, gap=0.0005):
    # n：规则网格每边节点数；或者字符串 = 非均匀网格名（meshes/step1_5/<名>.npz，M3 筛选起使用），此时 h = 标称 h_local
    mesh = common.load_mesh(n) if isinstance(n, str) else None
    h = float(mesh["h_local"]) if mesh is not None else L / (n - 1)
    pts = contact.tip_patch_points(LEN, a, h, tip) if tip != "point" else np.array([[0.0, 0.0, -LEN]])
    zf, T = V.piecewise(L + gap, [(gap / v, -v), (depth / v, -v)])
    root = common.make_root(dt=dt, loop="free")
    traj = lambda t: (np.array([L / 2, L / 2, zf(t) + LEN]), V.Q)
    nd, base, body, tip_mo, local = common.add_kinematic_needle(root, length=LEN, n_body=V.NB, pose0=traj(0.0))
    pmo, pg = contact.add_tip_patch(nd, pts)
    root.addObject(common.NeedleDriver(base, traj, dt, children=[(body, local), (tip_mo, local[-1:]), (pmo, pts)],
                                       name="driver"))
    root.addObject("CollisionLoop")
    tis = contact.add_tissue_with_surface(root, n=n, L=L, E=V.E, nu=V.NU, linear="ldl") if mesh is None else \
        contact.add_tissue_with_surface(root, L=L, E=V.E, nu=V.NU, linear="ldl", mesh=mesh)
    _, bg = contact.add_needle_geometry(nd, body, tip_mo, V.NB)
    contact.add_tip_contact(root, pg, bg, tis, distance=0.02)
    Sofa.Simulation.init(root)
    tets = np.array(tis["node"].topo.tetrahedra.value)
    X0 = np.asarray(tis["X"])
    V0 = tet_volumes(X0, tets)
    cs = root.csolver
    rec = {k: [] for k in ["depth", "F", "M", "n_active", "gs_iter", "gs_err", "vr_min", "vr_r", "vr_z", "t_step"]}
    d01 = dinv = None
    stop_at = None
    finite = True
    nsteps = int(round(T / dt))
    for k in range(nsteps):
        a0 = time.perf_counter(); Sofa.Simulation.animate(root, dt); ts = time.perf_counter() - a0
        z_tip = zf((k + 1) * dt)
        d = L - z_tip
        lam = np.array(cs.constraintForces.value)
        x = tis["dofs"].position.array()
        if not (np.isfinite(x).all() and np.isfinite(lam).all()):
            finite = False
            break
        F, M, _ = common.needle_wrench([pmo], lam, dt, np.array([L / 2, L / 2, z_tip]))
        vr = tet_volumes(x, tets) / V0
        i = int(np.argmin(vr))
        cen = X0[tets[i]].mean(0)
        if d > 0:
            rec["depth"].append(float(d)); rec["F"].append(F.tolist()); rec["M"].append(M.tolist())
            rec["n_active"].append(int((lam > 1e-12 * max(lam.max(), 1e-30)).sum()) if len(lam) else 0)
            rec["gs_iter"].append(int(cs.currentIterations.value)); rec["gs_err"].append(float(cs.currentError.value))
            rec["vr_min"].append(float(vr[i])); rec["t_step"].append(ts)
            rec["vr_r"].append(float(np.hypot(cen[0] - L / 2, cen[1] - L / 2))); rec["vr_z"].append(float(L - cen[2]))
            if d01 is None and vr[i] < 0.1:
                d01 = float(d)
            if dinv is None and vr[i] <= 0.0:
                dinv = float(d); stop_at = k + 10
        if stop_at is not None and k >= stop_at:
            break
    Sofa.Simulation.unload(root)
    return dict(n=n if mesh is None else str(n), h_mm=h * 1e3, tip=tip, a_mm=a * 1e3, v_mm_s=v * 1e3, n_points=int(len(pts)),
                d01_mm=None if d01 is None else d01 * 1e3, dinv_mm=None if dinv is None else dinv * 1e3,
                finite=finite, depth_reached_mm=rec["depth"][-1] * 1e3 if rec["depth"] else 0.0, rec=rec)


def fmt(x):
    return "> 20" if x is None else f"{x:.2f}"


if __name__ == "__main__":
    g = sys.argv[1]
    out = []
    for n, tip, a, v, note in GROUPS[g]:
        t0 = time.perf_counter()
        r = run(n, tip, a, v)
        r["note"] = note
        out.append(r)
        rr = r["rec"]
        i = int(np.argmin(rr["vr_min"])) if rr["vr_min"] else 0
        h = r["h_mm"]
        ratio = lambda d, s: "—" if d is None or s == 0 else f"{d/s:.2f}"
        print(f"[{g}] h={h:4.1f} {tip:5s} a={r['a_mm']:4.1f} v={r['v_mm_s']:.0f}mm/s：d_0.1 = {fmt(r['d01_mm'])} mm"
              f"（d/h {ratio(r['d01_mm'], h)}，d/a {ratio(r['d01_mm'], r['a_mm'])}），d_inv = {fmt(r['dinv_mm'])} mm"
              f"（d/h {ratio(r['dinv_mm'], h)}，d/a {ratio(r['dinv_mm'], r['a_mm'])}）；最差 V/V0 {min(rr['vr_min']):.3f} "
              f"在 r = {rr['vr_r'][i]*1e3:.1f} mm、深 {rr['vr_z'][i]*1e3:.1f} mm；到达 {r['depth_reached_mm']:.1f} mm，"
              f"F_z(末) {rr['F'][-1][2]:.3f} N，GS 最多 {max(rr['gs_iter'])} 次，有限 {r['finite']}，"
              f"{np.median(rr['t_step'])*1e3:.0f} ms/步，{time.perf_counter()-t0:.0f} s {note}", flush=True)
        json.dump(out, open(os.path.join(OUT, f"v4_large_depth_{g}.json"), "w"))
