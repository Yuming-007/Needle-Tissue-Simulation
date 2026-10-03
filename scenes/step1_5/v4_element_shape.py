"""第 1.5 步 V4 补充：翻转前最差单元是"被竖直压扁"还是"被剪歪"。见 step1_5_results.md V4。

同 v4_large_depth.run 的测试台和加载；最差单元 V/V0 首次降到 0.2 以下时停下，对该单元：
  - 变形梯度 F = Dx · DX⁻¹（线性四面体内为常数），奇异值分解 F = U Σ Vᵀ：
    最小主伸长 σ_min 及其在变形后构形中的方向 u_min（U 的对应列）与竖直方向的夹角；
    竖直压扁 ⇒ u_min 接近竖直（夹角小）；剪切 ⇒ 夹角大，并伴随大的转动 / 歪斜；
  - 单元在 z 方向的高度比（变形后 / 初始）和水平尺寸比；
  - 位置（形心到针轴的距离、深度），属于哪一层；以及相对圆盘边缘 a 的位置。
运行：python3 v4_element_shape.py → results/step1_5/v4_element_shape.json
"""
import sys
sys.dont_write_bytecode = True
import os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import Sofa, Sofa.Simulation
import v4_large_depth as M

CASES = [(9, "point", 0.0), (9, "disk", 0.010), (9, "cap", 0.010), (17, "disk", 0.010), (17, "cap", 0.010)]


def analyse(n, tip, a, v=0.005, depth=0.020, gap=0.0005, vr_stop=0.2):
    L, LEN, dt, V = M.L, M.LEN, M.dt, M.V
    h = L / (n - 1)
    pts = M.contact.tip_patch_points(LEN, a, h, tip) if tip != "point" else np.array([[0.0, 0.0, -LEN]])
    zf, T = V.piecewise(L + gap, [(gap / v, -v), (depth / v, -v)])
    root = M.common.make_root(dt=dt, loop="free")
    traj = lambda t: (np.array([L / 2, L / 2, zf(t) + LEN]), V.Q)
    nd, base, body, tip_mo, local = M.common.add_kinematic_needle(root, length=LEN, n_body=V.NB, pose0=traj(0.0))
    pmo, pg = M.contact.add_tip_patch(nd, pts)
    root.addObject(M.common.NeedleDriver(base, traj, dt, children=[(body, local), (tip_mo, local[-1:]), (pmo, pts)],
                                         name="driver"))
    root.addObject("CollisionLoop")
    tis = M.contact.add_tissue_with_surface(root, n=n, L=L, E=V.E, nu=V.NU, linear="ldl")
    _, bg = M.contact.add_needle_geometry(nd, body, tip_mo, V.NB)
    M.contact.add_tip_contact(root, pg, bg, tis, distance=0.02)
    Sofa.Simulation.init(root)
    tets = np.array(tis["node"].topo.tetrahedra.value)
    X0 = np.asarray(tis["X"])
    V0 = M.tet_volumes(X0, tets)
    res = None
    for k in range(int(round(T / dt))):
        Sofa.Simulation.animate(root, dt)
        x = tis["dofs"].position.array()
        vr = M.tet_volumes(x, tets) / V0
        i = int(np.argmin(vr))
        if vr[i] < vr_stop:
            t = tets[i]
            DX = (X0[t[1:]] - X0[t[0]]).T
            Dx = (x[t[1:]] - x[t[0]]).T
            F = Dx @ np.linalg.inv(DX)
            U, S, Vt = np.linalg.svd(F)
            umin = U[:, 2]
            ang = float(np.degrees(np.arccos(min(1.0, abs(umin[2])))))
            ext0 = X0[t].max(0) - X0[t].min(0); ext = x[t].max(0) - x[t].min(0)
            cen = X0[t].mean(0)
            # 旋转部分（极分解 F = R S）的转角
            R = U @ Vt
            rot = float(np.degrees(np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1))))
            res = dict(n=n, h_mm=h * 1e3, tip=tip, a_mm=a * 1e3, depth_mm=(L - zf((k + 1) * dt)) * 1e3, vr=float(vr[i]),
                       r_mm=float(np.hypot(cen[0] - L / 2, cen[1] - L / 2) * 1e3), z_mm=float((L - cen[2]) * 1e3),
                       stretches=S.tolist(), angle_min_to_vertical_deg=ang, rotation_deg=rot,
                       height_ratio=float(ext[2] / ext0[2]), horiz_ratio=(ext[:2] / ext0[:2]).tolist(),
                       nodes_r_mm=(np.hypot(X0[t, 0] - L / 2, X0[t, 1] - L / 2) * 1e3).tolist(),
                       nodes_z_mm=((L - X0[t, 2]) * 1e3).tolist())
            break
    Sofa.Simulation.unload(root)
    return res


if __name__ == "__main__":
    out = []
    which = sys.argv[1:] or ["all"]
    for n, tip, a in CASES:
        if which != ["all"] and str(n) not in which:
            continue
        r = analyse(n, tip, a)
        out.append(r)
        print(f"h={r['h_mm']:4.1f} {tip:5s} a={r['a_mm']:4.1f}：压深 {r['depth_mm']:.2f} mm 时 V/V0 = {r['vr']:.3f}；单元位置 r = {r['r_mm']:.1f} mm、"
              f"深 {r['z_mm']:.1f} mm（节点 r {np.round(r['nodes_r_mm'],1).tolist()}，深 {np.round(r['nodes_z_mm'],1).tolist()}）；"
              f"主伸长 {np.round(r['stretches'],3).tolist()}，最小伸长方向与竖直夹角 {r['angle_min_to_vertical_deg']:.0f}°，"
              f"转角 {r['rotation_deg']:.0f}°；高度比 {r['height_ratio']:.2f}，水平尺寸比 {np.round(r['horiz_ratio'],2).tolist()}", flush=True)
        json.dump(out, open(os.path.join(M.OUT, f"v4_element_shape{'' if which == ['all'] else '_' + '_'.join(which)}.json"), "w"), indent=1)
