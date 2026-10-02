"""第 1.5 步 V0（B2a 可行性）：针尖一侧取样的有限半径针尖（圆盘 / 半球）+ 插件多点接触。

测试台：8 cm 立方体，n = 9（h = 10 mm）和 n = 17（h = 5 mm），swapping=True，FTCfix polar，ν = 0.45，E = 5 kPa，
R1（SparseLDLSolver），dt = 0.01。针尖面片：contact.tip_patch_points（同心圆环，点距约等于 h），挂在纯运动学针上；
接触：插件 InsertionAlgorithm（多点）+ ConstraintUnilateral + SecondDirection，不刺穿。
针尖最低点从表面上方 0.5 mm 出发，5 mm/s 压到 5 mm（a_eff = 10 mm）。
记录每步：接触状态的点数（λ > 0）、总轴向力、侧向力、GS 迭代次数和残差。
回答 V0 的问题：能否稳定产生非零反力；接触数量 / 位置 / 总力是否合理；总力能否可靠读出。
运行：python3 v0_b2a_patch.py → results/step1_5/v0_b2a_patch.json
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


def run(n, shape, a_eff, depth=0.005, gap=0.0005, v=0.005, linear="ldl"):
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
    tis = contact.add_tissue_with_surface(root, n=n, L=L, E=V.E, nu=V.NU, linear=linear)
    _, bg = contact.add_needle_geometry(nd, body, tip, V.NB)
    contact.add_tip_contact(root, pg, bg, tis, distance=0.02)
    Sofa.Simulation.init(root)
    cs = root.csolver
    rec = dict(depth=[], F=[], n_active=[], n_rows=[], gs_iter=[], gs_err=[], lam_max=[])
    ts = []
    for k in range(int(round(T / dt))):
        a = time.perf_counter(); Sofa.Simulation.animate(root, dt); ts.append(time.perf_counter() - a)
        lam = np.array(cs.constraintForces.value)
        F, _, _ = common.needle_wrench([pmo], lam, dt, np.asarray(base.free_position.value)[0][:3])
        rec["depth"].append(float(L - zf((k + 1) * dt)))
        rec["F"].append(F.tolist()); rec["n_active"].append(int((lam > 1e-14).sum())); rec["n_rows"].append(int(len(lam)))
        rec["lam_max"].append(float(lam.max()) if len(lam) else 0.0)
        rec["gs_iter"].append(int(cs.currentIterations.value) if hasattr(cs, "currentIterations") else -1)
        rec["gs_err"].append(float(cs.currentError.value) if hasattr(cs, "currentError") else -1)
    Sofa.Simulation.unload(root)
    rec.update(n=n, h_mm=h * 1e3, shape=shape, a_eff_mm=a_eff * 1e3, n_points=int(len(pts)), t_step_ms=1e3 * float(np.median(ts)))
    return rec


if __name__ == "__main__":
    out = []
    for n in [9, 17]:
        for shape in ["disk", "cap"]:
            r = run(n, shape, 0.010)
            out.append(r)
            d = np.array(r["depth"]); F = np.array(r["F"])
            for dd in [1.0, 2.5, 5.0]:
                i = int(np.argmin(np.abs(d - dd * 1e-3)))
                print(f"h={r['h_mm']:4.1f}mm {shape:4s} a_eff=10mm 点数 {r['n_points']:2d}：压深 {dd:3.1f} mm 时 F_ax = {F[i,2]:.4f} N，"
                      f"接触点 {r['n_active'][i]:2d}/{r['n_rows'][i]:2d}，侧向 / 轴向 = {np.abs(F[i,:2]).max()/max(abs(F[i,2]),1e-30):.1e}，"
                      f"GS 迭代 {r['gs_iter'][i]}，残差 {r['gs_err'][i]:.1e}", flush=True)
            print(f"   λ 最大值 {max(r['lam_max']):.3e}（有无异常大值），每步 {r['t_step_ms']:.1f} ms", flush=True)
            json.dump(out, open(os.path.join(OUT, "v0_b2a_patch.json"), "w"))
