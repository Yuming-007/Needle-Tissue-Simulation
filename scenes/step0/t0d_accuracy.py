"""0d 第二阶段：R3（AsyncSparseLDLSolver，异步分解）相对 R1（SparseLDLSolver，每步分解）的精度检查。

A 准静态压头加载（0a-3 的场景：实时网格 n = 9，ν = 0.45，压头 a = 10 mm 压到 2 mm 再保持 1 s），
  v = 5 / 10 mm/s，dt = 0.01 / 0.02：比较压头力 F(t)，以及各自相对静力解的偏差。
B W2 负载（t0d_timing.build：21 行约束，针沿轴向 5 mm/s），针运动 0.5 s 后突然停止、再保持 0.5 s，
  n = 9（以及 n = 11），dt = 0.01 / 0.02：比较每步的约束力 λ 和嵌入点位置。
运行：python3 t0d_accuracy.py → results/step0/t0d_accuracy.json
"""
import sys
sys.dont_write_bytecode = True
import os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
import Sofa, Sofa.Simulation
import common
import t0a3_dynamics as A
import t0d_timing as T

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../results/step0/t0d_accuracy.json")


def test_A(delta_tab, F_tab):
    out = []
    for v in [5e-3, 10e-3]:
        for dt in [0.01, 0.02]:
            r = {lin: A.quasi_static(v, dt, 2e-3, 1.0, delta_tab, F_tab, linear=lin) for lin in ["ldl", "async"]}
            F1, F3, Fs = r["ldl"]["F"], r["async"]["F"], r["ldl"]["F_static"]
            nl = r["ldl"]["n_load"] - 1
            de = r["ldl"]["delta"]
            m = (np.arange(len(F1)) < nl) & (de >= 1e-3)
            d = dict(v=v, dt=dt,
                     max_diff_R3_R1_over_Fend=float(np.abs(F3 - F1).max() / F_tab[-1]),
                     seg_1_2mm_R1=float(np.abs(F1[m] / Fs[m] - 1).max()),
                     seg_1_2mm_R3=float(np.abs(F3[m] / Fs[m] - 1).max()),
                     final_R1=float(F1[-1] / F_tab[-1]), final_R3=float(F3[-1] / F_tab[-1]))
            out.append(d)
            print(f"A v={v*1e3:.0f}mm/s dt={dt}: max|F_R3-F_R1|/F(2mm) = {d['max_diff_R3_R1_over_Fend']:.2e}；"
                  f"δ=1-2mm 相对静力解：R1 {d['seg_1_2mm_R1']:.4f}，R3 {d['seg_1_2mm_R3']:.4f}；"
                  f"保持结束 F/F_static：R1 {d['final_R1']:.5f}，R3 {d['final_R3']:.5f}", flush=True)
    return out


def run_B(n, route, dt, t_move=0.5, t_hold=0.5):
    T.DT = dt
    root, n_pts = T.build(n, route, "W2")
    drv = root.driver
    v_cmd = 0.005
    base0 = np.asarray(root.Needle.base.position.value)[0][:3].copy()
    drv.traj = lambda tt: (base0 + np.array([0, 0, -v_cmd * min(tt, t_move)]), np.array([0, 0, 0, 1.0]))
    Sofa.Simulation.init(root)
    emo = root.TissueRoot.Tissue.Embedded.dofs
    lam, pos = [], []
    for k in range(int(round((t_move + t_hold) / dt))):
        Sofa.Simulation.animate(root, dt)
        lam.append(np.array(root.csolver.constraintForces.value) / dt)   # 力 = λ/dt
        pos.append(emo.position.array().copy())
    Sofa.Simulation.unload(root)
    T.DT = 0.01
    return np.array(lam), np.array(pos)


def test_B():
    out = []
    for n in [9, 11]:
        for dt in [0.01, 0.02]:
            l1, p1 = run_B(n, "R1", dt)
            l3, p3 = run_B(n, "R3", dt)
            scale = np.abs(l1).max()
            d = dict(n=n, dt=dt, max_force=float(scale),
                     max_lam_diff_rel=float(np.abs(l3 - l1).max() / scale),
                     final_lam_diff_rel=float(np.abs(l3[-1] - l1[-1]).max() / scale),
                     max_pos_diff=float(np.abs(p3 - p1).max()),
                     total_force_R1_end=float(l1[-1].reshape(-1, 3).sum(0)[2]),
                     total_force_R3_end=float(l3[-1].reshape(-1, 3).sum(0)[2]))
            out.append(d)
            print(f"B n={n} dt={dt}: 最大约束力 {scale:.4f} N；max|λ_R3-λ_R1|/max|λ| = {d['max_lam_diff_rel']:.2e}，"
                  f"保持结束时 {d['final_lam_diff_rel']:.2e}；嵌入点位置最大差 {d['max_pos_diff']:.2e} m；"
                  f"保持结束时 z 向合力 R1 {d['total_force_R1_end']:.5f} N，R3 {d['total_force_R3_end']:.5f} N", flush=True)
    return out


if __name__ == "__main__":
    dmax = 2e-3
    delta_tab = np.linspace(0, dmax, 9)
    F_tab = np.array([0.0] + [A.static_punch(d)[1] for d in delta_tab[1:]])
    out = dict(A=test_A(delta_tab, F_tab), B=test_B())
    json.dump(out, open(RES, "w"), indent=1)
