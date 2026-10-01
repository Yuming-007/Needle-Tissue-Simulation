"""第 1 步 P8：顶点接触在实时配置（R3，AsyncSparseLDLSolver）下的每步耗时；附 R3 对 R1 的一致性检查。
针以 5 mm/s 往复（压深 0–10 mm），连续 2000 步，去掉前 20 步。
运行：python3 t1_timing.py → results/step1/t1_timing.json
"""
import sys
sys.dont_write_bytecode = True
import os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
import t1_vertex as V

if __name__ == "__main__":
    L, gap, v, dt = V.L, 0.001, 0.005, 0.01
    segs = [(gap / v, -v)] + [(0.010 / v, -v), (0.010 / v, +v)] * 5
    zf, T = V.piecewise(L + gap, segs)
    out = {}
    for lin in ["async", "ldl"]:
        r = V.run(dt, zf, T, linear=lin)
        ts = 1e3 * np.array(r["t_step"][20:])
        F = np.array(r["F_needle"])[:, 2]
        out[lin] = dict(median=float(np.median(ts)), p95=float(np.percentile(ts, 95)), max=float(ts.max()),
                        frac_over_dt=float((ts > 1e3 * dt).mean()), F=F.tolist())
        print(f"{lin}: 中位数 {np.median(ts):.2f} ms，p95 {np.percentile(ts,95):.2f} ms，最大 {ts.max():.2f} ms，"
              f"超过 dt 的比例 {100*(ts > 1e3*dt).mean():.1f}%", flush=True)
    Fa, Fl = np.array(out["async"]["F"]), np.array(out["ldl"]["F"])
    print(f"R3 vs R1 针根力：最大差 / 最大力 = {np.abs(Fa-Fl).max()/np.abs(Fl).max():.2e}", flush=True)
    out["R3_vs_R1"] = float(np.abs(Fa - Fl).max() / np.abs(Fl).max())
    json.dump(out, open(os.path.join(V.OUT, "t1_timing.json"), "w"))
