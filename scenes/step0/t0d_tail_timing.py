"""0d 补充：冻结配置的长时间计时，记录中位数、p95 和最大值（无界面，只有物理计算）。

冻结配置：n = 9，R3（AsyncSparseLDLSolver），dt = 0.01，W2（21 行约束）。针以 5 mm/s 往复运动（振幅 10 mm，
和 GUI 场景相同），连续 2000 步，去掉前 20 步。另外用 R1 跑同样的工况作对照。
运行：python3 t0d_tail_timing.py → results/step0/t0d_tail_timing.json
"""
import sys
sys.dont_write_bytecode = True
import os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np
import Sofa, Sofa.Simulation
import t0d_timing as T

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../results/step0/t0d_tail_timing.json")


def run(route, nsteps=2000, skip=20):
    root, n_pts = T.build(9, route, "W2")
    base0 = np.asarray(root.Needle.base.position.value)[0][:3].copy()
    amp, v = 0.010, 0.005

    def traj(t):
        period = 2 * amp / v; ph = t % period
        s = v * ph if ph < period / 2 else amp - v * (ph - period / 2)
        return base0 + np.array([0, 0, -s]), np.array([0, 0, 0, 1.0])
    root.driver.traj = traj
    Sofa.Simulation.init(root)
    ts = []
    for k in range(nsteps):
        a = time.perf_counter(); Sofa.Simulation.animate(root, T.DT); ts.append(time.perf_counter() - a)
    Sofa.Simulation.unload(root)
    ts = 1e3 * np.array(ts[skip:])
    return dict(route=route, steps=len(ts), median_ms=float(np.median(ts)), p95_ms=float(np.percentile(ts, 95)),
                p99_ms=float(np.percentile(ts, 99)), max_ms=float(ts.max()),
                frac_over_dt=float((ts > 1e3 * T.DT).mean()))


if __name__ == "__main__":
    out = []
    for route in ["R3", "R1"]:
        r = run(route); out.append(r)
        print(f"{route}: {r['steps']} 步，中位数 {r['median_ms']:.2f} ms，p95 {r['p95_ms']:.2f} ms，p99 {r['p99_ms']:.2f} ms，"
              f"最大 {r['max_ms']:.2f} ms，超过 dt 的比例 {100*r['frac_over_dt']:.1f}%", flush=True)
    json.dump(out, open(RES, "w"), indent=1)
