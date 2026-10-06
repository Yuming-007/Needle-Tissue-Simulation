"""第 1.5 步 M3 筛选（A+ / B 系列）：加密区是否够大、落点依赖。见 docs/plan/step1_5_results.md "工作网格 M3 筛选"。

工具：lcp_tool.MeshArrays（任意四面体网格、全模型、线性小压深 LCP；SOFA 静力切线刚度）。
  已验证：规则网格数组输入时和原工具完全相同（163.9330 N/m）；A+ 上和动态接触（M1）差 4.7e-4，受力点数和压力中心一致。
(a) 加密区敏感性：圆盘对准中心，a_eff = 10 mm，点距 s = 标称 h_local。
    A 系列（h_local = 5 mm）：A（R = D = 15）、A+（20）、A++（25），同尺度均匀参照 U5；
    B 系列（h_local = 4 mm）：B（15）、B+（20）、B++（25），同尺度均匀参照 U4。
    "相对均匀参照的额外偏硬" = k / k_U − 1：非均匀网格相对同一标称尺寸的均匀网格引入的额外偏差（不是总误差；
    U5、U4 本身相对 Sneddon 仍偏硬约 30%）。
(b) 落点依赖：A+、B、B+；圆盘中心在加密区内随机偏移 8 个位置（|偏移| ≤ h_local，固定随机种子）加中心，
    记录刚度的最大 / 最小 − 1 和压力中心偏移。判据（初定）：≤ 约 10%（V2 规则测试台为 2–4%）。
运行：python3 m3_screen.py → results/step1_5/m3_screen.json
"""
import sys
sys.dont_write_bytecode = True
import os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import common
import lcp_tool as T

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../results/step1_5")
A = 0.010
K_S = 2 * A * T.E / (1 - T.NU ** 2)
C = np.array([T.L / 2, T.L / 2])
SERIES = {"A": (["A", "A+", "A++"], "U5"), "B": (["B", "B+", "B++"], "U4")}
PLACEMENT = ["A+", "B", "B+"]


def mesh_of(name):
    m = common.load_mesh(name)
    return T.MeshArrays(m["X"], m["tets"], h=float(m["h_local"])), float(m["h_local"]), int(len(m["X"]))


if __name__ == "__main__":
    out = dict(region=[], placement=[])
    cache = {}
    for ser, (names, ref) in SERIES.items():
        rows = []
        for name in names + [ref]:
            a0 = time.perf_counter(); m, h, nn = mesh_of(name); t_m = time.perf_counter() - a0
            r = T.disk_stiffness_at(m, A, h, C)
            r.update(name=name, series=ser, nodes=nn, h_local_mm=h * 1e3, k_over_kS=r["k"] / K_S, t_mesh=t_m)
            rows.append(r)
            if name in PLACEMENT:
                cache[name] = (m, h)
        kref = [r for r in rows if r["name"] == ref][0]["k"]
        for r in rows:
            r["extra_vs_uniform"] = r["k"] / kref - 1
            out["region"].append(r)
            print(f"[区域] {r['name']:4s}（{r['nodes']:5d} 节点，h_local {r['h_local_mm']:.0f} mm）：k/k_S = {r['k_over_kS']:.4f}，"
                  f"相对 {ref} {r['extra_vs_uniform']:+.2%}，受力点 {r['n_active']}/{r['n_points']}", flush=True)
    rng = np.random.default_rng(20261003)
    for name in PLACEMENT:
        m, h = cache[name]
        offs = [np.zeros(2)]
        while len(offs) < 9:
            o = rng.uniform(-h, h, 2)
            if np.hypot(*o) <= h:
                offs.append(o)
        ks, cops = [], []
        for o in offs:
            r = T.disk_stiffness_at(m, A, h, C + o)
            ks.append(r["k"]); cops.append(np.hypot(*r["cop"]))
        ks = np.array(ks)
        d = dict(name=name, offsets_mm=(np.array(offs) * 1e3).tolist(), k=ks.tolist(), spread=float(ks.max() / ks.min() - 1),
                 k_center=float(ks[0]), cop_max_mm=float(max(cops) * 1e3))
        out["placement"].append(d)
        print(f"[落点] {name:4s}：9 个位置 k = {ks.min():.2f}–{ks.max():.2f} N/m，最大/最小 − 1 = {d['spread']:+.2%}，"
              f"压力中心最大偏移 {d['cop_max_mm']:.2f} mm", flush=True)
    json.dump(out, open(os.path.join(OUT, "m3_screen.json"), "w"), indent=1)
