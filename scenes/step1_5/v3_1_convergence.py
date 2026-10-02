"""第 1.5 步 V3-1：平底圆盘（a_eff = 10 mm）的网格 / 布点收敛。见 step1_5_tip_survey.md §7、step1_5_results.md V3。

工具：lcp_tool.py（V3-0 已验证：和 V1 参照一致到 2.5e-6，四分之一模型 = 全模型到 1e-12）。线性、小压深。
网格：四分之一模型，h = 10、5、2.5（可选 1.25）mm。
布点：点距 s = 固定 10 mm（"固定点数"规则）、h（"点距约等于 h"规则）、h/2、h/4。
  - 固定 h 时，s 越细越接近"该网格上精确的有限元平底压头解"；看 h、h/2、h/4 是否形成平台，平台时最后一档作为
    有限元极限的近似 k_FE(h)。
  - 布点误差 = k(h, s) / k_FE(h) − 1；有限元误差 = k_FE(h) / k_S − 1。
参照：Sneddon 刚性平底压头 k_S = 2aE/(1−ν²)（0a-2：同一立方体的连续体结果约等于 k_S）。只作比较，不作约束；
拟合 k = k∞ + c·h（允许 k∞ ≠ k_S），4 个网格点时另试 k = k∞ + c·h^p。
运行：python3 v3_1_convergence.py [--fine] → results/step1_5/v3_1_convergence.json
"""
import sys
sys.dont_write_bytecode = True
import os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import lcp_tool as T

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../results/step1_5")
A = 0.010
K_S = 2 * A * T.E / (1 - T.NU ** 2)

if __name__ == "__main__":
    hs = [0.010, 0.005, 0.0025] + ([0.00125] if "--fine" in sys.argv else [])
    out = []
    for h in hs:
        a = time.perf_counter(); m = T.Mesh(h, quarter=True); t_mesh = time.perf_counter() - a
        print(f"== h = {h*1e3:.2f} mm（h/a = {h/A:.3f}），四分之一模型自由度 {m.ndof}，建网格 + 分解 {t_mesh:.1f} s", flush=True)
        for s in sorted({0.010, h, h / 2, h / 4}, reverse=True):
            a = time.perf_counter(); r = T.disk_stiffness(m, A, s); t = time.perf_counter() - a
            r.update(h_mm=h * 1e3, s_mm=s * 1e3, h_over_a=h / A, s_over_h=s / h, k_over_kS=r["k"] / K_S, t=t, ndof=m.ndof)
            out.append(r)
            print(f"  s = {s*1e3:6.3f} mm（s/h = {s/h:.2f}）：k = {r['k']:9.4f} N/m，k/k_S = {r['k']/K_S:.4f}，"
                  f"受力点 {r['n_active']}/{r['n_points']}，W 条件数 {r['W_cond']:.1e}，互补残差 {r['complementarity']:.1e}，{t:.1f} s", flush=True)
            json.dump(dict(k_S=K_S, a=A, rows=out), open(os.path.join(OUT, "v3_1_convergence.json"), "w"), indent=1)
        del m
