"""第 1.5 步 V3-2：换 a_eff = 5 mm，检验误差是否只取决于 h/a_eff。见 step1_5_results.md V3。

同 V3-1（lcp_tool.py，四分之一模型，线性小压深）。h = 10、5、2.5 mm（h/a = 2、1、0.5），s = h、h/2、h/4。
和 V3-1（a = 10 mm）在相同 h/a 下比较 k/k_S：如果一致，误差可以用一条 h/a 曲线描述。
注意：a = 5 mm 时立方体（8 cm）相对压头更大，有限尺寸修正比 a = 10 mm 小（0a-2），k_S 的可比性略有不同。
运行：python3 v3_2_aeff.py → results/step1_5/v3_2_aeff.json
"""
import sys
sys.dont_write_bytecode = True
import os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import lcp_tool as T

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../results/step1_5")
A = 0.005
K_S = 2 * A * T.E / (1 - T.NU ** 2)

if __name__ == "__main__":
    out = []
    for h in [0.010, 0.005, 0.0025]:
        a = time.perf_counter(); m = T.Mesh(h, quarter=True); t_mesh = time.perf_counter() - a
        print(f"== a = 5 mm，h = {h*1e3:.2f} mm（h/a = {h/A:.2f}），自由度 {m.ndof}，{t_mesh:.1f} s", flush=True)
        for s in [h, h / 2, h / 4]:
            r = T.disk_stiffness(m, A, s)
            r.update(h_mm=h * 1e3, s_mm=s * 1e3, h_over_a=h / A, s_over_h=s / h, k_over_kS=r["k"] / K_S)
            out.append(r)
            print(f"  s = {s*1e3:6.3f} mm（s/h = {s/h:.2f}）：k = {r['k']:9.4f} N/m，k/k_S = {r['k']/K_S:.4f}，"
                  f"受力点 {r['n_active']}/{r['n_points']}，互补残差 {r['complementarity']:.1e}", flush=True)
            json.dump(dict(k_S=K_S, a=A, rows=out), open(os.path.join(OUT, "v3_2_aeff.json"), "w"), indent=1)
