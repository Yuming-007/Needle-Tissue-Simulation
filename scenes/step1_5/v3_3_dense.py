"""第 1.5 步 V3-3：更密的 h 和更多的 a，检验"误差只取决于 h/a"并给出实测曲线。见 step1_5_results.md V3。

背景：V3-1 / V3-2 每个 a 只有 2–3 个网格点，不足以确定误差和 h/a 的函数形式（用户 2026-10-02）。
工具同 V3-1（lcp_tool.py，四分之一模型，线性小压深 LCP；V3-0 已验证）。
网格：h = 80 mm / N，N = 8、10、12、14、16、20、24、28、32（h = 10 … 2.5 mm；N 为偶数，四分之一模型的中心正好是节点）。
圆盘：a = 5、7.5、10、12.5 mm；点距 s = h（推荐规则）、h/2、h/4（平台检查，最后一档作为有限元极限的近似）。
输出：每个 (h, a, s) 的 k/k_S、受力点比例。只记录和画曲线，不预设函数形式。
运行：python3 v3_3_dense.py → results/step1_5/v3_3_dense.json
"""
import sys
sys.dont_write_bytecode = True
import os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import lcp_tool as T

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../results/step1_5")
AS = [0.005, 0.0075, 0.010, 0.0125]

if __name__ == "__main__":
    out = []
    for N in [8, 10, 12, 14, 16, 20, 24, 28, 32]:
        h = T.L / N
        a0 = time.perf_counter(); m = T.Mesh(h, quarter=True); t_mesh = time.perf_counter() - a0
        print(f"== h = {h*1e3:.3f} mm（N = {N}），自由度 {m.ndof}，{t_mesh:.1f} s", flush=True)
        for a in AS:
            kS = 2 * a * T.E / (1 - T.NU ** 2)
            line = []
            for s in [h, h / 2, h / 4]:
                r = T.disk_stiffness(m, a, s)
                r.update(h_mm=h * 1e3, a_mm=a * 1e3, s_mm=s * 1e3, h_over_a=h / a, s_over_h=s / h, k_S=kS, k_over_kS=r["k"] / kS)
                out.append(r); line.append(r)
            print(f"  a = {a*1e3:4.1f} mm（h/a = {h/a:.3f}）：k/k_S  s=h {line[0]['k_over_kS']:.4f}（{line[0]['n_active']}/{line[0]['n_points']}），"
                  f"s=h/2 {line[1]['k_over_kS']:.4f}，s=h/4 {line[2]['k_over_kS']:.4f}", flush=True)
            json.dump(out, open(os.path.join(OUT, "v3_3_dense.json"), "w"), indent=1)
        del m
